import openpathsampling as paths
import openmmtools
from simtk import unit
from simtk import openmm as mm
import mdtraj as md
import numpy as np
from paths_cli.commands.visit_all import visit_all_main
from paths_cli.commands.equilibrate import equilibrate_main

### USE SIMSTORE (NEW STORAGE) #############################################
from openpathsampling.experimental.storage import Storage, monkey_patch_all
from openpathsampling.experimental.storage.collective_variables import \
        MDTrajFunctionCV, CoordinateFunctionCV
paths = monkey_patch_all(paths)

### OPENMM SETUP ###########################################################
print("Setting up the OpenMM engine...")
INPUT_PDB = 'c7ax_input.pdb'

pdb = mm.app.PDBFile(INPUT_PDB)
forcefield = mm.app.ForceField('amber99sbildn.xml', 'tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=mm.app.PME,
    nonbondedCutoff=1.0*unit.nanometers, constraints=mm.app.HBonds,
    rigidWater=True, ewaldErrorTolerance=0.0005)

integrator = openmmtools.integrators.VVVRIntegrator(
    temperature=300.0*unit.kelvin,
    collision_rate=1.0/unit.picosecond,
    timestep=2.0*unit.femtosecond
)
integrator.setConstraintTolerance(0.00001)


### OPENPATHSAMPLING SETUP #################################################
print("Setting up the path sampling...")

##### engine ###############################################################
topology = paths.engines.MDTrajTopology(md.load(INPUT_PDB).topology)
init_snap = paths.engines.openmm.tools.ops_load_trajectory(INPUT_PDB)[0]

engine = paths.engines.openmm.Engine(
    topology=topology,
    system=system,
    integrator=integrator,
    options={'n_steps_per_frame': 10,
             'n_frames_max': 10000}
).named('300K')


##### collective variables #################################################
# here we use a trick that we calculate *all* dihedrals in one go, which is,
# in principle, faster that doing each one separately. This is particularly
# useful when calculating multiple CVs.
dihedral_indices = [
    [4, 6, 8, 14],  # phi
    [6, 8, 14, 16],  # phi
    [8, 14, 16, 18],  # omega
    [1, 4, 6, 8],  # theta
]
angle_indices = [
    [6, 8, 10], # NCaR
    [1, 8, 16], # end_Ca_end
]
distance_indices = [
    [5, 15],  # dOO
]

# when writing a function to wrap as an MDTrajFunctionCV, write it as if it
# were taking an mdtraj.Trajectory and returning the appropriate numpy array
def all_cvs(traj):
    import mdtraj as md
    import numpy as np
    dihedrals = md.compute_dihedrals(traj, indices=dihedral_indices)
    angles = md.compute_angles(traj, angle_indices=angle_indices)
    dists = md.compute_distances(traj, atom_pairs=distance_indices)
    stacked = np.stack([dihedrals, angles, dists])
    return stacked

dihedrals = MDTrajFunctionCV(md.compute_dihedrals,
                             topology=topology,
                             indices=dihedral_indices).named('dihedrals')

angles = MDTrajFunctionCV(md.compute_angles,
                          topology=topology,
                          angle_indices=angle_indices).named('angles')

helper_cvs = [dihedrals, angles]

def getitem(snapshot, cv, num):
    return cv(snapshot)[num]

phi = CoordinateFunctionCV(getitem, cv=dihedrals, num=0).named('phi')
psi = CoordinateFunctionCV(getitem, cv=dihedrals, num=1).named('psi')
omega = CoordinateFunctionCV(getitem, cv=dihedrals, num=2).named('omega')
theta = CoordinateFunctionCV(getitem, cv=dihedrals, num=3).named('theta')

NCaR = CoordinateFunctionCV(getitem, cv=angles, num=0).named("NCaR")
end_Ca_end = CoordinateFunctionCV(getitem, cv=angles,
                                  num=1).named("end_Ca_end")

dOO = MDTrajFunctionCV(md.compute_distances, topology=topology,
                       atom_pairs=distance_indices).named('dOO')

all_cvs = [phi, psi, omega, theta, NCaR, end_Ca_end, dOO] + helper_cvs


# we don't need to save these to disk, because they're already saved as part
# of `dihedrals` or `angles`
for cv in [phi, psi, omega, theta, NCaR, end_Ca_end]:
    cv.diskcache_enabled = False


##### TPS details: states, networks, move scheme ###########################
C7eq = (
    paths.CVDefinedVolume(phi, lambda_min=-np.pi, lambda_max=-1)
    & paths.CVDefinedVolume(psi, lambda_min=0, lambda_max=np.pi)
).named("C7eq")

C7ax = (
    paths.CVDefinedVolume(phi, lambda_min=0.75, lambda_max=1.75)
    & paths.CVDefinedVolume(psi, lambda_min=-np.pi, lambda_max=0.5)
).named("C7ax")

network = paths.TPSNetwork(C7eq, C7ax).named("tps-network")
scheme = paths.OneWayShootingMoveScheme(network,
                                        engine=engine).named("one-way")


### INITIAL TRAJECTORY ###################################################
print("Getting an initial trajectory from high temperature MD")
hi_temp_integrator = openmmtools.integrators.VVVRIntegrator(
    temperature=500.0*unit.kelvin,
    collision_rate=1.0/unit.picosecond,
    timestep=2.0*unit.femtosecond
)
hi_temp_integrator.setConstraintTolerance(0.00001)

hi_temp_engine = paths.engines.openmm.Engine(
    topology=topology,
    system=system,
    integrator=hi_temp_integrator,
    options={'n_steps_per_frame': 10,
             'n_frames_max': 100_000}
).named('500K')

# turn off all caching during trajectory generation
old_modes = {}
for cv in all_cvs:
    old_modes[cv] = cv.mode
    cv.mode = 'no-caching'

# we'll use the convenience functions from the CLI and chain them together
traj, _ = visit_all_main(
    output_storage=None,
    states=[C7eq, C7ax],
    engine=hi_temp_engine,
    initial_frame=init_snap
)

print("Equilibrating the high-temperature MD")
initial_conditions, _ = equilibrate_main(
    output_storage=None,
    scheme=scheme,
    init_conds=traj,
    multiplier=1,
    extra_steps=0
)

# reset modes
for cv in all_cvs:
    cv.mode = old_modes[cv]
    print(cv.name, cv.mode)

### SAVE EVERYTHING TO DISK ##############################################
print("Saving everything to the file ad_setup.db")
storage = Storage('ad_setup.db', mode='w')
storage.save(scheme)
storage.tags['initial_conditions'] = initial_conditions
# phi and psi were already saved in the scheme, but it doesn't hurt to save
# them again
storage.save([phi, psi, omega, theta, NCaR, end_Ca_end, dOO])
storage.close()
print("DONE!")
