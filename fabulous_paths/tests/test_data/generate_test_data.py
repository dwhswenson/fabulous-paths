"""
Script to generate test data for fabulous-paths!

Usage:
    python generate_test_data.py

Creates a file called test_data.db with the brief tests (fake TPS with
6 steps of one-way shooting in a fixed length 5-frame ensemble -- all paths
will be accepted).

This file is used in the test suite.
"""

import mdtraj as md
import openmmtools
import openpathsampling as paths


### SIMSTORE SETUP #########################################################
from openpathsampling.experimental.storage import Storage, monkey_patch_all
from openpathsampling.experimental.storage.collective_variables import \
        MDTrajFunctionCV
paths = monkey_patch_all(paths)
# don't need to worry about interface sets here


### MAIN SETUP #############################################################
print("Setting up the simulation")
testsystem = openmmtools.testsystems.AlanineDipeptideVacuum()
integrator = openmmtools.integrators.VVVRIntegrator()
topology = paths.engines.openmm.topology.MDTrajTopology(testsystem.mdtraj_topology)

engine = paths.engines.openmm.Engine(topology=topology,
                                     system=testsystem.system,
                                     integrator=integrator,
                                     options={
                                         'n_frames_max': 100,
                                         'n_steps_per_frame': 10,
                                     }).named('fake_engine')

ensemble=paths.LengthEnsemble(5).named("length 5")
mover = paths.OneWayShootingMover(ensemble=ensemble,
                                  selector=paths.UniformSelector(),
                                  engine=engine)

scheme = paths.LockedMoveScheme(root_mover=mover)


### INITIAL SNAPSHOT #######################################################
print("Getting energy-minimized initial snapshot")
engine.simulation.context.setPositions(testsystem.positions)
engine.simulation.minimizeEnergy()
initial_snapshot = engine.current_snapshot


### SETUP INITIAL CONDITIONS ###############################################
print("Setting up initial conditions")
initial_trajectory = engine.generate(initial_snapshot, ensemble.can_append)
sample = paths.Sample(replica=0,
                      ensemble=ensemble,
                      trajectory=initial_trajectory)
initial_conditions = paths.SampleSet([sample])

ref_frame = initial_trajectory.to_mdtraj()[0]
ref_frame.save("ref.pdb")

### RUN PATH SAMPLING ######################################################
storage = Storage('test_data.db', mode='w')
# add CVs
phi = MDTrajFunctionCV(md.compute_dihedrals,
                       topology=engine.topology,
                       indices=[[4, 6, 8, 14]]).named("phi")

psi = MDTrajFunctionCV(md.compute_dihedrals,
                       topology=engine.topology,
                       indices=[[6, 8, 14, 16]]).named("psi")

storage.save([phi, psi])

sim = paths.PathSampling(
    storage=storage,
    move_scheme=scheme,
    sample_set=initial_conditions
)
sim.run(6)
storage.close()
