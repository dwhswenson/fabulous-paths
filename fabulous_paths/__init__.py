from tqdm.auto import tqdm
import pandas as pd
import numpy as np

def extract_CV(traj, cvs):
    """Create a CV dataframe for FABULOUS from an OPS trajectory and CVs.

    Returns
    -------
    pandas.DataFrame :
        each row is a snapshot
    """
    columns = [cv.name for cv in cvs]
    results = np.array([cv(traj) for cv in cvs]).T
    return pd.DataFrame(results, columns=columns)


def extract_MD(trajectory, ref, keep_atoms):
    """Create an MD dataframe for FABULOUS from an OPS trajectory.

    Parameters
    ----------
    trajectory : :class:`openpathsampling.Trajectory`
    ref : :class:`mdtraj.Trajectory`
    keep_atoms : List[int]

    Returns
    -------
    pandas.DataFrame :
        each row is a snapshot, all coordinates in the columns
    """
    traj = trajectory.to_mdtraj().atom_slice(keep_atoms)
    traj.superpose(ref.atom_slice(keep_atoms))  # TODO: move this out

    xyz = traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3)
    return pd.DataFrame(xyz)


def _get_keep_atoms(topology, keep):
    if isinstance(keep, str):
        return topology.select(keep)
    else:
        return [int(i) for i in keep]


def extract_OPS(steps, ref, keep, cvs):
    """Full process of extracting CV and trajectory dataframes fom OPS.

    Parameters
    ----------
    steps: List[:class:`openpathsampling.MCStep`]
        input steps
    """
    from openpathsampling.analysis.tis.core import \
        steps_to_weighted_trajectories
    ensembles = [steps[0].active[0].ensemble]  # TODO: remove later
    keep_atoms = None
    desc1 = "Identifying trajectories"
    weighted_trajs = steps_to_weighted_trajectories(tqdm(steps, desc=desc1),
                                                    ensembles)
    if len(weighted_trajs) != 1:
        # NOTE: this is currently impossible, but might become possible when
        # OPS fixes steps_to_weighted to not require ensembles as input
        raise RuntimeError("FABULOUS requires all steps from the same",
                           f"ensemble. Found {len(weighted_trajs)} ",
                           "ensembles.")

    traj_counter = list(weighted_trajs.values())[0]

    desc2 = "Extracting data for FABULOUS"
    for traj in tqdm(traj_counter, desc=desc2):
        if keep_atoms is None:
            # this only runs on the first trajectory
            topology = traj[0].topology.mdtraj
            keep_atoms = _get_keep_atoms(topology, keep)

        trj_df = extract_MD(traj, ref, keep_atoms)
        cv_df = extract_CV(traj, cvs)
        yield trj_df, cv_df

def _concat_dfs(dfs):
    traj_dfs, cv_dfs = zip(*dfs)
    trajs = pd.concat(traj_dfs, ignore_index=True)
    cvs = pd.concat(cv_dfs, ignore_index=True)
    return trajs, cvs

def main(steps, cvs, ref, keep, yaml_file, n_gen, results_dir, label):
    dfs = extract_OPS(steps, ref, keep, cvs)
    trajs, cvs = _concat_dfs(dfs)
    print("Trajectories:")
    print(trajs.head())
    print("CVs:")
    print(cvs.head())
    ... # from here there should be something in FABULOUS that does the work
    pass
