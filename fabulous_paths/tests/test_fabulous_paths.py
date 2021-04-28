import pytest
from unittest import mock
from openpathsampling.tests.test_helpers import make_1d_traj
from openpathsampling.analysis.tis.core import \
    steps_to_weighted_trajectories
import pandas as pd
import numpy as np
import os

import pkg_resources

from fabulous_paths import *
from fabulous_paths import _get_keep_atoms, _concat_dfs

BACKBONE = [4, 5, 6, 8, 14, 15, 16, 18]

@pytest.fixture(scope="session")
def datafile():
    import openpathsampling as paths
    from openpathsampling.experimental import storage
    paths = storage.monkey_patch_all(paths)
    filename = pkg_resources.resource_filename(
        'fabulous_paths',
        os.path.join('tests', 'test_data', 'test_data.db')
    )
    if not os.path.exists(filename):  # pragma: no cover
        pytest.skip("Test file missing")
    datafile = storage.Storage(filename, mode='r')
    return datafile

@pytest.fixture
def trajectory_data(datafile):
    ensemble = datafile.ensembles['length 5']
    weighted_trajs = steps_to_weighted_trajectories(datafile.steps,
                                                    ensembles=[ensemble])
    traj_counter = list(weighted_trajs.values())[0]
    # we run 6 trajectories, plus initial
    assert sum(traj_counter.values()) == 7
    trajectories = list(traj_counter.keys())
    for traj in trajectories:
        assert len(traj) == 5  # each trajectory is 5 frames
    return trajectories

def test_extract_CV(datafile, trajectory_data):
    phi = datafile.cvs['phi']
    psi = datafile.cvs['psi']
    cvs = [phi, psi]
    for traj in trajectory_data:
        extracted = extract_CV(traj, cvs)
        assert extracted.shape == (len(traj), len(cvs))
        assert all(extracted.columns == ['phi', 'psi'])

def test_extract_MD(trajectory_data):
    for traj in trajectory_data:
        # use the first frame as ref; gives us a no-change to test
        ref = traj.to_mdtraj()
        extracted = extract_MD(traj, ref[0], BACKBONE)
        n_cols_expected = len(BACKBONE) * 3
        assert extracted.shape == (len(traj), n_cols_expected)
        first_frame = extracted.loc[0,:]
        sliced_ref = ref.atom_slice(BACKBONE)
        np.testing.assert_allclose(first_frame, sliced_ref.xyz[0].flatten(),
                                   rtol=1e-6)  # seems we need extra room

@pytest.mark.parametrize('type_', ['str', 'list'])
def test_get_keep_atoms(trajectory_data, type_):
    keep = {'str': 'backbone', 'list': BACKBONE}[type_]
    topology = trajectory_data[0][0].topology.mdtraj
    np.testing.assert_array_equal(_get_keep_atoms(topology, keep), BACKBONE)

@mock.patch('fabulous_paths.tqdm', new=lambda x, desc: x)
def test_extract_OPS(datafile):
    ref = datafile.steps[0].active[0].trajectory.to_mdtraj()[0]
    keep = BACKBONE
    cvs = [datafile.cvs['phi'], datafile.cvs['psi']]
    dfs = list(extract_OPS(datafile.steps, ref, keep, cvs))
    assert len(dfs) == len(datafile.steps)
    for trj_df, cv_df in dfs:
        # 5 == len(traj)
        # 2 == len(cvs)
        assert trj_df.shape == (5, len(BACKBONE) * 3)
        assert cv_df.shape == (5, 2)

def test_extract_OPS_bad_ensembles():
    # this can't be tested yet, because it is (temporarily) impossible
    pytest.skip()

def test_concat_dfs():
    sizes = [(4, 2), (4, 3)]
    dfs = [tuple(pd.DataFrame(np.random.random(size)) for size in sizes)
           for _ in range(2)]
    # sanity check
    for df_pair in dfs:
        for df, shape in zip(df_pair, sizes):
            assert df.shape == shape

    trajs, cvs = _concat_dfs(dfs)
    assert trajs.shape == (8, 2)
    assert cvs.shape == (8, 3)

def test_main_integration(datafile, tmpdir):
    # smoke test to ensure that the integration with FABULOUS works
    ref = datafile.steps[0].active[0].trajectory.to_mdtraj()[0]
    keep = BACKBONE
    cvs = [datafile.cvs['phi'], datafile.cvs['psi']]
    results_dir=tmpdir / "results"
    label='1'
    yaml_file=... # TODO
    main(steps=datafile.steps, cvs=cvs, ref=ref, keep=keep,
         yaml_file=yaml_file, n_gen=0, results_dir=results_dir, label=label)
