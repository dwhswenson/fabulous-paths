import pytest
from unittest import mock
from openpathsampling.tests.test_helpers import make_1d_traj
import pandas as pd

import pkg_resources

from paths_fabulous import *

BACKBONE = []

@pytest.fixture(scope="session")
def datafile():
    import openpathsampling as paths
    from openpathsampling.experimental import storage
    paths = storage.monkey_patch_all(paths)
    filename = pkg_resources.resource_filename(
        'paths_fabulous',
        os.path.join('tests', 'test_data', 'test_data.db')
    )
    datafile = storage.Storage(filename, mode='r')
    return datafile

@pytest.fixture
def trajectory_data(datafile):
    weighted_tajs = steps_to_weighted_trajectories(datafile.steps)
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
        assert extracted.shape == (len(cvs), len(traj))
        assert extracted.columns = ['phi', 'psi']
        # TODO: more asserts

    pytest.skip()

def test_extract_MD(trajectory_data):
    for traj in trajectory_data:
        ref = traj.to_mdtraj().atom_slice(BACKBONE)[0]
        extracted = extract_MD(traj, ref, BACKBONE)
        n_cols_expected = len(BACKBONE) * 3
        assert extracted.shape == (len(traj), n_cols_expected)
        np.testing.assert_all_close(extracted[0], np.zeros(n_cols_expected))

@pytest.mark.parametrize('type_', ['str', 'list'])
def test_get_keep_atoms(trajectory_data, type_):
    keep = {'str': 'backbone', 'list': BACKBONE}[type_]
    topology = trajectory_data[0][0].topology.mdtraj_topology
    assert _get_keep_atoms(topology, keep) == backbone

def test_extract_OPS():
    pytest.skip()

def test_main_integration():
    pytest.skip()
