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
from fabulous_paths import _get_keep_atoms

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
        extracted = extract_MD(traj, ref, BACKBONE)
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

def test_extract_OPS():
    pytest.skip()

def test_extract_OPS_bad_ensembles():
    pytest.skip()

def test_main_integration():
    pytest.skip()
