"""Regressions for bounded HDF reads and separately stored snapshots.

Run with the RUR checkout on PYTHONPATH:
    python -m unittest discover -s tests -p 'test_hdf_compatibility.py'
"""
import os
import tempfile
import unittest
from unittest import mock

import h5py
import numpy as np

from rur import uri


class BoundedArray:
    """Reject accidental whole-dataset reads in a partial-chunk request."""
    def __init__(self, data):
        self.data = data
        self.dtype = data.dtype
        self.shape = data.shape
        self.size = data.size

    def __getitem__(self, key):
        if isinstance(key, slice) and key.start is None and key.stop is None:
            raise AssertionError("Partial request attempted to read every row")
        return self.data[key]


class HdfCompatibilityTests(unittest.TestCase):
    def setUp(self):
        uri.timer.verbose = 0
        self.rows = np.zeros(10, dtype=[('id', 'i4'), ('x', 'f8')])
        self.rows['id'] = np.arange(10)
        self.boundaries = np.array([0, 2, 10])

    def test_first_chunk_reads_only_selected_rows(self):
        result = uri.domain_slice(BoundedArray(self.rows), np.array([0]), self.boundaries)
        np.testing.assert_array_equal(result, self.rows[:2])

    def test_one_chunk_parallel_request_stays_bounded(self):
        with tempfile.TemporaryDirectory() as root:
            path = os.path.join(root, 'part.h5')
            with h5py.File(path, 'w') as handle:
                handle.create_dataset('star/data', data=self.rows)
            original = h5py.Dataset.__getitem__
            def checked(dataset, key, *args, **kwargs):
                if (dataset.name == '/star/data' and isinstance(key, slice)
                        and key.start is None and key.stop is None):
                    raise AssertionError('Parallel fallback attempted a whole-file read')
                return original(dataset, key, *args, **kwargs)
            with mock.patch.object(h5py.Dataset, '__getitem__', checked):
                result = uri.domain_slice(None, np.array([0]), self.boundaries,
                                          nthread=4, hdf_path=path, hdf_group='star')
            np.testing.assert_array_equal(result, self.rows[:2])

    def test_empty_chunk_returns_typed_empty_result(self):
        result = uri.domain_slice(self.rows, np.array([0]), np.array([0, 0, 10]),
                                  nthread=4, target_fields=['id'])
        self.assertEqual(result.size, 0)
        self.assertEqual(result.dtype.names, ('id',))

    def test_full_volume_request_keeps_all_rows(self):
        result = uri.domain_slice(self.rows, np.array([0, 1]), self.boundaries)
        np.testing.assert_array_equal(result, self.rows)

    def test_native_fields_compound_keys_and_separate_storage(self):
        with tempfile.TemporaryDirectory() as root:
            repo = os.path.join(root, 'catalogue_repository')
            storage = os.path.join(root, 'snapshot_storage')
            os.makedirs(repo)
            os.makedirs(storage)
            data = np.zeros(4, dtype=[('position_x', 'f8'), ('position_y', 'f8'),
                                     ('position_z', 'f8'), ('mass', 'f4'),
                                     ('initial_mass', 'f4'), ('identity', 'i4'),
                                     ('birth_time', 'f4'), ('metallicity', 'f4'),
                                     ('levelp', 'i1')])
            for axis in 'xyz':
                data['position_' + axis] = [.1, .1001, .7, .9]
            data['mass'] = 1e-12
            data['initial_mass'] = 2e-12
            data['identity'] = [11, 12, 13, 14]
            data['birth_time'] = -2
            data['metallicity'] = .01
            data['levelp'] = 8
            coords = np.column_stack([data['position_' + axis] for axis in 'xyz'])
            keys = uri.utool.hilbert3d_map(coords, 21)
            order = np.argsort(keys)
            data = data[order]
            keys = keys[order]
            compound = np.zeros(3, dtype=[('hi', 'u8'), ('lo', 'u8')])
            compound['lo'] = [0, int(keys[2]), 2 ** 63]
            attrs = dict(H0=70., omega_m=.3, omega_l=.7, omega_k=0., omega_b=.045,
                         unit_l=3e24, unit_d=1e-27, unit_t=1e17, unit_m=2.7e46, boxlen=1.,
                         time=-1., aexp=.5, z=1., icoarse=1, ncpu=1, ndim=3,
                         levelmin=8, levelmax=21, iout=1)
            # Part-only metadata is deliberately supported: no raw files or cell file.
            for iout in [1, 2]:
                with h5py.File(os.path.join(storage, 'part_%05d.h5' % iout), 'w') as handle:
                    for key, value in attrs.items():
                        handle.attrs[key] = value
                    handle.attrs['iout'] = iout
                    group = handle.create_group('star')
                    group.attrs['levelmax'] = 21
                    group.create_dataset('data', data=data)
                    group.create_dataset('chunk_boundary', data=np.array([0, 2, 4]))
                    group.create_dataset('hilbert_boundary', data=compound)
            snapshot = uri.RamsesSnapshot(repo, 1, mode='nc', hdf=True,
                                          hdf_path=storage, verbose=0)
            self.assertEqual(snapshot.repo, repo)
            np.testing.assert_array_equal(snapshot.get_iout_avail(), [1, 2])
            box = np.array([[.1 - 1e-8, .1 + 1e-8]] * 3)
            for nthread in [1, 4]:
                stars = snapshot.get_part(box=box, pname='star', nthread=nthread,
                                          target_fields=['id', 'x', 'y', 'z', 'm', 'm0', 'epoch'])
                self.assertEqual(stars['id'].tolist(), [11])
                self.assertTrue(np.isfinite(stars['age', 'Gyr']).all())
                self.assertIs(snapshot.part, stars)
                self.assertIs(snapshot.part_data, stars.table)
            switched = snapshot.switch_iout(2)
            self.assertTrue(switched.hdf)
            self.assertEqual(switched.hdf_path, storage)
            self.assertEqual(switched.repo, repo)

            # Julia passes NumPy arrays for field selection, including gas reads.
            extra = [('level', 'i4'), ('density', 'f4'), ('pressure', 'f4'),
                     ('chem_H', 'f4'), ('dust_bin01', 'f4'), ('dust_bin02', 'f4'),
                     ('dust_bin03', 'f4'), ('dust_bin04', 'f4')]
            cells = np.zeros(len(data), dtype=data.dtype.descr + extra)
            for field in data.dtype.names:
                cells[field] = data[field]
            cells['level'] = 8
            cells['density'] = 2
            cells['pressure'] = 3
            with h5py.File(os.path.join(storage, 'cell_00001.h5'), 'w') as handle:
                for key, value in attrs.items():
                    handle.attrs[key] = value
                group = handle.create_group('leaf')
                group.attrs['levelmax'] = 21
                group.create_dataset('data', data=cells)
                group.create_dataset('chunk_boundary', data=np.array([0, 2, 4]))
                group.create_dataset('hilbert_boundary', data=compound)
            gas = snapshot.get_cell(box=box, target_fields=np.array(
                ['x', 'y', 'z', 'rho', 'P', 'level', 'H', 'd1', 'd2', 'd3', 'd4']))
            self.assertIs(snapshot.cell, gas)
            self.assertIs(snapshot.cell_data, gas.table)
            self.assertGreater(len(gas), 0)
            self.assertTrue(np.isfinite(gas['T', 'K']).all())


if __name__ == '__main__':
    unittest.main()
