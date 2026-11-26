"""Compare bounded real RAMSES/HDF regions without writing simulation storage."""
import argparse
import json
import os
import time

import h5py
import numpy as np

from rur import uri


def check(raw_repo, hdf_dir, iout):
    uri.timer.verbose = 0
    results = []
    for group, kind in [('star', 'part'), ('leaf', 'cell')]:
        with h5py.File(os.path.join(hdf_dir, '%s_%05d.h5' % (kind, iout)), 'r') as handle:
            first = handle[group]['data'][0]
            center = np.array([first['position_' + axis] for axis in 'xyz'])
        box = np.stack([center - 1e-8, center + 1e-8], axis=1)
        raw = uri.RamsesSnapshot(raw_repo, iout, mode='nc', box=box, verbose=0)
        cpus = raw.get_involved_cpu()
        types = ['part'] if kind == 'part' else ['amr', 'hydro']
        size = sum(os.path.getsize(raw.get_path(typ, int(cpu))) for typ in types for cpu in cpus)
        if len(cpus) > 4 or size > 2 * 1024 ** 3:
            raise RuntimeError('Probe would exceed bounded raw-I/O limit: %d files bytes' % size)
        fields = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'metal', 'H', 'O', 'Fe', 'Mg', 'C', 'N', 'Si', 'S', 'D']
        fields += ['id', 'm', 'm0', 'epoch', 'rho0'] if kind == 'part' else ['rho', 'P', 'level', 'd1', 'd2', 'd3', 'd4']
        start = time.monotonic()
        method = raw.get_part if kind == 'part' else raw.get_cell
        options = dict(pname='star') if kind == 'part' else {}
        before = method(target_fields=fields, nthread=1, hdf=False, **options)
        reference = before.table.copy()
        raw.clear(verbose=0)
        hdf = uri.RamsesSnapshot(raw_repo, iout, mode='nc', box=box,
                               hdf=True, hdf_path=hdf_dir, verbose=0)
        method = hdf.get_part if kind == 'part' else hdf.get_cell
        after = method(target_fields=fields, nthread=1, **options)
        actual = after.table
        assert len(reference) == len(actual) and len(actual) > 0, (len(reference), len(actual))
        order_fields = ['id'] if kind == 'part' else ['x', 'y', 'z', 'level']
        reference = np.sort(reference, order=order_fields)
        actual = np.sort(actual, order=order_fields)
        relative_errors = {}
        for field in fields:
            left, right = reference[field], actual[field]
            if np.issubdtype(left.dtype, np.integer):
                np.testing.assert_array_equal(left, right, err_msg=field)
            else:
                np.testing.assert_allclose(left, right, rtol=2e-6, atol=0, err_msg=field)
                nonzero = left != 0
                relative_errors[field] = float(np.max(np.abs((left[nonzero] - right[nonzero]) / left[nonzero]))) if nonzero.any() else 0.
        record = dict(group=group, rows=len(actual), selected_raw_bytes=size,
                      cpus=cpus.tolist(), elapsed_seconds=time.monotonic() - start,
                      maximum_relative_errors=relative_errors)
        results.append(record)
        print(json.dumps(record), flush=True)
        hdf.clear(verbose=0)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-repo', default='/storage8/NewCluster')
    parser.add_argument('--hdf-dir', default='/storage8/NewCluster/hdf')
    parser.add_argument('--iout', type=int, default=1010)
    parser.add_argument('--output', help='Optional JSON validation report outside simulation storage')
    args = parser.parse_args()
    results = check(args.raw_repo, args.hdf_dir, args.iout)
    if args.output:
        with open(args.output, 'w') as handle:
            json.dump(results, handle, indent=2)
            handle.write('\n')
