# ---------------------------------------------------------------------------
# Copyright (C) 2026 Huanchen Zhai. All Rights Reserved.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ---------------------------------------------------------------------------


import itertools
import numpy as np, sys
from pyscf import lib
from pyscf.tools.fcidump import to_scf

fid = 0
iiid = int(sys.argv[1])

lib.param.TMPDIR = '/tmp'

mf = to_scf('../Fe4S4.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf().newton()

n_sites = umf.mol.nao

fe_idxs = [
    [2, 3, 4, 5, 6], [7, 8, 9, 10, 11], [24, 25, 26, 27, 28], [29, 30, 31, 32, 33],
    [0, 1, *range(12, 24), 34, 35]
]

occs = []

for ix, x in enumerate([x for x in itertools.combinations(range(4), 2) if x[0] == 0]):
    y = tuple(z for z in range(4) if z not in x)
    xoccs = []
    for kx, ky in [(kx, ky) for kx in range(2) for ky in range(2)]:
        occ = [None] * 4 + ['2' * 16]
        for i in x:
            occ[i] = 'aaaaa'
        for i in y:
            occ[i] = 'bbbbb'
        for fx in ['2aaaa', 'a2aaa', 'aa2aa', 'aaa2a', 'aaaa2']:
            for fy in ['2bbbb', 'b2bbb', 'bb2bb', 'bbb2b', 'bbbb2']:
                occ[x[kx]], occ[y[ky]] = fx, fy
                xidets = [(ix, xx) for ixx, xxx in zip(fe_idxs, occ) for ix, xx in zip(ixx, xxx)]
                xdstr = ''.join([z for _, z in sorted(xidets, key=lambda x: x[0])])
                xoccs.append(xdstr)
    occs.append(xoccs)

print('INIT LEN =', [len(x) for x in occs])

data = []

for ix, xdstr in enumerate(occs[iiid]):

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    umf.verbose = 4
    umf.max_cycle = 200
    umf.conv_tol = 1E-13
    umf.ah_conv_tol = 1E-16
    umf.kernel(dm0=dm0)
    dm1 = np.array(umf.make_rdm1())

    found, ifx = -1, len(data)
    for jx, (kx, xener, xdm1) in enumerate(data):
        if abs(xener - umf.e_tot) < 1E-8 and np.linalg.norm(dm1 - xdm1) < 1E-4:
            found, ifx = kx, jx
            break
    if found == -1:
        data.append((ix, umf.e_tot, dm1))

    print(":: %5d %s %20.12f %5d %5d" % (ix, xdstr, umf.e_tot, found, ifx))
