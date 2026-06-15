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


import zipfile
import itertools
import numpy as np, sys
from pyscf import lib, cc
from pyscf.tools.fcidump import to_scf

fid = 4
iiid = int(sys.argv[1])

lib.param.TMPDIR = '/tmp'

with zipfile.ZipFile('../FEMOCO.FCIDUMP.zip', 'r') as zip_ref:
    zip_ref.extract('FEMOCO.FCIDUMP', '/tmp/')

mf = to_scf('/tmp/FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf().newton()

n_sites = umf.mol.nao

for ix in range(iiid * 24, (iiid + 1) * 24):
    l = open('../00-UHF/02-fm-uhf-24.out', 'r').readlines()[ix].split()
    xdstr, ref_e = l[3], float(l[4])

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    umf.verbose = 3
    umf.max_cycle = 200
    umf.conv_tol = 1E-12
    umf.kernel(dm0=dm0)
    if abs(umf.e_tot - ref_e) >= 1E-8:
        print('!!WARNING', umf.e_tot, ref_e)
    umf.verbose = 4

    mmc = cc.CCSD(umf)
    mmc.diis_space = 36
    mmc.max_cycle = 200
    mmc.kernel()
    e_ccsd_t = mmc.e_tot + mmc.ccsd_t()

    print(":: %5d %s %20.12f %20.12f %d" % (ix, xdstr, mmc.e_tot, e_ccsd_t, int(mmc.converged)))
