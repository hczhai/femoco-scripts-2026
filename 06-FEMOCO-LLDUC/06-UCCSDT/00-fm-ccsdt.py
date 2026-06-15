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
import numpy as np, sys
from pyscf import lib
from pyscf.tools.fcidump import to_scf

sys.path[:0] = ['/mnt/home/hzhai/Projects/hast']

import pyhast

fid = 0
iiid = int(sys.argv[1])
jjjd = int(sys.argv[2])

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

# for ix in range(iiid * 24, (iiid + 1) * 24):
for ix in [iiid * 24 + jjjd]:
    l = open('02-fm-uhf-24.out', 'r').readlines()[ix].split()
    xdstr, ref_e = l[3], float(l[4])

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    umf.verbose = 4
    umf.max_cycle = 400
    umf.conv_tol = 1E-13
    umf.ah_conv_tol = 1E-16
    umf.kernel(dm0=dm0)
    if abs(umf.e_tot - ref_e) >= 1E-8:
        print('!!WARNING', umf.e_tot, ref_e)

    print('INIT =', xdstr)
    print('EHF  =', umf.e_tot)

    cc_opts = dict(symm_amps=True, label_symm=True, lift=True, sum_lift=True, split=True, eval_t='hastctr')

    mcc = pyhast.sr.ucc.UCC(umf, t_order=2, verbose=4, diis=True, **cc_opts)
    mcc.diis_space = 36
    mcc.kernel(max_cycle=1600)

    mcc2 = pyhast.sr.ucc.UCC(umf, t_order=3, verbose=4, diis=True, **cc_opts)
    mcc2.diis_space = 36
    tamps = mcc2.init_amps()
    tamps[:2] = mcc.tamps
    mcc2.kernel(tamps=tamps, max_cycle=1600)

    print(":: %5d %s %20.12f %20.12f %d" % (ix, xdstr, mcc.e_tot, mcc2.e_tot, int(mcc2.converged)))
