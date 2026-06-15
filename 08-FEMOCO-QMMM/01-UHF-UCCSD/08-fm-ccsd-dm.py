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
import numpy as np, sys, tempfile
from pyscf import lib, cc
from pyscf.tools.fcidump import to_scf

fid = 8
iiid = int(sys.argv[1])

lib.param.TMPDIR = '/tmp'

tmpdir = tempfile.mkdtemp()

with zipfile.ZipFile('data/FEMOCO.QMMM.FCIDUMP.zip', 'r') as zip_ref:
    zip_ref.extract('FEMOCO.QMMM.FCIDUMP', tmpdir)

mf = to_scf(tmpdir + '/FEMOCO.QMMM.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf().newton()

n_sites = umf.mol.nao

mo_coeff = np.load('data/qmmm-lmo.npy')
ncore = 517
ncas = 76

idx = np.array([
    4, 3, 45, 73, 70, 64, 36, 17, 12, 23,
    15, 19, 16, 9, 18, 13, 38, 57, 69, 66,
    52, 60, 53, 50, 65, 39, 42, 59, 54, 68,
    41, 5, 6, 11, 10, 33, 35, 34, 14, 8,
    7, 24, 31, 20, 63, 43, 46, 48, 49, 37,
    61, 71, 58, 51, 44, 40, 72, 47, 55, 25,
    26, 22, 21, 32, 29, 56, 28, 27, 30, 62,
    67, 74, 75, 0, 1, 2,
], dtype=int)
assert sorted(idx) == list(range(76))

mo_coeff[:, ncore:ncore + ncas] = mo_coeff[:, ncore:ncore + ncas][:, idx]
xmo = mo_coeff[:, ncore:ncore + ncas]

ovlp = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (18, 10))
orig_mo_coeff = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-orig-mo-coeff.npy' % (18, 11))

act_76 = np.array([
    3, 4,
    45, 73, 71, 65, 36,
    18, 13, 21, 11, 16, 17, 8, 20, 15,
    37, 57, 69, 67, 53,
    59, 51, 49, 64, 39,
    44, 60, 54, 68, 41,
    5, 6, 12, 9, 33, 35, 34, 14, 7, 10, 22, 29, 19,
    62, 42, 46, 48, 50,
    38, 61, 70, 58, 52,
    43, 40, 72, 47, 55,
    25, 26, 24, 23, 32, 31,
    56,
    27, 28, 30,
    63, 66, 74, 75,
    0, 1, 2], dtype=int)
assert sorted(act_76) == list(range(76))

orig_mo_coeff[:, ncore:ncore + ncas] = orig_mo_coeff[:, ncore:ncore + ncas][:, act_76]
zmo = orig_mo_coeff[:, ncore:ncore + ncas]

proja = np.einsum('ai,ab,bj->ij', zmo, ovlp, xmo, optimize=True)

with zipfile.ZipFile('../FEMOCO.FCIDUMP.zip', 'r') as zip_ref:
    zip_ref.extract('FEMOCO.FCIDUMP', tmpdir)

pmf = to_scf(tmpdir + '/FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None)
intor_symmetric = pmf.mol.intor_symmetric
pmf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(pmf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)
pmf.mol.groupname = 'c1'
pmf.mol.symmetry = True
pmf.mol.verbose = 4
pmf = pmf.to_uhf().newton()

xdstrs = open('06-collect-uccsd.out', 'r').readlines()

for ix in [iiid]:
    x = xdstrs[ix].strip()
    ix, jx, xdstr = int(x.split()[0]), int(x.split()[2]), x.split()[3]
    e_ccsd = float(x.split()[4])
    ttx = x.split()[-1]
    print('--- %5d : %5d %s ---' % (ix, jx, ttx))

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    if ttx == 'OPT':
        pmf.verbose = 4
        pmf.max_cycle = 800
        pmf.conv_tol = 1E-13
        pmf.ah_conv_tol = 1E-16
        pmf.kernel(dm0=dm0)
        dm1 = np.array(pmf.make_rdm1())
        dm0 = np.einsum('ia,jb,xij->xab', proja, proja, dm1, optimize=True)
    elif ttx == 'DIA':
        pass
    else:
        assert False

    umf.verbose = 4
    umf.max_cycle = 800
    umf.conv_tol = 1E-13
    umf.ah_conv_tol = 1E-16
    umf.kernel(dm0=dm0)
    dm1 = np.array(umf.make_rdm1())

    np.save('data/uccsd-mo-coeff-%02d-%04d.npy' % (fid, iiid), umf.mo_coeff)

    mmc = cc.CCSD(umf)
    mmc.diis_space = 36
    mmc.max_cycle = 400
    mmc.kernel()
    e_ccsd_t = mmc.e_tot + mmc.ccsd_t()

    if abs(mmc.e_tot - e_ccsd) >= 1E-6:
        print('!!WARNING', mmc.e_tot, e_ccsd)

    dm = mmc.make_rdm1()
    np.save('data/uccsd-dm1-%02d-%04d.npy' % (fid, iiid), dm)

    print("%5d :: %5d %s %20.12f %20.12f %d" % (ix, jx, xdstr, mmc.e_tot, e_ccsd_t, int(mmc.converged)))
