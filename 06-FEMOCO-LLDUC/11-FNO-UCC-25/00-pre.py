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
from pyscf import gto, scf, lib
import numpy as np, sys, os
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.ao2mo import integrals as itg

fid = 0
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

# xdstr = open('../FM-INIT', 'r').readlines()[iiid].strip().split()[-1]
xdstr = "22aaaaa2222222222aaaa2bbbbbbbbb2222222222222bbbbbaaa2aaaaaa222222b222ba00222"

dm0 = np.zeros((2, n_sites, n_sites))
occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]
print(n_sites, dm0.shape, np.trace(dm0[0]), np.trace(dm0[1]))

umf.verbose = 4
umf.max_cycle = 800
umf.conv_tol = 1E-13
umf.ah_conv_tol = 1E-16
umf.kernel(dm0=dm0)

print('EHF =', umf.e_tot)

from pyscf import cc

mmc = cc.CCSD(umf)
mmc.diis_space = 12
mmc.max_cycle = 800

mmc.kernel()

print('ECCSD = ', mmc.e_tot)

dm1 = mmc.make_rdm1(with_mf=False)

noa, nob = umf.mol.nelec
wa, ua = np.linalg.eigh(dm1[0][:noa, :noa])
wb, ub = np.linalg.eigh(dm1[1][:nob, :nob])

wa, wb = wa[::-1], wb[::-1]
nat_coeff_occ_a = np.dot(umf.mo_coeff[0][:, :noa], ua)[..., ::-1]
nat_coeff_occ_b = np.dot(umf.mo_coeff[1][:, :nob], ub)[..., ::-1]

print('alpha nat occ (corr)')
print(wa)
print('beta nat occ (corr)')
print(wb)

mo_coeff_vir_a = umf.mo_coeff[0][:, noa:]
mo_coeff_vir_b = umf.mo_coeff[1][:, nob:]

moa = np.hstack([nat_coeff_occ_a, mo_coeff_vir_a])
mob = np.hstack([nat_coeff_occ_b, mo_coeff_vir_b])

umf.mo_coeff = [moa, mob]

np.save("data/%02d-mo-coeff-%02d.npy" % (fid, iiid), umf.mo_coeff)
