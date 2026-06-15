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


from pyscf.tools.fcidump import to_scf
import numpy as np, sys, os

fid = 3
pfid = 1
iiid = int(sys.argv[1])

npdm_dir = './data'

mf = to_scf('../Fe2S2.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.nelectron = 31
mf.mol.spin = 1
mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf().newton()

n_sites = umf.mol.nao

xdstr = '22bbbbb222222a2aaa22'

dm0 = np.zeros((2, n_sites, n_sites))
occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]
print(n_sites, dm0.shape, np.trace(dm0[0]), np.trace(dm0[1]))

umf.verbose = 4
umf.max_cycle = 800
umf.conv_tol = 1E-14
umf.ah_conv_tol = 1E-17
umf.kernel(dm0=dm0)

print('INIT =', xdstr)
print('EHF  =', umf.e_tot)

umf.mo_coeff = np.load(npdm_dir + '/%02d-mo-coeff-%02d.npy' % (pfid, iiid))

from pyscf import cc

mmc = cc.CCSD(umf)
mmc.diis_space = 12
mmc.max_cycle = 1000
mmc.kernel()

print('ECCSD = ', mmc.e_tot)

e_ccsd_t = mmc.e_tot + mmc.ccsd_t()
print('ECCSD(T) = ', e_ccsd_t)

from pyscf.cc import uccsd_t_lambda, uccsd_t_rdm
from pyscf.fci.spin_op import spin_square_general

conv, mmc.l1, mmc.l2 = uccsd_t_lambda.kernel(mmc)
dm1 = uccsd_t_rdm.make_rdm1(mmc, mmc.t1, mmc.t2, mmc.l1, mmc.l2)
dm2 = uccsd_t_rdm.make_rdm2(mmc, mmc.t1, mmc.t2, mmc.l1, mmc.l2)

ss, mult = spin_square_general(*dm1, *dm2, umf.mo_coeff, np.identity(umf.mo_coeff[0].shape[1]))
print('<S^2> = %.6f S = %.6f' % (ss, (mult - 1)))

dm1ao = np.einsum('xpi,xij,xqj->xpq', umf.mo_coeff, dm1, umf.mo_coeff, optimize=True)
scr = np.zeros((2, ))
for ii, ix in enumerate([[2, 3, 4, 5, 6], [13, 14, 15, 16, 17]]):
    for ki in ix:
        scr[ii] += dm1ao[0][ki, ki] - dm1ao[1][ki, ki]
print("SZ ::", "%s" % (' '.join('%.6f' % x for x in scr)))
