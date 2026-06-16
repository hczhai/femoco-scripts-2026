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
import numpy as np, sys, sys

sys.path[:0] = ['/mnt/home/hzhai/Projects/hast']

import pyhast

fid = 24
iiid = int(sys.argv[1])

mf = to_scf('../Fe4S4.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf().newton()

n_sites = umf.mol.nao

xdstr = open('04-fe4-uccsd-min.out', 'r').readlines()[iiid].strip().split()[3]
print('xdstr =', xdstr)

dm0 = np.zeros((2, n_sites, n_sites))
occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]
print(n_sites, dm0.shape, np.trace(dm0[0]), np.trace(dm0[1]))

umf.verbose = 4
umf.max_cycle = 400
umf.verbose = 0
umf.conv_tol = 1E-13
umf.ah_conv_tol = 1E-16
umf.kernel(dm0=dm0)

np.save("data/%02d-uhf-mo-coeff-%02d.npy" % (fid, iiid), umf.mo_coeff)
np.save("data/%02d-uhf-1pdm-%02d.npy" % (fid, iiid), umf.make_rdm1())

print('INIT =', xdstr)
print('EHF  =', umf.e_tot)

cc_opts = dict(symm_amps=True, label_symm=True, lift=True, sum_lift=True, split=True, eval_t='hastctr',
    level_shift=0.15, gen_lamb_eq=True, gen_npdm_eq=True, npdm_order=2)

mcc = pyhast.sr.ucc.UCC(umf, t_order=2, verbose=4, diis=True, **cc_opts)
mcc.diis_space = 48
mcc.kernel(max_cycle=2400)

mcc2 = pyhast.sr.ucc.UCC(umf, t_order=4, verbose=4, diis=True, **cc_opts)
mcc2.diis_space = 48
tamps = mcc2.init_amps()
tamps[:2] = mcc.tamps
mcc2.kernel(tamps=tamps, max_cycle=2400)
mcc2.solve_lambda(max_cycle=2400)

dm1, dm2 = mcc2.make_npdms(order=2)
np.save('data/%02d-cc-%02d-%02d-pdm1.npy' % (fid, mcc2.order, iiid), dm1)
np.save('data/%02d-cc-%02d-%02d-pdm2.npy' % (fid, mcc2.order, iiid), dm2)
