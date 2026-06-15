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


from pyscf import gto, scf, cc, lib
import numpy as np
from pyscf.data.elements import chemcore
from pyscf.cc import dfuccsd

fid = 5
pfid = 4

fm_basis = {'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH', 'H': 'def2-SVP'}

mol = gto.M(atom="../MoFe3.xyz", basis=fm_basis, charge=-3, spin=3, symmetry=False, verbose=3)
mol.max_memory = 1400000

frozen = chemcore(mol)
print('frozen =', frozen)

mf = scf.sfx2c(scf.UHF(mol))
mf = mf.density_fit()
mf.max_cycle = 200
mf.verbose = 4

dm0 = np.load('%02d-uhf-dm1.npy' % pfid)

mf = mf.newton()
mf.max_cycle = 1200
mf.conv_tol = 1E-10
mf.ah_conv_tol = 1E-14
mf.kernel(dm0=dm0)
dm1 = np.array(mf.make_rdm1())

np.save('%02d-mo-coeff.npy' % fid, mf.mo_coeff)
np.save('%02d-mo-occ.npy' % fid, mf.mo_occ)
np.save('%02d-mo-energy.npy' % fid, mf.mo_energy)
np.save('%02d-uhf-dm1.npy' % fid, dm1)

mmc = dfuccsd.UCCSD(mf, frozen=frozen)
mmc.verbose = 5
mmc.diis_space = 12
mmc.max_cycle = 800
mmc.incore_complete = True
mmc.kernel()

np.save("/mnt/ceph/users/hzhai/tmp/t1a.npy", mmc.t1[0])
np.save("/mnt/ceph/users/hzhai/tmp/t1b.npy", mmc.t1[1])
np.save("/mnt/ceph/users/hzhai/tmp/t2aa.npy", mmc.t2[0])
np.save("/mnt/ceph/users/hzhai/tmp/t2ab.npy", mmc.t2[1])
np.save("/mnt/ceph/users/hzhai/tmp/t2bb.npy", mmc.t2[2])

dm = mmc.make_rdm1()

np.save('%02d-uccsd-dm1.npy' % fid, dm)
