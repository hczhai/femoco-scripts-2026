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


from pyscf import gto, scf
import numpy as np

fid = 13
ofid = 12
xfid = 8
pfid = 5

ncas, n_elec, spin = 994, 423, 35
ecore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/ecore.npy' % ofid)
h1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/h1e.npy' % ofid)
g2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/g2e.npy' % ofid)

print('load done')

mol = gto.M()
mol.nelectron = n_elec
mol.spin = spin
mol.nao = ncas
mol.energy_nuc = lambda *args: ecore
mol.incore_anyway = True
mol.max_memory = 970000
mol.symmetry = False

mf = scf.UHF(mol)
mf.energy_nuc = lambda *args: ecore
mf.get_hcore = lambda *args: h1e
mf.get_ovlp = lambda *args: np.eye(mol.nao)
mf._eri = g2e
mf.max_memory = mol.max_memory

mf = mf.newton()

dm0 = np.load('%02d-proj-dm1.npy' % xfid)

mf.verbose = 4
mf.max_cycle = 200
mf.conv_tol = 1E-13
mf.ah_conv_tol = 1E-16
mf.kernel(dm0=dm0)

dm1 = np.array(mf.make_rdm1())

np.save('%02d-mo-coeff.npy' % fid, mf.mo_coeff)
np.save('%02d-mo-occ.npy' % fid, mf.mo_occ)
np.save('%02d-mo-energy.npy' % fid, mf.mo_energy)
np.save('%02d-uhf-dm1.npy' % fid, dm1)
