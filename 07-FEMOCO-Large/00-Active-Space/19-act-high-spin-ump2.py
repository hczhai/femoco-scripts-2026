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


from pyscf import gto, scf, mp
import numpy as np

fid = 19
ofid = 18

im = 3
ncas = [117, 180, 285, 404][im]
n_elec = [105, 143, 217, 277][im]
nfrozen = [159, 140, 103, 73][im]
ext_nfrozen = 362
spin = 35

assert nfrozen * 2 + n_elec == 423

ecore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, im))
h1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, im))
g2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, im))
mo_coeff = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-mo-coeff.npy' % (ofid, 10))
dmao = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-hs-dmao.npy' % (ofid, 10))
ovlp = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (ofid, 10))

print('load done')

mol = gto.M()
mol.nelectron = n_elec
mol.spin = spin
mol.nao = ncas
mol.energy_nuc = lambda *args: ecore
mol.incore_anyway = True
mol.max_memory = 1400000
mol.symmetry = False

mf = scf.UHF(mol)
mf.energy_nuc = lambda *args: ecore
mf.get_hcore = lambda *args: h1e
mf.get_ovlp = lambda *args: np.eye(mol.nao)
mf._eri = g2e
mf.max_memory = mol.max_memory

mf = mf.newton()

xmo = mo_coeff[:, ext_nfrozen + nfrozen:ext_nfrozen + nfrozen + ncas]
print(xmo.shape)
print(dmao.shape)
print(ovlp.shape)

dm0 = np.einsum('xab,ac,ci,bd,dj->xij', dmao, ovlp, xmo, ovlp, xmo, optimize=True)
print(np.trace(dm0[0]), np.trace(dm0[1]), n_elec)
print(np.trace(dm0[0]) + np.trace(dm0[1]))
print(np.trace(dm0[0]) - np.trace(dm0[1]))

mf.verbose = 4
mf.max_cycle = 200
mf.conv_tol = 1E-12
mf.ah_conv_tol = 1E-16
mf.kernel(dm0=dm0)

mmp = mp.MP2(mf).run()
