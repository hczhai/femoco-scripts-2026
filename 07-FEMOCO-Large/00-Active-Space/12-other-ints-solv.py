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


from pyscf import gto, scf, solvent, ao2mo, lib
import numpy as np, sys

print(lib.param.TMPDIR)

fid = 12
pfid = 5
xfid = 11

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

mol = gto.M(atom="../../FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=3)
mol.max_memory = 1200000

mf = scf.sfx2c(scf.RHF(mol))
mf = solvent.ddcosmo.ddcosmo_for_scf(mf)
mf.with_solvent.eps = 4.0
mf.max_cycle = 1200
mf.verbose = 4

mo_core = np.load('%02d-sel-core-mo-coeff.npy' % pfid)
mo_cas = np.load('%02d-sel-mo-coeff.npy' % pfid)
sel_occ = np.load('%02d-sel-mo-occ.npy' % pfid)

ecore = mol.energy_nuc()
core_dmao = 2 * mo_core @ mo_core.T.conj()
vj, vk = mf.get_jk(mol, core_dmao)
hveff_ao = vj - 0.5 * vk
hcore_ao = mf.get_hcore()
dmao = np.load('01-uhf-dm1.npy')
_, solvent_v = mf.with_solvent.kernel(dmao[0] + dmao[1])
hcore_ao += solvent_v
ecore += np.einsum("ij,ji->", core_dmao, hcore_ao + 0.5 * hveff_ao, optimize=True)
print('ecore =', ecore)

n_elec = mol.nelectron - mo_core.shape[1] * 2
ncas = mo_cas.shape[1]

print('Ncas =', ncas, 'Nelec =', n_elec)

h1e = mo_cas.T.conj() @ (hcore_ao + hveff_ao) @ mo_cas
print(h1e.shape)

np.save('/mnt/ceph/users/hzhai/%02d-fm-active/ecore.npy' % fid, ecore)
np.save('/mnt/ceph/users/hzhai/%02d-fm-active/h1e.npy' % fid, h1e)
