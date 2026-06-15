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


from pyscf import gto, scf, solvent
import numpy as np

fid = 8
ofid = 5
pfid = 1

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

mol = gto.M(atom="../../FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=3)
mol.max_memory = 500000

ovlp = mol.intor('int1e_ovlp')

mo_coeff = np.load('%02d-mo-coeff.npy' % pfid)
mo_cas = np.load('%02d-sel-mo-coeff.npy' % ofid)

dm0 = np.load('%02d-uhf-dm1.npy' % pfid)

print(dm0.shape)
print(np.trace(dm0[0]), np.trace(dm0[1]))

dm1 = np.einsum('xab,ac,ci,bd,dj->xij', dm0, ovlp, mo_cas, ovlp, mo_cas, optimize=True)

print(dm1.shape)
print(np.trace(dm1[0]), np.trace(dm1[1]), np.trace(dm1[0]) + np.trace(dm1[1]), np.trace(dm1[0]) - np.trace(dm1[1]))

np.save('%02d-proj-dm1.npy' % fid, dm1)
