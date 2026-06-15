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


import numpy as np
from pyscf import lib, gto
from pyblock2._pyscf.scf import transform_rdm1_to_ao

lib.param.TMPDIR = '/tmp'

fid = 12
pfid = 7
xfid = 4

ncore = 517
ncas = 76

mo_coeff = np.load('data/qmmm-lmo.npy')

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

xmol = gto.M(atom="../FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=0)
xmol.max_memory = 1200000

idx = np.array([
    4, 3, 45, 73, 70, 64, 36, 17, 12, 23,
    15, 19, 16, 9, 18, 13, 38, 57, 69, 66,
    52, 60, 53, 50, 65, 39, 42, 59, 54, 68,
    41, 5, 6, 11, 10, 33, 35, 34, 14, 8,
    7, 24, 31, 20, 63, 43, 46, 48, 49, 37,
    61, 71, 58, 51, 44, 40, 72, 47, 55, 25,
    26, 22, 21, 32, 29, 28, 27, 30, 56, 62,
    67, 74, 75, 0, 1, 2,
], dtype=int)

assert sorted(idx) == list(range(76))

mo_coeff[:, ncore:ncore + ncas] = mo_coeff[:, ncore:ncore + ncas][:, idx]
xmo = mo_coeff[:, :ncore + ncas]


def nuclear_density(mol, dmao):
    rhos = []
    for ia in range(8):
        coord = mol.atom_coord(ia)[None, :]
        ao_values = mol.eval_gto('GTOval_sph', coord)
        rho_alpha = np.einsum('ij,pi,pj->p', dmao[0], ao_values, ao_values, optimize=True)
        rho_beta = np.einsum('ij,pi,pj->p', dmao[1], ao_values, ao_values, optimize=True)
        rho = rho_alpha + rho_beta
        rhos.append(rho[0])
    return rhos

for ix in range(35):
    xdstr = open('%02d-uhf-min.out' % xfid, 'r').readlines()[ix].split()[3]
    jx = int(open('%02d-uhf-min.out' % xfid, 'r').readlines()[ix].split()[2])
    dm1 = np.load('data/uhf-dm1-%02d-%04d.npy' % (pfid, ix))

    dm1_full = np.zeros((2, ncore + ncas, ncore + ncas))
    dm1_full[:, ncore:, ncore:] = dm1
    dm1_full[:, :ncore, :ncore] = np.identity(ncore)[None, :, :]
    dmao = transform_rdm1_to_ao(dm1_full, xmo[:, :ncore + ncas])
    rho = nuclear_density(xmol, dmao)
    print('%5d :: %5d %s' % (ix, jx, xdstr), ' '.join('%20.12f' % x for x in rho))
