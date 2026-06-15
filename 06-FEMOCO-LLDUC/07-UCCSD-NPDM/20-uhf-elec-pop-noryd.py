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
from pyscf import lib, gto, scf, lo
from pyblock2._pyscf.scf import transform_rdm1_to_ao, transform_rdm1_to_mo

lib.param.TMPDIR = '/tmp'

fid = 20
pfid = 6

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

xmol = gto.M(atom="../FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=0)
xmol.max_memory = 1200000
xmo = np.load('../llduc-lmo.npy')

ncore = 517
ncas = 76

idx = np.array([
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

assert sorted(idx) == list(range(76))

xmo[:, ncore:ncore + ncas] = xmo[:, ncore:ncore + ncas][:, idx]

def mulliken_pop_dmao(mol, dmao):
    rmf = scf.RHF(mol)
    ld = lo.orth_ao(mol, "meta_lowdin")
    dm = transform_rdm1_to_mo(dmao, ld, rmf.get_ovlp())
    s = np.eye(dm.shape[-1])
    if isinstance(dm, np.ndarray) and dm.ndim == 2:
        dm = np.array((dm * 0.5, dm * 0.5))
    pop_a = np.einsum("ij,ji->i", dm[0], s).real
    pop_b = np.einsum("ij,ji->i", dm[1], s).real
    rydbglst = lo.nao._core_val_ryd_list(mol)[2]

    nelec_a = np.zeros(mol.natm)
    nelec_b = np.zeros(mol.natm)
    for i, s in enumerate(mol.ao_labels(fmt=None)):
        if i not in rydbglst:
            nelec_a[s[0]] += pop_a[i]
            nelec_b[s[0]] += pop_b[i]
    return (nelec_a + nelec_b)[:8], (nelec_a - nelec_b)[:8]

for ix in range(35):
    for jx in range(24):
        xdstr = open('../00-UHF/02-fm-uhf-24.out', 'r').readlines()[ix * 24 + jx].split()[3]
        mo_occ = [np.array([1] * 58 + [0] * 18), np.array([1] * 55 + [0] * 21)]
        mo_coeff = np.load('data/mo-coeff-%02d-%04d.npy' % (pfid, ix * 24 + jx))
        dm1 = np.einsum('xai,xi,xbi->xab', mo_coeff, mo_occ, mo_coeff, optimize='optimal')

        dm1_full = np.zeros((2, ncore + ncas, ncore + ncas))
        dm1_full[:, ncore:, ncore:] = dm1
        dm1_full[:, :ncore, :ncore] = np.identity(ncore)[None, :, :]
        dmao = transform_rdm1_to_ao(dm1_full, xmo[:, :ncore + ncas])
        pop, spn = mulliken_pop_dmao(xmol, dmao)
        print('%5d :: %5d %s' % (ix, jx, xdstr), ' '.join('%8.4f' % x for x in pop))
