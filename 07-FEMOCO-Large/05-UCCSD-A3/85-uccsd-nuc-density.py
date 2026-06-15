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

fid = 85
ofid = 18
pfid = 67
xfid = 9

im = 3
ncas = [117, 180, 285, 404][im]
n_elec = [105, 143, 217, 277][im]
nfrozen = [159, 140, 103, 73][im]
ext_nfrozen = 362

mo_coeff = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-mo-coeff.npy' % (ofid, 10))
ovlp = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (ofid, 10))

xmo = mo_coeff[:, :ext_nfrozen + nfrozen + ncas]

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

xmol = gto.M(atom="../FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=0)
xmol.max_memory = 1200000


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

for ix in [6, 7, 8, 9, 10]:
    xdstr = {int(l.split()[0]): l for l in open('%02d-uccsd-min-merge.out' % xfid, 'r').readlines()}[ix].split()[3]
    jx = int({int(l.split()[0]): l for l in open('%02d-uccsd-min-merge.out' % xfid, 'r').readlines()}[ix].split()[2])
    mo_coeff = np.load('data/uccsd-mo-coeff-%02d-%04d.npy' % (pfid, ix))
    dm1_mo = np.load('data/uccsd-dm1-%02d-%04d.npy' % (pfid, ix))
    dm1 = np.einsum('xpi,xij,xqj->xpq', mo_coeff, dm1_mo, mo_coeff, optimize=True)

    ncore = ext_nfrozen + nfrozen
    dm1_full = np.zeros((2, ncore + ncas, ncore + ncas))
    dm1_full[:, ncore:, ncore:] = dm1
    dm1_full[:, :ncore, :ncore] = np.identity(ncore)[None, :, :]
    dmao = transform_rdm1_to_ao(dm1_full, xmo[:, :ncore + ncas])
    rho = nuclear_density(xmol, dmao)
    print('%5d :: %5d %s' % (ix, jx, xdstr), ' '.join('%20.12f' % x for x in rho))
