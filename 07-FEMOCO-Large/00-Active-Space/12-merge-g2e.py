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

def compute_ijs(ncas, nb):
    r = []
    for i in range((ncas + nb - 1) // nb):
        for j in range(i + 1):
            r.append((i, j))
    return r

ncas = mo_cas.shape[1]
nb = 200
ijs = compute_ijs(ncas, nb)

npair = ncas * (ncas + 1) // 2
g2e = np.empty((npair * (npair + 1) // 2, ))
for ix, (i, j) in enumerate(ijs):
    print(ix, i, j)
    frag = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/frag-%03d.npy' % (xfid, ix))
    ishape = min(ncas, (i + 1) * nb) - i * nb
    jshape = min(ncas, (j + 1) * nb) - j * nb
    if i == j:
        xtril = np.tril_indices(ishape)
        a, b  = xtril[0] + i * nb, xtril[1] + i * nb
        lidx = a * (a + 1) // 2 + b
    else:
        a, b  = np.mgrid[:ishape] + i * nb, np.mgrid[:jshape] + j * nb
        lidx = (a[:, None] * (a[:, None] + 1) // 2 + b[None, :]).ravel()
    a, b = np.tril_indices(ncas)
    ridx = a * (a + 1) // 2 + b
    idx = lidx[:, None] * (lidx[:, None] + 1) // 2 + ridx[None, :]
    mask = lidx[:, None] >= ridx[None, :]
    g2e[idx[mask]] = frag[mask]
    idx = ridx[None, :] * (ridx[None, :] + 1) // 2 + lidx[:, None]
    mask = lidx[:, None] < ridx[None, :]
    g2e[idx[mask]] = frag[mask]

np.save('/mnt/ceph/users/hzhai/%02d-fm-active/g2e.npy' % fid, g2e)
g2e = None
