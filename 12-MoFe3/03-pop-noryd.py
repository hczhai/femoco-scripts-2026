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
from pyblock2._pyscf.scf import transform_rdm1_to_mo

fid = 3
pfid = 2

fm_basis = {'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH', 'H': 'def2-SVP'}

mol = gto.M(atom="../MoFe3.xyz", basis=fm_basis, charge=-3, spin=3, symmetry=False, verbose=3)
mol.max_memory = 500000

mf = scf.sfx2c(scf.UKS(mol))
mf.xc = 'b3lyp'
mf = solvent.ddcosmo.ddcosmo_for_scf(mf)
mf.with_solvent.eps = 4.0
mf.max_cycle = 200
mf.verbose = 4

dm0 = np.load('%02d-uhf-dm1.npy' % pfid)

mf.kernel(dm0=dm0)

def mulliken_pop_dmao(mol, dmao):
    from pyscf import scf, lo

    rmf = scf.RHF(mol)
    ld = lo.orth_ao(mol, "meta_lowdin")
    dm = transform_rdm1_to_mo(dmao, ld, rmf.get_ovlp())
    s = np.eye(dm.shape[-1])
    if isinstance(dm, np.ndarray) and dm.ndim == 2:
        dm = np.array((dm * 0.5, dm * 0.5))
    pop_a = np.einsum("ij,ji->i", dm[0], s).real
    pop_b = np.einsum("ij,ji->i", dm[1], s).real
    corelst, vallst, rydbglst = lo.nao._core_val_ryd_list(mol)

    nelec_a = np.zeros(mol.natm)
    nelec_b = np.zeros(mol.natm)
    for i, s in enumerate(mol.ao_labels(fmt=None)):
        if i in rydbglst:
            continue
        nelec_a[s[0]] += pop_a[i]
        nelec_b[s[0]] += pop_b[i]
    chg = (nelec_a + nelec_b)
    return chg[:4], (nelec_a - nelec_b)[:4]

chg, spn = mulliken_pop_dmao(mol, mf.make_rdm1())

print('CHG ', ' '.join('%8.4f' % x for x in chg))
print('SPN ', ' '.join('%8.4f' % x for x in spn))
