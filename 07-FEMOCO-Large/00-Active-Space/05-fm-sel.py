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
from pyblock2._pyscf import mcscf as b2mcscf
from pyscf import tools

fid = 5
pfid = 4

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

mol = gto.M(atom="../../FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=3)
mol.max_memory = 500000

lo_coeff = np.load('%02d-lo-coeff.npy' % pfid)
lo_occ = np.load('%02d-lo-occ.npy' % pfid)
lo_energy = np.load('%02d-lo-energy.npy' % pfid)

norb = lo_coeff.shape[1]
ova = mol.intor_symmetric("cint1e_ovlp_sph")
s12 = b2mcscf.uno.sqrtm(ova)
lcoeff = s12.dot(lo_coeff)

diff = lcoeff.conj().T @ lcoeff - np.identity(norb)
assert np.linalg.norm(diff) < 1e-7

ncore, nact = 556, 35

fm_atom_idxs = [
    [0, 1, 2, 3, 4, 5, 6], # Fe
    [7], # Mo
    [8], # C
    [9, 10, 11, 12, 13, 14, 15, 16, 17, 18], # S
    [90, 92, 75], # O, O, N (Mo-end)
]

il = 0
nelec = 0
sel_core_mo = []
sel_virt_mo = []
sel_mo = []
sel_occ = []

print("\nLowdin population for LMOs:")
labels = mol.ao_labels(None)
for iorb in range(norb):
    vec = lcoeff[:, iorb] ** 2
    ivs = np.argsort(vec)
    if iorb < ncore:
        text = "[C %3d] occ = %8.5f" % (iorb, lo_occ[iorb])
        ftext = " fii = %10.3f" % lo_energy[iorb]
    elif iorb >= ncore and iorb < ncore + nact:
        text = "[A %3d] occ = %8.5f" % (iorb, lo_occ[iorb])
        ftext = " fii = %10.3f" % lo_energy[iorb]
    else:
        text = "[V %3d] occ = %8.5f" % (iorb, lo_occ[iorb])
        ftext = " fii = %10.3f" % lo_energy[iorb]
    gtext = ""
    for iao in ivs[::-1][:3]:
        gtext += "(%3d-%2s-%7s = %5.3f) " % (
            labels[iao][0],
            labels[iao][1],
            labels[iao][2] + labels[iao][3],
            vec[iao],
        )
    ok = any(labels[iao][0] in [x for xx in fm_atom_idxs for x in xx] and vec[iao] > 0.05 for iao in ivs[::-1][:16])
    kt = 'SEL%4d' % il if ok else '       '
    il += ok
    nelec += lo_occ[iorb] if ok else 0
    if ok:
        sel_mo.append(lo_coeff[:, iorb])
        sel_occ.append(lo_occ[iorb])
    elif iorb < ncore:
        sel_core_mo.append(lo_coeff[:, iorb])
    else:
        sel_virt_mo.append(lo_coeff[:, iorb])
    print(kt + " " + text + ftext + " " + gtext)

print('Selected =', il, 'Nelec =', nelec, 'Ncore =', len(sel_core_mo))

np.save('%02d-sel-core-mo-coeff.npy' % fid, np.array(sel_core_mo).T)
np.save('%02d-sel-virt-mo-coeff.npy' % fid, np.array(sel_virt_mo).T)
np.save('%02d-sel-mo-coeff.npy' % fid, np.array(sel_mo).T)
np.save('%02d-sel-mo-occ.npy' % fid, np.array(sel_occ))
