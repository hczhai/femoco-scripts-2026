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

fid = 1

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

mol = gto.M(atom="../../FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=3)
mol.max_memory = 500000

mf = scf.sfx2c(scf.UKS(mol))
mf.xc = 'b3lyp'
mf = solvent.ddcosmo.ddcosmo_for_scf(mf)
mf.with_solvent.eps = 4.0
mf.max_cycle = 1200
mf.verbose = 4

dm0 = mf.from_chk('00-ref-hs_bp86.chk')
mf.kernel(dm0=dm0)
dm1 = np.array(mf.make_rdm1())

np.save('%02d-mo-coeff.npy' % fid, mf.mo_coeff)
np.save('%02d-mo-occ.npy' % fid, mf.mo_occ)
np.save('%02d-mo-energy.npy' % fid, mf.mo_energy)
np.save('%02d-uhf-dm1.npy' % fid, dm1)
