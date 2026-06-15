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

fid = 4

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

mf.mo_coeff = np.load('01-mo-coeff.npy')
mf.mo_energy = np.load('01-mo-energy.npy')
mf.mo_occ = np.load('01-mo-occ.npy')

print(mf.mo_occ, np.sum(mf.mo_occ, axis=1))

texts = []
lo_coeff, lo_occ, lo_energy = b2mcscf.get_uno(mf, loc_texts=texts)

print(lo_occ, np.sum(lo_occ))

np.save('%02d-lo-coeff.npy' % fid, lo_coeff)
np.save('%02d-lo-occ.npy' % fid, lo_occ)
np.save('%02d-lo-energy.npy' % fid, lo_energy)

coeff = lo_coeff
fname = 'orbs'
idx_a = 0
idx_b = 2818

tools.molden.from_mo(mol, fname + '.molden', coeff)
jmol_script = 'orbs.spt'
fspt = open(jmol_script, 'w')
fspt.write('''
initialize;
set background [xffffff];
set frank off
set autoBond true;
set bondRadiusMilliAngstroms 66;
set bondTolerance 0.5;
set forceAutoBond false;
''')
fspt.write('load %s.molden;\n' % fname)
fspt.write('rotate -30 y;\n')
fspt.write('rotate 60 z;\n')
fspt.write('zoom 130;\n')
fspt.write('translate x 2;\n')
for i in range(idx_a, idx_b):
    fspt.write('isoSurface MO %d fill noMesh noDots;\n' % (i+1))
    fspt.write('color isoSurface translucent 0.6 [x3300cc];\n')
    fspt.write('color isoSurface phase [x3E92CC] [xFFB238];\n')
    fspt.write('set echo top left; echo "%s"; color echo gray;\n' % texts[i])
    fspt.write('font echo 14 monospaced;\n')
    fspt.write('write JPG 90 "orbs/%s-%03d.jpg";\n' % (jmol_script, i))
