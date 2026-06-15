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


import numpy as np, sys, os, shutil
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes

fid = 16
pfid = {26: 5, 25: 10}[int(sys.argv[1])]
ofid = 22
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(
    scratch="/tmp/fm-m0-%02d-%02d" % (fid, iiid),
    symm_type=SymmetryTypes.SZ, stack_mem=400 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-14, n_threads=96, n_mkl_threads=1, mpi=False)

mf = to_scf('../FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None,)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf()

n_sites = umf.mol.nao

umf.mo_coeff = np.load("ref/%02d-mo-coeff-%02d.npy" % (ofid, iiid))

ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_uhf_integrals(umf, ncore=0, ncas=None, g2e_symm=1, pg_symm=False)

print('h1e.shape', [x.shape for x in h1e])
print('g2e.shape', [x.shape for x in g2e])

idx = np.load("ref/%02d-reorder-%02d.npy" % (ofid, iiid))
print('reorder = ', idx)

h1e = [x[idx][:, idx] for x in h1e]
g2e = [x[idx][:, idx][:, :, idx][:, :, :, idx] for x in g2e]

driver.n_sites = ncas
driver.n_elec = n_elec
driver.ecore = ecore
driver.h1e = h1e
driver.g2e = g2e
driver.orb_sym = orb_sym

driver.spin = 3
twos = 3

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

kets = []

restart_dir = "/mnt/ceph/users/hzhai/fm-m0-%02d-%02d-mps" % (pfid, iiid)

for ir in range(nroots):
    for k in os.listdir(restart_dir):
        if ('.KET-%d.' % ir) in k or k == 'KET-%d-mps_info.bin' % ir:
            shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)
    ket = driver.load_mps("KET-%d" % ir)
    ket = driver.adjust_mps(ket, dot=1)[0]
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

driver.get_csf_coefficients(kets[0], cutoff=0.02, iprint=1)
