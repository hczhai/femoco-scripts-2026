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
from pyscf.fci.spin_op import spin_square_general

fid = 9
ofid = 5
xfid = 7
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/tmp",
    mps_dir="/mnt/ceph/users/hzhai/fe4-udmrg-%02d-%02d" % (fid, iiid),
    restart_dir="/mnt/ceph/users/hzhai/fe4-udmrg-%02d-%02d-mps" % (fid, iiid),
    symm_type=SymmetryTypes.SZ, stack_mem=400 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-16, n_threads=64, n_mkl_threads=1, mpi=True)

mf = to_scf('../Fe4S4.FCIDUMP', molpro_orbsym=True, mf=None,)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf().newton()

n_sites = umf.mol.nao

umf.mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (ofid, iiid))

ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_uhf_integrals(umf, ncore=0, ncas=None, g2e_symm=1, pg_symm=False)

print('h1e.shape', [x.shape for x in h1e])
print('g2e.shape', [x.shape for x in g2e])

idx = np.load("data/%02d-reorder-%02d.npy" % (ofid, iiid))
print('reorder = ', idx)

h1e = [x[idx][:, idx] for x in h1e]
g2e = [x[idx][:, idx][:, :, idx][:, :, :, idx] for x in g2e]

driver.n_sites = ncas
driver.n_elec = n_elec
driver.ecore = ecore
driver.h1e = h1e
driver.g2e = g2e
driver.orb_sym = orb_sym

driver.spin = 0
twos = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

kets = []

restart_dir = "/mnt/ceph/users/hzhai/fe4-udmrg-%02d-%02d-mps" % (xfid, iiid)

for ir in range(nroots):
    driver.mpi.barrier()
    if driver.mpi.rank == 0 and restart_dir != driver.scratch:
        for k in os.listdir(restart_dir):
            if ('.KET-%d.' % ir) in k or k == 'KET-%d-mps_info.bin' % ir:
                shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)
    driver.mpi.barrier()
    ket = driver.load_mps("KET-%d" % ir)
    driver.mpi.barrier()
    ket = driver.adjust_mps(ket, dot=2)[0]
    driver.mpi.barrier()
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

    print('===== STATE %d =====' % ir)

    pdm1 = driver.get_1pdm(kets[ir], iprint=2, max_bond_dim=3000)
    pdm2 = driver.get_2pdm(kets[ir], iprint=2, max_bond_dim=3000)
    pdm2 = [x.transpose(0, 3, 1, 2) for x in pdm2]

    if driver.mpi.rank == 0:
        np.save('data/%02d-%02d-1pdm.npy' % (fid, iiid), pdm1)
        np.save('data/%02d-%02d-2pdm.npy' % (fid, iiid), pdm2)

        mo_coeff = umf.mo_coeff[:, :, idx]
        ss, mult = spin_square_general(*pdm1, *pdm2, mo_coeff, np.identity(mo_coeff[0].shape[1]))
        print('SSQ =', ss, mult)
