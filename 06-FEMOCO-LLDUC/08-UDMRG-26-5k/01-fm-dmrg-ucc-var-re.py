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


import block2 as b
from pyscf import gto, scf, lib
import numpy as np, sys
import itertools, zipfile, os
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes, ParallelTypes

fid = 1
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(
    scratch="/tmp/fm-m1-%02d-%02d" % (fid, iiid),
    mps_dir="/mnt/ceph/users/hzhai/fm-m1-%02d-%02d" % (fid, iiid),
    restart_dir="/mnt/ceph/users/hzhai/fm-m1-%02d-%02d-mps" % (fid, iiid),
    symm_type=SymmetryTypes.SZ, stack_mem=400 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-16, n_threads=96, n_mkl_threads=1, mpi=True)
driver.prule = driver.bw.bs.ParallelRuleQC(driver.mpi)

bond_dims = [50] * 8 + [100] * 4 + [250] * 4 + [500] * 4 + [1000] * 4 + [1500] * 4 + [2000] * 4 + [3000] * 4 + [4000] * 4 + [5000] * 8
noises = [1E-5] * (len(bond_dims) - 4) + [0] * 4
thrds = [1E-7] * len(bond_dims)
n_sweeps = len(bond_dims)

lib.param.TMPDIR = '/tmp'

with zipfile.ZipFile('../FEMOCO.FCIDUMP.zip', 'r') as zip_ref:
    zip_ref.extract('FEMOCO.FCIDUMP', '/tmp/')

mf = to_scf('/tmp/FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf().newton()

n_sites = umf.mol.nao

print('EHF =', umf.e_tot)

driver.mpi.barrier()
umf.mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (fid, iiid))
idx = np.load("data/%02d-reorder-%02d.npy" % (fid, iiid))
driver.mpi.barrier()

ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_uhf_integrals(umf, ncore=0, ncas=None, g2e_symm=1, pg_symm=False)

print('h1e.shape', [x.shape for x in h1e])
print('g2e.shape', [x.shape for x in g2e])
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

mpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=True,
    algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-12, para_type=ParallelTypes.Nothing)
print("MPO = ", mpo.get_bond_dims())

driver.mpi.barrier()
ket = driver.load_mps("KET-%d" % 0)
driver.mpi.barrier()
print('bond dim = ', ket.info.bond_dim)

prev_data = open('01-fm-dmrg-ucc-var-%d.out' % iiid, 'r').readlines()
sweep_start = [int(x.split()[2]) for x in prev_data if x.startswith('Sweep =')][-1]
forward = sweep_start % 2 == 0
print(iiid, sweep_start, forward)

print('===== STATE 0 =====')

energy = driver.dmrg(mpo, ket, n_sweeps=n_sweeps,
    bond_dims=bond_dims, noises=noises, thrds=thrds,
    lowmem_noise=True, twosite_to_onesite=None, tol=1E-12, cutoff=1E-24, iprint=2,
    dav_max_iter=400, dav_def_max_size=50, sweep_start=sweep_start, forward=forward)
print('STATE %d DMRG energy = %20.15f' % (0, energy))

pdm1 = driver.get_1pdm(ket, iprint=2, max_bond_dim=5000)
pdm2 = driver.get_2pdm(ket, iprint=2, max_bond_dim=5000)

if driver.mpi.rank == 0:
    np.save('data/%02d-%02d-1pdm.npy' % (fid, iiid), pdm1)
    np.save('data/%02d-%02d-2pdm.npy' % (fid, iiid), pdm2)
