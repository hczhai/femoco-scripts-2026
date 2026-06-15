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
from pyscf import gto, scf
import numpy as np, sys, os, shutil
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes, ParallelTypes

fid = 10
pfid = 3
ofid = 22
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(
    scratch="/tmp/fm-m0-%02d-%02d" % (fid, iiid),
    scratch_quota=1700 << 30,
    alt_scratch="/mnt/ceph/users/hzhai/tmp/fm-m0-%02d-%02d" % (fid, iiid),
    mps_dir="/mnt/ceph/users/hzhai/fm-m0-%02d-%02d" % (fid, iiid),
    restart_dir="/mnt/ceph/users/hzhai/fm-m0-%02d-%02d-mps" % (fid, iiid),
    symm_type=SymmetryTypes.SZ, stack_mem=400 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-14, n_threads=96, n_mkl_threads=1, mpi=True)
driver.prule = driver.bw.bs.ParallelRuleQC(driver.mpi)

bond_dims = [10000] * 4 + [12000] * 3 + [15000] * 3 + [18000] * 2
noises = [1E-5] * len(bond_dims)
thrds = [1E-7] * len(bond_dims)
n_sweeps = len(bond_dims)

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

mpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=True,
    algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-12, para_type=ParallelTypes.Nothing)
print("MPO = ", mpo.get_bond_dims())

kets = []

for ir in range(nroots):
    driver.mpi.barrier()
    ket = driver.load_mps("KET-%d" % ir)
    driver.mpi.barrier()
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

def f(name, iprint):
    import resource
    if driver.mpi.rank == 0 and name == "DMRG::sweep::iter.start":
        print("%10s" % driver.bw.b.Parsing.to_size_string(
            1024 * resource.getrusage(resource.RUSAGE_SELF).ru_maxrss), end='', flush=True)

driver.set_callback(f)

for ir in range(nroots):

    print('===== STATE %d =====' % ir)

    energy = driver.dmrg(mpo, kets[ir], n_sweeps=n_sweeps,
        bond_dims=bond_dims, noises=noises, thrds=thrds,
        lowmem_noise=True, twosite_to_onesite=None, tol=1E-12, cutoff=1E-24, iprint=2,
        dav_max_iter=400, dav_def_max_size=20, proj_weights=[5.0] * ir, proj_mpss=kets[:ir],
        sweep_start=9, forward=False)
    print('STATE %d DMRG energy = %20.15f' % (ir, energy))
