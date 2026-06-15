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


import zipfile
from pyscf import lib
import numpy as np, sys, os, shutil
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes, ParallelTypes

fid = 13
ofid = 3
pfid = 0
iiid = int(sys.argv[1])
frozen = 38

driver = DMRGDriver(scratch="/tmp/fm-fno-udmrg-new-%02d-%02d" % (fid, iiid),
    mps_dir="/mnt/sdceph/users/hzhai/fm-fno-udmrg-new-%02d-%02d" % (fid, iiid),
    restart_dir="/mnt/sdceph/users/hzhai/fm-fno-udmrg-new-%02d-%02d-mps" % (fid, iiid),
    symm_type=SymmetryTypes.SZ, stack_mem=300 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-16, n_threads=48, n_mkl_threads=1, mpi=True)
driver.prule = driver.bw.bs.ParallelRuleQC(driver.mpi)

bond_dims = [7500] * 4 + [7000] * 4 + [6500] * 4 + [6000] * 4 + [5500] * 4 + [5000] * 4
noises = [0] * len(bond_dims)
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

if iiid == 25:
    xdstr = "22aaaaa2222222222aaaa2bbbbbbbbb2222222222222bbbbbaaa2aaaaaa222222b222ba00222"
elif iiid == 26:
    xdstr = "222aaaa2222222222aaaabbb2bbbbbb2222222222222aaaaabbbbbaaaaa222222b222ba00222"
elif iiid == 27:
    xdstr = "22aaaaa2222222222aaaa2bbbbbbbbb22222222222222aaaaaaaaabbbbb222222b222ba00222"

dm0 = np.zeros((2, n_sites, n_sites))
occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]
print(n_sites, dm0.shape, np.trace(dm0[0]), np.trace(dm0[1]))

if driver.mpi.rank != 0:
    umf.stdout = open(os.devnull, 'w')

umf.verbose = 4
umf.max_cycle = 800
umf.conv_tol = 1E-13
umf.ah_conv_tol = 1E-16
umf.kernel(dm0=dm0)

print('EHF =', umf.e_tot)

umf.mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (ofid, iiid))

ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_uhf_integrals(umf, ncore=frozen, ncas=None, g2e_symm=1, pg_symm=False)

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

driver.spin = 3
twos = 3

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

mpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=True,
    algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-12, para_type=ParallelTypes.Nothing)
print("MPO = ", mpo.get_bond_dims())

restart_dir = "/mnt/sdceph/users/hzhai/fm-fno-udmrg-new-%02d-%02d-mps" % (ofid, iiid)

driver.mpi.barrier()
if driver.mpi.rank == 0 and restart_dir != driver.scratch:
    for k in os.listdir(restart_dir):
        if '.KET-0.' in k or k == 'KET-0-mps_info.bin':
            shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)
driver.mpi.barrier()
ket = driver.load_mps("KET-0")
driver.mpi.barrier()
ket = driver.adjust_mps(ket, dot=2)[0]
driver.mpi.barrier()
print('bond dim = ', ket.info.bond_dim)

print('===== STATE 0 =====')

energy = driver.dmrg(mpo, ket, n_sweeps=n_sweeps,
    bond_dims=bond_dims, noises=noises, thrds=thrds,
    lowmem_noise=True, twosite_to_onesite=None, tol=0.0, cutoff=1E-24, iprint=2,
    dav_max_iter=400, dav_def_max_size=20)
print('STATE %d DMRG energy = %20.15f' % (0, energy))
