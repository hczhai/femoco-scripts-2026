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
from pyscf import gto, scf, lib
import numpy as np, sys, os
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.mcscf import pmloc
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes, ParallelTypes

fid = 4
pfid = 0
iiid = int(sys.argv[1])
frozen = 34

driver = DMRGDriver(scratch="/tmp/fm-fno-udmrg-new-%02d-%02d" % (fid, iiid),
    mps_dir="/mnt/sdceph/users/hzhai/fm-fno-udmrg-new-%02d-%02d" % (fid, iiid),
    restart_dir="/mnt/sdceph/users/hzhai/fm-fno-udmrg-new-%02d-%02d-mps" % (fid, iiid),
    symm_type=SymmetryTypes.SZ, stack_mem=300 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-16, n_threads=48, n_mkl_threads=1, mpi=True)
driver.prule = driver.bw.bs.ParallelRuleQC(driver.mpi)

bond_dims = [50] * 8 + [100] * 4 + [250] * 4 + [500] * 4 + [1000] * 4 + [1500] * 4 + [2000] * 4 + [3000] * 4 + [4500] * 4 + [6000] * 4 + [8000] * 8
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

umf.mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (pfid, iiid))

from pyscf import cc

mmc = cc.CCSD(umf, frozen=frozen)
mmc.diis_space = 12
mmc.max_cycle = 800

if driver.mpi.rank != 0:
    mmc.stdout = open(os.devnull, 'w')

mmc.kernel()

print('ECCSD = ', mmc.e_tot)

dm = mmc.make_rdm1()

dmao = np.array([umf.mo_coeff[i, :, frozen:] @ dm[i][frozen:, frozen:] @ umf.mo_coeff[i, :, frozen:].T for i in [0, 1]])

ma, mb = umf.mo_coeff
neleca, nelecb = umf.mol.nelec

print('nelec =', neleca, nelecb)

partition = [
    [0, 1], [73], [74], [75], # ends
    [2, 3, 4, 5, 6], [16, 17, 18, 19, 20], [21, 22, 23, 24, 25],
    [26, 27, 28, 29, 30], [44, 45, 46, 47, 48], [49, 50, 51, 52, 53],
    [54, 55, 56, 57, 58], # fe
    [65, 69, 70, 71, 72], # mo
    [59, 60, 63], [61, 62, 67], [64, 66, 68], # S right
    [7, 12, 14], [8, 13, 15], [9, 10, 11], # S left
    [31, 32, 34], [33, 39, 40], [41, 42, 43], # S mid
    [35, 36, 37, 38] # C
]
partition = [np.array(x, dtype=int) for x in partition]

_, u = pmloc(partition, ma[:, frozen:neleca], tol=1e-6, maxcycle=1000, iop=1, iprint=driver.mpi.rank == 0)
ma[:, frozen:neleca] = ma[:, frozen:neleca].dot(u)
_, u = pmloc(partition, ma[:, neleca:], tol=1e-6, maxcycle=1000, iop=1, iprint=driver.mpi.rank == 0)
ma[:, neleca:] = ma[:, neleca:].dot(u)

_, u = pmloc(partition, mb[:, frozen:nelecb], tol=1e-6, maxcycle=1000, iop=1, iprint=driver.mpi.rank == 0)
mb[:, frozen:nelecb] = mb[:, frozen:nelecb].dot(u)
_, u = pmloc(partition, mb[:, nelecb:], tol=1e-6, maxcycle=1000, iop=1, iprint=driver.mpi.rank == 0)
mb[:, nelecb:] = mb[:, nelecb:].dot(u)

# orbital transform rot[a, b]
orb_rot = np.linalg.pinv(ma[:, frozen:]) @ mb[:, frozen:]
print('old overlap =', np.sum(np.diag(orb_rot ** 2)))

import scipy.optimize
kmidx = scipy.optimize.linear_sum_assignment(1 - orb_rot ** 2)[1]
mb[:, frozen:] = mb[:, frozen:][:, kmidx]

orb_rot = np.linalg.pinv(ma[:, frozen:]) @ mb[:, frozen:]
print('new overlap =', np.sum(np.diag(orb_rot ** 2)))

print('pairing = ', kmidx)

umf.mo_coeff = np.array([ma, mb])

if driver.mpi.rank == 0:
    np.save("data/%02d-mo-coeff-%02d.npy" % (fid, iiid), umf.mo_coeff)

driver.mpi.barrier()
umf.mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (fid, iiid))
driver.mpi.barrier()

dmmo = np.array([umf.mo_coeff[i, :, frozen:].T @ dmao[i] @ umf.mo_coeff[i, :, frozen:] for i in [0, 1]])

occ_cc_a = np.copy(np.diag(dmmo[0]))
occ_cc_b = np.copy(np.diag(dmmo[1]))
occ_cc_a[occ_cc_a < 0] = 0
occ_cc_b[occ_cc_b < 0] = 0
occ_cc_a[occ_cc_a > 1] = 1
occ_cc_b[occ_cc_b > 1] = 1
occ_cc_a = occ_cc_a + 0.025 * (occ_cc_a < 0.5) - 0.025 * (occ_cc_a > 0.5)
occ_cc_b = occ_cc_b + 0.025 * (occ_cc_b < 0.5) - 0.025 * (occ_cc_b > 0.5)

print('occ cc a', occ_cc_a)
print('occ cc b', occ_cc_b)

ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_uhf_integrals(umf, ncore=frozen, ncas=None, g2e_symm=1, pg_symm=False)

print('h1e.shape', [x.shape for x in h1e])
print('g2e.shape', [x.shape for x in g2e])

import block2 as b2

xmat = np.abs(np.einsum("ijji->ij", g2e[0] + g2e[1] * 2 + g2e[2], optimize=True))
kmat = np.abs(h1e[0] + h1e[1]) * 1e-7 + xmat
kmat = b2.VectorDouble(kmat.ravel())
idx = b2.OrbitalOrdering.fiedler(ncas, kmat)

opts = dict(n_generations=1000, n_configs=ncas * 2,
    n_elite=8, clone_rate=0.1, mutate_rate=0.1
)
n_tasks = 64
idxs = []
for i_task in range(0, n_tasks):
    b2.Random.rand_seed(1234 + i_task)
    idx = b2.OrbitalOrdering.ga_opt(ncas, kmat, **opts)
    f = b2.OrbitalOrdering.evaluate(ncas, kmat, idx)
    idx = tuple(idx)
    idxs.append(idx)
idx = np.array(sorted(list(set(idxs)))[0], dtype=int)

if driver.mpi.rank == 0:
    np.save("data/%02d-reorder-%02d.npy" % (fid, iiid), idx)

print('reorder = ', idx)

h1e = [x[idx][:, idx] for x in h1e]
g2e = [x[idx][:, idx][:, :, idx][:, :, :, idx] for x in g2e]

occ_alpha = np.array(occ_cc_a)
occ_beta = np.array(occ_cc_b)

print('occ_alpha = ', occ_alpha)
print('occ_beta = ', occ_beta)

occs = np.array([x for xa, xb in zip(occ_alpha[idx], occ_beta[idx]) for x in [xa, xb]], dtype=float)

print('occs = ', occs)

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

ket = driver.get_random_mps(tag="KET-%d" % 0, dot=1, bond_dim=50, occs=occs)
driver.mpi.barrier()
driver.align_mps_center(ket, ref=0)
driver.mpi.barrier()
ket = driver.adjust_mps(ket, dot=2)[0]
driver.mpi.barrier()
print('bond dim = ', ket.info.bond_dim)

print('===== STATE 0 =====')

energy = driver.dmrg(mpo, ket, n_sweeps=n_sweeps,
    bond_dims=bond_dims, noises=noises, thrds=thrds,
    lowmem_noise=True, twosite_to_onesite=None, tol=1E-12, cutoff=1E-24, iprint=2,
    dav_max_iter=400, dav_def_max_size=20)
print('STATE %d DMRG energy = %20.15f' % (0, energy))
