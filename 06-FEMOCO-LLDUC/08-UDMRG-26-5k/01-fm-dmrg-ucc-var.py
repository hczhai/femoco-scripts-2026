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
from pyblock2._pyscf.mcscf import pmloc
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

if driver.mpi.rank != 0:
    umf.stdout = open(os.devnull, 'w')

n_sites = umf.mol.nao

xdstrs = """
   26 ::   628 22aa2aa222222222aaaaabbbbbbbbbb2222222222222a2aaabbb2baaaaa222222b222ba00222  -22140.306452593275  -22140.350446506614 1
   26 ::   627 22aaaaa2222222222aaaab2bbbbbbbb2222222222222a2aaabbbbbaaaaa222222b222ba00222  -22140.301345811393  -22140.348974068926 1
   26 ::   624 222aaaa222222222aaaaabbbb2bbbbb2222222222222aa2aabbbbbaaaaa222222b222ba00222  -22140.301272142111  -22140.349396899401 1
   26 ::   626 22aa2aa2222222222aaaabbbb2bbbbb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.293610397483  -22140.348079427644 1
   26 ::   625 22a2aaa2222222222aaaab2bbbbbbbb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.291955339446  -22140.342491973453 1
   26 ::   646 222aaaa2222222222aaaabbbbbbb2bb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.290567769043  -22140.355020085412 1
   26 ::   629 22aaaaa222222222aaaa2bbbbbbbb2b2222222222222aaaa2bbbbbaaaaa222222b222ba00222  -22140.286834711475  -22140.340654236712 1
   26 ::   631 22aaaaa222222222a2aaabb2bbbbbbb2222222222222aaaaabbbbbaaaa2222222b222ba00222  -22140.285356996927  -22140.346147141376 1
   26 ::   647 222aaaa2222222222aaaabbb2bbbbbb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.284755506444  -22140.347684014167 1
   26 ::   633 22aaa2a222222222a2aaabb2bbbbbbb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.279397359442  -22140.347944674286 1
   26 ::   632 22aaaaa2222222222aaaabbbbb2bbbb2222222222222aa2aabbbbbaaaaa222222b222ba00222  -22140.276448836965  -22140.343661423893 1
   26 ::   639 22aaaaa222222222aa2aabbbbbbbbb22222222222222aaa2abbbbbaaaaa222222b222ba00222  -22140.275055099130  -22140.357669742039 1
""".strip().split('\n')
assert len(xdstrs) == 12

xdstr = xdstrs[iiid].split()[3]
print('xdstr =', xdstr)

dm0 = np.zeros((2, n_sites, n_sites))
occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]
print(n_sites, dm0.shape, np.trace(dm0[0]), np.trace(dm0[1]))

umf.verbose = 4
umf.max_cycle = 800
umf.kernel(dm0=dm0)

print('EHF =', umf.e_tot)

if driver.mpi.rank == 0:
    np.save("data/%02d-uhf-mo-coeff-%02d.npy" % (fid, iiid), umf.mo_coeff)
    np.save("data/%02d-uhf-1pdm-%02d.npy" % (fid, iiid), umf.make_rdm1())

from pyscf import cc

mmc = cc.CCSD(umf)
mmc.diis_space = 12
mmc.max_cycle = 800

if driver.mpi.rank != 0:
    mmc.stdout = open(os.devnull, 'w')

mmc.kernel()

print('ECCSD = ', mmc.e_tot)

dm = mmc.make_rdm1()

nat_occ, u = np.linalg.eigh(dm)
nat_coeff = np.einsum('...pi,...ij->...pj', umf.mo_coeff, u, optimize=True)
nat_coeff = nat_coeff[..., ::-1]
nat_occ = nat_occ[..., ::-1]

if driver.mpi.rank == 0:
    np.save("data/%02d-uccsd-nat-coeff-%02d.npy" % (fid, iiid), nat_coeff)
    np.save("data/%02d-uccsd-nat-occ-%02d.npy" % (fid, iiid), nat_occ)
    np.save("data/%02d-uccsd-1pdm-%02d.npy" % (fid, iiid), dm)

e_ccsd_t = mmc.e_tot + mmc.ccsd_t()
print('ECCSD(T) = ', e_ccsd_t)

dmao = np.array([umf.mo_coeff[i] @ dm[i] @ umf.mo_coeff[i].T for i in [0, 1]])

ma, mb = nat_coeff
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

_, u = pmloc(partition, ma[:, :neleca], tol=1e-6, maxcycle=1000, iop=1, iprint=driver.mpi.rank == 0)
ma[:, :neleca] = ma[:, :neleca].dot(u)
_, u = pmloc(partition, ma[:, neleca:], tol=1e-6, maxcycle=1000, iop=1, iprint=driver.mpi.rank == 0)
ma[:, neleca:] = ma[:, neleca:].dot(u)

_, u = pmloc(partition, mb[:, :nelecb], tol=1e-6, maxcycle=1000, iop=1, iprint=driver.mpi.rank == 0)
mb[:, :nelecb] = mb[:, :nelecb].dot(u)
_, u = pmloc(partition, mb[:, nelecb:], tol=1e-6, maxcycle=1000, iop=1, iprint=driver.mpi.rank == 0)
mb[:, nelecb:] = mb[:, nelecb:].dot(u)

# orbital transform rot[a, b]
orb_rot = np.linalg.pinv(ma) @ mb
print('old overlap =', np.sum(np.diag(orb_rot ** 2)))

import scipy.optimize
kmidx = scipy.optimize.linear_sum_assignment(1 - orb_rot ** 2)[1]
mb = mb[:, kmidx]

orb_rot = np.linalg.pinv(ma) @ mb
print('new overlap =', np.sum(np.diag(orb_rot ** 2)))

print('pairing = ', kmidx)

umf.mo_coeff = np.array([ma, mb])

if driver.mpi.rank == 0:
    np.save("data/%02d-mo-coeff-%02d.npy" % (fid, iiid), umf.mo_coeff)

driver.mpi.barrier()
umf.mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (fid, iiid))
driver.mpi.barrier()

dmmo = np.array([umf.mo_coeff[i].T @ dmao[i] @ umf.mo_coeff[i] for i in [0, 1]])

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

ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_uhf_integrals(umf, ncore=0, ncas=None, g2e_symm=1, pg_symm=False)

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
    dav_max_iter=400, dav_def_max_size=50)
print('STATE %d DMRG energy = %20.15f' % (0, energy))
