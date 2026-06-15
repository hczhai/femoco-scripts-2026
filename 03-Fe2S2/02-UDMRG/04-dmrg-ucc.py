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
import numpy as np, sys, itertools
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.mcscf import pmloc
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes
from pyscf import lib

lib.param.TMPDIR = '/tmp'

fid = 4
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/tmp/fe2-%02d-%02d" % (fid, iiid), symm_type=SymmetryTypes.SZ,
    restart_dir="/mnt/ceph/users/hzhai/fe2-%02d-%02d" % (fid, iiid),
    stack_mem=300 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-16, n_threads=96, n_mkl_threads=1, mpi=False)

bond_dims = [50] * 8 + [100] * 4 + [250] * 4 + [500] * 4 + [1000] * 4 + [1500] * 4 + [2000] * 4 + [3000] * 8
noises = [1E-5] * (len(bond_dims) - 4) + [0] * 4
thrds = [1E-7] * len(bond_dims)
n_sweeps = len(bond_dims)

mf = to_scf('../Fe2S2.FCIDUMP', molpro_orbsym=True, mf=None,)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf()

n_sites = umf.mol.nao

fe_idxs = [[2, 3, 4, 5, 6], [13, 14, 15, 16, 17]]

umf = umf.newton()

xdstr = '22bbbbb222222aaaaa22'

dm0 = np.zeros((2, n_sites, n_sites))
occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]
print(n_sites, dm0.shape, np.trace(dm0[0]), np.trace(dm0[1]))

umf.verbose = 4
umf.max_cycle = 800
umf.kernel(dm0=dm0)

print('EHF =', umf.e_tot)

np.save("data/%02d-uhf-mo-coeff-%02d.npy" % (fid, iiid), umf.mo_coeff)
np.save("data/%02d-uhf-1pdm-%02d.npy" % (fid, iiid), umf.make_rdm1())

from pyscf import cc

mmc = cc.CCSD(umf)
mmc.diis_space = 12
mmc.max_cycle = 1000
mmc.kernel()

print('ECCSD = ', mmc.e_tot)

e_ccsd_t = mmc.e_tot + mmc.ccsd_t()
print('ECCSD(T) = ', e_ccsd_t)

dm = mmc.make_rdm1()

nat_occ, u = np.linalg.eigh(dm)
nat_coeff = np.einsum('...pi,...ij->...pj', umf.mo_coeff, u, optimize=True)
nat_coeff = nat_coeff[..., ::-1]
nat_occ = nat_occ[..., ::-1]

np.save("data/%02d-uccsd-nat-coeff-%02d.npy" % (fid, iiid), nat_coeff)
np.save("data/%02d-uccsd-nat-occ-%02d.npy" % (fid, iiid), nat_occ)
np.save("data/%02d-uccsd-1pdm-%02d.npy" % (fid, iiid), dm)

dmao = np.array([umf.mo_coeff[i] @ dm[i] @ umf.mo_coeff[i].T for i in [0, 1]])

ma, mb = nat_coeff
neleca, nelecb = umf.mol.nelec

print('nelec =', neleca, nelecb)

partition = [
    [0], [1], [18], [19], # ends
    [2, 3, 4, 5, 6], [13, 14, 15, 16, 17], # fe
    [8, 10, 12], [7, 9, 11], # s
]
partition = [np.array(x, dtype=int) for x in partition]

_, u = pmloc(partition, ma[:, :neleca], tol=1e-6, maxcycle=1000, iop=1, iprint=1)
ma[:, :neleca] = ma[:, :neleca].dot(u)
_, u = pmloc(partition, ma[:, neleca:], tol=1e-6, maxcycle=1000, iop=1, iprint=1)
ma[:, neleca:] = ma[:, neleca:].dot(u)

_, u = pmloc(partition, mb[:, :nelecb], tol=1e-6, maxcycle=1000, iop=1, iprint=1)
mb[:, :nelecb] = mb[:, :nelecb].dot(u)
_, u = pmloc(partition, mb[:, nelecb:], tol=1e-6, maxcycle=1000, iop=1, iprint=1)
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

np.save("data/%02d-mo-coeff-%02d.npy" % (fid, iiid), umf.mo_coeff)

dmmo = np.array([umf.mo_coeff[i].T @ dmao[i] @ umf.mo_coeff[i] for i in [0, 1]])

occ_cc_a = np.copy(np.diag(dmmo[0]))
occ_cc_b = np.copy(np.diag(dmmo[1]))
occ_cc_a[occ_cc_a < 0] = 0
occ_cc_b[occ_cc_b < 0] = 0
occ_cc_a[occ_cc_a > 1] = 1
occ_cc_b[occ_cc_b > 1] = 1
occ_cc_a = occ_cc_a + 0.25 * (occ_cc_a < 0.5) - 0.25 * (occ_cc_a > 0.5)
occ_cc_b = occ_cc_b + 0.25 * (occ_cc_b < 0.5) - 0.25 * (occ_cc_b > 0.5)

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

print('A = ', umf.mo_occ[0, idx])
print('B = ', umf.mo_occ[1, idx])

driver.n_sites = ncas
driver.n_elec = n_elec
driver.ecore = ecore
driver.h1e = h1e
driver.g2e = g2e
driver.orb_sym = orb_sym

driver.spin = 0
twos = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

mpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=False,
    algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-12)
print("MPO = ", mpo.get_bond_dims())

kets = []

for ir in range(nroots):
    ket = driver.get_random_mps(tag="KET-%d" % ir, dot=1, bond_dim=250, occs=occs)
    driver.align_mps_center(ket, ref=0)
    ket = driver.adjust_mps(ket, dot=2)[0]
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

for ir in range(nroots):

    print('===== STATE %d =====' % ir)

    energy = driver.dmrg(mpo, kets[ir], n_sweeps=n_sweeps,
        bond_dims=bond_dims, noises=noises, thrds=thrds,
        lowmem_noise=True, twosite_to_onesite=None, tol=1E-12, cutoff=1E-24, iprint=2,
        dav_max_iter=400, dav_def_max_size=50, proj_weights=[5.0] * ir, proj_mpss=kets[:ir])
    print('STATE %d DMRG energy = %20.15f' % (ir, energy))

    kets[ir] = driver.adjust_mps(kets[ir], dot=1)[0]
    driver.align_mps_center(kets[ir], ref=0)
    print(ir, 'bond dim = ', kets[ir].info.bond_dim)

    dets, dvals = driver.get_csf_coefficients(kets[ir], cutoff=0.00001, iprint=1)

    np.save('data/%02d-dets-%02d.npy' % (fid, iiid), dets)
    np.save('data/%02d-dvals-%02d.npy' % (fid, iiid), dvals)
