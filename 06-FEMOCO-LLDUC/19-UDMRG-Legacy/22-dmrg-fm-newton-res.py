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
import numpy as np, sys
import itertools
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.mcscf import get_uno
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes
from pyblock2.algebra.io import MPSTools

def sqrtm(s):
    e, v = np.linalg.eigh(s)
    return np.dot(v * np.sqrt(e), v.T.conj())

def pmloc(partition, mocoeff, tol=1e-6, maxcycle=1000, iop=1, iprint=1):
    if iprint:
        print()
        print("[pm_loc_kernel]")
        print(" mocoeff.shape=", mocoeff.shape)
        print(" tol=", tol)
        print(" maxcycle=", maxcycle)
        print(" partition=", len(partition), "\\n", partition)
    k = mocoeff.shape[0]
    n = mocoeff.shape[1]
    natom = len(partition)

    def genPaij(mocoeff, partition):
        s12c = mocoeff.copy()
        # Lowdin
        natom = len(partition)
        pija = np.zeros((natom, n, n))
        for iatom in range(natom):
            idx = partition[iatom]
            pija[iatom] = np.dot(s12c[idx, :].T, s12c[idx, :])
        # P[i,j,a]
        pija = pija.transpose(1, 2, 0).copy()
        return pija

    u = np.identity(n)
    pija = genPaij(mocoeff, partition)

    # Start
    def funval(pija):
        return np.einsum("iia,iia", pija, pija)

    fun = funval(pija)
    if iprint:
        print(" initial funval = ", fun)
    for icycle in range(maxcycle):
        delta = 0.0
        # i>j
        ijdx = []
        for i in range(n - 1):
            for j in range(i + 1, n):
                bij = abs(np.sum(pija[i, j] * (pija[i, i] - pija[j, j])))
                ijdx.append((i, j, bij))
        ijdx = sorted(ijdx, key=lambda x: x[2], reverse=True)
        for i, j, bij in ijdx:
            # determine angle
            vij = pija[i, i] - pija[j, j]
            aij = np.dot(pija[i, j], pija[i, j]) - 0.25 * np.dot(vij, vij)
            bij = np.dot(pija[i, j], vij)
            if abs(aij) < 1.0e-10 and abs(bij) < 1.0e-10:
                continue
            p1 = np.sqrt(aij ** 2 + bij ** 2)
            cos4a = -aij / p1
            sin4a = bij / p1
            cos2a = np.sqrt((1 + cos4a) * 0.5)
            sin2a = np.sqrt((1 - cos4a) * 0.5)
            cosa = np.sqrt((1 + cos2a) * 0.5)
            sina = np.sqrt((1 - cos2a) * 0.5)
            # Why? Because we require alpha in [0,pi/2]
            if sin4a < 0.0:
                cos2a = -cos2a
                sina, cosa = cosa, sina
            # stationary condition
            if abs(cosa - 1.0) < 1.0e-10:
                continue
            if abs(sina - 1.0) < 1.0e-10:
                continue
            # incremental value
            delta += p1 * (1 - cos4a)
            # Transformation
            # Urot
            ui = u[:, i] * cosa + u[:, j] * sina
            uj = -u[:, i] * sina + u[:, j] * cosa
            u[:, i] = ui.copy()
            u[:, j] = uj.copy()
            # Bra-transformation of Integrals
            tmp_ip = pija[i, :, :] * cosa + pija[j, :, :] * sina
            tmp_jp = -pija[i, :, :] * sina + pija[j, :, :] * cosa
            pija[i, :, :] = tmp_ip.copy()
            pija[j, :, :] = tmp_jp.copy()
            # Ket-transformation of Integrals
            tmp_ip = pija[:, i, :] * cosa + pija[:, j, :] * sina
            tmp_jp = -pija[:, i, :] * sina + pija[:, j, :] * cosa
            pija[:, i, :] = tmp_ip.copy()
            pija[:, j, :] = tmp_jp.copy()
        fun = fun + delta
        if iprint:
            print("icycle=", icycle, "delta=", delta, "fun=", fun)
        if delta < tol:
            break

    # Check
    ierr = 0
    if delta < tol:
        if iprint:
            print("CONG: PMloc converged!")
    else:
        ierr = 1
        if iprint:
            print("WARNING: PMloc not converged")
    return ierr, u

fid = 22
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/central/scratch/hczhai/h2d-fm8-%02d-%02d" % (fid, iiid),
    restart_dir="/central/scratch/hczhai/h2d-fm8-%02d-%02d-mps" % (fid, iiid),
    symm_type=SymmetryTypes.SZ, stack_mem=120 << 30, min_mpo_mem=False,
    fp_codec_cutoff=0.0,
    n_threads=64, n_mkl_threads=1, mpi=False)

bond_dims = [50] * 8 + [100] * 4 + [250] * 4 + [500] * 4 + [1000] * 4 + [1500] * 4 + [2000] * 4 + [3000] * 4
noises = [1E-5] * (len(bond_dims) - 4) + [0] * 4
thrds = [1E-7] * len(bond_dims)
n_sweeps = len(bond_dims)

mf = to_scf('FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None,)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf()

n_sites = umf.mol.nao

fe_idxs = [
    [2, 3, 4, 5, 6], [16, 17, 18, 19, 20], [21, 22, 23, 24, 25],
    [26, 27, 28, 29, 30], [44, 45, 46, 47, 48], [49, 50, 51, 52, 53],
    [54, 55, 56, 57, 58]
]
mo_idxs = [[65, 69, 70, 71, 72]]
end_idxs = [[0, 1], [73, 74, 75]]
s_idxs = [[7, 8, 9, 10, 11, 12, 13, 14, 15], [59, 60, 61, 62, 63, 64, 68, 66, 67]]
sc_idxs = [[31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43]]

idxs = list(fe_idxs) + list(mo_idxs) + list(end_idxs) + list(s_idxs) + list(sc_idxs)
assert max([max(x) for x in idxs]) == n_sites - 1
assert min([min(x) for x in idxs]) == 0

occs = []

for ix, x in enumerate(itertools.combinations(range(7), 3)):
    y = tuple(z for z in range(7) if z not in x)
    xoccs = []
    for jx in range(3 * len(list(itertools.combinations(range(4), 2)))):
        kajd, kbjd = jx // 6, jx % 6
        occ = ['2'] * n_sites
        for ii, ip in enumerate(mo_idxs[0]):
            occ[ip] = 'bba00'[ii]
        for iv, v in enumerate(x):
            for ii, ip in enumerate(fe_idxs[v]):
                pxocc = ['bbbbb', 'bbbbb', 'bbbbb']
                pxocc[kajd] = '2bbbb'
                occ[ip] = pxocc[iv][ii]
        for iv, v in enumerate(y):
            for ii, ip in enumerate(fe_idxs[v]):
                pxocc = ['aaaaa', 'aaaaa', 'aaaaa', 'aaaaa']
                pxopt = list(itertools.combinations(range(4), 2))[kbjd]
                pxocc[pxopt[0]] = '2aaaa'
                pxocc[pxopt[1]] = '2aaaa'
                occ[ip] = pxocc[iv][ii]
        xoccs.append(''.join(occ))
    occs.append(xoccs)

umf = umf.newton()

emin = 0.0
xd_min = -1
for xd in list(range(len(occs[iiid]))) + [-1]:
    if xd == -1:
        xd = xd_min
    xdstr = occs[iiid][xd]
    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]
    print(n_sites, dm0.shape, np.trace(dm0[0]), np.trace(dm0[1]))

    umf.verbose = 4
    umf.max_cycle = 800
    umf.kernel(dm0=dm0)
    if umf.e_tot < emin:
        xd_min = xd
        emin = umf.e_tot

print('XD =', xd_min)
print('EMIN =', emin)
print('EHF =', umf.e_tot)

np.save("%02d-uhf-mo-coeff-%02d.npy" % (fid, iiid), umf.mo_coeff)
np.save("%02d-uhf-1pdm-%02d.npy" % (fid, iiid), umf.make_rdm1())

from pyscf import cc

mmc = cc.CCSD(umf)
mmc.diis_space = 12
mmc.max_cycle = 800
mmc.kernel()

print('ECCSD = ', mmc.e_tot)

dm = mmc.make_rdm1()

nat_occ, u = np.linalg.eigh(dm)
nat_coeff = np.einsum('...pi,...ij->...pj', umf.mo_coeff, u, optimize=True)
nat_coeff = nat_coeff[..., ::-1]
nat_occ = nat_occ[..., ::-1]

np.save("%02d-uccsd-nat-coeff-%02d.npy" % (fid, iiid), nat_coeff)
np.save("%02d-uccsd-nat-occ-%02d.npy" % (fid, iiid), nat_occ)
np.save("%02d-uccsd-1pdm-%02d.npy" % (fid, iiid), dm)

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

np.save("%02d-mo-coeff-%02d.npy" % (fid, iiid), umf.mo_coeff)

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

np.save("%02d-reorder-%02d.npy" % (fid, iiid), idx)

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

mpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=False,
    algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-12)
print("MPO = ", mpo.get_bond_dims())

kets = []

for ir in range(nroots):
    # ket = driver.get_mps_from_csf_coefficients(xdstrs, dvals=dvals, tag="KET-%d" % ir, dot=1)
    ket = driver.get_random_mps(tag="KET-%d" % ir, dot=1, bond_dim=50, occs=occs)
    driver.align_mps_center(ket, ref=0)
    ket = driver.adjust_mps(ket, dot=2)[0]
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

for ir in range(nroots):

    print('===== STATE %d =====' % ir)

    energy = driver.dmrg(mpo, kets[ir], n_sweeps=n_sweeps,
        bond_dims=bond_dims, noises=noises, thrds=thrds,
        lowmem_noise=True, twosite_to_onesite=None, tol=1E-12, cutoff=1E-24, iprint=2,
        dav_max_iter=400, dav_def_max_size=20, proj_weights=[5.0] * ir, proj_mpss=kets[:ir])
    print('STATE %d DMRG energy = %20.15f' % (ir, energy))

    kets[ir] = driver.adjust_mps(kets[ir], dot=1)[0]
    driver.align_mps_center(kets[ir], ref=0)

    impo = driver.get_identity_mpo()
    hmpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=True,
        algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-12)

    energy = driver.expectation(kets[ir], hmpo, kets[ir], iprint=2)
    print('STATE %d Expt energy = %20.15f' % (ir, energy))

    dm = driver.get_2pdm(kets[ir], mask=(0, 1, 1, 0), iprint=2, max_bond_dim=3000)
    np.save("%02d-dm-orig-%02d-%02d.npy" % (fid, iiid, ir), dm)

    dm = driver.get_1pdm(kets[ir], iprint=2, max_bond_dim=3000)
    np.save("%02d-1pdm-orig-%02d-%02d.npy" % (fid, iiid, ir), dm)

    hmpo, impo = None, None
