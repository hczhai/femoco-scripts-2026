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


from pyscf import gto, scf
import numpy as np, sys, zipfile
from pyscf import lib
from pyscf.tools.fcidump import to_scf

fid = 11
ofid = 18
iiid = int(sys.argv[1])

im = 3
ncas = [117, 180, 285, 404][im]
n_elec = [105, 143, 217, 277][im]
nfrozen = [159, 140, 103, 73][im]
ext_nfrozen = 362

assert nfrozen * 2 + n_elec == 423

ecore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, im))
h1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, im))
g2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, im))

mol = gto.M()
mol.nelectron = n_elec
mol.spin = 3
mol.nao = ncas
mol.energy_nuc = lambda *args: ecore
mol.incore_anyway = True
mol.max_memory = 1400000
mol.symmetry = False

mf = scf.UHF(mol)
mf.energy_nuc = lambda *args: ecore
mf.get_hcore = lambda *args: h1e
mf.get_ovlp = lambda *args: np.eye(mol.nao)
mf._eri = g2e
mf.max_memory = mol.max_memory

mf = mf.newton()

lib.param.TMPDIR = '/tmp'

pm = 2
pncas = [117, 180, 285, 404][pm]
pn_elec = [105, 143, 217, 277][pm]
pnfrozen = [159, 140, 103, 73][pm]

assert pnfrozen * 2 + pn_elec == 423

pecore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, pm))
ph1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, pm))
pg2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, pm))
proja = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, pm))
projc = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, pm))

pmol = gto.M()
pmol.nelectron = pn_elec
pmol.spin = 3
pmol.nao = pncas
pmol.energy_nuc = lambda *args: pecore
pmol.incore_anyway = True
pmol.max_memory = 1400000
pmol.symmetry = False

pmf = scf.UHF(pmol)
pmf.energy_nuc = lambda *args: pecore
pmf.get_hcore = lambda *args: ph1e
pmf.get_ovlp = lambda *args: np.eye(pmol.nao)
pmf._eri = pg2e
pmf.max_memory = pmol.max_memory

pmf = pmf.newton()

qm = 1
qncas = [117, 180, 285, 404][qm]
qn_elec = [105, 143, 217, 277][qm]
qnfrozen = [159, 140, 103, 73][qm]

assert qnfrozen * 2 + qn_elec == 423

qecore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, qm))
qh1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, qm))
qg2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, qm))
qproja = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, qm))
qprojc = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, qm))

qmol = gto.M()
qmol.nelectron = qn_elec
qmol.spin = 3
qmol.nao = qncas
qmol.energy_nuc = lambda *args: qecore
qmol.incore_anyway = True
qmol.max_memory = 1400000
qmol.symmetry = False

qmf = scf.UHF(qmol)
qmf.energy_nuc = lambda *args: qecore
qmf.get_hcore = lambda *args: qh1e
qmf.get_ovlp = lambda *args: np.eye(qmol.nao)
qmf._eri = qg2e
qmf.max_memory = qmol.max_memory

qmf = qmf.newton()

rm = 0
rncas = [117, 180, 285, 404][rm]
rn_elec = [105, 143, 217, 277][rm]
rnfrozen = [159, 140, 103, 73][rm]

assert rnfrozen * 2 + rn_elec == 423

recore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, rm))
rh1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, rm))
rg2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, rm))
rproja = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, rm))
rprojc = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, rm))

rmol = gto.M()
rmol.nelectron = rn_elec
rmol.spin = 3
rmol.nao = rncas
rmol.energy_nuc = lambda *args: recore
rmol.incore_anyway = True
rmol.max_memory = 1400000
rmol.symmetry = False

rmf = scf.UHF(rmol)
rmf.energy_nuc = lambda *args: recore
rmf.get_hcore = lambda *args: rh1e
rmf.get_ovlp = lambda *args: np.eye(rmol.nao)
rmf._eri = rg2e
rmf.max_memory = rmol.max_memory

rmf = rmf.newton()

with zipfile.ZipFile('../FEMOCO.FCIDUMP.zip', 'r') as zip_ref:
    zip_ref.extract('FEMOCO.FCIDUMP', '/tmp/')

umf = to_scf('/tmp/FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None)
intor_symmetric = umf.mol.intor_symmetric
umf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(umf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)
umf.mol.groupname = 'c1'
umf.mol.symmetry = True
umf.mol.verbose = 4
umf = umf.to_uhf().newton()
n_sites = umf.mol.nao

xdstrs = open('02-uhf-min-24.out', 'r').readlines()

mo_coeff = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-mo-coeff.npy' % (ofid, 10))
ovlp = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (ofid, 10))

xmo = mo_coeff[:, ext_nfrozen + nfrozen:ext_nfrozen + nfrozen + ncas]
pmo_c = mo_coeff[:, :ext_nfrozen + pnfrozen]
pmo_a = mo_coeff[:, ext_nfrozen + pnfrozen:ext_nfrozen + pnfrozen + pncas]
pprojc = np.einsum('ai,ab,bj->ij', pmo_c, ovlp, xmo, optimize=True)
pproja = np.einsum('ai,ab,bj->ij', pmo_a, ovlp, xmo, optimize=True)

qmo_c = mo_coeff[:, :ext_nfrozen + qnfrozen]
qmo_a = mo_coeff[:, ext_nfrozen + qnfrozen:ext_nfrozen + qnfrozen + qncas]
qqprojc = np.einsum('ai,ab,bj->ij', qmo_c, ovlp, pmo_a, optimize=True)
qqproja = np.einsum('ai,ab,bj->ij', qmo_a, ovlp, pmo_a, optimize=True)

rmo_c = mo_coeff[:, :ext_nfrozen + rnfrozen]
rmo_a = mo_coeff[:, ext_nfrozen + rnfrozen:ext_nfrozen + rnfrozen + rncas]
rrprojc = np.einsum('ai,ab,bj->ij', rmo_c, ovlp, qmo_a, optimize=True)
rrproja = np.einsum('ai,ab,bj->ij', rmo_a, ovlp, qmo_a, optimize=True)

data = []

for x in xdstrs:
    ix, jx, xdstr = int(x.split()[0]), int(x.split()[2]), x.split()[3]
    e_ref = float(x.split()[4])
    ttx = x.split()[-1]
    if ix != iiid:
        continue
    print('--- %5d : %5d %s ---' % (ix, jx, ttx))

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    if ttx == 'OPT':
        umf.verbose = 4
        umf.max_cycle = 200
        umf.conv_tol = 1E-12
        umf.ah_conv_tol = 1E-16
        umf.kernel(dm0=dm0)
        dm1 = np.array(umf.make_rdm1())
    elif ttx == 'DIA':
        dm1 = dm0
    elif ttx in ['PJ1-OPT', 'PJ1-DIA']:
        if ttx == 'PJ1-OPT':
            umf.verbose = 4
            umf.max_cycle = 200
            umf.conv_tol = 1E-12
            umf.ah_conv_tol = 1E-16
            umf.kernel(dm0=dm0)
            dm1 = np.array(umf.make_rdm1())
        elif ttx == 'PJ1-DIA':
            dm1 = dm0
        
        dm0 = np.einsum('ia,jb,xij->xab', qproja, qproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', qprojc, qprojc, [1, 1], optimize=True)

        qmf.verbose = 4
        qmf.max_cycle = 200
        qmf.conv_tol = 1E-12
        qmf.ah_conv_tol = 1E-16
        qmf.kernel(dm0=dm0)
        dm1 = np.array(qmf.make_rdm1())
    elif ttx in ['PJ1-PJ0-OPT', 'PJ1-PJ0-DIA']:
        if ttx == 'PJ1-PJ0-OPT':
            umf.verbose = 4
            umf.max_cycle = 200
            umf.conv_tol = 1E-12
            umf.ah_conv_tol = 1E-16
            umf.kernel(dm0=dm0)
            dm1 = np.array(umf.make_rdm1())
        elif ttx == 'PJ1-PJ0-DIA':
            dm1 = dm0
        
        dm0 = np.einsum('ia,jb,xij->xab', rproja, rproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', rprojc, rprojc, [1, 1], optimize=True)

        rmf.verbose = 4
        rmf.max_cycle = 200
        rmf.conv_tol = 1E-12
        rmf.ah_conv_tol = 1E-16
        rmf.kernel(dm0=dm0)
        dm1 = np.array(rmf.make_rdm1())

        dm0 = np.einsum('ia,jb,xij->xab', rrproja, rrproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', rrprojc, rrprojc, [1, 1], optimize=True)

        qmf.verbose = 4
        qmf.max_cycle = 200
        qmf.conv_tol = 1E-12
        qmf.ah_conv_tol = 1E-16
        qmf.kernel(dm0=dm0)
        dm1 = np.array(qmf.make_rdm1())
    else:
        assert False

    if ttx in ['PJ1-PJ0-OPT', 'PJ1-PJ0-DIA', 'PJ1-OPT', 'PJ1-DIA']:
        dm0 = np.einsum('ia,jb,xij->xab', qqproja, qqproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', qqprojc, qqprojc, [1, 1], optimize=True)
    else:
        dm0 = np.einsum('ia,jb,xij->xab', proja, proja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', projc, projc, [1, 1], optimize=True)

    pmf.verbose = 4
    pmf.max_cycle = 200
    pmf.conv_tol = 1E-12
    pmf.ah_conv_tol = 1E-16
    pmf.kernel(dm0=dm0)
    dm1 = np.array(pmf.make_rdm1())

    if abs(pmf.e_tot - e_ref) >= 1E-6:
        print('!!WARNING', pmf.e_tot, e_ref)

    dm0 = np.einsum('ia,jb,xij->xab', pproja, pproja, dm1, optimize=True)
    dm0 += np.einsum('ia,ib,x->xab', pprojc, pprojc, [1, 1], optimize=True)

    mf.verbose = 4
    mf.max_cycle = 200
    mf.conv_tol = 1E-12
    mf.ah_conv_tol = 1E-16
    mf.kernel(dm0=dm0)
    dm1 = np.array(mf.make_rdm1())

    found, ifx = -1, len(data)
    for lx, (kx, xener, xdm1) in enumerate(data):
        if abs(xener - mf.e_tot) < 1E-6 and np.linalg.norm(dm1 - xdm1) < 1E-3:
            found, ifx = kx, lx
            break
    if found == -1:
        data.append((jx, mf.e_tot, dm1))

    print(":: %5d %s %20.12f %5d %5d" % (jx, xdstr, mf.e_tot, found, ifx))
