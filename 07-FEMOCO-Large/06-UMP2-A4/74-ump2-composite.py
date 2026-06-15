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


from pyscf import gto, scf, mp
import numpy as np, sys, os, zipfile
from pyscf import lib
from pyscf.tools.fcidump import to_scf

fid = 73
ofid = 18
gfid = 12
iiid = int(sys.argv[1])

im = 3
ncas = [117, 180, 285, 404][im]
n_elec = [105, 143, 217, 277][im]
nfrozen = [159, 140, 103, 73][im]
ext_nfrozen = 362

assert nfrozen * 2 + n_elec == 423

gncas, gn_elec, gspin = 994, 423, 3
gecore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/ecore.npy' % gfid)
gh1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/h1e.npy' % gfid)
gg2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/g2e.npy' % gfid)
guno = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/ump2-uno-coeff.npy' % gfid)

print('large load done')

gmol = gto.M()
gmol.nelectron = gn_elec
gmol.spin = gspin
gmol.nao = gncas
gmol.energy_nuc = lambda *args: gecore
gmol.incore_anyway = True
gmol.max_memory = 5700000
gmol.symmetry = False

gmf = scf.UHF(gmol)
gmf.energy_nuc = lambda *args: gecore
gmf.get_hcore = lambda *args: gh1e
gmf.get_ovlp = lambda *args: np.eye(gmol.nao)
gmf._eri = gg2e
gmf.max_memory = gmol.max_memory

gmf = gmf.newton()

ecore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, im))
h1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, im))
g2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, im))
proja = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, im))
projc = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, im))

mol = gto.M()
mol.nelectron = n_elec
mol.spin = 3
mol.nao = ncas
mol.energy_nuc = lambda *args: ecore
mol.incore_anyway = True
mol.max_memory = 1490000
mol.symmetry = False

mf = scf.UHF(mol)
mf.energy_nuc = lambda *args: ecore
mf.get_hcore = lambda *args: h1e
mf.get_ovlp = lambda *args: np.eye(mol.nao)
mf._eri = g2e
mf.max_memory = mol.max_memory

mf = mf.newton()

sm = 2
sncas = [117, 180, 285, 404][sm]
sn_elec = [105, 143, 217, 277][sm]
snfrozen = [159, 140, 103, 73][sm]

assert snfrozen * 2 + sn_elec == 423

secore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, sm))
sh1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, sm))
sg2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, sm))
sproja = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, sm))
sprojc = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, sm))

smol = gto.M()
smol.nelectron = sn_elec
smol.spin = 3
smol.nao = sncas
smol.energy_nuc = lambda *args: secore
smol.incore_anyway = True
smol.max_memory = mol.max_memory
smol.symmetry = False

smf = scf.UHF(smol)
smf.energy_nuc = lambda *args: secore
smf.get_hcore = lambda *args: sh1e
smf.get_ovlp = lambda *args: np.eye(smol.nao)
smf._eri = sg2e
smf.max_memory = smol.max_memory

smf = smf.newton()

pm = 1
pncas = [117, 180, 285, 404][pm]
pn_elec = [105, 143, 217, 277][pm]
pnfrozen = [159, 140, 103, 73][pm]

assert pnfrozen * 2 + pn_elec == 423

pecore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, pm))
ph1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, pm))
pg2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, pm))
rproja = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, pm))
rprojc = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, pm))

pmol = gto.M()
pmol.nelectron = pn_elec
pmol.spin = 3
pmol.nao = pncas
pmol.energy_nuc = lambda *args: pecore
pmol.incore_anyway = True
pmol.max_memory = mol.max_memory
pmol.symmetry = False

pmf = scf.UHF(pmol)
pmf.energy_nuc = lambda *args: pecore
pmf.get_hcore = lambda *args: ph1e
pmf.get_ovlp = lambda *args: np.eye(pmol.nao)
pmf._eri = pg2e
pmf.max_memory = pmol.max_memory

pmf = pmf.newton()

qm = 0
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
qmol.max_memory = mol.max_memory
qmol.symmetry = False

qmf = scf.UHF(qmol)
qmf.energy_nuc = lambda *args: qecore
qmf.get_hcore = lambda *args: qh1e
qmf.get_ovlp = lambda *args: np.eye(qmol.nao)
qmf._eri = qg2e
qmf.max_memory = qmol.max_memory

qmf = qmf.newton()

lib.param.TMPDIR = '/tmp'

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

xdstrs = open('03-uccsd-min.out', 'r').readlines()

mo_coeff = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-mo-coeff.npy' % (ofid, 10))
ovlp = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (ofid, 10))

gmo_c = mo_coeff[:, :ext_nfrozen]
gmo_a = mo_coeff[:, ext_nfrozen:ext_nfrozen + gncas]

xmo_c = mo_coeff[:, :ext_nfrozen + nfrozen]
xmo_a = mo_coeff[:, ext_nfrozen + nfrozen:ext_nfrozen + nfrozen + ncas]
gprojc = np.einsum('ai,ab,bj->ij', xmo_c, ovlp, gmo_a, optimize=True)
gproja = np.einsum('ai,ab,bj->ij', xmo_a, ovlp, gmo_a, optimize=True)

smo_c = mo_coeff[:, :ext_nfrozen + snfrozen]
smo_a = mo_coeff[:, ext_nfrozen + snfrozen:ext_nfrozen + snfrozen + sncas]
ssprojc = np.einsum('ai,ab,bj->ij', smo_c, ovlp, xmo_a, optimize=True)
ssproja = np.einsum('ai,ab,bj->ij', smo_a, ovlp, xmo_a, optimize=True)

pmo_c = mo_coeff[:, :ext_nfrozen + pnfrozen]
pmo_a = mo_coeff[:, ext_nfrozen + pnfrozen:ext_nfrozen + pnfrozen + pncas]
pprojc = np.einsum('ai,ab,bj->ij', pmo_c, ovlp, smo_a, optimize=True)
pproja = np.einsum('ai,ab,bj->ij', pmo_a, ovlp, smo_a, optimize=True)

qmo_c = mo_coeff[:, :ext_nfrozen + qnfrozen]
qmo_a = mo_coeff[:, ext_nfrozen + qnfrozen:ext_nfrozen + qnfrozen + qncas]
qqprojc = np.einsum('ai,ab,bj->ij', qmo_c, ovlp, pmo_a, optimize=True)
qqproja = np.einsum('ai,ab,bj->ij', qmo_a, ovlp, pmo_a, optimize=True)

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
    elif ttx in ['PJ2-OPT', 'PJ2-DIA']:
        if ttx == 'PJ2-OPT':
            umf.verbose = 4
            umf.max_cycle = 200
            umf.conv_tol = 1E-12
            umf.ah_conv_tol = 1E-16
            umf.kernel(dm0=dm0)
            dm1 = np.array(umf.make_rdm1())
        elif ttx == 'PJ2-DIA':
            dm1 = dm0
        
        dm0 = np.einsum('ia,jb,xij->xab', sproja, sproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', sprojc, sprojc, [1, 1], optimize=True)

        smf.verbose = 4
        smf.max_cycle = 200
        smf.conv_tol = 1E-12
        smf.ah_conv_tol = 1E-16
        smf.kernel(dm0=dm0)
        dm1 = np.array(smf.make_rdm1())
    elif ttx in ['PJ2-PJ1-OPT', 'PJ2-PJ1-DIA']:
        if ttx == 'PJ2-PJ1-OPT':
            umf.verbose = 4
            umf.max_cycle = 200
            umf.conv_tol = 1E-12
            umf.ah_conv_tol = 1E-16
            umf.kernel(dm0=dm0)
            dm1 = np.array(umf.make_rdm1())
        elif ttx == 'PJ2-PJ1-DIA':
            dm1 = dm0
        
        dm0 = np.einsum('ia,jb,xij->xab', rproja, rproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', rprojc, rprojc, [1, 1], optimize=True)

        pmf.verbose = 4
        pmf.max_cycle = 200
        pmf.conv_tol = 1E-12
        pmf.ah_conv_tol = 1E-16
        pmf.kernel(dm0=dm0)
        dm1 = np.array(pmf.make_rdm1())

        dm0 = np.einsum('ia,jb,xij->xab', pproja, pproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', pprojc, pprojc, [1, 1], optimize=True)

        smf.verbose = 4
        smf.max_cycle = 200
        smf.conv_tol = 1E-12
        smf.ah_conv_tol = 1E-16
        smf.kernel(dm0=dm0)
        dm1 = np.array(smf.make_rdm1())

    elif ttx in ['PJ2-PJ1-PJ0-OPT', 'PJ2-PJ1-PJ0-DIA']:
        if ttx == 'PJ2-PJ1-PJ0-OPT':
            umf.verbose = 4
            umf.max_cycle = 200
            umf.conv_tol = 1E-12
            umf.ah_conv_tol = 1E-16
            umf.kernel(dm0=dm0)
            dm1 = np.array(umf.make_rdm1())
        elif ttx == 'PJ2-PJ1-PJ0-DIA':
            dm1 = dm0
        
        dm0 = np.einsum('ia,jb,xij->xab', qproja, qproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', qprojc, qprojc, [1, 1], optimize=True)

        qmf.verbose = 4
        qmf.max_cycle = 200
        qmf.conv_tol = 1E-12
        qmf.ah_conv_tol = 1E-16
        qmf.kernel(dm0=dm0)
        dm1 = np.array(qmf.make_rdm1())

        dm0 = np.einsum('ia,jb,xij->xab', qqproja, qqproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', qqprojc, qqprojc, [1, 1], optimize=True)

        pmf.verbose = 4
        pmf.max_cycle = 200
        pmf.conv_tol = 1E-12
        pmf.ah_conv_tol = 1E-16
        pmf.kernel(dm0=dm0)
        dm1 = np.array(pmf.make_rdm1())

        dm0 = np.einsum('ia,jb,xij->xab', pproja, pproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', pprojc, pprojc, [1, 1], optimize=True)

        smf.verbose = 4
        smf.max_cycle = 200
        smf.conv_tol = 1E-12
        smf.ah_conv_tol = 1E-16
        smf.kernel(dm0=dm0)
        dm1 = np.array(smf.make_rdm1())

    else:
        assert False

    if ttx in ['PJ2-PJ1-PJ0-OPT', 'PJ2-PJ1-PJ0-DIA', 'PJ2-PJ1-OPT', 'PJ2-PJ1-DIA', 'PJ2-OPT', 'PJ2-DIA']:
        dm0 = np.einsum('ia,jb,xij->xab', ssproja, ssproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', ssprojc, ssprojc, [1, 1], optimize=True)
    else:
        dm0 = np.einsum('ia,jb,xij->xab', proja, proja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', projc, projc, [1, 1], optimize=True)

    mf.verbose = 4
    mf.max_cycle = 200
    mf.conv_tol = 1E-12
    mf.ah_conv_tol = 1E-16
    mf.kernel(dm0=dm0)
    dm1 = np.array(mf.make_rdm1())

    if abs(mf.e_tot - e_ref) >= 1E-6:
        print('!!WARNING', mf.e_tot, e_ref)

    mmp = mp.MP2(mf)
    mmp.kernel(with_t2=False)

    gdm0 = np.einsum('ia,jb,xij->xab', gproja, gproja, dm1, optimize=True)
    gdm0 += np.einsum('ia,ib,x->xab', gprojc, gprojc, [1, 1], optimize=True)
    gdm0 = np.einsum('ai,bj,xij->xab', guno, guno, gdm0, optimize=True)

    gmf.verbose = 4
    gmf.max_cycle = 200
    gmf.conv_tol = 1E-12
    gmf.ah_conv_tol = 1E-16
    gmf.kernel(dm0=gdm0)
    gdm1 = np.array(gmf.make_rdm1())

    gmp = mp.MP2(gmf)
    gmp.kernel(with_t2=False)

    print(":: %5d %s %20.12f %20.12f %20.12f %20.12f %s" % (jx, xdstr, mf.e_tot, mmp.e_tot, gmf.e_tot, gmp.e_tot, ttx))

    break
