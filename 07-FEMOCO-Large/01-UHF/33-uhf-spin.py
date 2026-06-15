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
import numpy as np, sys, zipfile, itertools, tempfile
from pyscf import lib
from pyscf.tools.fcidump import to_scf

fid = 33
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
proja = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, im))
projc = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, im))

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

tmpdir = tempfile.mkdtemp()

with zipfile.ZipFile('../FEMOCO.FCIDUMP.zip', 'r') as zip_ref:
    zip_ref.extract('FEMOCO.FCIDUMP', tmpdir)

umf = to_scf(tmpdir + '/FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None)
intor_symmetric = umf.mol.intor_symmetric
umf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(umf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)
umf.mol.groupname = 'c1'
umf.mol.symmetry = True
umf.mol.verbose = 4
umf = umf.to_uhf().newton()
n_sites = umf.mol.nao

xdstrs = open('23-uhf-collect-all.out', 'r').readlines()

mo_coeff = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-mo-coeff.npy' % (ofid, 10))
ovlp = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (ofid, 10))

xmo = mo_coeff[:, ext_nfrozen + nfrozen:ext_nfrozen + nfrozen + ncas]
smo_c = mo_coeff[:, :ext_nfrozen + snfrozen]
smo_a = mo_coeff[:, ext_nfrozen + snfrozen:ext_nfrozen + snfrozen + sncas]
ssprojc = np.einsum('ai,ab,bj->ij', smo_c, ovlp, xmo, optimize=True)
ssproja = np.einsum('ai,ab,bj->ij', smo_a, ovlp, xmo, optimize=True)

pmo_c = mo_coeff[:, :ext_nfrozen + pnfrozen]
pmo_a = mo_coeff[:, ext_nfrozen + pnfrozen:ext_nfrozen + pnfrozen + pncas]
pprojc = np.einsum('ai,ab,bj->ij', pmo_c, ovlp, smo_a, optimize=True)
pproja = np.einsum('ai,ab,bj->ij', pmo_a, ovlp, smo_a, optimize=True)

qmo_c = mo_coeff[:, :ext_nfrozen + qnfrozen]
qmo_a = mo_coeff[:, ext_nfrozen + qnfrozen:ext_nfrozen + qnfrozen + qncas]
qqprojc = np.einsum('ai,ab,bj->ij', qmo_c, ovlp, pmo_a, optimize=True)
qqproja = np.einsum('ai,ab,bj->ij', qmo_a, ovlp, pmo_a, optimize=True)

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

xmol = gto.M(atom="../FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=0)
xmol.max_memory = 1200000

from pyblock2._pyscf.scf import transform_rdm1_to_ao, transform_rdm1_to_mo

def mulliken_pop_dmao(mol, dmao):
    from pyscf import scf, lo
    rmf = scf.RHF(mol)
    ld = lo.orth_ao(mol, "meta_lowdin")
    dm = transform_rdm1_to_mo(dmao, ld, rmf.get_ovlp())
    s = np.eye(dm.shape[-1])
    if isinstance(dm, np.ndarray) and dm.ndim == 2:
        dm = np.array((dm * 0.5, dm * 0.5))
    pop_a = np.einsum("ij,ji->i", dm[0], s).real
    pop_b = np.einsum("ij,ji->i", dm[1], s).real
    nelec_a = np.zeros(mol.natm)
    nelec_b = np.zeros(mol.natm)
    for i, s in enumerate(mol.ao_labels(fmt=None)):
        nelec_a[s[0]] += pop_a[i]
        nelec_b[s[0]] += pop_b[i]
    return (nelec_a - nelec_b)[:8]

bs_map = {
    tuple(int(i not in ixx) for i in range(7)): jx
    for jx, ixx in enumerate(itertools.combinations(range(7), 3))
}

data = []

for x in xdstrs:
    if '::' not in x:
        continue
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

    dmao = transform_rdm1_to_ao(dm1, xmo)
    spn = mulliken_pop_dmao(xmol, dmao)
    spn_sign = tuple(int(xp >= 0) for xp in spn)[:7]
    print(spn.tolist(), spn_sign)
    print()
    print("%5d :: %5d %s %20.12f %5s %5s %5s BS %5d" % (ix, jx, xdstr, mf.e_tot,
        x.split()[5], x.split()[6], x.split()[7], bs_map.get(spn_sign, -1)))
