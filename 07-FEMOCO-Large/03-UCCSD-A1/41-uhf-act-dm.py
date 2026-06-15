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


from pyscf import gto, scf, cc
import numpy as np, sys, zipfile, tempfile
from pyscf import lib
from pyscf.tools.fcidump import to_scf

tmpdir = tempfile.mkdtemp()

fid = 41
ofid = 18

im = 1
ncas = [117, 180, 285, 404][im]
n_elec = [105, 143, 217, 277][im]
nfrozen = [159, 140, 103, 73][im]
ext_nfrozen = 362

assert nfrozen * 2 + n_elec == 423

ecore = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, im))
h1e = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, im))
g2e = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, im))
proja = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, im))
projc = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, im))

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

lib.param.TMPDIR = tmpdir

pm = 0
pncas = [117, 180, 285, 404][pm]
pn_elec = [105, 143, 217, 277][pm]
pnfrozen = [159, 140, 103, 73][pm]

assert pnfrozen * 2 + pn_elec == 423

pecore = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (ofid, pm))
ph1e = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (ofid, pm))
pg2e = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (ofid, pm))
rproja = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, pm))
rprojc = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, pm))

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

xdstrs = open('05-collected-uhf.out', 'r').readlines()

mo_coeff = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-full-mo-coeff.npy' % (ofid, 10))
ovlp = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (ofid, 10))

xmo = mo_coeff[:, ext_nfrozen + nfrozen:ext_nfrozen + nfrozen + ncas]
pmo_c = mo_coeff[:, :ext_nfrozen + pnfrozen]
pmo_a = mo_coeff[:, ext_nfrozen + pnfrozen:ext_nfrozen + pnfrozen + pncas]
pprojc = np.einsum('ai,ab,bj->ij', pmo_c, ovlp, xmo, optimize=True)
pproja = np.einsum('ai,ab,bj->ij', pmo_a, ovlp, xmo, optimize=True)


data = []

for x in xdstrs:
    ix, jx, xdstr = int(x.split()[0]), int(x.split()[2]), x.split()[3]
    e_ref = float(x.split()[4])
    ttx = x.split()[-1]
    print('--- %5d : %5d %s ---' % (ix, jx, ttx))
    print('xdstr =', xdstr)

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
    elif ttx in ['PJ0-OPT', 'PJ0-DIA']:
        if ttx == 'PJ0-OPT':
            umf.verbose = 4
            umf.max_cycle = 200
            umf.conv_tol = 1E-12
            umf.ah_conv_tol = 1E-16
            umf.kernel(dm0=dm0)
            dm1 = np.array(umf.make_rdm1())
        elif ttx == 'PJ0-DIA':
            dm1 = dm0
        
        dm0 = np.einsum('ia,jb,xij->xab', rproja, rproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', rprojc, rprojc, [1, 1], optimize=True)

        pmf.verbose = 4
        pmf.max_cycle = 200
        pmf.conv_tol = 1E-12
        pmf.ah_conv_tol = 1E-16
        pmf.kernel(dm0=dm0)
        dm1 = np.array(pmf.make_rdm1())

    else:
        assert False

    if ttx in ['PJ0-OPT', 'PJ0-DIA']:
        dm0 = np.einsum('ia,jb,xij->xab', pproja, pproja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', pprojc, pprojc, [1, 1], optimize=True)
    else:
        dm0 = np.einsum('ia,jb,xij->xab', proja, proja, dm1, optimize=True)
        dm0 += np.einsum('ia,ib,x->xab', projc, projc, [1, 1], optimize=True)

    mf.verbose = 4
    mf.max_cycle = 200
    mf.conv_tol = 1E-12
    mf.ah_conv_tol = 1E-16
    mf.kernel(dm0=dm0)
    dm1 = np.array(mf.make_rdm1())

    np.save('data/uhf-mo-coeff-%02d-%04d.npy' % (fid, ix), mf.mo_coeff)
    np.save('data/uhf-dm1-%02d-%04d.npy' % (fid, ix), dm1)

    if abs(mf.e_tot - e_ref) >= 1E-6:
        print('!!WARNING', mf.e_tot, e_ref)

    print("%5d :: %5d %s %20.12f %s" % (ix, jx, xdstr, mf.e_tot, ttx))
