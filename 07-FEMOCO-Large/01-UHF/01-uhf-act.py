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

fid = 1
ofid = 18
iiid = int(sys.argv[1])

im = 1
ncas = [117, 180, 285, 404][im]
n_elec = [105, 143, 217, 277][im]
nfrozen = [159, 140, 103, 73][im]

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

xdstrs = open('00-fm-init-76.txt', 'r').readlines()

data = []

for x in xdstrs:
    ix, jx, xdstr = int(x.split()[0]), int(x.split()[2]), x.split()[3]
    if ix != iiid:
        continue
    print('--- %5d : %5d ---' % (ix, jx))

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    umf.verbose = 4
    umf.max_cycle = 200
    umf.conv_tol = 1E-12
    umf.ah_conv_tol = 1E-16
    umf.kernel(dm0=dm0)
    dm1 = np.array(umf.make_rdm1())

    dm0 = np.einsum('ia,jb,xij->xab', proja, proja, dm1, optimize=True)
    dm0 += np.einsum('ia,ib,x->xab', projc, projc, [1, 1], optimize=True)

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
