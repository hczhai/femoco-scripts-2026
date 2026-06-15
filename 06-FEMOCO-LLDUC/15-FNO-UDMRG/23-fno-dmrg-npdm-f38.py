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
from pyblock2.driver.core import DMRGDriver, SymmetryTypes
from pyscf.fci.spin_op import spin_square_general

fid = 23
ofid = 3
xfid = 13
pfid = 0
iiid = int(sys.argv[1])
frozen = 38

driver = DMRGDriver(scratch="/tmp/fm-fno-udmrg-new-%02d-%02d" % (fid, iiid),
    mps_dir="/mnt/sdceph/users/hzhai/fm-fno-udmrg-new-%02d-%02d" % (fid, iiid),
    restart_dir="/mnt/sdceph/users/hzhai/fm-fno-udmrg-new-%02d-%02d-mps" % (fid, iiid),
    symm_type=SymmetryTypes.SZ, stack_mem=300 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-16, n_threads=64, n_mkl_threads=1, mpi=True)

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


driver.n_sites = ncas
driver.n_elec = n_elec
driver.orb_sym = orb_sym

driver.spin = 3
twos = 3

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

restart_dir = "/mnt/sdceph/users/hzhai/fm-fno-udmrg-new-%02d-%02d-mps" % (xfid, iiid)

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

pdm1 = driver.get_1pdm(ket, iprint=2, max_bond_dim=ket.info.bond_dim)
pdm2 = driver.get_2pdm(ket, iprint=2, max_bond_dim=ket.info.bond_dim)
pdm2 = [x.transpose(0, 3, 1, 2) for x in pdm2]

xpdm1 = [np.pad(x, (frozen, 0)) for x in pdm1]
xpdm2 = [np.pad(x, (frozen, 0)) for x in pdm2]

for i in range(frozen):
    xpdm2[0][i, i, :, :] += xpdm1[0]
    xpdm2[0][:, :, i, i] += xpdm1[0]
    xpdm2[0][:, i, i, :] -= xpdm1[0]
    xpdm2[0][i, :, :, i] -= xpdm1[0].T
    xpdm2[1][i, i, :, :] += xpdm1[1]
for i in range(frozen):
    xpdm2[2][i, i, :, :] += xpdm1[1]
    xpdm2[2][:, :, i, i] += xpdm1[1]
    xpdm2[2][:, i, i, :] -= xpdm1[1]
    xpdm2[2][i, :, :, i] -= xpdm1[1].T
    xpdm2[1][:, :, i, i] += xpdm1[0]

for i in range(frozen):
    for j in range(frozen):
        xpdm2[0][i, i, j, j] += 1
        xpdm2[0][i, j, j, i] -= 1
for i in range(frozen):
    for j in range(frozen):
        xpdm2[2][i, i, j, j] += 1
        xpdm2[2][i, j, j, i] -= 1
for i in range(frozen):
    for j in range(frozen):
        xpdm2[1][i, i, j, j] += 1

for i in range(frozen):
    xpdm1[0][i, i] += 1
    xpdm1[1][i, i] += 1

np.save('data/%02d-%02d-1pdm.npy' % (fid, iiid), xpdm1)
np.save('data/%02d-%02d-2pdm.npy' % (fid, iiid), xpdm2)

mo_coeff = np.copy(umf.mo_coeff)
mo_coeff[:, :, frozen:] = umf.mo_coeff[:, :, frozen:][:, :, idx]
ss, mult = spin_square_general(*xpdm1, *xpdm2, mo_coeff, np.identity(mo_coeff[0].shape[1]))
print('SSQ =', ss, mult)
