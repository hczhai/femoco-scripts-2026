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
import numpy as np, sys
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes
from pyscf.fci.spin_op import spin_square_general
from pyscf import lib

lib.param.TMPDIR = '/tmp'

fid = 9
ofid = 4
pfid = 7
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/tmp/fe2-%02d-%02d" % (fid, iiid), symm_type=SymmetryTypes.SZ,
    restart_dir="/mnt/ceph/users/hzhai/fe2-%02d-%02d" % (fid, iiid),
    stack_mem=300 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-16, n_threads=96, n_mkl_threads=1, mpi=False)

bond_dims = [2400, 2000, 1600, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]

mf = to_scf('../Fe2S2.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf()

umf.mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (ofid, iiid))

ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_uhf_integrals(umf, ncore=0, ncas=None, g2e_symm=1, pg_symm=False)

idx = np.load("data/%02d-reorder-%02d.npy" % (ofid, iiid))
print('reorder = ', idx)

driver.n_sites = ncas
driver.n_elec = n_elec
driver.ecore = ecore

driver.spin = 0
twos = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

for isw, bd in enumerate(bond_dims):
    kets = []
    restart_dir = "/mnt/ceph/users/hzhai/fe2-%02d-%02d.%d" % (pfid, iiid, isw * 4 + 3)
    for ir in range(nroots):
        import os, shutil
        if restart_dir != driver.scratch:
            for k in os.listdir(restart_dir):
                if ('.KET-%d.' % ir) in k or k == 'KET-%d-mps_info.bin' % ir:
                    shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)
        ket = driver.load_mps("KET-%d" % ir)
        ket = driver.adjust_mps(ket, dot=1)[0]
        print(ir, 'bond dim = ', ket.info.bond_dim)
        kets.append(ket)

        kets[ir] = driver.adjust_mps(kets[ir], dot=1)[0]
        dm1 = driver.get_1pdm(kets[ir], iprint=2)
        dm2 = driver.get_2pdm(kets[ir], iprint=2)
        dm2 = [x.transpose(0, 3, 1, 2) for x in dm2]

        mo_coeff = umf.mo_coeff[:, :, idx]
        ss, mult = spin_square_general(*dm1, *dm2, mo_coeff, np.identity(mo_coeff[0].shape[1]))
        print('D :: S ', bd, ss, mult)

        dmao = np.einsum('xpi,xij,xqj->xpq', mo_coeff, dm1, mo_coeff, optimize=True)

        scr = np.zeros((2, ))
        for ii, ix in enumerate([[2, 3, 4, 5, 6], [13, 14, 15, 16, 17]]):
            for ki in ix:
                scr[ii] += dmao[0][ki, ki] - dmao[1][ki, ki]
        print("D :: SZ ", bd, ' '.join('%.6f' % x for x in scr))
