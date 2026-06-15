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


import numpy as np, sys, os, shutil
from pyblock2.driver.core import DMRGDriver, SymmetryTypes
from pyscf import lib

lib.param.TMPDIR = '/tmp'

fid = 14
pfid = 2
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/tmp/fe2-%02d-%02d" % (fid, iiid), symm_type=SymmetryTypes.SU2,
    restart_dir="/mnt/ceph/users/hzhai/fe2-%02d-%02d" % (fid, iiid),
    stack_mem=500 << 30, min_mpo_mem=False, fp_codec_cutoff=1E-16, n_threads=96, n_mkl_threads=1, mpi=True)

driver.read_fcidump(filename='../Fe2S2.FCIDUMP', pg='d2h')
driver.n_elec = 30
driver.spin = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

kets = []

restart_dir = "/mnt/ceph/users/hzhai/fe2-%02d-%02d" % (pfid, iiid)

for ir in range(nroots):
    if restart_dir != driver.scratch:
        for k in os.listdir(restart_dir):
            if ('.KET-%d.' % ir) in k or k == 'KET-%d-mps_info.bin' % ir:
                shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)
    ket = driver.load_mps("KET-%d" % ir)
    ket = driver.adjust_mps(ket, dot=1)[0]
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

for ir in range(nroots):

    print('===== STATE %d =====' % ir)

    pdm1 = driver.get_1pdm(kets[ir], iprint=2, max_bond_dim=12000)
    np.save('%02d-fe2-1pdm.npy' % fid, pdm1)
    pdm2 = driver.get_2pdm(kets[ir], iprint=2, max_bond_dim=12000).transpose(0, 3, 1, 2)
    np.save('%02d-fe2-2pdm.npy' % fid, pdm2)
