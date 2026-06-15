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


import numpy as np, sys
from pyblock2.driver.core import DMRGDriver, SymmetryTypes
from pyscf import lib
import shutil, os

lib.param.TMPDIR = '/tmp'

fid = 16
pfid = 0
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/mnt/sdceph/users/hzhai/tmp/fe4-new-%02d-%02d" % (fid, iiid), symm_type=SymmetryTypes.SU2,
    mps_dir="/mnt/sdceph/users/hzhai/fe4-new-%02d-%02d" % (fid, iiid),
    restart_dir="/mnt/sdceph/users/hzhai/fe4-new-%02d-%02d-mps" % (fid, iiid),
    stack_mem=500 << 30, min_mpo_mem=False, fp_codec_cutoff=1E-16, n_threads=64, n_mkl_threads=1, mpi=True)

driver.read_fcidump(filename='../Fe4S4.FCIDUMP', pg='d2h')
driver.spin = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

kets = []

restart_dir = "/mnt/sdceph/users/hzhai/fe4-new-%02d-%02d-mps" % (pfid, iiid)

for ir in range(nroots):
    driver.mpi.barrier()
    if driver.mpi.rank == 0 and restart_dir != driver.scratch:
        for k in os.listdir(restart_dir):
            if ('.KET-%d.' % ir) in k or k == 'KET-%d-mps_info.bin' % ir:
                shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)
    driver.mpi.barrier()
    ket = driver.load_mps("KET-%d" % ir)
    driver.mpi.barrier()
    ket = driver.adjust_mps(ket, dot=2)[0]
    driver.mpi.barrier()
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

for ir in range(nroots):

    print('===== STATE %d =====' % ir)
    dm = driver.get_npdm(kets[ir], pdm_type=2, npdm_expr='((C+D)2+(C+D)2)0', mask=(0, 0, 1, 1), iprint=2, max_bond_dim=ket.info.bond_dim)
    dm = dm * (0.5 * -np.sqrt(3) / 2)
    if driver.mpi.rank == 0:
        np.save('data/%02d-%02d-spin-corr.npy' % (fid, iiid), dm)
