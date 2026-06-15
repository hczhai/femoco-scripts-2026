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
import numpy as np, sys, itertools
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.mcscf import pmloc
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes, ParallelTypes
from pyscf import lib
import shutil, os

lib.param.TMPDIR = '/tmp'

fid = 5
pfid = 0
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/tmp/fe4-%02d-%02d" % (fid, iiid), symm_type=SymmetryTypes.SU2,
    mps_dir="/mnt/sdceph/users/hzhai/fe4-new-%02d-%02d" % (fid, iiid),
    restart_dir="/mnt/sdceph/users/hzhai/fe4-new-%02d-%02d-mps" % (fid, iiid),
    stack_mem=500 << 30, min_mpo_mem=False, fp_codec_cutoff=1E-16, n_threads=64, n_mkl_threads=1, mpi=False)

driver.read_fcidump(filename='../Fe4S4.FCIDUMP', pg='d2h')
driver.spin = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

kets = []

restart_dir = "/mnt/sdceph/users/hzhai/fe4-new-%02d-%02d-mps" % (pfid, iiid)

for ir in range(nroots):
    if restart_dir != driver.scratch:
        for k in os.listdir(restart_dir):
            if ('.KET-%d.' % ir) in k or k == 'KET-%d-mps_info.bin' % ir:
                shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)
    ket = driver.load_mps("KET-%d" % ir)
    ket = driver.adjust_mps(ket, dot=2)[0]
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

for ir in range(nroots):

    print('===== STATE %d =====' % ir)
    pdm1 = driver.get_1pdm(kets[ir], iprint=2, max_bond_dim=ket.info.bond_dim)
    np.save('data/%02d-%02d-1pdm.npy' % (fid, iiid), pdm1)
