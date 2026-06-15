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
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes, ParallelTypes
from pyscf import lib
import shutil, os

lib.param.TMPDIR = '/tmp'

fid = 22
nroots = 1
iiid = int(sys.argv[1])

pfid = [5, 6, 5][iiid]

driver = DMRGDriver(scratch="/tmp/fe4-%02d-%02d" % (fid, iiid), symm_type=SymmetryTypes.SU2,
    mps_dir="/mnt/ceph/users/hzhai/fe4-newsp-%02d-%02d" % (fid, iiid),
    stack_mem=500 << 30, min_mpo_mem=False, fp_codec_cutoff=1E-16, n_threads=80, n_mkl_threads=1, mpi=True)

driver.read_fcidump(filename='../Fe4S4.FCIDUMP', pg='d2h')
driver.spin = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

restart_dir = "/mnt/ceph/users/hzhai/fe4-newsp-%02d-%02d-mps" % (pfid, iiid)

driver.mpi.barrier()
if driver.mpi.rank == 0 and restart_dir != driver.scratch:
    for k in os.listdir(restart_dir):
        if '.KETX.' in k or k == 'KETX-mps_info.bin':
            shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)
driver.mpi.barrier()
ket = driver.load_mps("KETX")
driver.mpi.barrier()
ket = driver.adjust_mps(ket, dot=1)[0]
driver.mpi.barrier()
print('bond dim = ', ket.info.bond_dim)

print('===== STATE %d =====' % 0)

fe_idxs = [[2, 3, 4, 5, 6], [7, 8, 9, 10, 11], [24, 25, 26, 27, 28], [29, 30, 31, 32, 33]]

dm = driver.get_npdm(ket, pdm_type=2, npdm_expr='((C+D)2+(C+D)2)0', mask=(0, 0, 1, 1), iprint=2, max_bond_dim=20000)
dm = dm * (0.5 * -np.sqrt(3) / 2)
fe_idxs = np.array([x for xx in fe_idxs for x in xx], dtype=int)
dm = np.einsum('ijkl->ik', dm[fe_idxs, :][:, fe_idxs].reshape((4, 5, 4, 5)))

if driver.mpi.rank == 0:
    scr = dm
    for ii in range(len(scr)):
        print('SPIN-CORR %d ::' % ii, ' '.join(['%12.6f' % x for x in scr[ii]]))
    import matplotlib.pyplot as plt
    plt.matshow(dm, cmap='ocean_r')
    plt.gcf().set_dpi(300)
    plt.savefig("22-%02d-bip-spin-corr.png" % iiid, dpi=300)
