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

fid = 4
pfid = 3
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/tmp/fe4-%02d-%02d" % (fid, iiid), symm_type=SymmetryTypes.SU2,
    mps_dir="/mnt/ceph/users/hzhai/fe4-newsp-%02d-%02d" % (fid, iiid),
    restart_dir="/mnt/ceph/users/hzhai/fe4-newsp-%02d-%02d-mps" % (fid, iiid),
    stack_mem=500 << 30, min_mpo_mem=False, fp_codec_cutoff=1E-16, n_threads=80, n_mkl_threads=1, mpi=True)
driver.prule = driver.bw.bs.ParallelRuleQC(driver.mpi)

bond_dims = [2000] * 8 + [4000] * 8 + [6000] * 4 + [8000] * 4 + [10000] * 4 + [14000] * 4 + [18000] * 4 + [22000] * 4
noises = [1E-5] * len(bond_dims)
thrds = [1E-7] * len(bond_dims)
n_sweeps = len(bond_dims)

driver.read_fcidump(filename='../Fe4S4.FCIDUMP', pg='d2h')
driver.spin = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

mpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=True,
    algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-16, para_type=ParallelTypes.Nothing)
print("MPO = ", mpo.get_bond_dims())

restart_dir = "/mnt/ceph/users/hzhai/fe4-newsp-%02d-%02d" % (pfid, iiid)

driver.mpi.barrier()
if driver.mpi.rank == 0 and restart_dir != driver.scratch:
    for k in os.listdir(restart_dir):
        if '.KETX.' in k or k == 'KETX-mps_info.bin':
            shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)
driver.mpi.barrier()
ket = driver.load_mps("KETX")
driver.mpi.barrier()
print('bond dim = ', ket.info.bond_dim)

def f(name, iprint):
    import resource
    if driver.mpi.rank == 0 and name == "DMRG::sweep::iter.start":
        print("%10s" % driver.bw.b.Parsing.to_size_string(
            1024 * resource.getrusage(resource.RUSAGE_SELF).ru_maxrss), end='', flush=True)

driver.set_callback(f)

print('===== STATE %d =====' % 0)

energy = driver.dmrg(mpo, ket, n_sweeps=n_sweeps,
    bond_dims=bond_dims, noises=noises, thrds=thrds,
    lowmem_noise=True, twosite_to_onesite=None, tol=1E-12, cutoff=1E-24, iprint=2,
    dav_max_iter=1000, dav_def_max_size=50, sweep_start=[33, 25, 32][iiid], forward=[False, False, True][iiid])
print('STATE %d DMRG energy = %20.15f' % (0, energy))
