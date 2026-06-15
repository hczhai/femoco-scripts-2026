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
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes
from pyscf import lib

lib.param.TMPDIR = '/tmp'

fid = 2
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/tmp/fe2-%02d-%02d" % (fid, iiid), symm_type=SymmetryTypes.SU2,
    restart_dir="/mnt/ceph/users/hzhai/fe2-%02d-%02d" % (fid, iiid),
    stack_mem=500 << 30, min_mpo_mem=False, fp_codec_cutoff=1E-16, n_threads=96, n_mkl_threads=1, mpi=True)

bond_dims = [250] * 8 + [500] * 8 + [1000] * 8 + [2000] * 4 + [3000] * 4 + [4000] * 4 + [5000] * 4 \
     + [6000] * 4 + [7000] * 4 + [8000] * 4 + [12000] * 8
noises = [1E-5] * (len(bond_dims) - 4) + [0] * 4
thrds = [1E-7] * len(bond_dims)
n_sweeps = len(bond_dims)

driver.read_fcidump(filename='../Fe2S2.FCIDUMP', pg='d2h')
driver.n_elec = 30
driver.spin = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

mpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=True,
    algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-16)
print("MPO = ", mpo.get_bond_dims())

kets = []

for ir in range(nroots):
    ket = driver.get_random_mps(tag="KET-%d" % ir, dot=1, bond_dim=250)
    driver.mpi.barrier()
    driver.align_mps_center(ket, ref=0)
    driver.mpi.barrier()
    ket = driver.adjust_mps(ket, dot=2)[0]
    driver.mpi.barrier()
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

for ir in range(nroots):

    print('===== STATE %d =====' % ir)

    energy = driver.dmrg(mpo, kets[ir], n_sweeps=n_sweeps,
        bond_dims=bond_dims, noises=noises, thrds=thrds,
        lowmem_noise=True, twosite_to_onesite=None, tol=1E-12, cutoff=1E-24, iprint=2,
        dav_max_iter=1000, dav_def_max_size=50, proj_weights=[5.0] * ir, proj_mpss=kets[:ir])
    print('STATE %d DMRG energy = %20.15f' % (ir, energy))
