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
from pyscf import lib

lib.param.TMPDIR = '/tmp'

fid = 7
pfid = 4
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/tmp/fe2-%02d-%02d" % (fid, iiid), symm_type=SymmetryTypes.SZ,
    restart_dir_per_sweep="/mnt/ceph/users/hzhai/fe2-%02d-%02d" % (fid, iiid),
    stack_mem=300 << 30, min_mpo_mem=False,
    fp_codec_cutoff=1E-16, n_threads=96, n_mkl_threads=1, mpi=False)

bond_dims = [2400] * 4 + [2000] * 4 + [1600] * 4 + [1200] * 4 + [1000] * 4 + \
    [900] * 4 + [800] * 4 + [700] * 4 + [600] * 4 + [500] * 4 + [400] * 4 + [300] * 4 + [200] * 4 + [100] * 4
noises = [0] * len(bond_dims)
thrds = [1E-7] * len(bond_dims)
n_sweeps = len(bond_dims)

mf = to_scf('../Fe2S2.FCIDUMP', molpro_orbsym=True, mf=None,)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf()

umf.mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (pfid, iiid))

ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_uhf_integrals(umf, ncore=0, ncas=None, g2e_symm=1, pg_symm=False)

idx = np.load("data/%02d-reorder-%02d.npy" % (pfid, iiid))
print('reorder = ', idx)

h1e = [x[idx][:, idx] for x in h1e]
g2e = [x[idx][:, idx][:, :, idx][:, :, :, idx] for x in g2e]

driver.n_sites = ncas
driver.n_elec = n_elec
driver.ecore = ecore
driver.h1e = h1e
driver.g2e = g2e
driver.orb_sym = orb_sym

driver.spin = 0
twos = 0

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

mpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=False,
    algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-12)
print("MPO = ", mpo.get_bond_dims())

kets = []

restart_dir = "/mnt/ceph/users/hzhai/fe2-%02d-%02d" % (pfid, iiid)

for ir in range(nroots):
    import os, shutil
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

    energy = driver.dmrg(mpo, kets[ir], n_sweeps=n_sweeps,
        bond_dims=bond_dims, noises=noises, thrds=thrds,
        lowmem_noise=True, twosite_to_onesite=None, tol=0.0, cutoff=1E-24, iprint=2,
        dav_max_iter=400, dav_def_max_size=50, proj_weights=[5.0] * ir, proj_mpss=kets[:ir])
    print('STATE %d DMRG energy = %20.15f' % (ir, energy))
