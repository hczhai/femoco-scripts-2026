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
from pyscf import gto, scf
import numpy as np, sys
import itertools
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.mcscf import get_uno
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes, MPOAlgorithmTypes
from pyblock2.algebra.io import MPSTools

fid = 22
nroots = 1
iiid = int(sys.argv[1])

driver = DMRGDriver(scratch="/central/scratch/hczhai/h2d-fm8-%02d-%02d" % (fid, iiid),
    restart_dir="/central/scratch/hczhai/h2d-fm8-%02d-%02d-mps" % (fid, iiid),
    symm_type=SymmetryTypes.SZ, stack_mem=120 << 30, min_mpo_mem=False,
    fp_codec_cutoff=0.0,
    n_threads=64, n_mkl_threads=1, mpi=False)

bond_dims = [50] * 8 + [100] * 4 + [250] * 4 + [500] * 4 + [1000] * 4 + [1500] * 4 + [2000] * 4 + [3000] * 4
noises = [1E-5] * (len(bond_dims) - 4) + [0] * 4
thrds = [1E-7] * len(bond_dims)
n_sweeps = len(bond_dims)

mf = to_scf('FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None,)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf()

umf.mo_coeff = np.load("%02d-mo-coeff-%02d.npy" % (fid, iiid))

ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_uhf_integrals(umf, ncore=0, ncas=None, g2e_symm=1, pg_symm=False)

print('h1e.shape', [x.shape for x in h1e])
print('g2e.shape', [x.shape for x in g2e])

import block2 as b2

idx = np.load("%02d-reorder-%02d.npy" % (fid, iiid))

print('reorder = ', idx)

h1e = [x[idx][:, idx] for x in h1e]
g2e = [x[idx][:, idx][:, :, idx][:, :, :, idx] for x in g2e]

driver.n_sites = ncas
driver.n_elec = n_elec
driver.ecore = ecore
driver.h1e = h1e
driver.g2e = g2e
driver.orb_sym = orb_sym

driver.spin = 3
twos = 3

driver.initialize_system(n_sites=driver.n_sites, n_elec=driver.n_elec, spin=driver.spin, orb_sym=driver.orb_sym)

mpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=False,
    algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-12)
print("MPO = ", mpo.get_bond_dims())

kets = []

for ir in range(nroots):

    import os, shutil
    restart_dir = "/central/scratch/hczhai/h2d-fm8-%02d-%02d-mps" % (fid, iiid)
    for k in os.listdir(restart_dir):
        if ('.KET-%d.' % ir) in k or k == 'KET-%d-mps_info.bin' % ir:
            shutil.copy(restart_dir + "/" + k, driver.scratch + "/" + k)

    ket = driver.load_mps(tag="KET-%d" % ir)
    ket = driver.adjust_mps(ket, dot=2)[0]
    print(ir, 'bond dim = ', ket.info.bond_dim)
    kets.append(ket)

for ir in range(nroots):

    print('===== STATE %d =====' % ir)

    fwd, sst = 9999, 9999
    with open('%d-dmrg-fm-newton-%d.out' % (fid, iiid), 'r') as ff:
        for ll in ff.readlines():
            if ll.startswith('Sweep ='):
                fwd = ll.split()[6] == 'forward'
                sst = int(ll.split()[2])

    print(fwd, sst)

    energy = driver.dmrg(mpo, kets[ir], n_sweeps=n_sweeps,
        bond_dims=bond_dims, noises=noises, thrds=thrds,
        lowmem_noise=True, twosite_to_onesite=None, tol=1E-12, cutoff=1E-24, iprint=2,
        dav_max_iter=400, dav_def_max_size=20, proj_weights=[5.0] * ir, proj_mpss=kets[:ir],
        forward=fwd, sweep_start=sst)
    print('STATE %d DMRG energy = %20.15f' % (ir, energy))

    kets[ir] = driver.adjust_mps(kets[ir], dot=1)[0]
    driver.align_mps_center(kets[ir], ref=0)

    impo = driver.get_identity_mpo()
    hmpo = driver.get_qc_mpo(h1e=driver.h1e, g2e=driver.g2e, ecore=driver.ecore, iprint=2, add_ident=True,
        algo_type=MPOAlgorithmTypes.Conventional, integral_cutoff=1E-12)

    energy = driver.expectation(kets[ir], hmpo, kets[ir], iprint=2)
    print('STATE %d Expt energy = %20.15f' % (ir, energy))

    dm = driver.get_2pdm(kets[ir], mask=(0, 1, 1, 0), iprint=2, max_bond_dim=3000)
    np.save("%02d-dm-orig-%02d-%02d.npy" % (fid, iiid, ir), dm)

    dm = driver.get_1pdm(kets[ir], iprint=2, max_bond_dim=3000)
    np.save("%02d-1pdm-orig-%02d-%02d.npy" % (fid, iiid, ir), dm)

    hmpo, impo = None, None
