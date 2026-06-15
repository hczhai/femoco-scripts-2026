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


import os, sys

sys.path[:0] = ['/mnt/home/hzhai/Projects/hast']

import pyhast

import zipfile
import itertools
import numpy as np, sys
from pyscf.tools.fcidump import to_scf
from pyscf import lib

lib.param.TMPDIR = '/tmp'

fid = 9
ofid = 7
pfid = 0
iiid = int(sys.argv[1])
frozen = 22

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

# xdstr = open('../FM-INIT', 'r').readlines()[iiid].strip().split()[-1]
xdstr = "22aaaaa2222222222aaaa2bbbbbbbbb2222222222222bbbbbaaa2aaaaaa222222b222ba00222"

dm0 = np.zeros((2, n_sites, n_sites))
occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]
print(n_sites, dm0.shape, np.trace(dm0[0]), np.trace(dm0[1]))

umf.verbose = 4
umf.max_cycle = 800
umf.kernel(dm0=dm0)

print('EHF =', umf.e_tot)

umf.mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (pfid, iiid))

cc_opts = dict(symm_amps=True, label_symm=True, lift=True, sum_lift=True, split=True, eval_t='hastctr',
    diis_scratch='/mnt/ceph/users/hzhai/fno-cc-36-%02d-%02d' % (fid, iiid))

if not os.path.exists(cc_opts['diis_scratch']):
    os.mkdir(cc_opts['diis_scratch'])

prev_dir = '/mnt/ceph/users/hzhai/fno-cc-36-%02d-%02d/DIIS-UCCSDTQ-X02.npy' % (ofid, iiid)

mcc3 = pyhast.sr.ucc.UCCSDTQ(umf, verbose=4, diis=True, frozen=frozen, **cc_opts)
mcc3.diis_space = 12
mcc3.kernel(tamps=mcc3.vector_to_amplitudes(np.load(prev_dir), mcc3.order), max_cycle=800)
