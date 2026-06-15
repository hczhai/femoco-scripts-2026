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


import zipfile
import itertools
import numpy as np, sys, tempfile
from pyscf import lib
from pyscf.tools.fcidump import to_scf
import scipy.linalg

fid = 0
iiid = int(sys.argv[1])

lib.param.TMPDIR = '/tmp'

tmpdir = tempfile.mkdtemp()

with zipfile.ZipFile('../FEMOCO.FCIDUMP.zip', 'r') as zip_ref:
    zip_ref.extract('FEMOCO.FCIDUMP', tmpdir)

mf = to_scf(tmpdir + '/FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

gmf = mf.to_ghf()

n_sites = mf.mol.nao
h1e = mf.get_hcore()
ovlp = mf.get_ovlp()

gmf.get_hcore = lambda *args: scipy.linalg.block_diag(h1e, h1e)
gmf.get_ovlp = lambda *args: scipy.linalg.block_diag(ovlp, ovlp)

gmf = gmf.newton()

fe_idxs = [[2, 3, 4, 5, 6], [16, 17, 18, 19, 20], [21, 22, 23, 24, 25], [26, 27, 28, 29, 30],
    [44, 45, 46, 47, 48], [49, 50, 51, 52, 53], [54, 55, 56, 57, 58]]
mo_idxs = [[65, 69, 70, 71, 72]]
end_idxs = [[0, 1], [73, 74, 75]]
s_idxs = [[7, 8, 9, 10, 11, 12, 13, 14, 15], [59, 60, 61, 62, 63, 64, 68, 66, 67]]
sc_idxs = [[31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43]]

idxs = list(fe_idxs) + list(mo_idxs) + list(end_idxs) + list(s_idxs) + list(sc_idxs)
assert max([max(x) for x in idxs]) == n_sites - 1
assert min([min(x) for x in idxs]) == 0

occs = []

for ix, x in enumerate(itertools.combinations(range(7), 3)):
    y = tuple(z for z in range(7) if z not in x)
    xoccs = []
    for jx in range(3 * len(list(itertools.combinations(range(4), 2)))):
        kajd, kbjd = jx // 6, jx % 6
        occ = ['2'] * n_sites
        for ii, ip in enumerate(mo_idxs[0]):
            occ[ip] = 'bba00'[ii]
        for iv, v in enumerate(x):
            for ii, ip in enumerate(fe_idxs[v]):
                pxocc = ['bbbbb', 'bbbbb', 'bbbbb']
                pxocc[kajd] = 'BBBBB'
                occ[ip] = pxocc[iv][ii]
        for iv, v in enumerate(y):
            for ii, ip in enumerate(fe_idxs[v]):
                pxocc = ['aaaaa', 'aaaaa', 'aaaaa', 'aaaaa']
                pxopt = list(itertools.combinations(range(4), 2))[kbjd]
                pxocc[pxopt[0]] = 'AAAAA'
                pxocc[pxopt[1]] = 'AAAAA'
                occ[ip] = pxocc[iv][ii]
        xoccs.append(''.join(occ))
    occs.append(xoccs)

fe2_a = ['2aaaa', 'a2aaa', 'aa2aa', 'aaa2a', 'aaaa2']
fe2_b = ['2bbbb', 'b2bbb', 'bb2bb', 'bbb2b', 'bbbb2']

xstrs = []
for ii in range(len(occs[iiid])):
    for i, ix in enumerate(fe2_a):
        for j, jx in enumerate(fe2_a):
            for k, kx in enumerate(fe2_b):
                base_xstr = occs[iiid][ii]
                xdstr = base_xstr.replace('AAAAA', ix, 1).replace('AAAAA', jx, 1).replace('BBBBB', kx, 1)
                xstrs.append(xdstr)

data = []

for ix, xdstr in enumerate(xstrs):

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    gdm = scipy.linalg.block_diag(dm0[0], dm0[1])

    np.random.seed(2025)
    gdm += (np.random.random(gdm.shape) - 0.5) * 0.1

    gmf.verbose = 4
    gmf.max_cycle = 800
    gmf.conv_tol = 1E-13
    gmf.ah_conv_tol = 1E-16
    gmf.kernel(dm0=gdm)

    dm1 = np.array(gmf.make_rdm1())

    found, ifx = -1, len(data)
    for jx, (kx, xener, xdm1) in enumerate(data):
        if abs(xener - gmf.e_tot) < 1E-8:
            found, ifx = kx, jx
            break
    if found == -1:
        data.append((ix, gmf.e_tot, dm1))

    print(":: %5d %s %20.12f %5d %5d" % (ix, xdstr, gmf.e_tot, found, ifx))
