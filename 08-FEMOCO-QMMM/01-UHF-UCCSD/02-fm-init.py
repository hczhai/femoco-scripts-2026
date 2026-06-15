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
from pyscf import lib, gto, scf, lo
from pyscf.tools.fcidump import to_scf
from pyblock2._pyscf.scf import transform_rdm1_to_ao, transform_rdm1_to_mo

fid = 2
iiid = int(sys.argv[1])

lib.param.TMPDIR = '/tmp'

tmpdir = tempfile.mkdtemp()

with zipfile.ZipFile('data/FEMOCO.QMMM.FCIDUMP.zip', 'r') as zip_ref:
    zip_ref.extract('FEMOCO.QMMM.FCIDUMP', tmpdir)

mf = to_scf(tmpdir + '/FEMOCO.QMMM.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf().newton()

n_sites = umf.mol.nao

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

mo_coeff = np.load('data/qmmm-lmo.npy')
ncore = 517
ncas = 76

idx = np.array([
    4, 3, 45, 73, 70, 64, 36, 17, 12, 23,
    15, 19, 16, 9, 18, 13, 38, 57, 69, 66,
    52, 60, 53, 50, 65, 39, 42, 59, 54, 68,
    41, 5, 6, 11, 10, 33, 35, 34, 14, 8,
    7, 24, 31, 20, 63, 43, 46, 48, 49, 37,
    61, 71, 58, 51, 44, 40, 72, 47, 55, 25,
    26, 22, 21, 32, 29, 56, 28, 27, 30, 62,
    67, 74, 75, 0, 1, 2,
], dtype=int)
assert sorted(idx) == list(range(76))

mo_coeff[:, ncore:ncore + ncas] = mo_coeff[:, ncore:ncore + ncas][:, idx]
xmo = mo_coeff[:, ncore:ncore + ncas]

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

xmol = gto.M(atom="../FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=0)
xmol.max_memory = 1200000

def mulliken_pop_dmao(mol, dmao):
    rmf = scf.RHF(mol)
    ld = lo.orth_ao(mol, "meta_lowdin")
    dm = transform_rdm1_to_mo(dmao, ld, rmf.get_ovlp())
    s = np.eye(dm.shape[-1])
    if isinstance(dm, np.ndarray) and dm.ndim == 2:
        dm = np.array((dm * 0.5, dm * 0.5))
    pop_a = np.einsum("ij,ji->i", dm[0], s).real
    pop_b = np.einsum("ij,ji->i", dm[1], s).real
    nelec_a = np.zeros(mol.natm)
    nelec_b = np.zeros(mol.natm)
    for i, s in enumerate(mol.ao_labels(fmt=None)):
        nelec_a[s[0]] += pop_a[i]
        nelec_b[s[0]] += pop_b[i]
    return (nelec_a - nelec_b)[:8]

bs_map = {
    tuple(int(i not in ixx) for i in range(7)): jx
    for jx, ixx in enumerate(itertools.combinations(range(7), 3))
}

data = []

for ix, xdstr in enumerate(xstrs):

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    umf.verbose = 4
    umf.max_cycle = 800
    umf.conv_tol = 1E-13
    umf.ah_conv_tol = 1E-16
    umf.kernel(dm0=dm0)
    dm1 = np.array(umf.make_rdm1())

    found, ifx = -1, len(data)
    for jx, (kx, xener, xdm1) in enumerate(data):
        if abs(xener - umf.e_tot) < 1E-8 and np.linalg.norm(dm1 - xdm1) < 1E-4:
            found, ifx = kx, jx
            break
    if found == -1:
        data.append((ix, umf.e_tot, dm1))

    dmao = transform_rdm1_to_ao(dm1, xmo)
    spn = mulliken_pop_dmao(xmol, dmao)
    spn_sign = tuple(int(xp >= 0) for xp in spn)[:7]
    print(spn.tolist(), spn_sign)
    print()

    print("%5d :: %5d %s %20.12f %5d %5d" % (bs_map.get(spn_sign, -1), ix, xdstr, umf.e_tot, found, ifx))
