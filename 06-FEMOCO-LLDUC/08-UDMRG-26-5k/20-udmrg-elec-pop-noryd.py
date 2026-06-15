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


import numpy as np
from pyscf import lib, gto, scf, lo
from pyblock2._pyscf.scf import transform_rdm1_to_ao, transform_rdm1_to_mo

lib.param.TMPDIR = '/tmp'

fid = 20
pfid = 1
qfid = 2

fm_basis = {
    'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH',
    'C': 'def2-SVP', 'H': 'def2-SVP', 'O': 'def2-SVP',
    'N':'def2-SVP'
}

xmol = gto.M(atom="../../22-FM/FM.xyz", basis=fm_basis, charge=-3, spin=35, symmetry=False, verbose=0)
xmol.max_memory = 1200000
xmo = np.load('../../22-FM/llduc-lmo.npy')

ncore = 517
ncas = 76

idx = np.array([
    3, 4,
    45, 73, 71, 65, 36,
    18, 13, 21, 11, 16, 17, 8, 20, 15,
    37, 57, 69, 67, 53,
    59, 51, 49, 64, 39,
    44, 60, 54, 68, 41,
    5, 6, 12, 9, 33, 35, 34, 14, 7, 10, 22, 29, 19,
    62, 42, 46, 48, 50,
    38, 61, 70, 58, 52,
    43, 40, 72, 47, 55,
    25, 26, 24, 23, 32, 31,
    56,
    27, 28, 30,
    63, 66, 74, 75,
    0, 1, 2], dtype=int)

assert sorted(idx) == list(range(76))

xmo[:, ncore:ncore + ncas] = xmo[:, ncore:ncore + ncas][:, idx]

def mulliken_pop_dmao(mol, dmao):
    rmf = scf.RHF(mol)
    ld = lo.orth_ao(mol, "meta_lowdin")
    dm = transform_rdm1_to_mo(dmao, ld, rmf.get_ovlp())
    s = np.eye(dm.shape[-1])
    if isinstance(dm, np.ndarray) and dm.ndim == 2:
        dm = np.array((dm * 0.5, dm * 0.5))
    pop_a = np.einsum("ij,ji->i", dm[0], s).real
    pop_b = np.einsum("ij,ji->i", dm[1], s).real
    rydbglst = lo.nao._core_val_ryd_list(mol)[2]

    nelec_a = np.zeros(mol.natm)
    nelec_b = np.zeros(mol.natm)
    for i, s in enumerate(mol.ao_labels(fmt=None)):
        if i not in rydbglst:
            nelec_a[s[0]] += pop_a[i]
            nelec_b[s[0]] += pop_b[i]
    return (nelec_a + nelec_b)[:8], (nelec_a - nelec_b)[:8]

xdstrs = """
   26 ::   628 22aa2aa222222222aaaaabbbbbbbbbb2222222222222a2aaabbb2baaaaa222222b222ba00222  -22140.306452593275  -22140.350446506614 1
   26 ::   627 22aaaaa2222222222aaaab2bbbbbbbb2222222222222a2aaabbbbbaaaaa222222b222ba00222  -22140.301345811393  -22140.348974068926 1
   26 ::   624 222aaaa222222222aaaaabbbb2bbbbb2222222222222aa2aabbbbbaaaaa222222b222ba00222  -22140.301272142111  -22140.349396899401 1
   26 ::   626 22aa2aa2222222222aaaabbbb2bbbbb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.293610397483  -22140.348079427644 1
   26 ::   625 22a2aaa2222222222aaaab2bbbbbbbb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.291955339446  -22140.342491973453 1
   26 ::   646 222aaaa2222222222aaaabbbbbbb2bb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.290567769043  -22140.355020085412 1
   26 ::   629 22aaaaa222222222aaaa2bbbbbbbb2b2222222222222aaaa2bbbbbaaaaa222222b222ba00222  -22140.286834711475  -22140.340654236712 1
   26 ::   631 22aaaaa222222222a2aaabb2bbbbbbb2222222222222aaaaabbbbbaaaa2222222b222ba00222  -22140.285356996927  -22140.346147141376 1
   26 ::   647 222aaaa2222222222aaaabbb2bbbbbb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.284755506444  -22140.347684014167 1
   26 ::   633 22aaa2a222222222a2aaabb2bbbbbbb2222222222222aaaaabbbbbaaaaa222222b222ba00222  -22140.279397359442  -22140.347944674286 1
   26 ::   632 22aaaaa2222222222aaaabbbbb2bbbb2222222222222aa2aabbbbbaaaaa222222b222ba00222  -22140.276448836965  -22140.343661423893 1
   26 ::   639 22aaaaa222222222aa2aabbbbbbbbb22222222222222aaa2abbbbbaaaaa222222b222ba00222  -22140.275055099130  -22140.357669742039 1
""".strip().split('\n')
assert len(xdstrs) == 12

for ix in range(12):
    xdstr = xdstrs[ix].split()[3]

    mox = np.load('data/%02d-mo-coeff-%02d.npy' % (pfid, ix))
    idx = np.load('data/%02d-reorder-%02d.npy' % (pfid, ix))
    dm1_mo = np.load('data/%02d-%02d-1pdm.npy' % (qfid, ix))
    ridx = np.argsort(idx)
    dm1 = np.einsum('xpi,xij,xqj->xpq', mox, dm1_mo[:, ridx, :][:, :, ridx], mox, optimize='optimal')

    dm1_full = np.zeros((2, ncore + ncas, ncore + ncas))
    dm1_full[:, ncore:, ncore:] = dm1
    dm1_full[:, :ncore, :ncore] = np.identity(ncore)[None, :, :]
    dmao = transform_rdm1_to_ao(dm1_full, xmo[:, :ncore + ncas])
    pop, spn = mulliken_pop_dmao(xmol, dmao)
    print('%5d :: %5d %s' % (26, ix, xdstr), ' '.join('%8.4f' % x for x in pop))
