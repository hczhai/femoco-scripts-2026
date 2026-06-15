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


from pyscf import gto, scf, mp
import numpy as np

fid = 19
ofid = 18

mo_coeff = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-mo-coeff.npy' % (ofid, 10))
ovlp = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (ofid, 10))
orig_mo_coeff = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-orig-mo-coeff.npy' % (ofid, 11))

print('load done')

act_76 = np.array([
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

act_76_frozen = 517
assert sorted(act_76) == list(range(76))

orig_cmo = orig_mo_coeff[:, :act_76_frozen]
orig_amo = orig_mo_coeff[:, act_76_frozen + act_76]

for im in range(4):

    ncas = [117, 180, 285, 404][im]
    n_elec = [105, 143, 217, 277][im]
    nfrozen = [159, 140, 103, 73][im]
    ext_nfrozen = 362
    spin = 35
    assert nfrozen * 2 + n_elec == 423

    xmo = mo_coeff[:, ext_nfrozen + nfrozen:ext_nfrozen + nfrozen + ncas]

    projc = np.einsum('ai,ab,bj->ij', orig_cmo, ovlp, xmo, optimize=True)
    proja = np.einsum('ai,ab,bj->ij', orig_amo, ovlp, xmo, optimize=True)
    print(projc.shape, proja.shape)

    np.save('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-projc.npy' % (ofid, im), projc)
    np.save('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-proja.npy' % (ofid, im), proja)
