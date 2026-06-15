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

ncore = 517
ncas = 76

mo_coeff = np.load('data/qmmm-lmo.npy')

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

ovlp = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (18, 10))
orig_mo_coeff = np.load('/mnt/sdceph/users/hzhai/%02d-fm-active/%02d-orig-mo-coeff.npy' % (18, 11))

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

assert sorted(act_76) == list(range(76))

orig_mo_coeff[:, ncore:ncore + ncas] = orig_mo_coeff[:, ncore:ncore + ncas][:, act_76]
zmo = orig_mo_coeff[:, ncore:ncore + ncas]

print(np.einsum('ai,ab,bj->ij', xmo, ovlp, zmo, optimize=True))
print(np.argmax(np.abs(np.einsum('ai,ab,bj->ij', xmo, ovlp, zmo, optimize=True)), axis=0))
assert np.argmax(np.abs(np.einsum('ai,ab,bj->ij', xmo, ovlp, zmo, optimize=True)), axis=0).tolist() == list(range(0, 76))
assert all(x >= 0.99 for x in np.max(np.abs(np.einsum('ai,ab,bj->ij', xmo, ovlp, zmo, optimize=True)), axis=0))
