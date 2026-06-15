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

pfid = 16
iiid = 0

fe_ix = [[2, 3, 4, 5, 6], [7, 8, 9, 10, 11], [24, 25, 26, 27, 28], [29, 30, 31, 32, 33]]
fe_idxs = np.array([x for xx in fe_ix for x in xx], dtype=int)
dm = np.load('data/%02d-%02d-spin-corr.npy' % (pfid, iiid))

scr = np.einsum('ijkl->ik', dm[fe_idxs, :][:, fe_idxs].reshape((4, 5, 4, 5)))
for ii in range(len(scr)):
    print('SPIN-CORR %d ::' % ii, ' '.join(['%12.6f' % x for x in scr[ii]]))
