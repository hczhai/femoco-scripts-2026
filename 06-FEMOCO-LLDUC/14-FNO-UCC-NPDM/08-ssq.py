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
from pyscf.fci.spin_op import spin_square_general

pfid = 0

for fid in [1, 2, 3, 4, 5, 6, 7]:
    for order in [[2, 3], [2, 3, 4]][fid == 1]:
        for iiid in [25, 26, 27]:
            mo_coeff = np.load("data/%02d-mo-coeff-%02d.npy" % (pfid, iiid))
            dm1 = np.load('data/%02d-cc-%02d-%02d-pdm1.npy' % (fid, order, iiid))
            dm2 = np.load('data/%02d-cc-%02d-%02d-pdm2.npy' % (fid, order, iiid))
            ss, mult = spin_square_general(*dm1, *dm2, mo_coeff, np.identity(mo_coeff[0].shape[1]))
            print('fid %d order = %d iiid = %d <S^2> = %.6f 2S+1 = %.6f' % (fid, order, iiid, ss, mult))
