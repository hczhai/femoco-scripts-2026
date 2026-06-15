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

pfid = 1
iiid = 0

npdm_dir = './data'

for order in [2, 3, 4, 5]:
    mo_coeff = np.load(npdm_dir + '/%02d-mo-coeff-%02d.npy' % (pfid, iiid))
    dm1 = np.load(npdm_dir + '/%02d-cc-%02d-%02d-pdm1.npy' % (pfid, order, iiid))
    dm2 = np.load(npdm_dir + '/%02d-cc-%02d-%02d-pdm2.npy' % (pfid, order, iiid))
    ss, mult = spin_square_general(*dm1, *dm2, mo_coeff, np.identity(mo_coeff[0].shape[1]))
    print('order = %d <S^2> = %.6f S = %.6f' % (order, ss, (mult - 1)))

for order in [2, 3, 4, 5]:
    mo_coeff = np.load(npdm_dir + '/%02d-mo-coeff-%02d.npy' % (pfid, iiid))
    dm1 = np.load(npdm_dir + '/%02d-cc-%02d-%02d-pdm1.npy' % (pfid, order, iiid))
    dmao = np.einsum('xpi,xij,xqj->xpq', mo_coeff, dm1, mo_coeff, optimize=True)
    scr = np.zeros((2, ))
    for ii, ix in enumerate([[2, 3, 4, 5, 6], [13, 14, 15, 16, 17]]):
        for ki in ix:
            scr[ii] += dmao[0][ki, ki] - dmao[1][ki, ki]
    print("SZ :: order = %5d %s" % (order, ' '.join('%.6f' % x for x in scr)))
