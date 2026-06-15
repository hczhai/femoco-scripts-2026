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


import sys, numpy as np

for iiid in [27]:

    alpha_occ = [None] * 76
    beta_occ = [None] * 76

    oka, okb = False, False
    ia, ib = 0, 0

    for l in open('00-pre-%d.out' % iiid, 'r').readlines():
        if l.strip() == 'alpha nat occ (corr)':
            oka = True
        elif l.strip() == 'beta nat occ (corr)':
            oka, okb = False, True
        elif l.startswith('real'):
            oka, okb = False, False
        elif oka:
            x = [float(v) for v in l.replace('[', '').replace(']', '').split()]
            alpha_occ[ia:ia + len(x)] = x
            ia += len(x)
        elif okb:
            x = [float(v) for v in l.replace('[', '').replace(']', '').split()]
            beta_occ[ib:ib + len(x)] = x
            ib += len(x)

    assert ia == 58 and ib == 55

    occs = []
    for x in range(22, 46 + 1, 4)[::-1]:
        occs.append(2.0 + (alpha_occ[x] + beta_occ[x]))

    print(', '.join(['%10.6f' % x for x in occs]))
