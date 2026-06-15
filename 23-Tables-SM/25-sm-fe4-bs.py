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


import numpy as np, json
import scipy.stats

data = json.load(open("../20-Data/19-fe4s4-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/20-fe4s4-uccsd.json", "r", encoding='utf-8'))

for ix in range(3):
    e_uhf = np.array([x[0] for x in sorted([(x['uhf'], ) for x in data['fe4s4-uccsd'][str(ix)]])])
    e_uccsd = np.array([x[1] for x in sorted([(x['uhf'], x['uccsd']) for x in data['fe4s4-uccsd'][str(ix)]])])
    e_uccsd_t = np.array([x[1] for x in sorted([(x['uhf'], x['uccsd(t)']) for x in data['fe4s4-uccsd'][str(ix)]])])
    if ix == 0:
        print('\\multirow{11}{*}{Fe1$\\uparrow$, Fe2$\\uparrow$, Fe3$\\downarrow$, Fe4$\\downarrow$}')
    elif ix == 1:
        print('\\multirow{9}{*}{Fe1$\\uparrow$, Fe2$\\downarrow$, Fe3$\\uparrow$, Fe4$\\downarrow$}')
    elif ix == 2:
        print('\\multirow{6}{*}{Fe1$\\uparrow$, Fe2$\\downarrow$, Fe3$\\downarrow$, Fe4$\\uparrow$}')
    mr = np.argmin(e_uccsd)
    for ir in range(len(e_uhf)):
        if ir == mr:
            print('& \\textbf{%4s} & \\textbf{%8.6f} & \\textbf{%8.6f} & \\textbf{%8.6f} \\\\' % ('BS%d-%d' % (ix + 1, ir),
                e_uhf[ir], e_uccsd[ir], e_uccsd_t[ir], ))
        else:
            print('& %10s & %12.6f & %12.6f & %12.6f \\\\' % ('BS%d-%d' % (ix + 1, ir),
                e_uhf[ir], e_uccsd[ir], e_uccsd_t[ir], ))
    if ix != 2:
        print('\\hline')
