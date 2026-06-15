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

data = json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))

bs_map = data['name_x']

e_uccsdt = np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_uccsd = np.array([min([(x['uccsdt'], x['uccsdt_conv'], x['uccsd']) for x in data['uccsdt'][str(ix)]])[2] for ix in range(35)])

noccs = np.array([data['fno_ucc'][str(ix)][0]['occ'] for ix in [25, 26, 27]])
ccsd_act = np.array([[data['fno_ucc'][str(ix)][0][str(ir)]['uccsd'] for ir in range(1, 8)] for ix in [25, 26, 27]])
ccsdt_act = np.array([[data['fno_ucc'][str(ix)][0][str(ir)]['uccsdt'] for ir in range(1, 8)] for ix in [25, 26, 27]])
ccsdtq_act = np.array([[data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq'] for ir in range(1, 8)] for ix in [25, 26, 27]])

reg = [scipy.stats.linregress(noccs[ii][-3:], (e_uccsdt[ix] + ccsdtq_act[ii] - ccsdt_act[ii])[-3:])
    for ii, ix in enumerate([25, 26, 27])]
intercept_cc = [reg[ii].intercept + reg[ii].slope * 2.00 for ii, _ in enumerate([25, 26, 27])]

for ii, ix in enumerate([25, 26, 27]):
    for ir in range(1, 8):
        print('%10s & %10s & %10.6f & %10.6f & %10.6f & %10.6f \\\\' % (
            bs_map[str(ix)] if ir == 1 else '',
            '(%2do, %2de)' % (26 + ir * 4, 13 + ir * 8),
            data['fno_ucc'][str(ix)][0][str(ir)]['uccsd'] + 22140,
            data['fno_ucc'][str(ix)][0][str(ir)]['uccsdt'] + 22140,
            data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq'] + 22140,
            data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq'] - data['fno_ucc'][str(ix)][0][str(ir)]['uccsdt']
            + e_uccsdt[ix] + 22140,
        ))
    print('%10s & %10s & %10.6f & %10.6f & \\textbf{extrapolated} & \\textbf{%.6f} \\\\' % (
        '', '(76o, 113e)', e_uccsd[ix] + 22140, e_uccsdt[ix] + 22140, intercept_cc[ii] + 22140))
    if ii != 2:
        print("     \\hline")
