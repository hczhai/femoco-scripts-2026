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
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))

bs_map = data['name_x']

sci_notation = lambda x: "$%.2f\\times 10^{%.0f}$" % tuple([float(x) for x in ("%.2E" % x).split('E+')])

noccs_26 = [58 - data['fno_ucc']['26'][0][str(ir)]['nfrozen'] for ir in range(1, 8)]
times_26 = [data['fno_ucc']['26'][0][str(ir)]['uccsdtq_time'] * (96 if ir >= 5 else 48)
    / 3600 / data['fno_ucc']['26'][0][str(ir)]['uccsdtq_niter'] for ir in range(1, 8)]
flops_26 = data['fno_ucc']['26'][0]['flops']

reg = scipy.stats.linregress(np.log(noccs_26[-2:]), np.log(times_26[-2:]))
cpu_hours_58 = np.exp(reg.intercept + reg.slope * np.log(58))

reg = scipy.stats.linregress(np.log(noccs_26[-2:]), np.log(flops_26[-2:]))
flop_58 = np.exp(reg.intercept + reg.slope * np.log(58))

for ii, ix in enumerate([25, 26, 27]):
    for ir in range(1, 8):
        print('%10s & %10s & %4d & %5d & %20s & %6.2f hours & %9.1f \\\\' % (
            bs_map[str(ix)] if ir == 1 else '',
            '(%2do, %2de)' % (26 + ir * 4, 13 + ir * 8), 1,
            data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq_niter'],
            sci_notation(data['fno_ucc'][str(ix)][0]['flops'][ir - 1]),
            data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq_time']
                / 3600 / data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq_niter'],
            data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq_time'] * (96 if ir >= 5 else 48)
                / 3600 / data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq_niter'],
        ))
    print("     \\hline")
print('     \\textbf{estimated} & (76o, 113e) & & & %20s & & %.0f \\\\' % (
    sci_notation(flop_58), cpu_hours_58
))
