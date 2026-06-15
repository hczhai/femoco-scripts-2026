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


import json, numpy as np
import scipy.stats

data = json.load(open("../20-Data/19-fe4s4-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/20-fe4s4-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/21-fe4s4-uccsdtq.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/18-fe4s4-sa-dmrg.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/10-fe2s2.json", "r", encoding='utf-8'))

fe2_mix_corr = np.array([data['fe2s2']['mix'][x] - data['fe2s2']['mix'][y]
    for x, y in zip(['uccsd', 'uccsdt', 'uccsdtq', 'sa_dmrg_12k'],
        ['uhf', 'uccsd', 'uccsdt', 'uccsdtq'])])

for ix in range(3):
    e_uhf = np.array([x[0] for x in sorted([(x['uhf'], ) for x in data['fe4s4-uccsd'][str(ix)]])])
    e_uccsd = np.array([x[1] for x in sorted([(x['uhf'], x['uccsd']) for x in data['fe4s4-uccsd'][str(ix)]])])
    e_uccsd_t = np.array([x[1] for x in sorted([(x['uhf'], x['uccsd(t)']) for x in data['fe4s4-uccsd'][str(ix)]])])
    mr = np.argmin(e_uccsd)
    e_uccsdt = data['fe4s4-uccsdtq'][str(ix)][0]['uccsdt']
    e_uccsdtq = data['fe4s4-uccsdtq'][str(ix)][0]['uccsdtq']
    fe4_corr = np.array([e_uccsd[mr] - e_uhf[mr], e_uccsdt - e_uccsd[mr], e_uccsdtq - e_uccsdt])
    print('%10s & %10s & %12.6f & %10s & %12.6f & %10s & %8s \\\\' % (
        'BS%d' % (ix + 1), "UHF",
        data['fe2s2']['mix']['uhf'] + 116, '',
        e_uhf[mr] + 327, '', ''
    ))
    print('%10s & %10s & %12.6f & %10.6f & %12.6f & %10.6f & %8.2f \\\\' % (
        '', 'UCCSD', data['fe2s2']['mix']['uccsd'] + 116, fe2_mix_corr[0],
        e_uccsd[mr] + 327, fe4_corr[0],
        fe4_corr[0] / fe2_mix_corr[0]
    ))
    print('%10s & %10s & %12.6f & %10s & %12.6f & %10s & %8s \\\\' % (
        '', 'UCCSD(T)', data['fe2s2']['mix']['uccsd(t)'] + 116, '',
        e_uccsd_t[mr] + 327, '', '',
    ))
    print('%10s & %10s & %12.6f & %10.6f & %12.6f & %10.6f & %8.2f \\\\' % (
        '', 'UCCSDT', data['fe2s2']['mix']['uccsdt'] + 116, fe2_mix_corr[1],
        e_uccsdt + 327, fe4_corr[1],
        fe4_corr[1] / fe2_mix_corr[1]
    ))
    print('%10s & %10s & %12.6f & %10.6f & %12.6f & %10.6f & %8.2f \\\\' % (
        '', 'UCCSDTQ', data['fe2s2']['mix']['uccsdtq'] + 116, fe2_mix_corr[2],
        e_uccsdtq + 327, fe4_corr[2],
        fe4_corr[2] / fe2_mix_corr[2]
    ))
    print('%10s & %10s & %12.6f & %10.6f & \\textbf{%.6f} & \\textbf{%.6f} & \\textbf{%.2f} \\\\' % (
        '', 'exact',
        data['fe2s2']['mix']['sa_dmrg_12k'] + 116,
        fe2_mix_corr[-1],
        e_uccsdtq + 327 + fe2_mix_corr[-1] * np.average(fe4_corr / fe2_mix_corr[:3]),
        fe2_mix_corr[-1] * np.average(fe4_corr / fe2_mix_corr[:3]),
        np.average(fe4_corr / fe2_mix_corr[:3])
    ))
    if ix != 2:
        print("     \\hline")
