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

data = json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/10-fe2s2.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))

bs_map = data['name_x']

sci_notation = lambda x: "$%.2f\\times 10^{-%.0f}$" % tuple([float(x) for x in ("%.2E" % x).split('E-')])

e_uccsdt = np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_uccsd = np.array([min([(x['uccsd'], x['uccsdt_conv'], x['uccsd']) for x in data['uccsdt'][str(ix)]])[2] for ix in range(35)])
e_uhf = np.array([min([(x['uhf'], x['uccsdt_conv'], x['uhf']) for x in data['uccsdt'][str(ix)]])[2] for ix in range(35)])

noccs = np.array([data['fno_ucc'][str(ix)][0]['occ'] for ix in [25, 26, 27]])
ccsd_act = np.array([[data['fno_ucc'][str(ix)][0][str(ir)]['uccsd'] for ir in range(1, 8)] for ix in [25, 26, 27]])
ccsdt_act = np.array([[data['fno_ucc'][str(ix)][0][str(ir)]['uccsdt'] for ir in range(1, 8)] for ix in [25, 26, 27]])
ccsdtq_act = np.array([[data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq'] for ir in range(1, 8)] for ix in [25, 26, 27]])

reg = [scipy.stats.linregress(noccs[ii][-3:], (e_uccsdt[ix] + ccsdtq_act[ii] - ccsdt_act[ii])[-3:])
    for ii, ix in enumerate([25, 26, 27])]
intercept_cc = [reg[ii].intercept + reg[ii].slope * 2.00 for ii, _ in enumerate([25, 26, 27])]

fe2_mix_corr = np.array([data['fe2s2']['mix'][x] - data['fe2s2']['mix'][y]
    for x, y in zip(['uccsd', 'uccsdt', 'uccsdtq', 'sa_dmrg_12k'],
        ['uhf', 'uccsd', 'uccsdt', 'uccsdtq'])])

femoco_corr = np.array([[
    float(e_uccsd[ix] - e_uhf[ix]),
    float(e_uccsdt[ix] - e_uccsd[ix]),
    float(intercept_cc[ii] - e_uccsdt[ix])] for ii, ix in enumerate([25, 26, 27])])

for ii, ix in enumerate([25, 26, 27]):
    print('%10s & %10s & %12.6f & %10s & %12.6f & %10s & %8s \\\\' % (
        bs_map[str(ix)], "UHF",
        data['fe2s2']['mix']['uhf'] + 116, '',
        e_uhf[ix] + 22140, '', ''
    ))
    for jj, key in enumerate(['uccsd', 'uccsdt', 'uccsdtq']):
        print('%10s & %10s & %12.6f & %10.6f & %12.6f & %10.6f & %8.2f \\\\' % (
            '', key.upper(),
            data['fe2s2']['mix'][key] + 116,
            fe2_mix_corr[jj],
            {'uccsd': e_uccsd[ix], 'uccsdt': e_uccsdt[ix], 'uccsdtq': intercept_cc[ii]}[key]  + 22140, femoco_corr[ii][jj],
            femoco_corr[ii][jj] / fe2_mix_corr[jj]
        ))
    print('%10s & %10s & %12.6f & %10.6f & \\textbf{estimated} & \\textbf{%.6f} & \\textbf{%.2f} \\\\' % (
        '', 'exact',
        data['fe2s2']['mix']['sa_dmrg_12k'] + 116,
        fe2_mix_corr[-1],
        fe2_mix_corr[-1] * np.average(femoco_corr[ii] / fe2_mix_corr[:3]),
        np.average(femoco_corr[ii] / fe2_mix_corr[:3])
    ))
    if ii != 2:
        print("     \\hline")
