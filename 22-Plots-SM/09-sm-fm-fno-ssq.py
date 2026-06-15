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
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(8, 6))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))

bs_map = data['name_x']
colors = data['colors']

# 2S = (1 + 4 * <S^2>) ** 0.5 - 1
fsexpt = lambda x: ((1 + 4 * x) ** 0.5 - 1) / 2

noccs = np.array([data['fno_ucc'][str(ix)][0]['occ'] for ix in [25, 26, 27]])
ccsd_act = np.array([[data['fno_ucc'][str(ix)][0][str(ir)]['uccsd_ssq'][0] for ir in range(1, 8)] for ix in [25, 26, 27]])
ccsdt_act = np.array([[data['fno_ucc'][str(ix)][0][str(ir)]['uccsdt_ssq'][0] for ir in range(1, 8)] for ix in [25, 26, 27]])

ccsd_act = fsexpt(ccsd_act)
ccsdt_act = fsexpt(ccsdt_act)

# nocc_25 = np.array([1.939955, 1.964367, 1.97616, 1.984085, 1.989804, 1.992379, 1.994594])
# nocc_26 = np.array([1.940123, 1.965564, 1.977781, 1.985933, 1.990291, 1.992982, 1.994971])
# nocc_27 = np.array([1.940128, 1.966361, 1.977811, 1.98679, 1.989957, 1.993155, 1.995021])

# ccsd_act_25 = np.array([16.70, 16.48, 16.12, 15.76, 15.29, 14.78])
# ccsd_act_26 = np.array([16.85, 16.59, 16.22, 15.90, 15.58, 15.13])
# ccsd_act_27 = np.array([16.85, 16.55, 16.26, 15.92, 15.55, 15.11])

# ccsdt_act_25 = np.array([16.67, 16.43, 16.01, 15.59, 15.00, 14.28])
# ccsdt_act_26 = np.array([16.87, 16.51, 16.13, 15.74, 15.35, 14.69])
# ccsdt_act_27 = np.array([16.85, 16.47, 16.16, 15.77, 15.35, 14.73])

# ccsd_act_25 = ((1 + 4 * ccsd_act_25) ** 0.5 - 1) / 2
# ccsd_act_26 = ((1 + 4 * ccsd_act_26) ** 0.5 - 1) / 2
# ccsd_act_27 = ((1 + 4 * ccsd_act_27) ** 0.5 - 1) / 2

# ccsdt_act_25 = ((1 + 4 * ccsdt_act_25) ** 0.5 - 1) / 2
# ccsdt_act_26 = ((1 + 4 * ccsdt_act_26) ** 0.5 - 1) / 2
# ccsdt_act_27 = ((1 + 4 * ccsdt_act_27) ** 0.5 - 1) / 2

plt.plot(noccs[2], ccsd_act[2], '--', marker='D', color=colors[10], linewidth=2, markersize=4.5, mfc='white', mew=2, label='%s (UCCSD)' % bs_map[str(27)], alpha=0.7)
plt.plot(noccs[1], ccsd_act[1], '--', marker='p', color=colors[6], linewidth=2, markersize=6, mfc='white', mew=2, label='%s (UCCSD)' % bs_map[str(26)], alpha=0.7)
plt.plot(noccs[0], ccsd_act[0], '--', marker='H', color=colors[4], linewidth=2, markersize=6, mfc='white', mew=2, label='%s (UCCSD)' % bs_map[str(25)], alpha=0.7)

plt.plot(noccs[2], ccsdt_act[2], marker='D', color=colors[10], linewidth=2, markersize=4.5, mfc='white', mew=2, label='%s (UCCSDT)' % bs_map[str(27)])
plt.plot(noccs[1], ccsdt_act[1], marker='p', color=colors[6], linewidth=2, markersize=6, mfc='white', mew=2, label='%s (UCCSDT)' % bs_map[str(26)])
plt.plot(noccs[0], ccsdt_act[0], marker='H', color=colors[4], linewidth=2, markersize=6, mfc='white', mew=2, label='%s (UCCSDT)' % bs_map[str(25)])


# for i in range(7):
#     plt.text(nocc_25[i], (ccsdt_25 + ccsdtq_act_25 - ccsdt_act_25)[i] - 0.0003, '%do' % (30 + i * 4),
#         ha=['center', 'left'][i != 0], va=['top', 'center'][i != 0], color=colors[0], fontsize=12)

plt.xlim((2.000, 1.928))
# plt.ylim((1.5, 3.65))
# plt.ylabel('E[CCSDTQ-act] - E[CCSDT-act] + E[CCSDT]')
plt.ylabel('Total spin expectation $\\langle S \\rangle$')
plt.xlabel('NO occupation cutoff')
plt.legend(loc='lower right', frameon=False, fontsize=14, ncols=1)
plt.subplots_adjust(left=0.13, bottom=0.1, right=0.97, top=0.95)
plt.savefig('09-sm-fm-fno-ssq.png', dpi=300)
