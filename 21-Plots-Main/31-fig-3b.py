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


import matplotlib.pyplot as plt
import numpy as np, json

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(8, 6))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))

bs_map = data['name_ijk']
colors = data['colors']

nocc_25 = np.array(data['fno_ucc']['25'][0]['occ'])
nocc_26 = np.array(data['fno_ucc']['26'][0]['occ'])

e_uccsdt = 22140 + np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])

ccsdt_act_25 = np.array([data['fno_ucc']['25'][0][str(ir)]['uccsdt'] for ir in range(1, 8)])
ccsdt_act_26 = np.array([data['fno_ucc']['26'][0][str(ir)]['uccsdt'] for ir in range(1, 8)])
ccsdtq_act_25 = np.array([data['fno_ucc']['25'][0][str(ir)]['uccsdtq'] for ir in range(1, 8)])
ccsdtq_act_26 = np.array([data['fno_ucc']['26'][0][str(ir)]['uccsdtq'] for ir in range(1, 8)])

import scipy.stats
reg_25 = scipy.stats.linregress(nocc_25[-3:], (e_uccsdt[25] + ccsdtq_act_25 - ccsdt_act_25)[-3:])
reg_26 = scipy.stats.linregress(nocc_26[-3:], (e_uccsdt[26] + ccsdtq_act_26 - ccsdt_act_26)[-3:])
reg_x_25 = np.linspace(nocc_25[-3], 2.00, 10)
reg_x_26 = np.linspace(nocc_26[-3], 2.00, 10)

print('%10s E[54o] = %25.12f' % (bs_map['25'], (e_uccsdt[25] + ccsdtq_act_25 - ccsdt_act_25)[-1]))
print('%10s E[54o] = %25.12f' % (bs_map['26'], (e_uccsdt[26] + ccsdtq_act_26 - ccsdt_act_26)[-1]))
print('%10s E[fit] = %25.12f' % (bs_map['25'], reg_25.intercept + reg_25.slope * 2.00))
print('%10s E[fit] = %25.12f' % (bs_map['26'], reg_26.intercept + reg_26.slope * 2.00))

plt.plot(reg_x_25, reg_25.intercept + reg_25.slope * reg_x_25, '--', linewidth=1.5, color=colors[4])
plt.plot(reg_x_26, reg_26.intercept + reg_26.slope * reg_x_26, '--', linewidth=1.5, color=colors[6])
plt.plot(nocc_25, e_uccsdt[25] + ccsdtq_act_25 - ccsdt_act_25, marker='H', color=colors[4], linewidth=2, markersize=6, mfc='white', mew=2,
    label=bs_map['25'] + " ($E^{\\text{estimated}}_{\\text{UCCSDTQ}} = %.6f$)" % (reg_25.intercept + reg_25.slope * 2.00))
plt.plot(nocc_26, e_uccsdt[26] + ccsdtq_act_26 - ccsdt_act_26, marker='p', color=colors[6], linewidth=2, markersize=6, mfc='white', mew=2,
    label=bs_map['26'] + " ($E^{\\text{estimated}}_{\\text{UCCSDTQ}} = %.6f$)" % (reg_26.intercept + reg_26.slope * 2.00))

for i in range(7):
    plt.text(nocc_25[i], (e_uccsdt[25] + ccsdtq_act_25 - ccsdt_act_25)[i] - 0.0004, '%do' % (30 + i * 4),
        ha=['center', 'left'][i != 0], va=['top', 'center'][i != 0], color=colors[0], fontsize=12)

plt.text(1.987, -0.4051, '$E^{\\text{estimated}}_{\\text{UCCSDTQ}} = E_{\\text{UCCSDT}}+E_{\\text{FNO-UCCSDTQ}}-E_{\\text{FNO-UCCSDT}}$',
    ha='left', va='center', color=colors[0], fontsize=15)

plt.xlim((2.000, 1.93))
plt.ylabel('Estimated CCSDTQ energy + 22140 (Hartree)')
plt.xlabel('NO occupation cutoff')
plt.legend(loc='lower right', frameon=False, fontsize=14, ncols=1)
plt.subplots_adjust(left=0.17, bottom=0.1, right=0.95, top=0.95)
plt.savefig('31-fig-3b.png', dpi=300)
