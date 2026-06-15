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
import scipy.stats

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(8, 6))

ha_to_kcal_mol = 627.5094740631

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/03-llduc-udmrg-5k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/05-llduc-udmrg-8k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/26-llduc-composite.json", "r", encoding='utf-8'))

data |= json.load(open("../20-Data/10-fe2s2.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))

bs_map = data['name_ijk']
colors = data['colors']

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

fe2_mix_corr = [data['fe2s2']['mix'][x] - data['fe2s2']['mix'][y]
    for x, y in zip(['uccsd', 'uccsdt', 'uccsdtq', 'sa_dmrg_12k'],
        ['uhf', 'uccsd', 'uccsdt', 'uccsdtq'])]

femoco_corr = [[
    float(e_uccsd[ix] - e_uhf[ix]),
    float(e_uccsdt[ix] - e_uccsd[ix]),
    float(intercept_cc[ii] - e_uccsdt[ix])] for ii, ix in enumerate([25, 26, 27])]

femoco_ratio = [[x / y for x, y in zip(xs, fe2_mix_corr)] for xs in femoco_corr]

plt.plot([1, 2, 3], femoco_ratio[0], '-', linewidth=2, marker='H', mec=colors[4], markersize=7, color=colors[4], mfc='white', mew=2, label=bs_map['25'], clip_on=False)
plt.plot([1, 2, 3], femoco_ratio[1], '-', linewidth=2, marker='p', mec=colors[6], markersize=7, color=colors[6], mfc='white', mew=2, label=bs_map['26'], clip_on=False)
plt.plot([4], [np.average(femoco_ratio[0])], '-', linewidth=0, marker='X', mec=colors[4], markersize=8, color=colors[4], mfc='white', mew=1.7, label=bs_map['25'] + ' (estimated)', clip_on=False)
plt.plot([4], [np.average(femoco_ratio[1])], '-', linewidth=0, marker='X', mec=colors[6], markersize=8, color=colors[6], mfc='white', mew=1.7, label=bs_map['26'] + ' (estimated)', clip_on=False)

# ik = 0
# tks = []
# for k in inits:
#     plt.plot(np.arange(ik, ik + 1), [(e_uccsd[k] - e_ref) * ha_to_kcal_mol], '-', linewidth=2, marker='^', mec=colors[1], markersize=6, color=colors[1], mfc=colors[1], mew=2, label='UCCSD' if k == 25 else None, clip_on=False)
#     plt.plot(np.arange(ik, ik + 1), [(e_uccsd_t[k] - e_ref) * ha_to_kcal_mol], '-', linewidth=2, marker='s', mec=colors[2], markersize=5, color=colors[2], mfc=colors[2], mew=2, label='UCCSD(T)' if k == 25 else None, clip_on=False)
#     tks.append(ik)
#     ik += 2.4

# plt.axvline(x=ik - 0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)

# xr = list(np.arange(ik, ik + 9)) + list(np.arange(ik + 9, ik + 9 + 7 * 1.5, 1.5))

# plt.plot(xr, (e_gccsd - e_ref) * ha_to_kcal_mol, '-', linewidth=2, marker='^', mec=colors[1], markersize=6, color=colors[1], mfc='white', mew=2, label='GCCSD', clip_on=False)
# plt.plot(xr, (e_gccsd_t - e_ref) * ha_to_kcal_mol, '-', linewidth=2, marker='s', mec=colors[2], markersize=5, color=colors[2], mfc='white', mew=2, label='GCCSD(T)', clip_on=False)

# tks += xr

plt.ylim((0, 8.9))
plt.xlim((0.3, 4.7))

plt.xticks([1, 2, 3, 4], ['D', 'T', 'Q', 'post Q'])
# plt.yticks(np.arange(7, 9.1, 0.5), ["%.1f" % x for x in np.arange(7, 9.1, 0.5)])

plt.text(1.034, 4.25, '$E^{\\text{estimated}}_{\\text{UCC}} = E^{\\text{estimated}}_{\\text{UCCSDTQ}}+\\Delta E_{\\text{post-Q}}^{\\text{Fe(III)-Fe(II)}}' +
    '\\cdot \\left[ \\frac{\\Delta E_{\\text{post-Q}}^{\\text{FeMo-co}}}{\\Delta E_{\\text{post-Q}}^{\\text{Fe(III)-Fe(II)}}} \\right]$',
    ha='left', va='center', color=colors[0], fontsize=15)
plt.text(3.98, 4.472, '${}^{\\text{estimated}}$',
    ha='left', va='center', color=colors[0], fontsize=15)

for ii in range(1, 4):
    plt.text(ii, femoco_ratio[0][ii - 1] + [0.144, -0.144][ii == 1], '%.2f' % femoco_ratio[0][ii - 1],
        ha='center', va=['bottom', 'top'][ii == 1], color=colors[0], fontsize=14, alpha=0.7)
    plt.text(ii, femoco_ratio[1][ii - 1] + [0.144, -0.144][ii != 1], '%.2f' % femoco_ratio[1][ii - 1],
        ha='center', va=['bottom', 'top'][ii != 1], color=colors[0], fontsize=14, alpha=0.7)

plt.text(4.08, np.average(femoco_ratio[0]) - 0.12, '%.2f' % np.average(femoco_ratio[0]),
    ha='left', va='center', color=colors[0], fontsize=14, alpha=0.7)
plt.text(4.08, np.average(femoco_ratio[1]) + 0.12, '%.2f' % np.average(femoco_ratio[1]),
    ha='left', va='center', color=colors[0], fontsize=14, alpha=0.7)
    
plt.ylabel('Correlation energy increment ratio')
plt.xlabel('UCC excitation level')
# plt.xticks(tks, [bs_map[str(ix)] for ix in inits] + ['%d' % x for x in ghf_ids])
plt.legend(loc='lower right', frameon=False, fontsize=16, ncols=2)
plt.subplots_adjust(left=0.11, right=0.97, top=0.95, bottom=0.1)
plt.savefig('32-fig-3c.png', dpi=300)
