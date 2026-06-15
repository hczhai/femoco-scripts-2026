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
plt.figure(figsize=(9.5, 6))
plt.grid(which='major', axis='y', linestyle='--')

ha_to_kcal_mol = 627.5094740631
kj_mol_to_kcal_mol = 0.239005736

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/00-llduc-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/72-dft-10.1002-qua.25627.json", "r", encoding='utf-8'))

bs_init = data['bs']
colors = data['colors']

e_uhf = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['uhf'][str(ix)]])[0] for ix in range(35)])
e_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[0] for ix in range(35)])
e_uccsdt = np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_b3lyp = np.array(data['B3LYP/def2-TZVPD']) * kj_mol_to_kcal_mol / ha_to_kcal_mol
e_tpss = np.array(data['TPSS-D3/def2-TZVPD']) * kj_mol_to_kcal_mol / ha_to_kcal_mol

e_uhf = (e_uhf - np.min(e_uhf)) * ha_to_kcal_mol
e_uccsd = (e_uccsd - np.min(e_uccsd)) * ha_to_kcal_mol
e_uccsdt = (e_uccsdt - np.min(e_uccsdt)) * ha_to_kcal_mol
e_b3lyp = (e_b3lyp - np.min(e_b3lyp)) * ha_to_kcal_mol
e_tpss = (e_tpss - np.min(e_tpss)) * ha_to_kcal_mol

ik = 0
for k in bs_init:
    plt.plot([ik - 0.20] * len(k), e_uhf[np.array(k)], '-', linewidth=2, marker='.', mec=colors[0], markersize=6, color=colors[0], mfc=colors[0], clip_on=False)
    plt.plot([ik - 0.10] * len(k), e_uccsd[np.array(k)], '-', linewidth=2, marker='.', mec=colors[1], markersize=6, color=colors[1], mfc=colors[1], clip_on=False)
    plt.plot([ik + 0.00] * len(k), e_uccsdt[np.array(k)], '-', linewidth=2, marker='.', mec=colors[3], markersize=6, color=colors[3], mfc=colors[3], clip_on=False)
    plt.plot([ik + 0.10] * len(k), e_b3lyp[np.array(k)], '-', linewidth=2, marker='.', mec=colors[7], markersize=6, color=colors[7], mfc=colors[7], clip_on=False)
    plt.plot([ik + 0.20] * len(k), e_tpss[np.array(k)], '-', linewidth=2, marker='.', mec=colors[8], markersize=6, color=colors[8], mfc=colors[8], clip_on=False)
    if k != bs_init[-1]:
        plt.axvline(x=ik + 0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
    ik += 1

ik = 0
for k in bs_init:
    plt.plot([ik - 0.20], np.min(e_uhf[np.array(k)]), '-', linewidth=2, marker='o', mec=colors[0], markersize=6, color=colors[0], mfc='white', mew=1.5, label='UHF (LLDUC)' if k == [34] else None, clip_on=False)
    plt.plot([ik - 0.10], np.min(e_uccsd[np.array(k)]), '-', linewidth=2, marker='^', mec=colors[1], markersize=6, color=colors[1], mfc='white', mew=1.5, label='UCCSD (LLDUC)' if k == [34] else None, clip_on=False)
    plt.plot([ik + 0.00], np.min(e_uccsdt[np.array(k)]), '-', linewidth=2, marker='p', mec=colors[3], markersize=6, color=colors[3], mfc='white', mew=1.5, label='UCCSDT (LLDUC)' if k == [34] else None, clip_on=False)
    plt.plot([ik + 0.10], np.min(e_b3lyp[np.array(k)]), '-', linewidth=2, marker='X', mec=colors[7], markersize=6, color=colors[7], mfc='white', mew=1.5, label='DFT/B3LYP (Ref. 22)' if k == [34] else None, clip_on=False)
    plt.plot([ik + 0.20], np.min(e_tpss[np.array(k)]), '-', linewidth=2, marker='D', mec=colors[8], markersize=5, color=colors[8], mfc='white', mew=1.5, label='DFT/TPSS-D3 (Ref. 22)' if k == [34] else None, clip_on=False)
    ik += 1

plt.ylabel('Relative energy (kcal/mol)')
plt.ylim((0, 40))
plt.xticks(list(range(len(bs_init))), ['BS%d' % x for x in range(1, 11)])
plt.legend(loc='upper right', frameon=False, fontsize=16, ncols=1)
plt.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.075)
plt.savefig('22-fig-2c.png', dpi=300)
