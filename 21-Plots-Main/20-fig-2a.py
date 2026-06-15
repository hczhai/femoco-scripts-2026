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

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/00-llduc-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/03-llduc-udmrg-5k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/26-llduc-composite.json", "r", encoding='utf-8'))

bs_init = data['bs']
colors = data['colors']

e_ref = min([v['energy'] for _, v in data['composite'].items() if 'energy' in v])

e_uhf = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['uhf'][str(ix)]])[0] for ix in range(35)])
e_uhf_conv = [min([(x['uhf'], x['uhf_conv']) for x in data['uhf'][str(ix)]])[1] for ix in range(35)]
assert all(e_uhf_conv)

e_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[0] for ix in range(35)])
e_uccsd_t = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])
e_uccsd_conv = [min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[2] for ix in range(35)]
assert all(e_uccsd_conv)

e_uccsdt = np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_uccsdt_conv = [min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[1] for ix in range(35)]
# assert all(e_uccsdt_conv)

e_dmrg_5k = np.array([data['udmrg_5k'][str(ix)][0]['udmrg'] for ix in range(35)])

ik = 0
tks = []
for k in bs_init:
    plt.plot(np.arange(ik, ik + len(k)), (e_uhf[np.array(k)] - e_ref) * ha_to_kcal_mol, '-', linewidth=2, marker='o', mec=colors[0], markersize=6, color=colors[0], mfc='white', mew=1.5, label='UHF' if k == [34] else None)
    plt.plot(np.arange(ik, ik + len(k)), (e_uccsd[np.array(k)] - e_ref) * ha_to_kcal_mol, '-', linewidth=2, marker='^', mec=colors[1], markersize=6, color=colors[1], mfc='white', mew=1.5, label='UCCSD' if k == [34] else None)
    plt.plot(np.arange(ik, ik + len(k)), (e_uccsd_t[np.array(k)] - e_ref) * ha_to_kcal_mol, '-', linewidth=2, marker='s', mec=colors[2], markersize=5, color=colors[2], mfc='white', mew=1.5, label='UCCSD(T)' if k == [34] else None)
    plt.plot(np.arange(ik, ik + len(k)), (e_uccsdt[np.array(k)] - e_ref) * ha_to_kcal_mol, '-', linewidth=2, marker='p', mec=colors[3], markersize=6, color=colors[3], mfc='white', mew=1.5, label='UCCSDT' if k == [34] else None)
    plt.plot(np.arange(ik, ik + len(k)), (e_dmrg_5k[np.array(k)] - e_ref) * ha_to_kcal_mol, '-', linewidth=2, marker='h', mec=colors[4], markersize=6, color=colors[4], mfc='white', mew=1.5, label='UDMRG ($D$=5k)' if k == [34] else None)
    tks.append((ik + ik + len(k) - 1) / 2)
    plt.axvline(x=ik + len(k) - 0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
    ik += len(k)

plt.ylabel('Estimated total energy error (kcal/mol)')
plt.ylim((0, 300))
plt.xticks(tks, [('BS%d' % x) if x > 2 else [None, 'BS1  ', ' 2'][x] for x in range(1, 11)])
plt.xlim((-0.5, 34.5))
plt.legend(loc='center left', frameon=False, fontsize=16, ncols=3, bbox_to_anchor=(0.0, 0.59))
plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.075)
plt.savefig('20-fig-2a.png', dpi=300)
