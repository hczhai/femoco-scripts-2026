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
plt.figure(figsize=(5, 6))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/14-large-uccsd.json", "r", encoding='utf-8'))

bs_init = data['bs']
bs_map = data['name_x']
colors = data['colors']

prop_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd-elec-pop-noryd'][:7])
    for x in data['uccsd'][str(ix)]])[2] for ix in range(35)] )

prop_uccsd_ls = np.array([[min([(x['uccsd'], x['uccsd(t)'], x.get('uccsd-elec-pop-noryd', [np.nan] * 8)[:7])
    for x in data['large_uccsd'][str(ia)][str(ix)] if 'uccsd' in x]
    + [(0, 0, [0] * 7)])[2] for ix in range(35)] for ia in range(4)])

ibs = [6, 7, 9][1]
inits = bs_init[ibs] if ibs != 6 else bs_init[6] + bs_init[8]

for ix, i in enumerate(inits):
    for p in range(7):
        plt.plot(prop_uccsd[i, p], ix + 0.24, '', linewidth=0, marker='o', mec=colors[p], markersize=4, color=colors[p], mfc='none', label='UCCSD (76o)' if ix == 0 and p == 0 else None)
        plt.plot(prop_uccsd_ls[0][i, p], ix + 0.12, '', linewidth=0, marker='s', mec=colors[p], markersize=4, color=colors[p], mfc='none', label='UCCSD (117o)' if ix == 0 and p == 0 else None)
        plt.plot(prop_uccsd_ls[1][i, p], ix - 0.00, '', linewidth=0, marker='^', mec=colors[p], markersize=4, color=colors[p], mfc='none', label='UCCSD (180o)' if ix == 0 and p == 0 else None)
        plt.plot(prop_uccsd_ls[2][i, p], ix - 0.12, '', linewidth=0, marker='D', mec=colors[p], markersize=4, color=colors[p], mfc='none', label='UCCSD (285o)' if ix == 0 and p == 0 else None)
        plt.plot(prop_uccsd_ls[3][i, p], ix - 0.24, '', linewidth=0, marker='p', mec=colors[p], markersize=5, color=colors[p], mfc='none', label='UCCSD (404o)' if ix == 0 and p == 0 else None)
plt.yticks([0, 1, 2, 3, 4, 5], [bs_map[str(x)] for x in inits])
plt.legend(bbox_to_anchor=(0.5, 1.13), loc='center', fontsize=12, ncols=2)
plt.xlabel('Nelec on Fe (no ryd)')
plt.tight_layout()
plt.savefig('16-sm-large-pop-e-elec-bs8.png', dpi=300)
