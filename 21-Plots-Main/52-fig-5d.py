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
plt.figure(figsize=(9, 6))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/14-large-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/25-fe-oxidation.json", "r", encoding='utf-8'))

bs_init = data['bs']
colors = data['colors']
grays = data['grays']

ia = 2
prop_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x.get('uccsd-elec-pop-noryd', [np.nan] * 8)[:7])
    for x in data['large_uccsd'][str(ia)][str(ix)] if 'uccsd' in x]
    + [(0, 0, [0] * 7)])[2] for ix in range(35)])

ik = 0
tks = []
for k in bs_init:
    for ix in range(7):
        plt.plot(np.arange(ik, ik + len(k)), prop_uccsd[np.array(k), ix], ' ', linewidth=2, marker='o', mec=colors[ix],
            markersize=6, color=colors[ix], mfc='white', mew=2)
    tks.append((ik + ik + len(k) - 1) / 2)
    plt.axvline(x=ik + len(k) - 0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
    ik += len(k)
plt.axhline(y=data['fe-oxidation']['fe2+'], color=grays[0], linewidth=2, linestyle='--', label='$\\mathrm{Fe^{2+}}$')
plt.axhline(y=data['fe-oxidation']['fe2.5+'], color=grays[0], linewidth=2, linestyle='--', label='$\\mathrm{Fe^{2.5+}}$', alpha=0.5)
plt.axhline(y=data['fe-oxidation']['fe3+'], color=grays[0], linewidth=2, linestyle='--', label='$\\mathrm{Fe^{3+}}$', alpha=0.25)

plt.ylabel('Number of electrons on Fe')
plt.xticks(tks, [('BS%d' % x) if x > 2 else [None, 'BS1  ', ' 2'][x] for x in range(1, 11)])
plt.xlim((-0.5, 34.5))
plt.legend(loc='upper right', frameon=False, fontsize=16, ncols=3)
plt.tight_layout()
plt.savefig('52-fig-5d.png', dpi=300)
