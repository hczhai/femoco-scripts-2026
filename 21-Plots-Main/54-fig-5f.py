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
plt.figure(figsize=(7, 5))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/14-large-uccsd.json", "r", encoding='utf-8'))

bs_init = data['bs']
bs_map = data['name_ijk']
colors = data['colors']

ia = 3
prop_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x.get('uccsd-nuc-density', [0.0] * 8))
    for x in data['large_uccsd'][str(ia)][str(ix)] if 'uccsd' in x] + [(0, 0, [0.0] * 8)])[2] for ix in range(35)])
prop_order = [np.argsort(x) for x in prop_uccsd]

inits = [18, 25, 26, 27, 8]
femap = [0, 1, 4, 3, 2, 5, 7, 6]

for ix, init in enumerate(inits):
    for jx in range(7):
        plt.barh(ix, 1, left=6 - jx + 0.5, color=colors[prop_order[init][jx]], height=0.5, alpha=0.25)

for ix, init in enumerate(inits):
    for jx in range(7):
        plt.text(6 - jx + 1, ix, 'Fe%d' % femap[prop_order[init][jx] + 1], color=colors[prop_order[init][jx]], fontsize=18,
            ha='center', va='center', zorder=10)

plt.yticks([0, 1, 2, 3, 4], [bs_map[str(ix)] for ix in inits])
plt.xlim(0, 8)
plt.xticks(list(range(1, 8)))

plt.gca().invert_yaxis()
plt.xlabel('Fe center reverse ordering \nby electron density at the nucleus')
plt.tight_layout()
plt.savefig('54-fig-5f.png', dpi=300)
