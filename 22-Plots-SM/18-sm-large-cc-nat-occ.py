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
import numpy as np, sys, json

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(7.5, 6))

iiid = int(sys.argv[1])

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/14-large-uccsd.json", "r", encoding='utf-8'))

bs_map = data['name_x']
bs_init = data['bs']
colors = data['colors']
grays = data['grays']

ix = iiid
name = '%s (#%d)' % (bs_map[str(ix)], ix)

uccsd_no = np.array([min([(x['uccsd'], x.get('uccsd_nat-occ', [np.nan] * 76)) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])

n_sites = [117, 180, 285, 404]
n_froz = [54 - 39, 91 - 39, 121 - 39]

large_uccsd = [np.array([min([(x['uccsd'], x.get('uccsd_nat-occ', [np.nan] * n_sites[ia]))
    for x in data['large_uccsd'][str(ia)][str(jx)] if 'uccsd' in x] + [(0, [0.0] * n_sites[ia])])[1] for jx in range(35)])
    for ia in [1, 2, 3]]

plt.plot(np.arange(len(uccsd_no[ix])), uccsd_no[ix], '-', linewidth=1.25, marker='s',
    mec=colors[1], color=colors[1], mfc='white', mew=1.5, label='UCCSD (LLDUC)', clip_on=False)
plt.plot(np.arange(len(uccsd_no[ix])), large_uccsd[1][ix][n_froz[1]:n_froz[1] + 76], '-', linewidth=1.25, marker='^',
    mec=grays[0], color=grays[0], mfc='white', mew=1.5, label='UCCSD (285o)', clip_on=False)
plt.plot(np.arange(len(uccsd_no[ix])), large_uccsd[2][ix][n_froz[2]:n_froz[2] + 76], '-', linewidth=1.25, marker='h',
    mec=grays[2], color=grays[2], mfc='white', mew=1.5, label='UCCSD (404o)', clip_on=False)

plt.text(2, 1.75, name, fontsize=18)

plt.ylabel('Natural orbital occupancy')
plt.xlabel('Natural orbital index')
plt.ylim((0, 2))
plt.xlim((0, 76))
plt.legend(loc='upper right', frameon=False, fontsize=16, ncols=1)
plt.subplots_adjust(left=0.12, right=0.97, top=0.95, bottom=0.11)
plt.savefig('18-sm-large-cc-nat-occ-%02d.png' % iiid, dpi=300)
