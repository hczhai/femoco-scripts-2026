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
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/03-llduc-udmrg-5k.json", "r", encoding='utf-8'))

bs_map = data['name_x']
bs_init = data['bs']
colors = data['colors']

ix = iiid
name = '%s (#%d)' % (bs_map[str(ix)], ix)

uccsd_no = np.array([min([(x['uccsd'], x.get('uccsd_nat-occ', [np.nan] * 76)) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])
uccsdt_no = np.array([min([(x['uccsdt'], x.get('uccsdt_nat-occ', [np.nan] * 76)) for x in data['uccsdt'][str(ix)]])[1] for ix in range(35)])
udmrg_5k_no = np.array([data['udmrg_5k'][str(ix)][0].get('nat-occ', [np.nan] * 76) for ix in range(35)])

plt.plot(np.arange(len(uccsd_no[ix])), uccsd_no[ix], '-', linewidth=1.25, marker='s',
    mec=colors[1], color=colors[1], mfc='white', mew=1.5, label='UCCSD', clip_on=False)
plt.plot(np.arange(len(uccsdt_no[ix])), uccsdt_no[ix], '-', linewidth=1.25, marker='p',
    mec=colors[3], color=colors[3], mfc='white', mew=1.5, label='UCCSDT', clip_on=False)
plt.plot(np.arange(len(udmrg_5k_no[ix])), udmrg_5k_no[ix], '-', linewidth=1.25, marker='h',
    mec=colors[4], color=colors[4], mfc='white', mew=1.5, label='UDMRG ($D$=5k)', clip_on=False)

plt.text(2, 1.85, name, fontsize=18)

plt.ylabel('Natural orbital occupancy')
plt.xlabel('Natural orbital index')
plt.ylim((0, 2))
plt.xlim((0, 76))
plt.legend(loc='upper right', frameon=False, fontsize=16, ncols=1)
plt.subplots_adjust(left=0.12, right=0.97, top=0.95, bottom=0.11)
plt.savefig('13-sm-dmrg-cc-nat-occ-%02d.png' % iiid, dpi=300)
