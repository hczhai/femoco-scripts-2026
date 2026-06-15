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
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({'font.size': 16})
f, axs = plt.subplots(2, 3, figsize=(6, 4))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/22-fe4s4-udmrg.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/18-fe4s4-sa-dmrg.json", "r", encoding='utf-8'))

colors = data['colors']
grays = data['grays']
lights = data['lights']

scorr = [
    [data['fe4s4-sadmrg-20k'][str(ix)]['spin-corr'][0] for ix in range(3)],
    [data['fe4s4-udmrg'][str(ix)]['spin-corr']['12000'] for ix in range(3)]
]

cmap = LinearSegmentedColormap.from_list("my_cmap", [colors[4], 'white', colors[3]])

vmin, vmax = np.min(scorr), np.max(scorr)

for ix in range(2):
    for iy in range(3):
        ax = axs[ix, iy]
        ax.matshow(scorr[ix][iy], cmap=cmap, vmin=vmin, vmax=vmax)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])

for ix in range(0, 3):
    axs[0, ix].text(1.5, -1.2, "SA-DMRG (20k)\nBS%d" % (ix + 1), fontsize=10, clip_on=False, ha='center', va='center')
    axs[1, ix].text(1.5, -1.2, "UDMRG (12k)\nBS%d" % (ix + 1), fontsize=10, clip_on=False, ha='center', va='center')

plt.tight_layout()
plt.savefig('20-sm-fe4-spin-corr.png', dpi=300)
