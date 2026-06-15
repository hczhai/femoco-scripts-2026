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
f, axs = plt.subplots(2, 8, figsize=(10, 4))

data = {
    "colors": [
        "#2F2D2E",
        "#BB4430",
        "#F18F01",
        "#99C24D",
        "#1E77CA",
        "#23C9FF",
        "#7137E5",
        "#A347BA",
        "#84602E",
        "#430B9D",
        "#18AD8A"
    ],
    "grays": [
        "#000000",
        "#3D3D3D",
        "#767676",
        "#C5C5C5",
    ],
    "blues": [
        "#1E77CA",
        "#A9D2F8",
    ],
    "lights": [
        "#FCDAAD",
        "#E3FCAC",
        "#A9D2F8",
        "#CFACFD",
    ]
}

colors = data['colors']
grays = data['grays']
lights = data['lights']

cmap = LinearSegmentedColormap.from_list("my_cmap", [colors[4], 'white', colors[3]])

fsexpt = lambda x: ((1 + 4 * x) ** 0.5 - 1) / 2

scorrs = np.load('24-spin-corr.npy')
ener = np.load('24-energies.npy') * 1000
sexpt = fsexpt(np.load('24-ssqs.npy'))
ener -= np.min(ener)

idx = np.array([1, 4, 3, 2, 5, 7, 6, 8]) - 1

for ix in range(2):
    for iy in range(8):
        ax = axs[ix, iy]
        scorrs[ix, iy] = scorrs[ix, iy][idx, :][:, idx]
        xx = scorrs[ix, iy]
        xx[xx > 0] = 5
        xx[xx < 0] = -5
        ax.matshow(scorrs[ix, iy], cmap=cmap, vmin=-5, vmax=5)
        axs[ix, iy].text(-0.5, -2.5, "SZ%d/2-R%d\n$E=%.3f$\n$\\langle S\\rangle=%.1f$" %
            (ix * 2 + 1, iy, ener[ix, iy], sexpt[ix, iy]), fontsize=10, clip_on=False, ha='left', va='center')
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])

plt.tight_layout()
plt.savefig('27-plot-24-disc.png', dpi=300)
