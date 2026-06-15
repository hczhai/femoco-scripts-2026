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
plt.figure(figsize=(8, 6))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/10-fe2s2.json", "r", encoding='utf-8'))

colors = data['colors']
ds = [100, 200, 400, 600, 800, 1000, 1200, 1600, 2000, 2400]

dmrg_ssq = np.array([data['fe2s2']['std']['udmrg']['reverse_ssq'][str(d)][0] for d in ds])
dmrg_twos = np.array([data['fe2s2']['std']['udmrg']['reverse_ssq'][str(d)][1] - 1 for d in ds])
dmrg_twosz1 = np.array([data['fe2s2']['std']['udmrg']['reverse_sz'][str(d)][0] for d in ds])
dmrg_twosz2 = np.array([data['fe2s2']['std']['udmrg']['reverse_sz'][str(d)][1] for d in ds])

plt.plot(ds, dmrg_ssq, '-', linewidth=2, marker='p', mec=colors[0],  markersize=7, color=colors[0], mfc='white', mew=2, label='$\\langle S^2 \\rangle$', clip_on=False)
plt.plot(ds, dmrg_twos, '-', linewidth=2, marker='p', mec=colors[1],  markersize=7, color=colors[1], mfc='white', mew=2, label='$2\\langle S \\rangle$', clip_on=False)
plt.plot(ds, -dmrg_twosz1, '-', linewidth=2, marker='p', mec=colors[2], markersize=7, color=colors[2], mfc='white', mew=2, label='$-2S_z$ (Fe1)', clip_on=False)
plt.plot(ds, dmrg_twosz2, '--', linewidth=2, marker='p', mec=colors[4], markersize=7, color=colors[4], mfc='white', mew=2, label='$2S_z$ (Fe2)', clip_on=False)

plt.xlabel('MPS bond dimension $D$')
plt.ylim((0, 5))
plt.legend(loc='upper right', frameon=False, fontsize=14, ncols=2)
plt.subplots_adjust(left=0.05, right=0.93, top=0.95, bottom=0.12)
plt.savefig('07-sm-fe2-ssq-c-std-udmrg.png', dpi=300)
