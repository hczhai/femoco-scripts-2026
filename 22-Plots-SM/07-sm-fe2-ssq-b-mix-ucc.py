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
orders = np.array([0, 2, 2.5, 3, 4, 5])
xs = ['uhf', 'uccsd', 'uccsd(t)', 'uccsdt', 'uccsdtq', 'uccsdtqp']

cc_ssq = np.array([data['fe2s2']['mix']["%s_ssq" % x][0] for x in xs])
cc_twos = np.array([data['fe2s2']['mix']["%s_ssq" % x][1] - 1 for x in xs])
cc_twosz1 = np.array([data['fe2s2']['mix']["%s_sz" % x][0] for x in xs])
cc_twosz2 = np.array([data['fe2s2']['mix']["%s_sz" % x][1] for x in xs])

plt.plot(orders, cc_ssq, '-', linewidth=2, marker='o', mec=colors[0],  markersize=7, color=colors[0], mfc='white', mew=2, label='$\\langle S^2 \\rangle$', clip_on=False)
plt.plot(orders, cc_twos, '-', linewidth=2, marker='o', mec=colors[1],  markersize=7, color=colors[1], mfc='white', mew=2, label='$2\\langle S \\rangle$', clip_on=False)
plt.plot(orders, -cc_twosz1, '-', linewidth=2, marker='o', mec=colors[2], markersize=7, color=colors[2], mfc='white', mew=2, label='$-2S_z$ (Fe1)', clip_on=False)
plt.plot(orders, cc_twosz2, '--', linewidth=2, marker='o', mec=colors[4], markersize=7, color=colors[4], mfc='white', mew=2, label='$2S_z$ (Fe2)', clip_on=False)

plt.xlabel('UHF/CC excitation order')
plt.ylim((0, 5))
plt.xlim((-0.5, 5.5))
plt.legend(loc='upper right', frameon=False, fontsize=14, ncols=2)
plt.subplots_adjust(left=0.05, right=0.96, top=0.95, bottom=0.12)
plt.savefig('07-sm-fe2-ssq-b-mix-ucc.png', dpi=300)
