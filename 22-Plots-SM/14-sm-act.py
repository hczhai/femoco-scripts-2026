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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams.update({'font.size': 12})
fig, ax = plt.subplots(figsize=(12, 6))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/12-large-act.json", "r", encoding='utf-8'))

colors = data['colors']
grays = data['grays']

occ = np.array(data['large_act']['4']['uno_occ'])
ovlp = np.array(data['large_act']['4']['uno_llduc_overlap'])

x_single = np.mgrid[:len(occ)][(occ >= 0.1) & (occ <= 1.8)]
y_single = occ[(occ >= 0.1) & (occ <= 1.8)]

x_dbl0 = np.mgrid[:len(occ)][(occ <= 1.9700)]
x_dbl1 = np.mgrid[:len(occ)][(occ <= 1.9960)]
x_dbl2 = np.mgrid[:len(occ)][(occ <= 1.9840)]
x_dbl3 = np.mgrid[:len(occ)][(occ <= 1.9985)]
y_dbl0 = occ[(occ <= 1.970)]
y_dbl1 = occ[(occ <= 1.996)]

x_dbl = np.mgrid[:len(occ)][(occ >= 1.8)]
y_dbl = occ[(occ >= 1.8)]

x_vir = np.mgrid[:len(occ)][(occ <= 0.1)]
y_vir = occ[(occ <= 0.1)]

x_vir0 = np.mgrid[:len(occ)][(occ <= 0.020)]
x_vir1 = np.mgrid[:len(occ)][(occ <= 0.010)]
x_vir2 = np.mgrid[:len(occ)][(occ <= 0.005)]
x_vir3 = np.mgrid[:len(occ)][(occ <= 0.0025)]

y_vir0 = occ[(occ <= 0.020)]
y_vir1 = occ[(occ <= 0.010)]
y_vir2 = occ[(occ <= 0.005)]
y_vir3 = occ[(occ <= 0.0025)]

axins = inset_axes(ax, width=2.0, height=1.2, loc='upper left', bbox_to_anchor=(350, 1.9), bbox_transform=ax.transData)
axins.plot(x_single, y_single, color=colors[0])
axins.set_xticks([x_single[0], x_single[-1]])
axins.set_yticks([y_single[0], 1.0, y_single[-1]])
axins.grid(True, alpha=0.3)
axins.text(194, 1.01, '[0.1 < occ < 1.8] = %d' % (x_single[-1] + 1 - x_single[0]), ha='left', va='center', fontsize=12, zorder=100)

axins = inset_axes(ax, width=3.5, height=1.2, loc='upper left', bbox_to_anchor=(650, 1.9), bbox_transform=ax.transData)
axins.plot(x_dbl, y_dbl, color=colors[0])
axins.set_xticks([x_dbl[0], x_dbl1[0], x_dbl2[0], x_dbl0[0], x_dbl3[0], x_dbl[-1]])
axins.set_yticks([y_dbl[0], 1.970, 1.984, y_dbl[-1]])
# axins.set_yscale('log')
axins.grid(True, alpha=0.3)
axins.text(00, 1.973, '[1.8 < occ < 1.9985] = %d' % (1 + x_dbl[-1] - x_dbl3[0]), ha='left', va='center', fontsize=12, zorder=100)
axins.text(00, 1.958, '[1.8 < occ < 1.9960] = %d' % (1 + x_dbl[-1] - x_dbl1[0]), ha='left', va='center', fontsize=12, zorder=100)
axins.text(00, 1.943, '[1.8 < occ < 1.9840] = %d' % (1 + x_dbl[-1] - x_dbl2[0]), ha='left', va='center', fontsize=12, zorder=100)
axins.text(00, 1.928, '[1.8 < occ < 1.9700] = %d' % (1 + x_dbl[-1] - x_dbl0[0]), ha='left', va='center', fontsize=12, zorder=100)

axins = inset_axes(ax, width=6, height=1.8, loc='lower left', bbox_to_anchor=(350, 0.25), bbox_transform=ax.transData)
axins.plot(x_vir, y_vir, color=colors[0])
axins.set_xticks([x_vir[0], x_vir0[0]-1, x_vir1[0]-1, x_vir2[0]-1, x_vir3[0]-1, x_vir[-1]])
axins.set_yticks([y_vir[0], 0.020, 0.010, 0.005, 0.0025, y_vir[-1]], ['%.4f' % y_vir[0], '0.0200', '0.0100', '', '0.0025', ''])
axins.grid(True, alpha=0.3)
axins.text(420, 0.042, '[0.0200 < occ < 0.1] = %d' % (x_vir0[0] - x_vir[0]), ha='left', va='center', fontsize=12, zorder=100)
axins.text(420, 0.035, '[0.0100 < occ < 0.1] = %d' % (x_vir1[0] - x_vir[0]), ha='left', va='center', fontsize=12, zorder=100)
axins.text(420, 0.028, '[0.0050 < occ < 0.1] = %d' % (x_vir2[0] - x_vir[0]), ha='left', va='center', fontsize=12, zorder=100)
axins.text(420, 0.021, '[0.0025 < occ < 0.1] = %d' % (x_vir3[0] - x_vir[0]), ha='left', va='center', fontsize=12, zorder=100)

ax.axvline(x_single[0] - 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)
ax.axvline(x_single[-1] + 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)
ax.axvline(x_dbl0[0] - 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)
ax.axvline(x_dbl3[0] - 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)
ax.axvline(x_dbl1[0] - 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)
ax.axvline(x_dbl2[0] - 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)

ax.axvline(x_vir0[0] - 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)
ax.axvline(x_vir1[0] - 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)
ax.axvline(x_vir2[0] - 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)
ax.axvline(x_vir3[0] - 0.5, linestyle='--', color=grays[0], alpha=0.5, linewidth=0.75)

ax.text(x_single[0], 2.00, '%s' % x_single[0], ha='center', va='bottom', fontsize=8)
ax.text(x_single[-1], 2.00, '%s' % x_single[-1], ha='center', va='bottom', fontsize=8)
ax.text(x_dbl0[0], 2.00, '%s' % x_dbl0[0], ha='center', va='bottom', fontsize=8)
ax.text(x_dbl3[0], 2.00, '%s' % x_dbl3[0], ha='center', va='bottom', fontsize=8)
ax.text(x_dbl1[0], 2.00, '%s' % x_dbl1[0], ha='center', va='bottom', fontsize=8)
ax.text(x_dbl2[0], 2.00, '%s' % x_dbl2[0], ha='center', va='bottom', fontsize=8)

ax.text(x_vir0[0] - 1, 2.00, '%s' % (x_vir0[0] - 1), ha='center', va='bottom', fontsize=8)
ax.text(x_vir1[0] - 1, 2.00, '%s' % (x_vir1[0] - 1), ha='center', va='bottom', fontsize=8)
ax.text(x_vir2[0] - 1, 2.00, '%s' % (x_vir2[0] - 1), ha='center', va='bottom', fontsize=8)
ax.text(x_vir3[0] - 1, 2.00, '%s' % (x_vir3[0] - 1), ha='center', va='bottom', fontsize=8)

ax.plot(np.arange(len(occ)), occ, color=colors[0], linewidth=2)
ax.plot(np.arange(len(ovlp)), ovlp, color=colors[4], linewidth=2)
ax.set_xlabel('UMP2 natural orbital index')
ax.set_ylabel('UMP2 natural orbital occ')
ax.set_xlim((0, 1000))
ax.set_ylim((0, 2))

plt.subplots_adjust(left=0.07, right=0.97, top=0.92, bottom=0.1)
plt.savefig('14-sm-act.png', dpi=300)