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
plt.figure(figsize=(10, 6))

iiid = 26

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/09-llduc-udmrg-5k-ref.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/73-expr-jacs-2000-122-4926.json", "r", encoding='utf-8'))

bs_map = data['name_x']
colors = data['colors']

expr = np.array(data['expr-isomer-shift'])
expr2 = expr.copy()
expr2[np.sum(expr[:, None] == expr[None, :], axis=1) != 2] = np.nan

def align_distri(x, ref):
    return (x - x.mean()) / x.std() * ref.std() + ref.mean()

nuc_ccsd = np.array([x[1] for x in sorted([(x['uhf'], x['uccsd-nuc-density'][:7]) for x in data['uccsd'][str(iiid)]])])

udmrg_map = {x['init']: x['udmrg-nuc-density'][:7] for x in data['udmrg_5k_ref'][str(iiid)] if 'udmrg-nuc-density' in x}
nuc_udmrg = np.array([udmrg_map.get(x[1], [np.nan] * 7) for x in sorted([(x['uhf'], x['init']) for x in data['uccsd'][str(iiid)]])])

nuc_ccsd = align_distri(-nuc_ccsd, expr)
nuc_udmrg[~np.isnan(nuc_udmrg)] = align_distri(-nuc_udmrg[~np.isnan(nuc_udmrg)], expr)

idxs_a, idxs_b = [0, 1, 2, 3, 4, 5], [7, 8, 9, 15, 22, 23]

plt.subplot(1, 2, 1)
for ix, i in enumerate(idxs_a):
    for p in range(7):
        plt.plot(nuc_ccsd[i, p], ix + 0.1, '', linewidth=0, marker='s', mec=colors[p], markersize=4, color=colors[p], mew=1.5, mfc='none', label='UCCSD (76o)' if ix == 0 and p == 0 else None)
        plt.plot(nuc_udmrg[i, p], ix - 0.1, '', linewidth=0, marker='h', mec=colors[p], markersize=4, color=colors[p], mew=1.5, mfc='none', label='UDMRG (D=5k,76o)' if ix == 0 and p == 0 else None)
plt.plot(expr, [6 + 0.1] * 7, '', linewidth=0, marker='o', mec=colors[0], markersize=4, color=colors[0], mew=1.5, mfc='none')
plt.plot(expr2, [6 - 0.1] * 7, '', linewidth=0, marker='o', mec=colors[0], markersize=4, color=colors[0], mew=1.5, mfc='none')
plt.yticks([0, 1, 2, 3, 4, 5, 6], ['%s-%d' % (bs_map[str(iiid)], x) for x in idxs_a] + ['EXP'])
plt.xlim(0.30, 0.51)
plt.legend(bbox_to_anchor=(0.5, 1.08), loc='center', fontsize=12, ncols=2)

plt.subplot(1, 2, 2)
for ix, i in enumerate(idxs_b):
    for p in range(7):
        plt.plot(nuc_ccsd[i, p], ix + 0.1, '', linewidth=0, marker='s', mec=colors[p], markersize=4, color=colors[p], mew=1.5, mfc='none', label='UCCSD (76o)' if ix == 0 and p == 0 else None)
        plt.plot(nuc_udmrg[i, p], ix - 0.1, '', linewidth=0, marker='h', mec=colors[p], markersize=4, color=colors[p], mew=1.5, mfc='none', label='UDMRG (D=5k,76o)' if ix == 0 and p == 0 else None)
plt.plot(expr, [6 + 0.1] * 7, '', linewidth=0, marker='o', mec=colors[0], markersize=4, color=colors[0], mew=1.5, mfc='none')
plt.plot(expr2, [6 - 0.1] * 7, '', linewidth=0, marker='o', mec=colors[0], markersize=4, color=colors[0], mew=1.5, mfc='none')
plt.yticks([0, 1, 2, 3, 4, 5, 6], ['%s-%d' % (bs_map[str(iiid)], x) for x in idxs_b] + ['EXP'])
plt.xlim(0.30, 0.51)
plt.xlabel('Rescaled isomer shift $\\alpha[\\rho(\\text{Fe}) - C] + \\beta$' + ' ' * 80)
plt.legend(bbox_to_anchor=(0.5, 1.08), loc='center', fontsize=12, ncols=2)

plt.tight_layout()
plt.savefig('05-sm-ref-dep-prop-a-nuc.png', dpi=300)
