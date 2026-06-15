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
plt.figure(figsize=(8, 6))

iiid = int(sys.argv[1])

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/00-llduc-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/03-llduc-udmrg-5k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/09-llduc-udmrg-5k-ref.json", "r", encoding='utf-8'))

bs_map = data['name_x']
bs_init = data['bs']
colors = data['colors']

ix = iiid
name = '%s (#%d)' % (bs_map[str(ix)], ix)

e_uhf = 22140 + np.array([x[0] for x in sorted([(x['uhf'], ) for x in data['uccsd'][str(ix)]])])
e_uccsd = 22140 + np.array([x[1] for x in sorted([(x['uhf'], x['uccsd']) for x in data['uccsd'][str(ix)]])])
e_uccsd_t = 22140 + np.array([x[1] for x in sorted([(x['uhf'], x['uccsd(t)']) for x in data['uccsd'][str(ix)]])])
e_uccsd_conv = np.array([x[1] for x in sorted([(x['uhf'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])])
t1diag_uccsd = np.array([x[1] for x in sorted([(x['uhf'], x['t1-diagnostic']) for x in data['uccsd'][str(ix)]])])

uccsdt_map = {x['init']: (x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]}
e_uccsdt = 22140 + np.array([uccsdt_map.get(x[1], [np.nan])[0] for x in sorted([(x['uhf'], x['init']) for x in data['uccsd'][str(ix)]])])
e_uccsdt_conv = np.array([uccsdt_map.get(x[1], [0, np.nan])[1] for x in sorted([(x['uhf'], x['init']) for x in data['uccsd'][str(ix)]])])

print('Min UHF init:', name + '-', min([(x, ix) for ix, x in enumerate(e_uhf) if not np.isnan(x)])[1])
print('Min UCCSD init:', name + '-', min([(x, ix) for ix, x in enumerate(e_uccsd) if not np.isnan(x)])[1])
print('Min UCCSDT init:', name + '-', min([(x, ix) for ix, x in enumerate(e_uccsdt) if not np.isnan(x)])[1])

udmrg_map = {x['init']: x['udmrg'] for x in data['udmrg_5k_ref'][str(ix)]}
e_udmrg = 22140 + np.array([udmrg_map.get(x[1], np.nan) for x in sorted([(x['uhf'], x['init']) for x in data['uccsd'][str(ix)]])])

plt.plot(np.arange(24), e_uhf, '-', linewidth=2, marker='o', mec=colors[0], markersize=5, color=colors[0], mfc='white', mew=2, label='UHF')
plt.plot(np.arange(24), e_uccsd + np.where(e_uccsd_conv, 0, np.nan), '-', linewidth=2, marker='^', mec=colors[1], color=colors[1], mfc='white', mew=2, label='UCCSD (converged)')
plt.plot(np.arange(24), e_uccsd_t + np.where(e_uccsd_conv, 0, np.nan), '-', linewidth=2, marker='s', mec=colors[2], color=colors[2], mfc='white', mew=2, label='UCCSD(T) (converged)')
plt.plot(np.arange(24), e_uccsdt + np.where(e_uccsdt_conv, 0, np.nan), '-', linewidth=2, marker='p', mec=colors[3], color=colors[3], mfc='white', mew=2, label='UCCSDT (converged)')
if not all(np.isnan(e_udmrg)):
    plt.plot(np.arange(24), e_udmrg, '-', linewidth=2, marker='h', mec=colors[4], color=colors[4], mfc='white', mew=2, label='UDMRG ($D$=5k)')
plt.plot(np.arange(24), e_uccsd + np.where(e_uccsd_conv, np.nan, 0), '-', linewidth=2, marker='X', mec=colors[1], color=colors[1], mfc='white', mew=1.5, label='UCCSD (not converged)')
plt.plot(np.arange(24), e_uccsd_t + np.where(e_uccsd_conv, np.nan, 0), '-', linewidth=2, marker='X', mec=colors[2], color=colors[2], mfc='white', mew=1.5, label='UCCSD(T) (not converged)')
plt.plot(np.arange(24), e_uccsdt + np.where(e_uccsdt_conv, np.nan, 0), '-', linewidth=2, marker='X', mec=colors[3], color=colors[3], mfc='white', mew=1.5, label='UCCSDT (not converged)')

plt.plot(np.arange(24), t1diag_uccsd, ':', linewidth=2, marker='^', mec=colors[0], color=colors[0], mfc='white', label='UCCSD $T_1$ diagnostic', mew=2, alpha=0.75)

minx = lambda x: min(x[~np.isnan(x)]) if len(x[~np.isnan(x)]) != 0 else np.nan

plt.axhline(minx(e_uhf), linestyle='--', linewidth=1, color=colors[0], zorder=-2)
plt.axhline(minx(e_uccsd + np.where(e_uccsd_conv, 0, np.nan)), linestyle='--', linewidth=1, color=colors[1], zorder=-2)
plt.axhline(minx(e_uccsd_t + np.where(e_uccsd_conv, 0, np.nan)), linestyle='--', linewidth=1, color=colors[2], zorder=-2)
plt.axhline(minx(e_uccsdt + np.where(e_uccsdt_conv, 0, np.nan)), linestyle='--', linewidth=1, color=colors[3], zorder=-2)

plt.text(0.5, 0.195, name, fontsize=18)

plt.ylabel('Total energy - -22140 (Hartree)')
plt.xlabel('UHF reference index')
plt.ylim((-0.4, 0.24))
plt.xlim((-0.5, 23.5))
plt.legend(loc='center left', frameon=False, fontsize=14, ncols=2, bbox_to_anchor=(0.00, 0.47))
plt.subplots_adjust(left=0.14, right=0.96, top=0.95, bottom=0.12)
plt.savefig('04-sm-ref-dep-cc-%02d.png' % iiid, dpi=300)
