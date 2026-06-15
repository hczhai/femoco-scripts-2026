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
data |= json.load(open("../20-Data/00-llduc-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/26-llduc-composite.json", "r", encoding='utf-8'))

colors = data['colors']

# 2S = (1 + 4 * <S^2>) ** 0.5 - 1
fsexpt = lambda x: ((1 + 4 * x) ** 0.5 - 1) / 2

ix = 26
ds = sorted([int(x) for x in data['udmrg_18k'][str(ix)]['reverse_ssq'][0]])
e_dmrg = {x['bond']: x['energy'] for _, x in
    sorted(data['udmrg_18k'][str(ix)]['reverse'][0].items(), key=lambda x: int(x[0]))}

dmrg_sexpt = np.array([fsexpt(data['udmrg_18k'][str(ix)]['reverse_ssq'][0][str(d)][0]) for d in ds])
dmrg_ener = 22140 + np.array([e_dmrg[d] for d in ds])

e_uhf = np.array([min([(x['uhf'], x['uhf_ssq'][0]) for x in data['uhf'][str(ix)]])[0] for ix in range(35)])
e_uhf_ssq = np.array([min([(x['uhf'], x['uhf_ssq'][0]) for x in data['uhf'][str(ix)]])[1] for ix in range(35)])
e_uccsd = np.array([min([(x['uccsd'], x.get('uccsd_ssq', [np.nan, np.nan])[0]) for x in data['uccsd'][str(ix)]])[0] for ix in range(35)])
e_uccsd_ssq = np.array([min([(x['uccsd'], x.get('uccsd_ssq', [np.nan, np.nan])[0]) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])
e_uccsdt = np.array([min([(x['uccsdt'], x.get('uccsdt_ssq', [np.nan, np.nan])[0]) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_uccsdt_ssq = np.array([min([(x['uccsdt'], x.get('uccsdt_ssq', [np.nan, np.nan])[0]) for x in data['uccsdt'][str(ix)]])[1] for ix in range(35)])

cc_sexpt = np.array([fsexpt(e_uhf_ssq[ix]), fsexpt(e_uccsd_ssq[ix]), fsexpt(e_uccsdt_ssq[ix])])
cc_ener = 22140 + np.array([e_uhf[ix], e_uccsd[ix], e_uccsdt[ix]])

e_ref = 22140 + min([v['energy'] for _, v in data['composite'].items() if 'energy' in v])

plt.plot(cc_sexpt, cc_ener, '-', linewidth=2, marker='o', mec=colors[6], markersize=7, color=colors[6], mfc='white', mew=2, clip_on=False, label='UHF/CC')
plt.plot(dmrg_sexpt, dmrg_ener, '--', linewidth=2, marker='p', mec=colors[4], markersize=7, color=colors[4], mfc='white', mew=2, clip_on=False, label='UDMRG')
plt.plot([1.5], [e_ref], '-.', color=colors[0], linewidth=2, marker='X', mfc='white', markersize=7, mew=1.5, clip_on=False, label='exact (estimated)')

plt.plot([1.5, 3.8], [e_ref, e_ref], '-.', color=colors[0])

plt.ylabel('Total energy $- -22140.0$ (Hartree)')
plt.xlabel('Total spin expectation $\\langle S \\rangle$')
plt.xlim((1.5, 3.8))
plt.legend(loc='upper left', frameon=False, fontsize=14, ncols=1)
plt.subplots_adjust(left=0.15, right=0.96, top=0.95, bottom=0.12)
plt.savefig('08-sm-fm-ssq.png', dpi=300)
