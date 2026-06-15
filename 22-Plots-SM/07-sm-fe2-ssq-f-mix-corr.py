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
xs = ['uhf', 'uccsd', 'uccsd(t)', 'uccsdt', 'uccsdtq', 'uccsdtqp']
ds = [100, 200, 400, 600, 800, 1000, 1200, 1600, 2000, 2400]

cc_ener = np.array([data['fe2s2']['mix'][x] for x in xs])
cc_sexpt = np.array([(data['fe2s2']['mix']["%s_ssq" % x][1] - 1) / 2 for x in xs])

e_udmrg = {x['bond']: x['energy'] for _, x in sorted(data['fe2s2']['mix']['udmrg']['reverse'].items(), key=lambda x: int(x[0]))}
dmrg_ener = np.array([e_udmrg[d] for d in ds])
dmrg_sexpt = np.array([(data['fe2s2']['mix']['udmrg']['reverse_ssq'][str(d)][1] - 1) / 2 for d in ds])

e_ref = data['fe2s2']['mix']['sa_dmrg_12k']

plt.plot(cc_sexpt, cc_ener - e_ref, '-', linewidth=2, marker='o', mec=colors[6],  markersize=7, color=colors[6], mfc='white', clip_on=False, label='UCC')
plt.plot(dmrg_sexpt, dmrg_ener - e_ref, '--', linewidth=2, marker='p', mec=colors[4],  markersize=7, color=colors[4], mfc='white', clip_on=False, label='UDMRG')

plt.ylabel('Total energy error (Hartree)')
plt.xlabel('Total spin expectation $\\langle S \\rangle$')
plt.yscale('log')
plt.legend(loc='upper left', frameon=False, fontsize=14, ncols=1)
plt.subplots_adjust(left=0.15, right=0.96, top=0.95, bottom=0.12)
plt.savefig('07-sm-fe2-ssq-f-mix-corr.png', dpi=300)
