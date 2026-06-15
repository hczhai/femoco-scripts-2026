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
import scipy.stats

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(8, 6))
plt.grid(which='major', axis='y', linestyle='--')

ha_to_kcal_mol = 627.5094740631

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/19-fe4s4-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/20-fe4s4-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/21-fe4s4-uccsdtq.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/22-fe4s4-udmrg.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/18-fe4s4-sa-dmrg.json", "r", encoding='utf-8'))

colors = data['colors']

rev_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-sadmrg-20k']['2']['reverse'][0].items(), key=lambda x: int(x[0]))}

bonds = sorted(rev_dmrg)
dmrg_rev_dw = [rev_dmrg[bond]["dw"] for bond in bonds]
dmrg_rev_energy = [rev_dmrg[bond]["energy"] for bond in bonds]
reg = scipy.stats.linregress(dmrg_rev_dw, dmrg_rev_energy)
e_ref = reg.intercept

e_uhf = np.array([[x[1] for x in sorted([(x['uccsd'], x['uhf']) for x in data['fe4s4-uccsd'][str(ix)]])][0] for ix in range(3)])
e_uccsd = np.array([[x[0] for x in sorted([(x['uccsd'], x['uhf']) for x in data['fe4s4-uccsd'][str(ix)]])][0] for ix in range(3)])
e_uccsd_t = np.array([[x[1] for x in sorted([(x['uccsd'], x['uccsd(t)']) for x in data['fe4s4-uccsd'][str(ix)]])][0] for ix in range(3)])
e_uccsdt = np.array([data['fe4s4-uccsdtq'][str(ix)][0]['uccsdt'] for ix in range(3)])
e_uccsdtq = np.array([data['fe4s4-uccsdtq'][str(ix)][0]['uccsdtq'] for ix in range(3)])

print('UCCSDTQ diff: %10.2f %10.2f %10.2f' % tuple((e_uccsdtq - e_ref) * ha_to_kcal_mol))

e_dmrg_8k = np.array([[x['energy'] for _, x in sorted(data['fe4s4-udmrg'][str(ix)]['reverse'][0].items(),
    key=lambda x: int(x[0])) if x['bond'] == 8000][-1] for ix in range(3)])
e_dmrg_12k = np.array([[x['energy'] for _, x in sorted(data['fe4s4-udmrg'][str(ix)]['forward'][0].items(),
    key=lambda x: int(x[0])) if x['bond'] == 12000][-1] for ix in range(3)])

e_uhf = (e_uhf - e_ref) * ha_to_kcal_mol
e_uccsd = (e_uccsd - e_ref) * ha_to_kcal_mol
e_uccsd_t = (e_uccsd_t - e_ref) * ha_to_kcal_mol
e_uccsdt = (e_uccsdt - e_ref) * ha_to_kcal_mol
e_uccsdtq = (e_uccsdtq - e_ref) * ha_to_kcal_mol
e_dmrg_8k = (e_dmrg_8k - e_ref) * ha_to_kcal_mol
e_dmrg_12k = (e_dmrg_12k - e_ref) * ha_to_kcal_mol

pp, dx = 0, 0.3

def separate(x):
    y = list(x)
    for i in range(len(x)):
        for j in range(len(x)):
            for k in range(len(x)):
                if y[k] > 5:
                    continue
                orig_i = y[i]
                if y[i] < y[j] and y[j] < y[k] and abs(y[i] - y[k]) < 0.4:
                    y[i], y[j], y[k] = (y[i] + y[j] + y[k]) / 3 - 0.2, (y[i] + y[j] + y[k]) / 3, (y[i] + y[j] + y[k]) / 3 + 0.2
                    if y[k] > 4.84:
                        dd = y[k] - 4.84
                        y[i], y[j], y[k] = y[i] - dd, y[j] - dd, y[k] - dd
                if y[i] > orig_i:
                    y[i] = orig_i
    for i in range(len(x)):
        for j in range(len(x)):
            if y[i] < y[j] and abs(y[i] - y[j]) < 0.2:
                y[i], y[j] = (y[i] + y[j]) / 2 - 0.1, (y[i] + y[j]) / 2 + 0.1
    return np.array(y)

for eners, cl in list(zip([e_uhf, e_uccsd, e_uccsd_t, e_uccsdt, e_uccsdtq, e_dmrg_12k], colors))[1:]:
    xeners = separate(eners / 5.5) * 5.5
    for kk in range(len(eners)):
        plt.plot([pp - dx - 0.45, pp + dx - 0.45], [eners[kk], eners[kk]], '-', linewidth=2.25, color=cl, clip_on=False)
        plt.text(pp + dx - 0.12 - 0.25, xeners[kk], 'BS%d' % (kk + 1), va='center', fontsize=14)
    pp += 2

plt.xticks([x for x in range(0, 10, 2)], ['UCCSD', 'UCCSD(T)', 'UCCSDT', 'UCCSDTQ', 'UDMRG\n($D=\\mathrm{12k}$)  '])

plt.ylabel('$E - E_{\\text{extrapolated SA-DMRG [BS3]}}$ (kcal/mol)')
plt.xlim((-1, 8.5))
plt.ylim((0, 32))
plt.subplots_adjust(left=0.10, right=0.985, top=0.95, bottom=0.1)
plt.savefig('22-sm-fe4-ener.png', dpi=300)
