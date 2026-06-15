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
import matplotlib
import scipy.stats

matplotlib.rc('font', size=14)
plt.figure(figsize=(12 / 2, 9 / 2), dpi=150)
f, axs = plt.subplots(1, 1, sharey=True)
f.set_size_inches(w=12 / 2, h=9 / 2)

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))

colors = data['colors']

sci_notation = lambda x: "%.2f\\times 10^{%.0f}" % tuple([float(x) for x in ("%.2E" % x).split('E+')])

rev_dmrg_26 = {x['bond']: x for _, x in
    sorted(data['udmrg_18k'][str(26)]['reverse'][0].items(), key=lambda x: int(x[0])) if x['bond'] >= 5000}
bonds_26 = [x[0] for x in sorted(rev_dmrg_26.items())]
times_26 = [x[1]['time'] / 3600 * x[1]["nnodes"] * 96 for x in sorted(rev_dmrg_26.items())]

reg = scipy.stats.linregress(np.log(bonds_26), np.log(times_26))
pred_bdim = 393000
cpu_hours_pred = np.exp(reg.intercept + reg.slope * np.log(pred_bdim))

x_reg = np.array([4500, 16500])

marker = ['o', '^', 's', 'p']

plt.plot(x_reg, np.exp(reg.intercept + reg.slope * np.log(x_reg)) / 1E6, '--', linewidth=2, color=colors[0],
    label='Fit $T = (D\\ /\\ %.2f)^{%.2f} $' % ((1 / np.exp(reg.intercept)) ** (1/reg.slope), reg.slope))
plt.plot(bonds_26, np.array(times_26) / 1E6, marker='p', mfc='white', mec=colors[4], color=colors[4], linewidth=0, mew=2, markersize=7, label='UDMRG cost')
plt.xlim(4500, 16500)
plt.xscale('log')
plt.yscale('log')
plt.minorticks_off()
plt.xticks([5000, 7000, 9000, 12000, 16000], list(map(str, [5000, 7000, 9000, 12000, 16000])))
plt.yticks([0.005, 0.01, 0.03, 0.05, 0.10], ['%.3f' % x for x in [0.005, 0.01, 0.03, 0.05, 0.10]])
axs.set_ylabel("CPU core hours per sweep $T (\\times 10^6)$")
axs.set_xlabel("MPS bond dimension $D$")

plt.text(5600, 0.003, 'Extrapolated $D=%s$ cost: $T=%s$' % (pred_bdim, sci_notation(cpu_hours_pred)), fontsize=12)

axs.legend(loc='upper left', fontsize=12)
plt.subplots_adjust(left=0.17, right=0.95, top=0.95, bottom=0.14, wspace=0.0, hspace=0.0)
plt.savefig('10-sm-dmrg-cc-cost-a-dmrg.png', dpi=300)
