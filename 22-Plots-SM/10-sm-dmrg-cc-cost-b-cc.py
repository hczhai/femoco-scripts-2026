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
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))

colors = data['colors']

sci_notation = lambda x: "%.2f\\times 10^{%.0f}" % tuple([float(x) for x in ("%.2E" % x).split('E+')])

noccs_26 = [58 - data['fno_ucc']['26'][0][str(ir)]['nfrozen'] for ir in range(1, 8)]
times_26 = [data['fno_ucc']['26'][0][str(ir)]['uccsdtq_time'] * (96 if ir >= 5 else 48)
    / 3600 / data['fno_ucc']['26'][0][str(ir)]['uccsdtq_niter'] for ir in range(1, 8)]

reg = scipy.stats.linregress(np.log(noccs_26[-2:]), np.log(times_26[-2:]))
cpu_hours_58 = np.exp(reg.intercept + reg.slope * np.log(58))

x_reg = np.array([10, 60])

plt.plot(x_reg, np.exp(reg.intercept + reg.slope * np.log(x_reg)), '--', linewidth=1.5, color=colors[0],
    label='Fit $T = (N_{\\text{occ}}\\ /\\ %.2f)^{%.2f} $' % ((1 / np.exp(reg.intercept)) ** (1/reg.slope), reg.slope))
plt.plot(noccs_26, np.array(times_26), marker='H', mfc='white', mec=colors[6], color=colors[6], linewidth=0, mew=1.5, markersize=7, label='UCCSDTQ cost')
plt.plot(noccs_26[-2:], np.array(times_26)[-2:], ' ', marker='H', markersize=7, mec='none', color='none', mfc=colors[0], mew=1.5, alpha=0.34)
plt.xlim(10, 60)
plt.yscale('log')
plt.xscale('log')
from matplotlib.ticker import LogLocator
axs.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 2), numticks=1))
axs.set_xticks([10] + noccs_26 + [58], [''] + list(map(str, noccs_26 + [58])))
axs.set_ylabel("CPU core hours per iteration $T$")
axs.set_xlabel("Number of occupied orbitals $N_{\\text{occ}}$")

plt.text(15.2, 0.1, 'Extrapolated $N_{\\text{occ}}=58$ cost: $T=%s$' % sci_notation(cpu_hours_58), fontsize=12)

axs.legend(loc='upper left', fontsize=12)
plt.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.14, wspace=0.0, hspace=0.0)
plt.savefig('10-sm-dmrg-cc-cost-b-cc.png', dpi=300)
