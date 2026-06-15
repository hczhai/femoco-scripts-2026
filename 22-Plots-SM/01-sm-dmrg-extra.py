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


import json
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('font', size=14)
plt.figure(figsize=(8, 6))

ha_to_mha = 1000

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/71-dmrg-arxiv-2208-02199.json", "r", encoding='utf-8'))

colors = data['colors']

extra_data = {}

for key in ['2feii-2feiii', 'pn', 'femoco']:
    bdim = np.array(sorted([int(k) for k in data['%s-energy' % key]]))
    ener = np.array([data['%s-energy' % key][str(k)] for k in bdim])
    dw = np.array([data['%s-dw' % key][str(k)] for k in bdim])
    reg = scipy.stats.linregress(dw[:-1], ener[:-1])
    extra_data[key] = [bdim, ener, dw, reg]
    print('%16s E[extra] = %15.6f error = %12.6f | Ref: E[extra] = %15.6f error = %12.6f' % (key,
        reg.intercept, (min(ener) - reg.intercept) / 5,
        data['%s-energy-extra' % key], data['%s-energy-extra-error' % key]))

extra_data['2feii-2feiii'].append('$\\mathrm{2Fe^{II}}$-$\\mathrm{2Fe^{III}}$')
extra_data['pn'].append('$\\mathrm{P^N}\\ $cluster')
extra_data['femoco'].append('FeMo-co')

x = np.array([0, 3E-4])

plt.grid(which='both', alpha=0.5)
for ik, key in enumerate(['2feii-2feiii', 'pn', 'femoco']):
    bdim, ener, dw, reg, label = extra_data[key]
    plt.plot(x, reg.slope * x * ha_to_mha, '--', linewidth=2, color=colors[ik])
    plt.plot(dw, (ener - reg.intercept) * ha_to_mha, ' ', marker='o', mfc='white', mec=colors[ik], color=colors[ik],
        markersize=7, mew=2, label=label)
    if ik == 0:
        plt.text(dw[-1], (ener[-1] - reg.intercept) * ha_to_mha - 5, "$D=%d$" % bdim[-1],
            verticalalignment='center', horizontalalignment='center', fontsize=14, color=colors[ik])
    else:
        plt.text(dw[-1] - 1E-5, (ener[-1] - reg.intercept) * ha_to_mha, "$D=%d$" % bdim[-1],
            verticalalignment='center', horizontalalignment='right', fontsize=14, color=colors[ik])

plt.legend(loc='upper left')
plt.xlim([0, 3E-4])
plt.ylim((-10, 100))
plt.xlabel("Largest discarded weight")
plt.ylabel("Energy error $E_{\\mathrm{DMRG}} - E_{\\mathrm{extrap}}\\ (\\mathrm{mHa})$")
plt.tight_layout()
plt.savefig('01-sm-dmrg-extra.png', dpi=300)
