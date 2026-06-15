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


import numpy as np, sys, os

fid = 44
pfid = 14
n_sweeps = 24
n_sites = 26 + (pfid % 10) * 4
iiid = int(sys.argv[1])
froz = 76 - n_sites

sw_eners = [[] for _ in range(n_sweeps)]
sw_dws = [[] for _ in range(n_sweeps)]
sw_bdims = [0 for  _ in range(n_sweeps)]
isw = -1
for l in open('%d-fno-dmrg-rev-f%d-%02d.out' % (pfid, froz, iiid), 'r').readlines():
    if l.startswith('Sweep ='):
        isw = int(l.split()[2])
        sw_bdims[isw] = int(l.split('Bond dimension = ')[1].split('|')[0].strip())
    elif l.startswith('Time elapsed'):
        isw = -1
    elif isw != -1 and ((' --> Site =' in l) or (' <-- Site =' in l)):
        ee, dw = [float(x.strip()) for x in l.split('E =')[1].split('FLOPS')[0].split('Error =')]
        sw_eners[isw].append(ee)
        sw_dws[isw].append(dw)
for x, y in zip(sw_eners, sw_dws):
    assert len(x) == n_sites - 1 and len(y) == n_sites - 1

eners, dws, bds = np.array([min(x) for x in sw_eners]), np.array([max(x) for x in sw_dws]), sw_bdims
eners = eners - -22140
eners, dws, bds = eners[3::4], dws[3::4], bds[3::4]

eners = np.array(eners)
bds = np.array(bds)

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import r2_score

def polynomial_fit(x, y, degree):
    poly_features = PolynomialFeatures(degree=degree)
    linear_reg = LinearRegression()

    poly_model = Pipeline([
        ('poly', poly_features),
        ('linear', linear_reg)
    ])

    poly_model.fit(x.reshape(-1, 1), y)
    coefficients = poly_model.named_steps['linear'].coef_
    intercept = poly_model.named_steps['linear'].intercept_

    y_pred = poly_model.predict(x.reshape(-1, 1)).flatten()
    r2 = r2_score(y, y_pred)

    return poly_model, coefficients, intercept, r2

import matplotlib.pyplot as plt

model, coefficients, intercept, r2 = polynomial_fit(1 / bds, eners, 2)
colors = ['#2F2D2E', '#BB4430', '#F18F01', '#99C24D', '#048BA8', '#23C9FF']

print('EXTRA = %24.15f' % intercept)

x_reg_x = np.linspace(0, 1 / (min(bds) - 500), 100)
plt.plot(x_reg_x, model.predict(x_reg_x.reshape(-1, 1)).flatten(), '--',
    linewidth=1, color=colors[0], label='$E(D=\\infty) = %.6f$ [$R^2 = %.6f$]' % (intercept, r2))
plt.plot(1 / bds, eners, linewidth=0, marker='o', mec=colors[0], color=colors[0], mfc='white')

plt.legend(loc='lower right')
plt.xlim((0, 1 / (min(bds) - 500)))
plt.ylabel("Sweep energy $E - (-22140) / \\mathrm{Ha}$")
plt.xticks(1 / bds[::3], ['1/%s' % x for x in bds[::3]])
plt.xlabel('1/M')
plt.tight_layout()
plt.savefig('%d-invm-f%d-%02d.png' % (fid, froz, iiid), dpi=400)
