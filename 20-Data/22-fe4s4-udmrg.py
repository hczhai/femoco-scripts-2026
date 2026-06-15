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


import json, numpy as np

data = {'fe4s4-udmrg': {ix: {'forward': [{}], 'reverse': [{}]} for ix in range(3)}}

def trans_flops(x):
    a, b = x.split()
    if b == 'MFLOP/SWP':
        return float(a) * 1024 ** 2
    elif b == 'GFLOP/SWP':
        return float(a) * 1024 ** 3
    elif b == 'TFLOP/SWP':
        return float(a) * 1024 ** 4
    elif b == 'PFLOP/SWP':
        return float(a) * 1024 ** 5
    elif b == 'EFLOP/SWP':
        return float(a) * 1024 ** 6
    else:
        assert False

def parse_dmrg_output(filename, res, n_sites, isw_start=0):
    raw_data = {}
    for xre in res:
        with open(filename % xre, "r") as f:
            isw, nnodes = -1, -1
            for l in f.readlines():
                if l.startswith('MPI INIT:'):
                    nnodes = int(l.split()[-1])
                elif l.startswith('Expectation |'):
                    break
                elif l.startswith('Sweep ='):
                    isw = int(l.split()[2]) + isw_start
                    if isw not in raw_data:
                        raw_data[isw] = {}
                        raw_data[isw]['sweep'] = {}
                        raw_data[isw]['files'] = []
                        raw_data[isw]['nnodes'] = []
                    else:
                        assert int(l.split()[11]) == raw_data[isw]['bond']
                    raw_data[isw]['bond'] = int(l.split()[11])
                    raw_data[isw]['files'].append(xre)
                    raw_data[isw]['nnodes'].append(nnodes)
                elif 'FLOPS =' in l:
                    assert isw != -1
                    sites = int(l.split('Site =')[1].split('-')[0])
                    ener = float(l.split(' E =')[1].split()[0])
                    dw = float(l.split(' Error =')[1].split()[0])
                    flops = float(l.split(' FLOPS =')[1].split()[0])
                    tdav = float(l.split(' Tdav =')[1].split()[0])
                    t = float(l.split(' T =')[1].split()[0])
                    raw_data[isw]['sweep'][sites] = [ener, dw, flops, tdav, t]
                elif l.startswith('Time elapsed ='):
                    assert isw != -1
                    raw_data[isw]['energy'] = float(l.split('| E =')[1].split()[0])
                    raw_data[isw]['dw'] = float(l.split('| DW =')[1].split()[0])
                elif l.startswith('Time sweep ='):
                    assert isw != -1
                    raw_data[isw]['time'] = float(l.split()[3])
                    raw_data[isw]['flop'] = trans_flops(l.split('|')[-1])

    for isw in sorted(raw_data):
        if len(raw_data[isw]['files']) != 1:
            assert len(raw_data[isw]['sweep']) == n_sites - 1
            raw_data[isw]['energy'] = min([x[0] for x in raw_data[isw]['sweep'].values()])
            raw_data[isw]['dw'] = max([x[1] for x in raw_data[isw]['sweep'].values()])
            raw_data[isw]['time'] = sum([x[4] for x in raw_data[isw]['sweep'].values()])
            raw_data[isw]['flop'] = sum([x[2] * x[3] for x in raw_data[isw]['sweep'].values()])
        del raw_data[isw]['sweep'], raw_data[isw]['files']
        # assert min(raw_data[isw]['nnodes']) == max(raw_data[isw]['nnodes'])
        raw_data[isw]['nnodes'] = max(raw_data[isw]['nnodes'])
        if 'energy' not in raw_data[isw] and isw - isw_start == len(raw_data) - 1:
            del raw_data[isw]
    return raw_data

for ix in range(3):
    data['fe4s4-udmrg'][ix]['forward'][0] |= parse_dmrg_output(
        "../04-Fe4S4/01-UHF-UCC-UDMRG/05-fe4-udmrg-%d.out", [ix], n_sites=36, isw_start=0)
    data['fe4s4-udmrg'][ix]['reverse'][0] |= parse_dmrg_output(
        "../04-Fe4S4/01-UHF-UCC-UDMRG/07-fe4-udmrg-rev-%d.out", [ix], n_sites=36, isw_start=0)

def polynomial_fit(x, y, degree):
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.linear_model import LinearRegression
    from sklearn.pipeline import Pipeline
    from sklearn.metrics import r2_score
    import numpy as np
    x = np.array(x)
    y = np.array(y)

    poly_features = PolynomialFeatures(degree=degree)
    linear_reg = LinearRegression()
    poly_model = Pipeline([('poly', poly_features), ('linear', linear_reg)])
    poly_model.fit(x.reshape(-1, 1), y)
    coefficients = poly_model.named_steps['linear'].coef_
    intercept = poly_model.named_steps['linear'].intercept_
    y_pred = poly_model.predict(x.reshape(-1, 1)).flatten()
    r2 = r2_score(y, y_pred)
    return poly_model, coefficients, intercept, r2

def exp_log_fit(x, y):
    from scipy.optimize import curve_fit
    M2, E = np.array(x), np.array(y) + 327
    E_model = lambda M, E_inf, A, kappa: E_inf + A * np.exp(-kappa * M2)

    E_inf_init, A_init, kappa_init = -0.25, 0.5, 0.01
    (E_inf, A, kappa), pcov = curve_fit(E_model, M2, E, xtol=1E-14, gtol=1e-14, ftol=1e-14, 
        p0=[E_inf_init, A_init, kappa_init], maxfev=20000)
    # perr = np.sqrt(np.diag(pcov))
    # poly_model = lambda xx, E_inf=E_inf, A=A, kappa=kappa: E_model(xx, E_inf, A, kappa)
    return E_inf - 327

for ix in range(3):
    data_x = {x['bond']: x['energy'] for _, x in
        sorted(data['fe4s4-udmrg'][ix]['reverse'][0].items(), key=lambda x: int(x[0]))}
    inv_bond_x = [1 / x for x, _ in sorted(data_x.items())]
    energy_x = [y for _, y in sorted(data_x.items())]
    fits = {i: polynomial_fit(inv_bond_x, energy_x, i) for i in range(1, 5)}
    data['fe4s4-udmrg'][ix]['fit'] = {i: fits[i][2] for i in range(1, 5)}
    data['fe4s4-udmrg'][ix]['exp-log-fit'] = exp_log_fit([np.log(x) ** 2 for x, _ in sorted(data_x.items())], energy_x)

for ix in range(3):
    data['fe4s4-udmrg'][ix]['spin-corr'] = {}
    for ib, bond in enumerate([12000, 8000, 3000]):
        data['fe4s4-udmrg'][ix]['spin-corr'][bond] = []
        with open("../04-Fe4S4/01-UHF-UCC-UDMRG/%02d-fe4-spin-corr-%dk-%d.out" % (12 + ib, bond // 1000, ix), "r") as f:
            for l in f.readlines():
                if l.startswith('SPIN-CORR '):
                    data['fe4s4-udmrg'][ix]['spin-corr'][bond].append([float(x) for x in l.split()[3:]])

with open("22-fe4s4-udmrg.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Fe4S4: Number of UDMRG (12k) calculations: ', sum([len(data['fe4s4-udmrg'][ix]['forward']) for ix in range(3)]))
