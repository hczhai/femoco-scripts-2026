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

data = {'udmrg_18k': {25: {'forward': [{}], 'reverse': [{}]}, 26: {'forward': [{}], 'reverse': [{}]}}}

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

### 25

data['udmrg_18k'][25]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/19-UDMRG-Legacy/22-dmrg-fm-newton%s-25.out",
    [''], n_sites=76, isw_start=0)

data['udmrg_18k'][25]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/19-UDMRG-Legacy/27-uhf-dmrg-5k-25%s.out",
    ['-bak', ''], n_sites=76, isw_start=len(data['udmrg_18k'][25]['forward'][0]))

data['udmrg_18k'][25]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/19-UDMRG-Legacy/28-fm-7k-25.out%s",
    ['.1', '.2'], n_sites=76, isw_start=len(data['udmrg_18k'][25]['forward'][0]))

data['udmrg_18k'][25]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/02-UDMRG-25/10-fm-dmrg-18k%s-25.out",
    ['', '-re', '-re2', '-re3', '-re4', '-re5', '-re6', '-re7', '-re8'],
    n_sites=76, isw_start=len(data['udmrg_18k'][25]['forward'][0]))

data['udmrg_18k'][25]['reverse'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/02-UDMRG-25/11-fm-dmrg-18k-rev%s-25.out",
    ['', '-re', '-re2', '-re3', '-re4', '-re5', '-re6', '-re7'], n_sites=76, isw_start=0)

data['udmrg_18k'][25]['reverse'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/02-UDMRG-25/17-fm-dmrg-12k-rev%s-25.out",
    ['', '-re', '-re2'], n_sites=76, isw_start=len(data['udmrg_18k'][25]['reverse'][0]))

### 26

data['udmrg_18k'][26]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/19-UDMRG-Legacy/22-dmrg-fm-newton%s-26.out",
    [''], n_sites=76, isw_start=0)

data['udmrg_18k'][26]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/19-UDMRG-Legacy/27-uhf-dmrg-5k%s-26.out",
    [''], n_sites=76, isw_start=len(data['udmrg_18k'][26]['forward'][0]))

data['udmrg_18k'][26]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/19-UDMRG-Legacy/28-fm-7k-26.out%s",
    ['.1', '.2'], n_sites=76, isw_start=len(data['udmrg_18k'][26]['forward'][0]))

data['udmrg_18k'][26]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/19-UDMRG-Legacy/33-uhf-dmrg-10k%s-26.out",
    [''], n_sites=76, isw_start=len(data['udmrg_18k'][26]['forward'][0]))

data['udmrg_18k'][26]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/19-UDMRG-Legacy/%d-uhf-dmrg-12k%s-26.out",
    [(34, ''), (35, '-re')], n_sites=76, isw_start=len(data['udmrg_18k'][26]['forward'][0]))

data['udmrg_18k'][26]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/03-UDMRG-26/02-fm-dmrg%s-26.out",
    ['', '-re', '-re2'], n_sites=76, isw_start=len(data['udmrg_18k'][26]['forward'][0]))

data['udmrg_18k'][26]['forward'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/03-UDMRG-26/05-fm-dmrg-18k%s-26.out",
    ['', '-re', '-re2'], n_sites=76, isw_start=len(data['udmrg_18k'][26]['forward'][0]))

data['udmrg_18k'][26]['reverse'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/03-UDMRG-26/06-fm-dmrg-18k-rev%s-26.out",
    ['', '-re', '-re2', '-re3', '-re4', '-re5', '-re6', '-re7'], n_sites=76, isw_start=0)

data['udmrg_18k'][26]['reverse'][0] |= parse_dmrg_output(
    "../06-FEMOCO-LLDUC/03-UDMRG-26/09-fm-dmrg-12k-rev%s-26.out",
    ['', '-re', '-re2'], n_sites=76, isw_start=len(data['udmrg_18k'][26]['reverse'][0]))

def exp_log_fit(x, y):
    from scipy.optimize import curve_fit
    M2, E = np.array(x), np.array(y) + 22140
    E_model = lambda M, E_inf, A, kappa: E_inf + A * np.exp(-kappa * M2)

    E_inf_init, A_init, kappa_init = -0.41, 0.28, 0.02
    (E_inf, A, kappa), pcov = curve_fit(E_model, M2, E, xtol=1E-14, gtol=1e-14, ftol=1e-14, 
        p0=[E_inf_init, A_init, kappa_init], maxfev=20000)
    # perr = np.sqrt(np.diag(pcov))
    # poly_model = lambda xx, E_inf=E_inf, A=A, kappa=kappa: E_model(xx, E_inf, A, kappa)
    return E_inf - 22140

for ix in [25, 26]:
    data_x = {x['bond']: x['energy'] for _, x in
        sorted(data['udmrg_18k'][ix]['reverse'][0].items(), key=lambda x: int(x[0]))
        if x['bond'] >= 5000}
    inv_bond_x = [1 / x for x, _ in sorted(data_x.items())]
    energy_x = [y for _, y in sorted(data_x.items())]
    fits = {i: np.polynomial.polynomial.polyfit(inv_bond_x, energy_x, i, full=True)[0] for i in range(1, 7)}
    data['udmrg_18k'][ix]['fit'] = {i: fits[i][0] for i in range(1, 7)}
    data['udmrg_18k'][ix]['exp-log-fit'] = exp_log_fit([np.log(x) ** 2 for x, _ in sorted(data_x.items())], energy_x)

for ix in [26]:
    data['udmrg_18k'][ix]['reverse_ssq'] = [{}]
    for jx, bond in enumerate([5, 6, 7, 8, 9, 11]):
        with open("../06-FEMOCO-LLDUC/17-UDMRG-NPDM/%02d-fm-dmrg-%dk%s-%d.out" % (jx + 1, bond, '-mpi' if bond >= 9 else '', ix), "r") as f:
            x = [l.split() for l in f.readlines() if len(l.split()) == 2][-1]
            data['udmrg_18k'][ix]['reverse_ssq'][0][bond * 1000] = [float(x[0]), float(x[1])]

for ix in [25, 26]:
    with open("../06-FEMOCO-LLDUC/05-UDMRG-DET/16-fm-dmrg-18k-det-%d.out" % (ix, ), "r") as f:
        lines = list(f.readlines())
        dets = [l.split()[2] for l in lines if l.startswith('DET ')]
        coeffs = [float(l.split()[4]) for l in lines if l.startswith('DET ')]
    data['udmrg_18k'][ix]['dets'] = [dets, coeffs]

with open("04-llduc-udmrg-18k.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of UDMRG (D=18k) calculations: ', sum([len(data['udmrg_18k'][ix]['forward']) for ix in [25, 26]]))
