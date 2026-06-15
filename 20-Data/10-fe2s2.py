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

data = {'fe2s2': {'std': {}, 'mix': {}}}
data['fe2s2']['std']['udmrg'] = {'forward': {}, 'reverse': {}, 'reverse_ssq': {}, 'reverse_sz': {}}
data['fe2s2']['mix']['udmrg'] = {'forward': {}, 'reverse': {}, 'reverse_ssq': {}, 'reverse_sz': {}}

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

for ik, key in enumerate(['std', 'mix']):
    with open("../03-Fe2S2/00-UCC/%02d-fe2-cc-%s-0.out" % (ik, key), "r") as f:
        for l in f.readlines():
            if l.startswith('INIT ='):
                data['fe2s2'][key]['init'] = l.split()[2]
            elif l.startswith('EHF  = '):
                data['fe2s2'][key]['uhf'] = float(l.split()[2])
            elif l.startswith('Init t2, MP2 energy ='):
                data['fe2s2'][key]['ump2'] = float(l.split()[-1]) + data['fe2s2'][key]['uhf']
            elif l.startswith('E(UCCSD) ='):
                data['fe2s2'][key]['uccsd'] = float(l.split()[2])
            elif l.startswith('E(UCCSDT) ='):
                data['fe2s2'][key]['uccsdt'] = float(l.split()[2])
            elif l.startswith('E(UCCSDTQ) ='):
                data['fe2s2'][key]['uccsdtq'] = float(l.split()[2])
            elif l.startswith('E(UCCSDTQP) ='):
                data['fe2s2'][key]['uccsdtqp'] = float(l.split()[2])
    with open("../03-Fe2S2/00-UCC/%02d-fe2-ccsd-t-%s-0.out" % (ik + 2, key), "r") as f:
        for l in f.readlines():
            if l.startswith('ECCSD(T) ='):
                data['fe2s2'][key]['uccsd(t)'] = float(l.split()[2])
            elif l.startswith('<S^2> ='):
                data['fe2s2'][key]['uccsd(t)_ssq'] = [float(l.split()[2]), float(l.split()[5]) + 1]
            elif l.startswith('SZ ::'):
                data['fe2s2'][key]['uccsd(t)_sz'] = [float(l.split()[2]), float(l.split()[3])]
    with open("../03-Fe2S2/00-UCC/%02d-fe2-cc-ssq-%s.out" % (ik + 4, key), "r") as f:
        for l in f.readlines():
            if '<S^2>' in l:
                order = int(l.split()[2])
                data['fe2s2'][key]['uccsdtqp'[:order + 3] + '_ssq'] = [float(l.split()[5]), float(l.split()[8]) + 1]
            elif l.startswith('SZ ::'):
                order = int(l.split()[4])
                data['fe2s2'][key]['uccsdtqp'[:order + 3] + '_sz'] = [float(l.split()[5]), float(l.split()[6])]
    with open("../03-Fe2S2/00-UCC/%02d-fe2-uhf-sz-%s.out" % (ik + 6, key), "r") as f:
        for l in f.readlines():
            if '<S^2> =' in l:
                data['fe2s2'][key]['uhf_ssq'] = [float(l.split()[-4]), float(l.split()[-1])]
            elif l.startswith('SZ ::'):
                data['fe2s2'][key]['uhf_sz'] = [float(l.split()[2]), float(l.split()[3])]
    with open("../03-Fe2S2/01-SA-DMRG/%02d-dmrg-0.out" % (ik + 2, ), "r") as f:
        for l in f.readlines():
            if l.startswith('STATE 0 DMRG energy ='):
                data['fe2s2'][key]['sa_dmrg_12k'] = float(l.split()[5])
    with open("../03-Fe2S2/02-UDMRG/%02d-dmrg-ucc%s-0.out" % (ik * 2 + 4, ['', '-mixed'][ik]), "r") as f:
        bond_dim = 0
        for l in f.readlines():
            if 'Bond dimension =' in l:
                bond_dim = int(l.split('Bond dimension =')[1].split()[0])
            elif l.startswith('STATE 0 DMRG energy ='):
                data['fe2s2'][key]['udmrg']['forward'][bond_dim] = {'udmrg': float(l.split()[5])}
    with open("../03-Fe2S2/02-UDMRG/%02d-dmrg-ucc%s-rev-0.out" % (ik + 7, ['', '-mix'][ik]), "r") as f:
        isw = -1
        raw_data = data['fe2s2'][key]['udmrg']['reverse']
        for l in f.readlines():
            if l.startswith('Sweep ='):
                isw = int(l.split()[2])
                if isw not in raw_data:
                    raw_data[isw] = {}
                else:
                    assert int(l.split()[11]) == raw_data[isw]['bond']
                raw_data[isw]['bond'] = int(l.split()[11])
            elif l.startswith('Time elapsed ='):
                assert isw != -1
                raw_data[isw]['energy'] = float(l.split('| E =')[1].split()[0])
                raw_data[isw]['dw'] = float(l.split('| DW =')[1].split()[0])
            elif l.startswith('Time sweep ='):
                assert isw != -1
                raw_data[isw]['time'] = float(l.split()[3])
                raw_data[isw]['flop'] = trans_flops(l.split('|')[-1])
    with open("../03-Fe2S2/02-UDMRG/%02d-dmrg-ucc%s-1pdm-0.out" % (ik + 9, ['', '-mix'][ik]), "r") as f:
        for l in f.readlines():
            if ' :: S ' in l:
                bdim = int(l.split()[3])
                data['fe2s2'][key]['udmrg']['reverse_ssq'][bdim] = [float(l.split()[4]), float(l.split()[5])]
            elif ' :: SZ ' in l:
                bdim = int(l.split()[3])
                data['fe2s2'][key]['udmrg']['reverse_sz'][bdim] = [float(l.split()[4]), float(l.split()[5])]

with open("10-fe2s2.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
