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

data = {'fe4s4-uccsdtq': {ix: [] for ix in range(3)}}

for ix in range(3):
    d = {}
    with open("../04-Fe4S4/01-UHF-UCC-UDMRG/06-ccsdtq-%d.out" % (ix, ), "r") as f:
        for l in f.readlines():
            if l.startswith('INIT ='):
                d['init'] = l.split()[2]
            elif l.startswith('EHF  ='):
                d['uhf'] = float(l.split()[2])
            elif l.startswith('E(UCCSD) ='):
                d['uccsd'] = float(l.split()[2])
            elif l.startswith('E(UCCSDT) ='):
                d['uccsdt'] = float(l.split()[2])
            elif l.startswith('E(UCCSDTQ) ='):
                d['uccsdtq'] = float(l.split()[2])
            elif l.startswith('UCCSD not converged'):
                d['uccsd_conv'] = False
            elif l.startswith('UCCSDT not converged'):
                d['uccsdt_conv'] = False
            elif l.startswith('UCCSDTQ not converged'):
                d['uccsdtq_conv'] = False
            elif l.startswith('UCCSD converged'):
                d['uccsd_conv'] = True
            elif l.startswith('UCCSDT converged'):
                d['uccsdt_conv'] = True
            elif l.startswith('UCCSDTQ converged'):
                d['uccsdtq_conv'] = True
    assert 'init' in d and 'uhf' in d and 'uccsd' in d and 'uccsd_conv' in d
    assert 'uccsdt' in d and 'uccsdt_conv' in d and 'uccsdtq' in d and 'uccsdtq_conv' in d
    data['fe4s4-uccsdtq'][ix].append(d)

for ix in range(3):
    data['fe4s4-uccsdtq'][ix][0]['uccsd-spin-corr'] = []
    with open("../04-Fe4S4/01-UHF-UCC-UDMRG/22-fe4-ccsd-spin-corr-%d.out" % (ix, ), "r") as f:
        for l in f.readlines():
            if l.startswith('SPIN-CORR '):
                data['fe4s4-uccsdtq'][ix][0]['uccsd-spin-corr'].append([float(x) for x in l.split()[3:]])
    data['fe4s4-uccsdtq'][ix][0]['uccsdt-spin-corr'] = []
    with open("../04-Fe4S4/01-UHF-UCC-UDMRG/23-fe4-ccsdt-spin-corr-%d.out" % (ix, ), "r") as f:
        for l in f.readlines():
            if l.startswith('SPIN-CORR '):
                data['fe4s4-uccsdtq'][ix][0]['uccsdt-spin-corr'].append([float(x) for x in l.split()[3:]])
    data['fe4s4-uccsdtq'][ix][0]['uccsdtq-spin-corr'] = []
    with open("../04-Fe4S4/01-UHF-UCC-UDMRG/28-fe4-ccsdtq-spin-corr-%d.out" % (ix, ), "r") as f:
        for l in f.readlines():
            if l.startswith('SPIN-CORR '):
                data['fe4s4-uccsdtq'][ix][0]['uccsdtq-spin-corr'].append([float(x) for x in l.split()[3:]])

with open("21-fe4s4-uccsdtq.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Fe4S4: Number of UCCSDTQ calculations: ', sum([len(data['fe4s4-uccsdtq'][ix]) for ix in range(3)]))
