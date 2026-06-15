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

data = {'heisenberg': {}}

for ix, key in enumerate(['uhf', 'uccsd', 'uccsd_t']):
    data['heisenberg']['fit-' + key] = {'param': {}, 'predicted': [0] * 35}
    with open('../06-FEMOCO-LLDUC/20-Heisenberg/%02d-fit-%s-c3.out' % (ix, key), 'r') as f:
        for l in f.readlines():
            if l.startswith('CONST ='):
                data['heisenberg']['fit-' + key]['param']['e0'] = float(l.split()[-1])
            elif l.startswith('NM = FeTerm-FeL'):
                data['heisenberg']['fit-' + key]['param']['jlt'] = float(l.split()[-1])
            elif l.startswith('NM = FeL-FeL'):
                data['heisenberg']['fit-' + key]['param']['jlc'] = float(l.split()[-1])
            elif l.startswith('NM = FeR-FeR'):
                data['heisenberg']['fit-' + key]['param']['jrc'] = float(l.split()[-1])
            elif l.startswith('NM = FeR-Mo'):
                data['heisenberg']['fit-' + key]['param']['jrt'] = float(l.split()[-1])
            elif l.startswith('NM = FeL-FeR (axial)'):
                data['heisenberg']['fit-' + key]['param']['jca'] = float(l.split()[-1])
            elif l.startswith('NM = FeL-FeR (cross)'):
                data['heisenberg']['fit-' + key]['param']['jcc'] = float(l.split()[-1])
            elif l.startswith('RMSE ='):
                data['heisenberg']['fit-' + key]['param']['rmse'] = float(l.split()[-1])
            elif '::' in l and len(l.split()) == 3:
                data['heisenberg']['fit-' + key]['predicted'][int(l.split()[0])] = float(l.split()[-1])

for fn, key in zip(['03-sa-dmrg', '06-nsa-dmrg-bs'], ['sa-dmrg', 'nsa-dmrg']):
    data['heisenberg']['fit-uccsd_t'][key] = {}
    with open('../06-FEMOCO-LLDUC/20-Heisenberg/%s.out' % fn, 'r') as f:
        rel_ener = False
        ssq = False
        sz = False
        for l in f.readlines():
            if l.startswith('Emin ='):
                data['heisenberg']['fit-uccsd_t'][key]['emin'] = float(l.split()[-1])
            elif l.startswith('Relative energies:'):
                data['heisenberg']['fit-uccsd_t'][key]['rel_energy'] = {}
                rel_ener = True
            elif rel_ener and ' E = ' in l:
                data['heisenberg']['fit-uccsd_t'][key]['rel_energy'][int(l.split()[2])] = \
                    [float(x) for x in l.split('E =')[1].split()]
            elif l.startswith('<S^2> for each state:'):
                data['heisenberg']['fit-uccsd_t'][key]['ssq'] = {}
                rel_ener = False
                ssq = True
            elif ssq and ' SSQ = ' in l:
                data['heisenberg']['fit-uccsd_t'][key]['ssq'][int(l.split()[2])] = \
                    [float(x) for x in l.split('SSQ =')[1].split()]
            elif l.startswith('Sz for each state:'):
                data['heisenberg']['fit-uccsd_t'][key]['sz'] = {}
                ssq = False
                sz = True
            elif sz and ' iroot = ' in l:
                if int(l.split()[2]) not in data['heisenberg']['fit-uccsd_t'][key]['sz']:
                    data['heisenberg']['fit-uccsd_t'][key]['sz'][int(l.split()[2])] = [None] * 8
                data['heisenberg']['fit-uccsd_t'][key]['sz'][int(l.split()[2])][int(l.split()[6])] = \
                    [float(x) for x in l.split()[7:]]
            elif l.startswith('Spin correlation for each state:'):
                sz = False

with open("24-heisenberg.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Heisenberg: Number of models: ', len([x for x in data['heisenberg']]))
