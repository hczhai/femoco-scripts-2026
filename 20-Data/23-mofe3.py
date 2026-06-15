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

data = {'mofe3': {0: {}}}

with open('../12-MoFe3/03-pop-noryd-0.out', 'r') as f:
    data['mofe3'][0]['uks_conv'] = True
    for l in f.readlines():
        if 'SCF energy =' in l:
            data['mofe3'][0]['uks'] = float(l.split('SCF energy =')[1].split()[0])
        elif l.startswith('SCF not converged'):
            data['mofe3'][0]['uks_conv'] = False
        elif l.startswith('XC functionals ='):
            data['mofe3'][0]['uks_xc'] = l.split()[3]
        elif '<S^2> =' in l:
            data['mofe3'][0]['uks_ssq'] = [float(l.split()[-4]), float(l.split()[-1])]
        elif l.startswith('CHG '):
            data['mofe3'][0]['uks-elec-pop-noryd'] = [float(x) for x in l.split()[1:5]]
        elif l.startswith('SPN '):
            data['mofe3'][0]['uks-spin'] = [float(x) for x in l.split()[1:5]]

with open('../12-MoFe3/05-df-uccsd-0.out', 'r') as f:
    data['mofe3'][0]['uhf_conv'] = True
    data['mofe3'][0]['uccsd_conv'] = False
    for l in f.readlines():
        if 'SCF energy =' in l:
            data['mofe3'][0]['uhf'] = float(l.split('SCF energy =')[1].split()[0])
        elif l.startswith('SCF not converged'):
            data['mofe3'][0]['uhf_conv'] = False
        elif '<S^2> =' in l:
            data['mofe3'][0]['uhf_ssq'] = [float(l.split()[-4]), float(l.split()[-1])]
        elif l.startswith('E(UCCSD) ='):
            data['mofe3'][0]['uccsd'] = float(l.split()[2])
        elif l.startswith('UCCSD converged'):
            data['mofe3'][0]['uccsd_conv'] = True

mode = None
with open('../12-MoFe3/06-pop-noryd.out', 'r') as f:
    for l in f.readlines():
        if l.startswith(':: UCCSD'):
            mode = 'uccsd'
        elif l.startswith(':: UHF'):
            mode = 'uhf'
        elif l.startswith('CHG '):
            data['mofe3'][0]['%s-elec-pop-noryd' % mode] = [float(x) for x in l.split()[1:5]]
        elif l.startswith('SPN '):
            data['mofe3'][0]['%s-spin' % mode] = [float(x) for x in l.split()[1:5]]

with open("23-mofe3.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('MoFe3: Number of calculations: ', len([x for x in data['mofe3'][0] if x.endswith('-elec-pop-noryd')]))
