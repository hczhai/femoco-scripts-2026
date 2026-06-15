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

data = {'fno_udmrg': {ix: [{}] for ix in [25, 26, 27]}}

for ix in [25, 26, 27]:
    for ir in range(1, 7):
        data['fno_udmrg'][ix][0][ir] = {'nfrozen': 50 - ir * 4}
        xre = '-re' if ir in [2, 3, 4] else ''
        with open('../06-FEMOCO-LLDUC/15-FNO-UDMRG/%02d-fno-dmrg-f%d%s-%d.out' % (ir, 50 - ir * 4, xre, ix), "r") as f:
            bond1 = 0
            for l in f.readlines():
                if l.startswith('EHF ='):
                    data['fno_udmrg'][ix][0][ir]['uhf'] = float(l.split()[2])
                elif 'Bond dimension =' in l:
                    bond1 = int(l.split('Bond dimension =')[1].split()[0])
                elif l.startswith('STATE 0 DMRG energy ='):
                    data['fno_udmrg'][ix][0][ir]['udmrg_%d' % bond1] = float(l.split()[5])
        with open('../06-FEMOCO-LLDUC/15-FNO-UDMRG/%02d-fno-dmrg-rev-f%d-%d.out' % (ir + 10, 50 - ir * 4, ix), "r") as f:
            bond2 = 0
            for l in f.readlines():
                if 'Bond dimension =' in l:
                    bond2 = int(l.split('Bond dimension =')[1].split()[0])
                elif l.startswith('STATE 0 DMRG energy ='):
                    data['fno_udmrg'][ix][0][ir]['udmrg_%d' % bond2] = float(l.split()[5])
        assert bond1 == 8000 and bond2 == 5000
        assert 'udmrg_5000' in data['fno_udmrg'][ix][0][ir] and 'udmrg_8000' in data['fno_udmrg'][ix][0][ir]
        with open('../06-FEMOCO-LLDUC/15-FNO-UDMRG/%02d-fno-dmrg-npdm-f%d-%d.out' % (ir + 20, 50 - ir * 4, ix), "r") as f:
            for l in f.readlines():
                if l.startswith('SSQ ='):
                    data['fno_udmrg'][ix][0][ir]['udmrg_%d_ssq' % bond2] = [float(l.split()[2]), float(l.split()[3])]
        with open('../06-FEMOCO-LLDUC/15-FNO-UDMRG/%02d-fno-dmrg-npdm-f%d-%d.out' % (ir + 30, 50 - ir * 4, ix), "r") as f:
            for l in f.readlines():
                if l.startswith('SSQ ='):
                    data['fno_udmrg'][ix][0][ir]['udmrg_%d_ssq' % bond1] = [float(l.split()[2]), float(l.split()[3])]
        with open('../06-FEMOCO-LLDUC/15-FNO-UDMRG/%02d-invm-%d.out' % (ir + 40, ix), "r") as f:
            for l in f.readlines():
                if l.startswith('EXTRA ='):
                    data['fno_udmrg'][ix][0][ir]['udmrg_fit_invm'] = float(l.split()[2]) - 22140
        with open('../06-FEMOCO-LLDUC/15-FNO-UDMRG/%02d-logm-%d.out' % (ir + 50, ix), "r") as f:
            for l in f.readlines():
                if l.startswith('EXTRA ='):
                    data['fno_udmrg'][ix][0][ir]['udmrg_fit_logm'] = float(l.split()[2]) - 22140

with open("11-llduc-fno-udmrg.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of FNO-UDMRG calculations: ', sum([len(data['fno_udmrg'][ix]) for ix in [25, 26, 27]]))
