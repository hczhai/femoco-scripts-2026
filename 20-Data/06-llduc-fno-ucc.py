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


import json, os

data = {'fno_ucc': {ix: [{}] for ix in [25, 26, 27]}}

for ix in [25, 26, 27]:
    with open('../06-FEMOCO-LLDUC/%d-FNO-UCC-%d/00-pre-occ.out' % (ix - 14, ix, ), "r") as f:
        data['fno_ucc'][ix][0]['occ'] = [float(x.strip()) for x in f.read().split(',')]
    data['fno_ucc'][ix][0]['flops'] = []
    for ir in range(1, 8):
        with open('../06-FEMOCO-LLDUC/18-FNO-UCC-FLOP/01-fno-nosymm-f%d.out' % (50 - ir * 4, ), "r") as f:
            data['fno_ucc'][ix][0]['flops'].append([int(l.split()[3]) for l in f.readlines() if l.startswith('TOTAL FLOP =')][-1])
    for ir in range(1, 8):
        data['fno_ucc'][ix][0][ir] = {'nfrozen': 50 - ir * 4}
        data['fno_ucc'][ix][0][ir]['uccsd_niter'] = 0
        data['fno_ucc'][ix][0][ir]['uccsdt_niter'] = 0
        data['fno_ucc'][ix][0][ir]['uccsdtq_niter'] = 0
        data['fno_ucc'][ix][0][ir]['uccsd_time'] = 0
        data['fno_ucc'][ix][0][ir]['uccsdt_time'] = 0
        data['fno_ucc'][ix][0][ir]['uccsdtq_time'] = 0
        with open('../06-FEMOCO-LLDUC/%d-FNO-UCC-%d/%02d-fno-f%d-%d.out' % (ix - 14, ix, ir, 50 - ir * 4, ix), "r") as f:
            for l in f.readlines():
                if l.startswith('EHF ='):
                    data['fno_ucc'][ix][0][ir]['uhf'] = float(l.split()[2])
                elif 'MP2 energy =' in l:
                    data['fno_ucc'][ix][0][ir]['ump2'] = float(l.split()[-1]) + data['fno_ucc'][ix][0][ir]['uhf']
                elif l.startswith('E(UCCSD) ='):
                    data['fno_ucc'][ix][0][ir]['uccsd'] = float(l.split()[2])
                elif l.startswith('E(UCCSDT) ='):
                    data['fno_ucc'][ix][0][ir]['uccsdt'] = float(l.split()[2])
                elif l.startswith('E(UCCSDTQ) ='):
                    data['fno_ucc'][ix][0][ir]['uccsdtq'] = float(l.split()[2])
                elif 'E_corr(UCCSD) =' in l and 'norm(d tamps)' in l:
                    data['fno_ucc'][ix][0][ir]['uccsd_niter'] += 1
                    data['fno_ucc'][ix][0][ir]['uccsd_time'] += float(l.split()[-1])
                elif 'E_corr(UCCSDT) =' in l and 'norm(d tamps)' in l:
                    data['fno_ucc'][ix][0][ir]['uccsdt_niter'] += 1
                    data['fno_ucc'][ix][0][ir]['uccsdt_time'] += float(l.split()[-1])
                elif 'E_corr(UCCSDTQ) =' in l and 'norm(d tamps)' in l:
                    data['fno_ucc'][ix][0][ir]['uccsdtq_niter'] += 1
                    data['fno_ucc'][ix][0][ir]['uccsdtq_time'] += float(l.split()[-1])
        assert 'uccsd' in data['fno_ucc'][ix][0][ir]
        assert 'uccsdt' in data['fno_ucc'][ix][0][ir]
    for ir, ifr, xre in zip([6, 7, 7], [8, 9, 10], ['-re', '-re', '-re2']):
        if not os.path.exists('../06-FEMOCO-LLDUC/%d-FNO-UCC-%d/%02d-fno-f%d%s-%d.out' % (ix - 14, ix, ifr, 50 - ir * 4, xre, ix)):
            continue
        with open('../06-FEMOCO-LLDUC/%d-FNO-UCC-%d/%02d-fno-f%d%s-%d.out' % (ix - 14, ix, ifr, 50 - ir * 4, xre, ix), "r") as f:
            for l in f.readlines():
                if l.startswith('E(UCCSDTQ) ='):
                    data['fno_ucc'][ix][0][ir]['uccsdtq'] = float(l.split()[2])
                elif 'E_corr(UCCSDTQ) =' in l and 'norm(d tamps)' in l:
                    data['fno_ucc'][ix][0][ir]['uccsdtq_niter'] += 1
                    data['fno_ucc'][ix][0][ir]['uccsdtq_time'] += float(l.split()[-1])
    for ir in range(1, 8):
        assert 'uccsdtq' in data['fno_ucc'][ix][0][ir]
    for ir, ifr in [[1, 11]]:
        data['fno_ucc'][ix][0][ir]['uccsdtqp_niter'] = 0
        data['fno_ucc'][ix][0][ir]['uccsdtqp_time'] = 0
        with open('../06-FEMOCO-LLDUC/%d-FNO-UCC-%d/%02d-fno-f%d-p-%d.out' % (ix - 14, ix, ifr, 50 - ir * 4, ix), "r") as f:
            for l in f.readlines():
                if l.startswith('E(UCCSDTQP) ='):
                    data['fno_ucc'][ix][0][ir]['uccsdtqp'] = float(l.split()[2])
                elif 'E_corr(UCCSDTQP) =' in l and 'norm(d tamps)' in l:
                    data['fno_ucc'][ix][0][ir]['uccsdtqp_niter'] += 1
                    data['fno_ucc'][ix][0][ir]['uccsdtqp_time'] += float(l.split()[-1])
        assert 'uccsdtqp' in data['fno_ucc'][ix][0][ir]
with open('../06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/08-ssq.out', "r") as f:
    for l in f.readlines():
        ir = int(l.split()[1])
        order = int(l.split()[4])
        ix = int(l.split()[7])
        ssq = float(l.split()[10])
        mult = float(l.split()[13])
        data['fno_ucc'][ix][0][ir]['ucc%s_ssq' % 'sdtqp'[:order]] = [ssq, mult]



with open("06-llduc-fno-ucc.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of FNO-UCC calculations: ', sum([len(data['fno_ucc'][ix]) for ix in [25, 26, 27]]))
