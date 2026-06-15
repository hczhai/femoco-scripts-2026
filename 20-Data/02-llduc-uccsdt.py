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


import json, os, numpy as np

data = {'uccsdt': {ix: [] for ix in range(35)}}

for ix in range(35):
    fid, nref = (0, 24) if ix in [25, 26, 27] else (10, 12)
    for ir in range(nref):
        d = {}
        with open("../06-FEMOCO-LLDUC/06-UCCSDT/%02d-fm-ccsdt-%d-%d.out" % (fid, ix, ir), "r") as f:
            for l in f.readlines():
                if l.startswith('INIT ='):
                    d['init'] = l.split()[2]
                elif l.startswith('EHF  ='):
                    d['uhf'] = float(l.split()[2])
                elif l.startswith('E(UCCSD) ='):
                    d['uccsd'] = float(l.split()[2])
                elif l.startswith('E(UCCSDT) ='):
                    d['uccsdt'] = float(l.split()[2])
                elif l.startswith('UCCSD not converged'):
                    d['uccsd_conv'] = False
                elif l.startswith('UCCSDT not converged'):
                    d['uccsdt_conv'] = False
                elif l.startswith('UCCSD converged'):
                    d['uccsd_conv'] = True
                elif l.startswith('UCCSDT converged'):
                    d['uccsdt_conv'] = True
        assert 'init' in d and 'uhf' in d and 'uccsd' in d and 'uccsd_conv' in d
        if 'uccsdt' not in d:
            d['uccsdt_conv'] = False
            with open("../06-FEMOCO-LLDUC/06-UCCSDT/%02d-fm-ccsdt-%d-%d.out" % (fid, ix, ir), "r") as f:
                for l in f.readlines():
                    if 'E_corr(UCCSDT) =' in l:
                        d['uccsdt'] = float(l.split('E_corr(UCCSDT) =')[1].split()[0]) + d['uhf']
        assert 'uccsdt' in d and 'uccsdt_conv' in d
        data['uccsdt'][ix].append(d)
    if os.path.exists("../06-FEMOCO-LLDUC/06-UCCSDT/11-fm-ccsdt-hard-%d.out" % (ix, )):
        d = {}
        with open("../06-FEMOCO-LLDUC/06-UCCSDT/11-fm-ccsdt-hard-%d.out" % (ix, ), "r") as f:
            for l in f.readlines():
                if l.startswith('INIT ='):
                    d['init'] = l.split()[2]
                elif l.startswith('EHF  ='):
                    d['uhf'] = float(l.split()[2])
                elif l.startswith('E(UCCSD) ='):
                    d['uccsd'] = float(l.split()[2])
                elif l.startswith('E(UCCSDT) ='):
                    d['uccsdt'] = float(l.split()[2])
                elif l.startswith('E(UCCSD) ='):
                    d['uccsd'] = float(l.split()[2])
                elif l.startswith('UCCSD not converged'):
                    d['uccsd_conv'] = False
                elif l.startswith('UCCSDT not converged'):
                    d['uccsdt_conv'] = False
                elif l.startswith('UCCSD converged'):
                    d['uccsd_conv'] = True
                elif l.startswith('UCCSDT converged'):
                    d['uccsdt_conv'] = True
        assert 'init' in d and 'uhf' in d and 'uccsd' in d and 'uccsd_conv' in d
        assert 'uccsdt' in d and 'uccsdt_conv' in d
        if d['init'] in [x['init'] for x in data['uccsdt'][ix]]:
            idx = [x['init'] for x in data['uccsdt'][ix]].index(d['init'])
            data['uccsdt'][ix][idx] = d
        else:
            data['uccsdt'][ix].append(d)

for ix in [25, 26, 27]:
    with open("../06-FEMOCO-LLDUC/16-UCC-NPDM/30-fm-cc-new-%d.out" % (ix, ), "r") as f:
        xinit = None
        for l in f.readlines():
            if l.startswith('INIT ='):
                xinit = l.split()[2]
            elif '<S^2> =' in l:
                idx = [x['init'] for x in data['uccsdt'][ix]].index(xinit)
                data['uccsdt'][ix][idx]['uhf_ssq'] = [float(l.split()[-4]), float(l.split()[-1])]
    with open("../06-FEMOCO-LLDUC/16-UCC-NPDM/31-ssq.out", "r") as f:
        for l in f.readlines():
            if int(l.split()[5]) == ix and int(l.split()[2]) == 3:
                idx = [x['init'] for x in data['uccsdt'][ix]].index(xinit)
                data['uccsdt'][ix][idx]['uccsdt_ssq'] = [float(l.split()[-4]), float(l.split()[-1])]
    dm1 = np.load('../06-FEMOCO-LLDUC/16-UCC-NPDM/data/30-cc-03-%02d-pdm1.npy' % ix)
    mo = np.load('../06-FEMOCO-LLDUC/16-UCC-NPDM/data/30-mo-coeff-%02d.npy' % ix)
    sfdm = mo[0] @ dm1[0] @ mo[0].T + mo[1] @ dm1[1] @ mo[1].T
    idx = [x['init'] for x in data['uccsdt'][ix]].index(xinit)
    data['uccsdt'][ix][idx]['uccsdt_nat-occ'] = np.linalg.eigh(sfdm)[0][::-1].tolist()

with open("02-llduc-uccsdt.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of UCCSDT calculations: ', sum([len(data['uccsdt'][ix]) for ix in range(35)]))
