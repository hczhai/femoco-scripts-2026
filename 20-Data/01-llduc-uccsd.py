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

data = {'uccsd': {ix: [] for ix in range(35)}}

for ix in range(35):
    with open("../06-FEMOCO-LLDUC/01-UCCSD/04-fm-ccsd-%d.out" % (ix, ), "r") as f:
        uhf, uhf_conv = 0.0, True
        for l in f.readlines():
            if '::' in l and int(l.strip().split()[-1]) <= 1:
                ll = l.strip().split()
                data['uccsd'][ix].append({
                    'init': ll[2],
                    'uccsd': float(ll[3]),
                    'uccsd(t)': float(ll[4]),
                    'uccsd_conv': int(ll[5]) == 1,
                    'uhf': uhf,
                    'uhf_conv': uhf_conv,
                })
                uhf, uhf_conv = 0.0, True
            elif l.startswith('SCF not converged'):
                uhf_conv = False
            elif 'SCF energy =' in l:
                uhf = float(l.split('SCF energy =')[1].split()[0])

ix = 31
with open("../06-FEMOCO-LLDUC/01-UCCSD/04-fm-ccsd-re-%d.out" % (ix, ), "r") as f:
    uhf, uhf_conv = 0.0, True
    for l in f.readlines():
        if '::' in l and int(l.strip().split()[-1]) <= 1:
            ll = l.strip().split()
            idx = [x['init'] for x in data['uccsd'][ix]].index(ll[2])
            data['uccsd'][ix][idx] = {
                'init': ll[2],
                'uccsd': float(ll[3]),
                'uccsd(t)': float(ll[4]),
                'uccsd_conv': int(ll[5]) == 1,
                'uhf': uhf,
                'uhf_conv': uhf_conv,
            }
        elif l.startswith('SCF not converged'):
            uhf_conv = False
        elif 'SCF energy =' in l:
            uhf = float(l.split('SCF energy =')[1].split()[0])

for ix in range(35):
    with open("../06-FEMOCO-LLDUC/07-UCCSD-NPDM/06-fm-ccsd-%d.out" % (ix, ), "r") as f:
        for l in f.readlines():
            if '::' in l and int(l.strip().split()[-1]) <= 1:
                ll = l.strip().split()
                idx = [x['init'] for x in data['uccsd'][ix]].index(ll[2])
                data['uccsd'][ix][idx]['t1-diagnostic'] = float(ll[-2])

for ik, key in enumerate(['uhf-elec-pop-noryd', 'uhf-spin', 'uhf-nuc-density',
    'uccsd-elec-pop-noryd', 'uccsd-spin', 'uccsd-nuc-density']):
    with open("../06-FEMOCO-LLDUC/07-UCCSD-NPDM/%d-%s.out" % (20 + ik, key), "r") as f:
        for l in f.readlines():
            ll = l.strip().split()
            ix = int(ll[0])
            idx = [x['init'] for x in data['uccsd'][ix]].index(ll[3])
            data['uccsd'][ix][idx][key] = [float(x) for x in ll[4:4 + 8]]

with open("../06-FEMOCO-LLDUC/07-UCCSD-NPDM/27-uhf-mo-pop.out", "r") as f:
    for l in f.readlines():
        ll = l.strip().split()
        ix = int(ll[0])
        idx = [x['init'] for x in data['uccsd'][ix]].index(ll[3])
        data['uccsd'][ix][idx]['uhf-mo-pop'] = [float(x) for x in ll[4:4 + 10]]

with open("../06-FEMOCO-LLDUC/07-UCCSD-NPDM/26-uccsd-mo-pop.out", "r") as f:
    for l in f.readlines():
        ll = l.strip().split()
        ix = int(ll[0])
        idx = [x['init'] for x in data['uccsd'][ix]].index(ll[3])
        data['uccsd'][ix][idx]['uccsd-mo-pop'] = [float(x) for x in ll[4:4 + 10]]

for ix in [25, 26, 27]:
    with open("../06-FEMOCO-LLDUC/16-UCC-NPDM/30-fm-cc-new-%d.out" % (ix, ), "r") as f:
        xinit = None
        for l in f.readlines():
            if l.startswith('INIT ='):
                xinit = l.split()[2]
            elif '<S^2> =' in l:
                idx = [x['init'] for x in data['uccsd'][ix]].index(xinit)
                data['uccsd'][ix][idx]['uhf_ssq'] = [float(l.split()[-4]), float(l.split()[-1])]
    with open("../06-FEMOCO-LLDUC/16-UCC-NPDM/31-ssq.out", "r") as f:
        for l in f.readlines():
            if int(l.split()[5]) == ix and int(l.split()[2]) == 2:
                idx = [x['init'] for x in data['uccsd'][ix]].index(xinit)
                data['uccsd'][ix][idx]['uccsd_ssq'] = [float(l.split()[-4]), float(l.split()[-1])]
    dm1 = np.load('../06-FEMOCO-LLDUC/16-UCC-NPDM/data/30-cc-02-%02d-pdm1.npy' % ix)
    mo = np.load('../06-FEMOCO-LLDUC/16-UCC-NPDM/data/30-mo-coeff-%02d.npy' % ix)
    sfdm = mo[0] @ dm1[0] @ mo[0].T + mo[1] @ dm1[1] @ mo[1].T
    idx = [x['init'] for x in data['uccsd'][ix]].index(xinit)
    data['uccsd'][ix][idx]['uccsd_nat-occ'] = np.linalg.eigh(sfdm)[0][::-1].tolist()
    with open("../06-FEMOCO-LLDUC/16-UCC-NPDM/32-fm-ccsd-new-%d.out" % (ix, ), "r") as f:
        xinit = None
        for l in f.readlines():
            if l.startswith('INIT ='):
                xinit = l.split()[2]
            elif '<S^2> =' in l:
                idx = [x['init'] for x in data['uccsd'][ix]].index(xinit)
                data['uccsd'][ix][idx]['uhf_ssq'] = [float(l.split()[-4]), float(l.split()[-1])]
    with open("../06-FEMOCO-LLDUC/16-UCC-NPDM/33-ssq.out", "r") as f:
        for l in f.readlines():
            if int(l.split()[5]) == ix and int(l.split()[2]) == 2:
                idx = [x['init'] for x in data['uccsd'][ix]].index(xinit)
                data['uccsd'][ix][idx]['uccsd_ssq'] = [float(l.split()[-4]), float(l.split()[-1])]
    dm1 = np.load('../06-FEMOCO-LLDUC/16-UCC-NPDM/data/32-cc-02-%02d-pdm1.npy' % ix)
    mo = np.load('../06-FEMOCO-LLDUC/16-UCC-NPDM/data/32-mo-coeff-%02d.npy' % ix)
    sfdm = mo[0] @ dm1[0] @ mo[0].T + mo[1] @ dm1[1] @ mo[1].T
    idx = [x['init'] for x in data['uccsd'][ix]].index(xinit)
    data['uccsd'][ix][idx]['uccsd_nat-occ'] = np.linalg.eigh(sfdm)[0][::-1].tolist()

with open("01-llduc-uccsd.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of UCCSD calculations: ', sum([len(data['uccsd'][ix]) for ix in range(35)]))
