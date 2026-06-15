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

data = {'qmmm_ccsd': {ix: [] for ix in range(35)} | {-1: []}}

for ix in range(35):
    with open("../08-FEMOCO-QMMM/01-UHF-UCCSD/05-fm-ccsd-%d.out" % (ix, ), "r") as f:
        uhf_conv, xtt = True, None
        for l in f.readlines():
            if '::' in l and int(l.strip().split()[-1]) <= 1:
                ll = l.strip().split()
                data['qmmm_ccsd'][int(ll[0])].append({
                    'init': ll[3],
                    'proj': xtt,
                    'uccsd': float(ll[4]),
                    'uccsd(t)': float(ll[5]),
                    'uccsd_conv': int(ll[7]) == 1,
                    't1-diagnostic': float(ll[6]),
                    'uhf': uhf,
                    'uhf_conv': uhf_conv,
                })
                uhf_conv = True
            elif l.startswith('--- '):
                xtt = l.split()[-2]
            elif l.startswith('SCF not converged'):
                uhf_conv = False
            elif 'SCF energy =' in l:
                uhf = float(l.split('SCF energy =')[1].split()[0])

for ik, key in enumerate(['uhf-elec-pop-noryd', 'uhf-spin', 'uhf-nuc-density',
    'uccsd-elec-pop-noryd', 'uccsd-spin', 'uccsd-nuc-density']):
    with open("../08-FEMOCO-QMMM/01-UHF-UCCSD/%d-%s.out" % (ik + 10, key), "r") as f:
        for l in f.readlines():
            ll = l.strip().split()
            ix = int(ll[0])
            idx = [x['init'] for x in data['qmmm_ccsd'][ix]].index(ll[3])
            data['qmmm_ccsd'][ix][idx][key] = [float(x) for x in ll[4:4 + 8]]

with open("16-qmmm-uccsd.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('QMMM: Number of UCCSD calculations: ', sum([len(data['qmmm_ccsd'][ix]) for ix in range(35)]))
