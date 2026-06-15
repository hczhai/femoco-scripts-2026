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

data = {'udmrg_5k_ref': {ix: [] for ix in range(35)}}

for ix in [26]:
    for jx in range(12):
        d = {}
        for xre in ['', '-re']:
            with open("../06-FEMOCO-LLDUC/08-UDMRG-%d-5k/01-fm-dmrg-ucc-var%s-%d.out" % (ix, xre, jx), "r") as f:
                for l in f.readlines():
                    if l.startswith('xdstr ='):
                        d['init'] = l.split()[2]
                    elif l.startswith('STATE 0 DMRG energy ='):
                        d['udmrg'] = float(l.split()[5])
        assert 'init' in d and 'udmrg' in d
        data['udmrg_5k_ref'][ix].append(d)
    for jx in range(12, 24):
        d = {}
        with open("../06-FEMOCO-LLDUC/08-UDMRG-%d-5k/03-fm-dmrg-ucc-other-%d.out" % (ix, jx), "r") as f:
            for l in f.readlines():
                if l.startswith('xdstr ='):
                    d['init'] = l.split()[2]
                elif l.startswith('STATE 0 DMRG energy ='):
                    d['udmrg'] = float(l.split()[5])
        assert 'init' in d and 'udmrg' in d
        data['udmrg_5k_ref'][ix].append(d)

    for ik, key in enumerate(['udmrg-elec-pop-noryd', 'udmrg-spin', 'udmrg-nuc-density']):
        with open("../06-FEMOCO-LLDUC/08-UDMRG-%d-5k/%d-%s.out" % (ix, 20 + ik, key), "r") as f:
            for l in f.readlines():
                ll = l.strip().split()
                idx = [x['init'] for x in data['udmrg_5k_ref'][ix]].index(ll[3])
                data['udmrg_5k_ref'][ix][idx][key] = [float(x) for x in ll[4:4 + 8]]

with open("09-llduc-udmrg-5k-ref.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of UDMRG (D=5k) references calculations: ', sum([len(data['udmrg_5k_ref'][ix]) for ix in range(35)]))
