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

data = {'udmrg_8k': {ix: [] for ix in range(35)}}

for ix in [4, 25, 26, 27]:
    d = {}
    with open("../06-FEMOCO-LLDUC/09-UDMRG-8k/21-fm-dmrg-uccsd-8k-%d.out" % (ix, ), "r") as f:
        for l in f.readlines():
            if l.startswith('xdstr ='):
                d['init'] = l.split()[2]
            elif l.startswith('STATE 0 DMRG energy ='):
                d['udmrg'] = float(l.split()[5])
        assert 'init' in d and 'udmrg' in d
        data['udmrg_8k'][ix].append(d)

with open("05-llduc-udmrg-8k.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of UDMRG (D=8k) calculations: ', sum([len(data['udmrg_8k'][ix]) for ix in range(35)]))
