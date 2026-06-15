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

data = json.load(open("../20-Data/07-llduc-ghf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/08-llduc-gccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))

gid_map = {x['init']: (x['id'], x['ghf']) for x in data['ghf']['unique']}
e_ccsd = sorted([x['gccsd'], x['gccsd(t)'], x['init']] for x in data['gccsd'])[:16]
ghf_ids = [gid_map.get(x[2], '?') for x in e_ccsd]

notes = {}
notes[61] = notes[21] = notes[44] = notes[153] = notes[135] = '%s (%s)' % (data['name_x']['25'], data['name_ijk']['25'].split('-')[1])
notes[62] = '%s (%s)' % (data['name_x']['26'], data['name_ijk']['26'].split('-')[1])
notes[66] = '%s (%s)' % (data['name_x']['27'], data['name_ijk']['27'].split('-')[1])

for (xid, ghf), (gccsd, gccsd_t, _) in zip(ghf_ids, e_ccsd):
    print('%10s & %15.6f & %15.6f & %15.6f & %12s \\\\' % (
        'GHF-%d' % xid, ghf, gccsd, gccsd_t, notes.get(xid, '')))

print('\\hline')

uccsd = [min([[x['uccsd'], x['uccsd(t)'], x['uhf']] for x in data['uccsd'][str(ix)]]) for ix in range(35)]

for ix in [25, 26, 27]:
    print('%10s & %15.6f & %15.6f & %15.6f & \\\\' % (
        data['name_x'][str(ix)], uccsd[ix][2], uccsd[ix][0], uccsd[ix][1]
    ))
