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

data = json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/14-large-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/23-mofe3.json", "r", encoding='utf-8'))

print('\\multicolumn{2}{c|}{\\multirow{3}{*}{$\\ce{[MoFe3S4(SH)6]^{3-}}$}}')
for ik, key in enumerate(['uks', 'uhf', 'uccsd']):
    if ik != 0:
       print('\\multicolumn{2}{c|}{}', end=' ')
    print('& %5s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2s & %8.2s & %8.2s & %8.2s \\\\' % (key.upper(),
        *data['mofe3']['0'][key + '-elec-pop-noryd'], '', '', '', ''))

print("    \\hline")
print("    \\multirow{3}{*}{FeMo-co}")

e_uccsd_pop = [min([(x['uccsd'], x.get('uccsd-elec-pop-noryd', [0] * 8)) for x in data['large_uccsd']['3'][str(ix)]])[1] for ix in [25, 26, 27]]

for ix in [25, 26, 27]:
    print('& %6s & %5s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\' % (
        data['name_x'][str(ix)],
        '\\multirow{3}{*}{UCCSD}' if ix == 25 else '',
        e_uccsd_pop[ix - 25][-1], *[e_uccsd_pop[ix - 25][x - 1] for x in [1, 4, 3, 2, 5, 7, 6]]))
