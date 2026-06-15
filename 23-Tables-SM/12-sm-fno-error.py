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
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/11-llduc-fno-udmrg.json", "r", encoding='utf-8'))

bs_map = data['name_x']

for ii, ix in enumerate([25, 26, 27]):
    for ir in range(1, 8):
        if ir < 7:
            print('%10s & %10s & %10.6f & %10.6f & %10.6f & %10.6f & %10.6f & %10.6f & %+10.6f \\\\' % (
                bs_map[str(ix)] if ir == 1 else '',
                '(%2do, %2de)' % (26 + ir * 4, 13 + ir * 8),
                data['fno_ucc'][str(ix)][0][str(ir)]['uhf'] + 22140,
                data['fno_ucc'][str(ix)][0][str(ir)]['uccsd'] + 22140,
                data['fno_ucc'][str(ix)][0][str(ir)]['uccsdt'] + 22140,
                data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq'] + 22140,
                data['fno_udmrg'][str(ix)][0][str(ir)]['udmrg_8000'] + 22140,
                data['fno_udmrg'][str(ix)][0][str(ir)]['udmrg_fit_logm'] + 22140,
                data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq'] - data['fno_udmrg'][str(ix)][0][str(ir)]['udmrg_fit_logm'],
            ))
        else:
            print('%10s & %10s & %10.6f & %10.6f & %10.6f & %10.6f & %10s & %10s & %10s \\\\' % (
                bs_map[str(ix)] if ir == 1 else '',
                '(%2do, %2de)' % (26 + ir * 4, 13 + ir * 8),
                data['fno_ucc'][str(ix)][0][str(ir)]['uhf'] + 22140,
                data['fno_ucc'][str(ix)][0][str(ir)]['uccsd'] + 22140,
                data['fno_ucc'][str(ix)][0][str(ir)]['uccsdt'] + 22140,
                data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq'] + 22140,
                '', '', ''
            ))
    if ii != 2:
        print("     \\hline")
