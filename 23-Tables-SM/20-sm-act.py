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

data = json.load(open("../20-Data/12-large-act.json", "r", encoding='utf-8'))

for ix in range(5):
    print('    (%3do, %3de) & %4d & %4d  & %4d  & %4d & %16.6f & %16.6f \\\\' % (
        data['large_act'][str(ix)]['ncas'],
        data['large_act'][str(ix)]['nelec'],
        data['large_act'][str(ix)]['occ > 1'],
        data['large_act'][str(ix)]['occ ~ 1'],
        data['large_act'][str(ix)]['occ < 1'],
        data['large_act'][str(ix)]['nfrozen'],
        data['large_act'][str(ix)]['uhf_high_spin'],
        data['large_act'][str(ix)]['ump2_high_spin'],
    ))
    if ix == 3:
        print('    \\hline')
