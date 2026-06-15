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


import numpy as np, json

data = json.load(open("../20-Data/24-heisenberg.json", "r", encoding='utf-8'))

for sz in [1, 3]:
    for ir in range(3):
        print('(3/2,%d/2) & %8.3f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\' % (
            sz,
            data['heisenberg']['fit-uccsd_t']['nsa-dmrg']['rel_energy'][str(sz)][ir] * 1000,
            *[data['heisenberg']['fit-uccsd_t']['nsa-dmrg']['sz'][str(sz)][ir][x] for x in [0, 3, 2, 1, 4, 6, 5, 7]]
        ))
    if sz == 1:
        print("\\hline")
