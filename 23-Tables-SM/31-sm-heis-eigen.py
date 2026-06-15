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

for k, v in sorted(data['heisenberg']['fit-uccsd_t']['sa-dmrg']['rel_energy'].items(), key=lambda x: int(x[0])):
    print('%4d/2 & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\' % (
        int(k), *[x * 1000 for x in v],
    ))
