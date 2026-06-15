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

for key in ['uhf', 'uccsd', 'uccsd_t']:
    name = 'UCCSD(T)' if key == 'uccsd_t' else key.upper()
    print('%10s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\' % (
        name,
        data['heisenberg']['fit-' + key]['param']['jlt'],
        data['heisenberg']['fit-' + key]['param']['jlc'],
        data['heisenberg']['fit-' + key]['param']['jrc'],
        data['heisenberg']['fit-' + key]['param']['jrt'],
        data['heisenberg']['fit-' + key]['param']['jca'],
        data['heisenberg']['fit-' + key]['param']['jcc'],
        (data['heisenberg']['fit-' + key]['param']['e0'] + 22140) * 1000,
        data['heisenberg']['fit-' + key]['param']['rmse'] * 1000,
    ))