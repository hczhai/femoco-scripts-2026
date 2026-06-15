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

data = json.load(open("../20-Data/10-fe2s2.json", "r", encoding='utf-8'))

ref = data['fe2s2']['std']

for key in ['uhf', 'uccsd', 'uccsd(t)', 'uccsdt', 'uccsdtq', 'uccsdtqp']:
    print('%16s & %12.6f & %10.6f & %10.6f & %10.6f & %10.6f \\\\' % (
        key.upper(), ref[key], ref[key + '_ssq'][0], ref[key + '_ssq'][1] - 1,
        ref[key + '_sz'][0], ref[key + '_sz'][1]
    ))

print('%16s & %12.6f & %10.6f & %10.6f & %10s & %10s \\\\' % (
    'SA-DMRG[$D$=12k]', ref['sa_dmrg_12k'], 0.0, 0.0, '', ''
    ))
