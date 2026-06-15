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

# data from: Int J Quantum Chem. 2018;e25627. https://doi.org/10.1002/qua.25627
# digitized from Fig. 3C and 3D

data = {}

data['TPSS-D3/def2-TZVPD'] = [
    92.1, 80.3, 64.9, 49.8, 47.9,
    82.0, 52.5, 62.1, 43.6, 48.1,
    45.0, 60.3, 26.0, 26.1, 27.5,
    28.7, 38.4, 37.3, 0.0, 38.3,
    2.6, 37.4, 67.8, 66.3, 43.3,
    0.0, 38.0, 37.1, 65.6, 36.8,
    63.1, 38.3, 70.2, 71.7, 137.0
]

data['B3LYP/def2-TZVPD'] = [
    124.3, 134.8, 51.1, 20.4, 15.6,
    133.8, 19.7, 64.7, 13.1, 12.6,
    11.9, 56.3, 49.1, 48.4, 51.3,
    65.8, 25.0, 28.6, 0.0, 26.5,
    4.7, 36.8, 96.7, 126.4, 99.4,
    5.4, 20.0, 27.4, 95.4, 94.1,
    118.1, 91.6, 110.8, 118.0, 268.0
]

with open("72-dft-10.1002-qua.25627.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
