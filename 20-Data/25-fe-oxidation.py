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

data = {'fe-oxidation': {
    'fe2+': 24.83 * 24.97 / 24.94,
    'fe2.5+': 24.66 * 24.97 / 24.94,
    'fe3+': 24.53 * 24.97 / 24.94,
}}

with open("25-fe-oxidation.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
