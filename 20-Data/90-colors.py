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

data = {
    "colors": [
        "#2F2D2E",
        "#BB4430",
        "#F18F01",
        "#99C24D",
        "#1E77CA",
        "#23C9FF",
        "#7137E5",
        "#A347BA",
        "#84602E",
        "#430B9D",
        "#18AD8A"
    ],
    "grays": [
        "#000000",
        "#3D3D3D",
        "#767676",
        "#C5C5C5",
    ],
    "blues": [
        "#1E77CA",
        "#A9D2F8",
    ],
    "lights": [
        "#FCDAAD",
        "#E3FCAC",
        "#A9D2F8",
        "#CFACFD",
    ]
}

with open("90-colors.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
