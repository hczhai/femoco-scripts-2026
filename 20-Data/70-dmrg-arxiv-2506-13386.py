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

# data from: https://arxiv.org/pdf/2506.13386

data = {}

data['2fe-2s-lmo'] = {
    '100': -116.602019,
    '500': -116.605538,
    '1000': -116.605602,
    '2000': -116.605609,
    '3000': -116.605609,
}

data['2fe-2s-emo'] = {
    '100': -116.603145,
    '500': -116.605528,
    '1000': -116.605599,
    '2000': -116.605608,
    '3000': -116.605609,
}

data['4fe-4s-lmo'] = {
    '100': -327.157422,
    '500': -327.214060,
    '1000': -327.226723,
    '2000': -327.234789,
    '3000': -327.237975,
    '4000': -327.239683,
    '5000': -327.240749,
    '6000': -327.241478,
    '7000': -327.242005,
    '8000': -327.242405,
}

data['4fe-4s-emo'] = {
    '100': -327.200321,
    '500': -327.227047,
    '1000': -327.234198,
    '2000': -327.238697,
    '3000': -327.240564,
    '4000': -327.241589,
    '5000': -327.242233,
    '6000': -327.242679,
    '7000': -327.243007,
    '8000': -327.243261,
}

data['p-sluter-lmo'] = {
    '100': -1075.372662,
    '500': -1075.459722,
    '1000': -1075.483101,
    '2000': -1075.499669,
    '3000': -1075.507308,
    '4000': -1075.511586,
    '5000': -1075.514333,
    '6000': -1075.516286,
    '7000': -1075.517731,
    '8000': -1075.518856,
}

data['p-sluter-emo'] = {
    '100': -1075.452082,
    '500': -1075.488247,
    '1000': -1075.498670,
    '2000': -1075.507117,
    '3000': -1075.511554,
    '4000': -1075.514302,
    '5000': -1075.516178,
    '6000': -1075.517544,
    '7000': -1075.518598,
    '8000': -1075.519439,
}

data['femoco-lmo'] = {
    '100': -1118.482507,
    '500': -1118.820009,
    '1000': -1118.920234,
    '2000': -1118.996538,
    '3000': -1119.033406,
    '4000': -1119.056419,
    '5000': -1119.072124,
    '6000': -1119.083773,
    '7000': -1119.092880,
    '8000': -1119.100234,
    '9000': -1119.106330,
    '10000': -1119.111498,
}

data['femoco-emo'] = {
    '100': -1118.868657,
    '500': -1119.005819,
    '1000': -1119.050973,
    '2000': -1119.087453,
    '3000': -1119.105478,
    '4000': -1119.116896,
    '5000': -1119.124937,
    '6000': -1119.131007,
    '7000': -1119.135803,
    '8000': -1119.139701,
    '9000': -1119.142944,
    '10000': -1119.145684,
}

for k in data['p-sluter-lmo']:
    data['p-sluter-lmo'][k] += -16416.706473
for k in data['p-sluter-emo']:
    data['p-sluter-emo'][k] += -16416.706473

for k in data['femoco-lmo']:
    data['femoco-lmo'][k] += -21021.214012
for k in data['femoco-emo']:
    data['femoco-emo'][k] += -21021.214012

with open("70-dmrg-arxiv-2506-13386.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
