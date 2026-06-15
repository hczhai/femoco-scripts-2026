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

# data from: https://arxiv.org/pdf/2208.02199

data = {}

data['2feii-2feiii-energy'] = {
    '1000': -8432.789966,
    '2000': -8432.798444,
    '3000': -8432.801780,
    '4000': -8432.803590,
    '5000': -8432.804722,
    '6000': -8432.805500,
    '7000': -8432.806069,
    '8000': -8432.806502,
}

data['2feii-2feiii-dw'] = {
    '1000': 3.73E-4,
    '2000': 2.01E-4,
    '3000': 1.34E-4,
    '4000': 9.76E-5,
    '5000': 7.43E-5,
    '6000': 5.84E-5,
    '7000': 4.65E-5,
    '8000': 3.50E-5,
}

data['2feii-2feiii-energy-extra'] = -8432.8084
data['2feii-2feiii-energy-extra-error'] = 0.0004

data['4feiii-energy'] = {
    '1000': -8432.606085,
    '2000': -8432.615484,
    '3000': -8432.619229,
    '4000': -8432.621279,
    '5000': -8432.622578,
    '6000': -8432.623476,
    '7000': -8432.624094,
}

data['4feiii-dw'] = {
    '1000': 3.03E-4,
    '2000': 1.78E-4,
    '3000': 1.29E-4,
    '4000': 1.01E-4,
    '5000': 8.18E-5,
    '6000': 6.70E-5,
    '7000': 5.44E-5,
}

data['4feiii-energy-extra'] = -8432.6287
data['4feiii-energy-extra-error'] = 0.0009

data['psyn-energy'] = {
    '1000': -17833.977588,
    '1500': -17833.989684,
    '2000': -17833.996979,
    '2500': -17834.001975,
    '3000': -17834.005658,
    '3500': -17834.008503,
    '4000': -17834.010764,
    '4500': -17834.012610,
    '5000': -17834.014051,
}

data['psyn-dw'] = {
    '1000': 5.08E-4,
    '1500': 4.41E-4,
    '2000': 3.88E-4,
    '2500': 3.50E-4,
    '3000': 3.21E-4,
    '3500': 2.95E-4,
    '4000': 2.74E-4,
    '4500': 2.53E-4,
    '5000': 2.27E-4,
}

data['psyn-energy-extra'] = -17834.049
data['psyn-energy-extra-error'] = 0.007

data['pn-energy'] = {
    '1000': -17492.190311,
    '1500': -17492.200594,
    '2000': -17492.206725,
    '2500': -17492.210894,
    '3000': -17492.213953,
    '3500': -17492.216294,
    '4000': -17492.218146,
    '4500': -17492.219644,
    '5000': -17492.220847,
}

data['pn-dw'] = {
    '1000': 4.17E-4,
    '1500': 3.40E-4,
    '2000': 2.89E-4,
    '2500': 2.51E-4,
    '3000': 2.23E-4,
    '3500': 2.00E-4,
    '4000': 1.82E-4,
    '4500': 1.63E-4,
    '5000': 1.36E-4,
}

data['pn-energy-extra'] = -17492.239
data['pn-energy-extra-error'] = 0.004

data['pox-energy'] = {
    '500': -17489.880043,
    '1000': -17489.895112,
    '1500': -17489.901003,
    '2000': -17489.904229,
    '2500': -17489.906285,
    '3000': -17489.907715,
    '3500': -17489.908781,
    '4000': -17489.909590,
}

data['pox-dw'] = {
    ' 500': 4.06E-4,
    '1000': 2.60E-4,
    '1500': 1.94E-4,
    '2000': 1.57E-4,
    '2500': 1.34E-4,
    '3000': 1.17E-4,
    '3500': 1.03E-4,
    '4000': 8.70E-5,
}

data['pox-energy-extra'] = -17489.919
data['pox-energy-extra-error'] = 0.002

data['femoco-energy'] = {
    '2000': -22140.21120,
    '3000': -22140.24039,
    '4000': -22140.25652,
    '4500': -22140.26141,
    '5000': -22140.26562,
    '5500': -22140.26995,
    '6000': -22140.27335,
    '6500': -22140.27633,
    '7000': -22140.27881,
}

data['femoco-dw'] = {
    '2000': 4.53E-4,
    '3000': 3.77E-4,
    '4000': 3.19E-4,
    '4500': 3.02E-4,
    '5000': 2.86E-4,
    '5500': 2.69E-4,
    '6000': 2.60E-4,
    '6500': 2.55E-4,
    '7000': 2.47E-4,
}

data['femoco-energy-extra'] = -22140.36
data['femoco-energy-extra-error'] = 0.02

with open("71-dmrg-arxiv-2208-02199.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
