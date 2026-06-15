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

data = json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/00-llduc-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))

bs_map = data['name_x']

e_uhf = 22140 + np.array([min([(x['uhf'], x['uhf_ssq'][0]) for x in data['uhf'][str(ix)]])[0] for ix in range(35)])
e_uhf_ssq = np.array([min([(x['uhf'], x['uhf_ssq'][0]) for x in data['uhf'][str(ix)]])[1] for ix in range(35)])

e_uccsd = 22140 + np.array([min([(x['uccsd'], x.get('uccsd_ssq', [np.nan, np.nan])[0]) for x in data['uccsd'][str(ix)]])[0] for ix in range(35)])
e_uccsd_ssq = np.array([min([(x['uccsd'], x.get('uccsd_ssq', [np.nan, np.nan])[0]) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])

e_uccsdt = 22140 + np.array([min([(x['uccsdt'], x.get('uccsdt_ssq', [np.nan, np.nan])[0]) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_uccsdt_ssq = np.array([min([(x['uccsdt'], x.get('uccsdt_ssq', [np.nan, np.nan])[0]) for x in data['uccsdt'][str(ix)]])[1] for ix in range(35)])

for ix in [25, 26, 27]:
    print('%10s & %10.6f & %8.2f & %10.6f & %8.2f & %10.6f & %8.2f \\\\' % (
        bs_map[str(ix)], e_uhf[ix], e_uhf_ssq[ix], e_uccsd[ix], e_uccsd_ssq[ix], e_uccsdt[ix], e_uccsdt_ssq[ix], 
    ))
