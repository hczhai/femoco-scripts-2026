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
data |= json.load(open("../20-Data/03-llduc-udmrg-5k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/05-llduc-udmrg-8k.json", "r", encoding='utf-8'))

names = [data['name_x'][str(ix)] for ix in range(35)]

e_uhf = 22140 + np.array([min([(x['uhf'], x['uhf_conv']) for x in data['uhf'][str(ix)]])[0] for ix in range(35)])
e_uccsd = 22140 + np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[0] for ix in range(35)])
e_uccsd_t = 22140 + np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])
e_uccsdt = 22140 + np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_dmrg_5k = 22140 + np.array([data['udmrg_5k'][str(ix)][0]['udmrg'] for ix in range(35)])
e_dmrg_8k = 22140 + np.array([data['udmrg_8k'][str(ix)][0]['udmrg'] if len(data['udmrg_8k'][str(ix)]) != 0 else np.nan for ix in range(35)])
t1diag_uccsd = np.array([min([(x['uccsd'], x['t1-diagnostic']) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])

for ix in range(35):
    print('%10s & %10.6f & %10.6f & %10.6f & %10.6f & %10.6f & %10s & %8.5f \\\\' % (
        data['name_x'][str(ix)],
        e_uhf[ix], e_uccsd[ix], e_uccsd_t[ix], e_uccsdt[ix], e_dmrg_5k[ix],
        '%10.6f' % e_dmrg_8k[ix] if not np.isnan(e_dmrg_8k[ix]) else '', t1diag_uccsd[ix]
    ))
