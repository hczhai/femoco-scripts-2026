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
data |= json.load(open("../20-Data/15-qmmm-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/16-qmmm-uccsd.json", "r", encoding='utf-8'))

e_uhf = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['qmmm_uhf'][str(ix)]])[0] for ix in range(35)])
e_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['qmmm_ccsd'][str(ix)]])[0] for ix in range(35)])
e_uccsd_t = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['qmmm_ccsd'][str(ix)]])[1] for ix in range(35)])
e_uccsd_conv = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['qmmm_ccsd'][str(ix)]])[2] for ix in range(35)])
t1diag_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x['t1-diagnostic']) for x in data['qmmm_ccsd'][str(ix)]])[2] for ix in range(35)])

for ix in range(35):
    print('%10s & %16.6f & %16.6f & %16.6f & %12.5f \\\\' % (
        data['name_x'][str(ix)], e_uhf[ix], e_uccsd[ix], e_uccsd_t[ix], t1diag_uccsd[ix]
    ))
