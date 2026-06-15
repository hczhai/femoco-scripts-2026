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
data |= json.load(open("../20-Data/13-large-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/14-large-uccsd.json", "r", encoding='utf-8'))

names = [data['name_x'][str(ix)] for ix in range(35)]
inits = [ix for ix in range(35) if int(names[ix].split('-')[0].split('BS')[1]) >= 7]

ia = 0
e_uhf_00 = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['large_uhf'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_00 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_t_00 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[1] for ix in inits])
e_uccsd_conv_00 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[2] for ix in inits])
assert all(e_uccsd_conv_00.tolist())

ia = 1
e_uhf_01 = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['large_uhf'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_01 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_t_01 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[1] for ix in inits])
e_uccsd_conv_01 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[2] for ix in inits])
assert all(e_uccsd_conv_01.tolist())

ia = 2
e_uhf_02 = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['large_uhf'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_02 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_t_02 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[1] for ix in inits])
e_uccsd_conv_02 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[2] for ix in inits])
assert all(e_uccsd_conv_02.tolist())

ia = 3
e_uhf_03 = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['large_uhf'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_03 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_t_03 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[1] for ix in inits])
e_uccsd_conv_03 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[2] for ix in inits])
assert all(e_uccsd_conv_03.tolist())

for ii, ix in enumerate(inits):
    print('%10s & %16.6f & %16.6f & %16.6f & %16.6f & %16.6f & %16.6f \\\\' % (
        data['name_x'][str(ix)],
        e_uhf_02[ii], e_uccsd_02[ii], e_uccsd_t_02[ii], e_uhf_03[ii], e_uccsd_03[ii], e_uccsd_t_03[ii]
    ))
