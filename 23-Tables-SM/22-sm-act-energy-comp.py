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
data |= json.load(open("../20-Data/17-large-ump2.json", "r", encoding='utf-8'))

names = [data['name_x'][str(ix)] for ix in range(35)]
inits = [ix for ix in range(35) if int(names[ix].split('-')[0].split('BS')[1]) >= 7]

ia = 3
e_uhf_03 = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['large_uhf'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_03 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[0] for ix in inits])
e_uccsd_t_03 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[1] for ix in inits])
e_uccsd_conv_03 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[2] for ix in inits])
assert all(e_uccsd_conv_03.tolist())

ump2_map = {(x['init'], x['proj']): x for ix in inits for x in data['large_ump2'][str(ix)]}
e_init_03 = [min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv'], x['init'], x['proj']) for x in data['large_uccsd'][str(ia)][str(ix)]])[3:] for ix in inits]
e_uhf_03 = [ump2_map[k]['uhf_small'] for k in e_init_03]
e_ump2_03 = [ump2_map[k]['ump2_small'] for k in e_init_03]
e_uhf_04 = [ump2_map[k]['uhf_large'] for k in e_init_03]
e_ump2_04 = [ump2_map[k]['ump2_large'] for k in e_init_03]

for ii, ix in enumerate(inits):
    print('%10s & %16.6f & %16.6f & %16.6f & %12.6f & %16.6f & %16.6f \\\\' % (
        data['name_x'][str(ix)],
        e_uhf_04[ii], e_ump2_04[ii], e_ump2_03[ii], e_ump2_04[ii] - e_ump2_03[ii],
        e_uccsd_03[ii] + e_ump2_04[ii] - e_ump2_03[ii],
        e_uccsd_t_03[ii] + e_ump2_04[ii] - e_ump2_03[ii],
    ))
