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


import json, numpy as np

data = {'large_act': {ix: {} for ix in range(5)}}

with open('../07-FEMOCO-Large/00-Active-Space/12-other-ints-solv-0.out', 'r') as f:
    for l in f.readlines():
        if l.startswith('Ncas ='):
            data['large_act'][4]['ncas'] = int(l.split()[2])
            data['large_act'][4]['nelec'] = int(l.split()[5])

with open('../07-FEMOCO-Large/00-Active-Space/05-fm-sel.out', 'r') as f:
    for l in f.readlines():
        if l.startswith('Selected ='):
            data['large_act'][4]['ext_nfrozen'] = int(l.split()[8])
            data['large_act'][4]['nfrozen'] = 0

with open('../07-FEMOCO-Large/00-Active-Space/14-ump2-0.out', 'r') as f:
    for l in f.readlines():
        if l.startswith('converged SCF energy ='):
            data['large_act'][4]['uhf_high_spin'] = float(l.split()[4])
            data['large_act'][4]['uhf_high_spin_ssq'] = [float(l.split()[7]), float(l.split()[10])]
        if l.startswith('E(UMP2) ='):
            data['large_act'][4]['ump2_high_spin'] = float(l.split()[2])

uno_occ = np.load('../07-FEMOCO-Large/00-Active-Space/17-uno-occ.npy')
data['large_act'][4]['uno_occ'] = uno_occ.tolist()
data['large_act'][4]['uno_llduc_overlap'] = np.load('../07-FEMOCO-Large/00-Active-Space/17-mo-ovlp.npy').tolist()
data['large_act'][4]['occ ~ 1'] = int(np.count_nonzero((0.1 < uno_occ) & (uno_occ < 1.8)))
data['large_act'][4]['occ > 1'] = int(np.count_nonzero(uno_occ > 1.8))
data['large_act'][4]['occ < 1'] = int(np.count_nonzero(uno_occ < 0.1))

for ix in range(4):
    with open('../07-FEMOCO-Large/00-Active-Space/19-act-high-spin-ump2-%d.out' % ix, 'r') as f:
        ls = list(f.readlines())
        data['large_act'][ix]['ncas'] = int(ls[1].split(', ')[1].split(')')[0])
        data['large_act'][ix]['nelec'] = int(ls[4].split()[2])
        data['large_act'][ix]['nfrozen'] = (data['large_act'][4]['nelec'] - data['large_act'][ix]['nelec']) // 2
        act_occ = uno_occ[data['large_act'][ix]['nfrozen']:data['large_act'][ix]['nfrozen'] + data['large_act'][ix]['ncas']]
        data['large_act'][ix]['occ ~ 1'] = int(np.count_nonzero((0.1 < act_occ) & (act_occ < 1.8)))
        data['large_act'][ix]['occ > 1'] = int(np.count_nonzero(act_occ > 1.8))
        data['large_act'][ix]['occ < 1'] = int(np.count_nonzero(act_occ < 0.1))
        for l in ls:
            if l.startswith('converged SCF energy ='):
                data['large_act'][ix]['uhf_high_spin'] = float(l.split()[4])
                data['large_act'][ix]['uhf_high_spin_ssq'] = [float(l.split()[7]), float(l.split()[10])]
            if l.startswith('E(UMP2) ='):
                data['large_act'][ix]['ump2_high_spin'] = float(l.split()[2])

with open("12-large-act.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
