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


import numpy as np, json, os

data = {'large_uccsd': {ia: {ix: [] for ix in range(35)} | {-1: []} for ia in range(4)}}

for ia in range(2):
    for ix in range(35):
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%02d-uccsd-act-%d.out" % (ia + 2, ia, ia + 4, ix), "r") as f:
            uhf_conv = True
            for l in f.readlines():
                if '::' in l and int(l.strip().split()[-2]) <= 1:
                    ll = l.strip().split()
                    data['large_uccsd'][ia][int(ll[0])].append({
                        'init': ll[3],
                        'proj': ll[8].strip(),
                        'uhf': float(ll[4]),
                        'uccsd': float(ll[5]),
                        'uccsd(t)': float(ll[6]),
                        'uccsd_conv': int(ll[7]) == 1,
                        'uhf_conv': uhf_conv,
                    })
                    uhf_conv = True
                elif l.startswith('SCF not converged'):
                    uhf_conv = False

for ia in [2]:
    for ix in range(35):
        bs_map = {(l.split()[3], l.split()[-3]): int(l.split('BS')[1].strip()) for l in open("../07-FEMOCO-Large/01-UHF/32-uhf-spin-%d.out" % (ix, ), "r").readlines() if '::' in l}
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/02-uccsd-act-%d.out" % (ia + 2, ia, ix), "r") as f:
            uhf_conv = True
            for l in f.readlines():
                if '::' in l and int(l.strip().split()[-2]) <= 1:
                    ll = l.strip().split()
                    data['large_uccsd'][ia][bs_map[(ll[2], ll[7].strip())]].append({
                        'init': ll[2],
                        'proj': ll[7].strip(),
                        'uhf': float(ll[3]),
                        'uccsd': float(ll[4]),
                        'uccsd(t)': float(ll[5]),
                        'uccsd_conv': int(ll[6]) == 1,
                        'uhf_conv': uhf_conv,
                    })
                    uhf_conv = True
                elif l.startswith('SCF not converged'):
                    uhf_conv = False

for ia in [3]:
    with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/07-uhf-act-dm.out" % (ia + 2, ia), "r") as f:
        uhf_conv = True
        for l in f.readlines():
            if '::' in l:
                ll = l.strip().split()
                idx = int(ll[0])
                if idx not in [2, 3, 4, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 25, 26, 27]:
                    data['large_uccsd'][ia][idx].append({
                        'init': ll[3],
                        'proj': ll[5].strip(),
                        'uhf': float(ll[4]),
                        'uhf_conv': uhf_conv,
                    })
                uhf_conv = True
            elif l.startswith('SCF not converged'):
                uhf_conv = False

for ia in [3]:
    for ix in [2, 3, 4, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 25, 26, 27]:
        bs_map = {(l.split()[3], l.split()[-3]): int(l.split('BS')[1].strip()) for l in open("../07-FEMOCO-Large/01-UHF/33-uhf-spin-%d.out" % (ix, ), "r").readlines() if '::' in l}
        init_map = {(int(l.split()[0]), int(l.split()[2]), l.split()[-1]): l.split()[3]
            for l in open("../07-FEMOCO-Large/%02d-UCCSD-A%d/03-uhf-min-24.out" % (ia + 2, ia), "r").readlines()}
        for jx in range(24):
            d = {}
            with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/03-uccsd-act-%d-%d.out" % (ia + 2, ia, ix, jx), "r") as f:
                uhf_conv = True
                for l in f.readlines():
                    if '::' in l and int(l.strip().split()[-2]) <= 1:
                        ll = l.strip().split()
                        d = {
                            'init': ll[2],
                            'proj': ll[7].strip(),
                            'uhf': float(ll[3]),
                            'uccsd': float(ll[4]),
                            'uccsd(t)': float(ll[5]),
                            'uccsd_conv': int(ll[6]) == 1,
                            'uhf_conv': uhf_conv,
                        }
                        uhf_conv = True
                    elif l.startswith('SCF not converged'):
                        uhf_conv = False
                        d['uhf_conv'] = False
                    elif l.startswith('--- '):
                        d['proj'] = l.split()[-2]
                        d['init'] = init_map[(int(l.split()[1]), int(l.split()[3]), l.split()[4])]
                    elif 'SCF energy =' in l:
                        d['uhf'] = float(l.split('SCF energy =')[1].split()[0])
                    elif 'UCCSD converged' in l:
                        d['uccsd_conv'] = True
                    elif l.startswith('E(UCCSD) ='):
                        d['uccsd'] = float(l.split()[2])
                    elif l.startswith('UCCSD(T) correction ='):
                        d['uccsd(t)'] = float(l.split()[3]) + d['uccsd']
                    elif 'E_corr(UCCSD) =' in l:
                        d['uccsd'] = float(l.split('E_corr(UCCSD) =')[1].split()[0]) + d['uhf']
            d['uhf_conv'] = d.get('uhf_conv', True)
            d['uccsd_conv'] = d.get('uccsd_conv', False)
            d['uccsd(t)'] = d.get('uccsd(t)', None)
            assert 'uccsd' in d and 'uhf' in d
            data['large_uccsd'][ia][bs_map[(d['init'], d['proj'])]].append(d)

for ia in [2]:
    for ix in range(29):
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/06-uccsd-act-add-%d.out" % (ia + 2, ia, ix), "r") as f:
            uhf_conv = True
            for l in f.readlines():
                if '::' in l and int(l.strip().split()[-2]) <= 1:
                    ll = l.strip().split()
                    assert (ll[3], ll[-1]) not in [(x['init'], x['proj']) for x in data['large_uccsd'][ia][int(ll[0])]]
                    data['large_uccsd'][ia][int(ll[0])].append({
                        'init': ll[3],
                        'proj': ll[-1].strip(),
                        'uhf': float(ll[4]),
                        'uccsd': float(ll[5]),
                        'uccsd(t)': float(ll[6]),
                        'uccsd_conv': int(ll[7]) == 1,
                        'uhf_conv': uhf_conv,
                    })
                    uhf_conv = True
                elif l.startswith('SCF not converged'):
                    uhf_conv = False

for ia in [3]:
    for ix in [*range(0, 11), 16, 17, 18, 19]:
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/07-uccsd-act-add-re-%d.out" % (ia + 2, ia, ix), "r") as f:
            uhf_conv = True
            for l in f.readlines():
                if '::' in l and int(l.strip().split()[-2]) <= 1:
                    ll = l.strip().split()
                    assert (ll[3], ll[-1]) not in [(x['init'], x['proj']) for x in data['large_uccsd'][ia][int(ll[0])]]
                    data['large_uccsd'][ia][int(ll[0])].append({
                        'init': ll[3],
                        'proj': ll[-1].strip(),
                        'uhf': float(ll[4]),
                        'uccsd': float(ll[5]),
                        'uccsd(t)': float(ll[6]),
                        'uccsd_conv': int(ll[7]) == 1,
                        'uhf_conv': uhf_conv,
                    })
                    uhf_conv = True
                elif l.startswith('SCF not converged'):
                    uhf_conv = False

for ia in range(0, 2):
    for ik, key in enumerate(['uhf-elec-pop-noryd', 'uhf-spin', 'uhf-nuc-density',
        'uccsd-elec-pop-noryd', 'uccsd-spin', 'uccsd-nuc-density']):
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-%s.out" % (ia + 2, ia, ik + [60, 70][ia], key), "r") as f:
            for l in f.readlines():
                ll = l.strip().split()
                ix = int(ll[0])
                idx = [x['init'] for x in data['large_uccsd'][ia][ix]].index(ll[3])
                data['large_uccsd'][ia][ix][idx][key] = [float(x) for x in ll[4:4 + 8]]
    for ix in range(35):
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-uccsd-act-dm-%d.out" % (ia + 2, ia, [50, 51][ia], ix), "r") as f:
            for l in f.readlines():
                if '::' in l:
                    ll = l.strip().split()
                    idx = [x['init'] for x in data['large_uccsd'][ia][ix]].index(ll[3])
                    dm1 = np.load('../07-FEMOCO-Large/%02d-UCCSD-A%d/data/uccsd-dm1-%d-%04d.npy' % (ia + 2, ia, [50, 51][ia], ix))
                    mo = np.load('../07-FEMOCO-Large/%02d-UCCSD-A%d/data/uccsd-mo-coeff-%d-%04d.npy' % (ia + 2, ia, [50, 51][ia], ix))
                    sfdm = mo[0] @ dm1[0] @ mo[0].T + mo[1] @ dm1[1] @ mo[1].T
                    data['large_uccsd'][ia][ix][idx]['uccsd_nat-occ'] = np.linalg.eigh(sfdm)[0][::-1].tolist()

for ia in range(2, 4):
    for ik, key in enumerate(['uhf-elec-pop-noryd', 'uhf-spin', 'uhf-nuc-density']):
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-%s.out" % (ia + 2, ia, ik + [10, 20][ia - 2], key), "r") as f:
            for l in f.readlines():
                ll = l.strip().split()
                ix = int(ll[0])
                idx = [x['init'] for x in data['large_uccsd'][ia][ix]].index(ll[3])
                data['large_uccsd'][ia][ix][idx][key] = [float(x) for x in ll[4:4 + 8]]

for ia in [2]:
    for ik, key in enumerate(['uccsd-elec-pop-noryd', 'uccsd-spin', 'uccsd-nuc-density']):
        init_map = {int(l.split()[0]): (l.split()[3], l.split()[8]) for l in open("../07-FEMOCO-Large/%02d-UCCSD-A%d/02-uccsd-min.out" % (ia + 2, ia), "r").readlines()}
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-%s.out" % (ia + 2, ia, ik + 83, key), "r") as f:
            for l in f.readlines():
                ll = l.strip().split()
                bs_map = {(l.split()[3], l.split()[-3]): int(l.split('BS')[1].strip())
                    for l in open("../07-FEMOCO-Large/01-UHF/32-uhf-spin-%d.out" % (int(ll[0]), ), "r").readlines() if '::' in l}
                ix = bs_map[init_map[int(ll[0])]]
                idx = [x['init'] for x in data['large_uccsd'][ia][ix]].index(ll[3])
                data['large_uccsd'][ia][ix][idx][key] = [float(x) for x in ll[4:4 + 8]]
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-%s.out" % (ia + 2, ia, ik + 86, key), "r") as f:
            for l in f.readlines():
                ll = l.strip().split()
                ix = int(ll[0])
                idx = [x['init'] for x in data['large_uccsd'][ia][ix]].index(ll[3])
                data['large_uccsd'][ia][ix][idx][key] = [float(x) for x in ll[4:4 + 8]]
        for fx in [52, 56]:
            for pix in range(35):
                if not os.path.exists("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-uccsd-act-dm-%d.out" % (ia + 2, ia, fx, pix)):
                    continue
                with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-uccsd-act-dm-%d.out" % (ia + 2, ia, fx, pix), "r") as f:
                    for l in f.readlines():
                        if '::' in l:
                            ll = l.strip().split()
                            bs_map = {(l.split()[3], l.split()[-3]): int(l.split('BS')[1].strip())
                                for l in open("../07-FEMOCO-Large/01-UHF/32-uhf-spin-%d.out" % (pix, ), "r").readlines() if '::' in l}
                            ix = bs_map[init_map[pix]] if fx == 52 else pix
                            idx = [x['init'] for x in data['large_uccsd'][ia][ix]].index(ll[2])
                            dm1 = np.load('../07-FEMOCO-Large/%02d-UCCSD-A%d/data/uccsd-dm1-%d-%04d.npy' % (ia + 2, ia, fx, pix))
                            mo = np.load('../07-FEMOCO-Large/%02d-UCCSD-A%d/data/uccsd-mo-coeff-%d-%04d.npy' % (ia + 2, ia, fx, pix))
                            sfdm = mo[0] @ dm1[0] @ mo[0].T + mo[1] @ dm1[1] @ mo[1].T
                            data['large_uccsd'][ia][ix][idx]['uccsd_nat-occ'] = np.linalg.eigh(sfdm)[0][::-1].tolist()

for ia in [3]:
    for ik, key in enumerate(['uccsd-elec-pop-noryd', 'uccsd-spin', 'uccsd-nuc-density']):
        init_map = {int(l.split()[0]): (l.split()[3], l.split()[8]) for l in open("../07-FEMOCO-Large/%02d-UCCSD-A%d/03-uccsd-min.out" % (ia + 2, ia), "r").readlines()}
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-%s.out" % (ia + 2, ia, ik + 80, key), "r") as f:
            for l in f.readlines():
                ll = l.strip().split()
                bs_map = {(l.split()[3], l.split()[-3]): int(l.split('BS')[1].strip())
                    for l in open("../07-FEMOCO-Large/01-UHF/33-uhf-spin-%d.out" % (int(ll[0]), ), "r").readlines() if '::' in l}
                ix = bs_map[init_map[int(ll[0])]]
                idx = [x['init'] for x in data['large_uccsd'][ia][ix]].index(ll[3])
                data['large_uccsd'][ia][ix][idx][key] = [float(x) for x in ll[4:4 + 8]]
        with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-%s.out" % (ia + 2, ia, ik + 83, key), "r") as f:
            for l in f.readlines():
                ll = l.strip().split()
                ix = int(ll[0])
                idx = [x['init'] for x in data['large_uccsd'][ia][ix]].index(ll[3])
                data['large_uccsd'][ia][ix][idx][key] = [float(x) for x in ll[4:4 + 8]]
        for fx in [63, 67]:
            for pix in range(35):
                if not os.path.exists("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-uccsd-act-dm-%d.out" % (ia + 2, ia, fx, pix)):
                    continue
                with open("../07-FEMOCO-Large/%02d-UCCSD-A%d/%d-uccsd-act-dm-%d.out" % (ia + 2, ia, fx, pix), "r") as f:
                    for l in f.readlines():
                        if '::' in l:
                            ll = l.strip().split()
                            bs_map = {(l.split()[3], l.split()[-3]): int(l.split('BS')[1].strip())
                                for l in open("../07-FEMOCO-Large/01-UHF/33-uhf-spin-%d.out" % (pix, ), "r").readlines() if '::' in l}
                            ix = bs_map[init_map[pix]] if fx == 63 else pix
                            idx = [x['init'] for x in data['large_uccsd'][ia][ix]].index(ll[2])
                            dm1 = np.load('../07-FEMOCO-Large/%02d-UCCSD-A%d/data/uccsd-dm1-%d-%04d.npy' % (ia + 2, ia, fx, pix))
                            mo = np.load('../07-FEMOCO-Large/%02d-UCCSD-A%d/data/uccsd-mo-coeff-%d-%04d.npy' % (ia + 2, ia, fx, pix))
                            sfdm = mo[0] @ dm1[0] @ mo[0].T + mo[1] @ dm1[1] @ mo[1].T
                            data['large_uccsd'][ia][ix][idx]['uccsd_nat-occ'] = np.linalg.eigh(sfdm)[0][::-1].tolist()

with open("14-large-uccsd.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

for ia in range(4):
    print('Number of UCCSD calculations: ', sum([len(data['large_uccsd'][ia][ix]) for ix in range(35)]))
