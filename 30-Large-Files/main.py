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


import subprocess, os, json, sys

tok = os.environ.get('TOK', None)
url = "https://api.github.com/repos/hczhai/femoco-scripts-2026/releases"
url2 = "https://uploads.github.com/repos/hczhai/femoco-scripts-2026/releases"

tok = """ -H "Authorization: Bearer %s" """ % tok if tok is not None else ""

rel_curl = """curl -L -H "Accept: application/vnd.github+json" \\
    %s -H "X-GitHub-Api-Version: 2022-11-28" "%s" """ % (tok, url)
post_curl = """curl -L -X POST -H "Accept: application/vnd.github+json" \\
    %s -H "X-GitHub-Api-Version: 2022-11-28" \\
    -H "Content-Type: application/octet-stream" "%s" -T "%s" --no-buffer"""
get_curl = """curl -L -H "Accept: application/octet-stream" %s \\
    -H "X-GitHub-Api-Version: 2022-11-28" "%s" -o "%s" """

process = subprocess.Popen(rel_curl, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = process.communicate()

if process.returncode == 0:
    try:
        jsd = json.loads(stdout.decode())
    except:
        print(stdout.decode())
        quit()
    rel_data = {xx["tag_name"]: xx["id"] for xx in jsd}
else:
    print('Error: ', stderr.decode())
    quit()

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        for a, b in sorted(rel_data.items(), key=lambda x: x[1]):
            print("%30s %15s" % (a, b))
    elif sys.argv[1] == "download":
        remote = sys.argv[2]
        name = sys.argv[3]
        local = sys.argv[4]
        rel = rel_data[remote]
        ia, ib = int(sys.argv[5]), int(sys.argv[6])
        for ix in range(ia, ib):
            print("=== %s : %5s / %5s - %5s" % (name, ix, ia, ib))
            xlocal = local if ix == -1 or '%' not in local else local % ix
            if xlocal.endswith('.npy') or xlocal.endswith('.chk'):
                xlocal = os.path.dirname(xlocal)
            if not os.path.exists(xlocal):
                os.makedirs(xlocal)
            assets = [(zz["name"], zz["id"]) for xx in jsd if xx["id"] == rel for zz in xx["assets"]]
            fname, fnid = [(nm, jx) for nm, jx in assets if nm.startswith("TOTAL-%s-%d-" % (name, ix))][0]
            xcurl = get_curl % (tok, url + "/assets/%s" % (fnid, ), fname)
            process = subprocess.Popen(xcurl, shell=True, stdout=sys.stdout, stderr=sys.stderr)
            process.wait()
            if process.returncode != 0:
                quit()
            print()
            total_n = int(fname[:-4].split('-')[-1])
            for ixf in range(total_n):
                print("--- %5s / %5s" % (ixf, total_n))
                zfn = "Z-%s-%d-%d.zip" % (name, ix, ixf)
                zfname, zfnid = [(nm, ix) for nm, ix in assets if nm == zfn][0]
                xcurl = get_curl % (tok, url + "/assets/%s" % (zfnid, ), zfname)
                process = subprocess.Popen(xcurl, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                process.wait()
                if process.returncode != 0:
                    quit()
                print()
                xcmd = "unzip %s -d %s" % (zfn, xlocal)
                process = subprocess.Popen(xcmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                process.wait()
                if process.returncode != 0:
                    quit()
                os.remove(zfn)
                print('unzip ok')
            print('checking')
            fns = list(os.listdir(xlocal))
            for l in open(fname, 'r').readlines():
                a, b = l.split()
                print(a, b)
                assert a.split('/')[-1] in fns
                assert os.path.getsize(os.path.join(xlocal, a.split('/')[-1])) == int(b)
            print('checking ok')
            print('merging')
            splits = sorted([fn for fn in os.listdir(xlocal) if '.SPLIT.' in fn])
            if len(splits) != 0:
                assert all([fn.split('.SPLIT.')[0] == splits[0].split('.SPLIT.')[0] for fn in splits])
                pre_name = splits[0].split('.SPLIT.')[0]
                xcmd = "cat %s/%s.SPLIT.* > %s/%s" % (xlocal, pre_name, xlocal, pre_name)
                process = subprocess.Popen(xcmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                process.wait()
                if process.returncode != 0:
                    quit()
                for fn in splits:
                    os.remove(xlocal + "/" + fn)
                print('merging ok')
            os.remove(fname)
    elif sys.argv[1] == "upload":
        remote = sys.argv[2]
        name = sys.argv[3]
        local = sys.argv[4]
        ia, ib = int(sys.argv[5]), int(sys.argv[6])
        for ix in range(ia, ib):
            print("=== %s : %5s / %5s - %5s" % (name, ix, ia, ib))
            cur_sz = 0
            fn_gp = [[]]
            xlocal = local if ix == -1 else local % ix
            if os.path.isdir(xlocal):
                for fn in sorted(os.listdir(xlocal), key=lambda x:
                    ['bin' not in x, 'INFO' not in x, 'LEFT' not in x,
                    None if x.split('.')[-1] == 'bin' or not x.split('.')[-1].lstrip('-').isdigit()
                    else int(x.split('.')[-1]), x]):
                    filepath = os.path.join(xlocal, fn)
                    if os.path.isfile(filepath):
                        sz = os.path.getsize(filepath)
                        assert sz <= 1.95 * 1024 ** 3
                        if cur_sz + sz > 1.9 * 1024 ** 3:
                            fn_gp.append([])
                            cur_sz = 0
                        fn_gp[-1].append((filepath, sz))
                        cur_sz += sz
            else:
                sz = os.path.getsize(xlocal)
                if sz <= 1.95 * 1024 ** 3:
                    fn_gp[-1].append((xlocal, sz))
                else:
                    fn_gp = []
                    xcmd = "split -b 2093796556 -a 4 --numeric-suffixes=1 %s %s.SPLIT. " % (xlocal, xlocal)
                    process = subprocess.Popen(xcmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                    process.wait()
                    if process.returncode != 0:
                        quit()
                    for i in range(1, 9999):
                        filepath = "%s.SPLIT.%04d" % (xlocal, i)
                        if os.path.isfile(filepath):
                            sz = os.path.getsize(filepath)
                            assert sz <= 1.95 * 1024 ** 3
                            fn_gp.append([(filepath, sz)])
            fname = "TOTAL-%s-%d-%d.txt" % (name, ix, len(fn_gp))
            xcurl = post_curl % (tok, url2 + "/%s/assets?name=%s" % (rel_data[remote], fname), fname)
            with open(fname, "w") as f:
                for x, z in [z for x in fn_gp for z in x]:
                    f.write("%s %s\n" % (x, z))
            process = subprocess.Popen(xcurl, shell=True, stdout=sys.stdout, stderr=sys.stderr)
            process.wait()
            if process.returncode != 0:
                quit()
            print()
            os.remove(fname)
            for ixf, fns in enumerate(fn_gp):
                print("--- %5s / %5s" % (ixf, len(fn_gp)))
                zfn = "Z-%s-%d-%d.zip" % (name, ix, ixf)
                xcmd = "zip -1 -j %s " % zfn + " ".join([x for x, _ in fns])
                process = subprocess.Popen(xcmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                process.wait()
                if process.returncode != 0:
                    quit()
                print('zip ok')
                xcurl = post_curl % (tok, url2 + "/%s/assets?name=%s" % (rel_data[remote], zfn), zfn)
                process = subprocess.Popen(xcurl, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                process.wait()
                if process.returncode != 0:
                    quit()
                print()
                os.remove(zfn)
                if len(fns) == 1 and '.SPLIT.' in fns[0][0]:
                    os.remove(fns[0][0])
