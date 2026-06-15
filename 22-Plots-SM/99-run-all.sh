echo "### 01-sm-dmrg-extra"
python -u 01-sm-dmrg-extra.py > 01-sm-dmrg-extra.out
echo "### 02-sm-ghf-ener"
python -u 02-sm-ghf-ener.py > 02-sm-ghf-ener.out
echo "### 04-sm-ref-dep-cc"
python -u 04-sm-ref-dep-cc.py 10 > 04-sm-ref-dep-cc-10.out
python -u 04-sm-ref-dep-cc.py 12 > 04-sm-ref-dep-cc-12.out
python -u 04-sm-ref-dep-cc.py 20 > 04-sm-ref-dep-cc-20.out
python -u 04-sm-ref-dep-cc.py 25 > 04-sm-ref-dep-cc-25.out
python -u 04-sm-ref-dep-cc.py 26 > 04-sm-ref-dep-cc-26.out
python -u 04-sm-ref-dep-cc.py 27 > 04-sm-ref-dep-cc-27.out
echo "### 05-sm-ref-dep-prop"
python -u 05-sm-ref-dep-prop-a-nuc.py
python -u 05-sm-ref-dep-prop-b-pop.py
echo "### 06-sm-llduc-bs-levels"
python -u 06-sm-llduc-bs-levels.py > 06-sm-llduc-bs-levels.out
echo "### 07-sm-fe2-ssq"
python -u 07-sm-fe2-ssq-a-std-ucc.py
python -u 07-sm-fe2-ssq-b-mix-ucc.py
python -u 07-sm-fe2-ssq-c-std-udmrg.py
python -u 07-sm-fe2-ssq-d-mix-udmrg.py
python -u 07-sm-fe2-ssq-e-std-corr.py
python -u 07-sm-fe2-ssq-f-mix-corr.py
echo "### 08-sm-fm-ssq"
python -u 08-sm-fm-ssq.py
echo "### 09-sm-fm-fno-ssq"
python -u 09-sm-fm-fno-ssq.py
echo "### 10-sm-dmrg-cc-cost"
python -u 10-sm-dmrg-cc-cost-a-dmrg.py
python -u 10-sm-dmrg-cc-cost-b-cc.py
echo "### 11-sm-fm-dmrg-energy"
python -u 11-sm-fm-dmrg-energy.py > 11-sm-fm-dmrg-energy.out
echo "### 12-sm-fm-dmrg-dw"
python -u 12-sm-fm-dmrg-dw-a-25.py > 12-sm-fm-dmrg-dw-a-25.out
python -u 12-sm-fm-dmrg-dw-b-26.py > 12-sm-fm-dmrg-dw-b-26.out
echo "### 13-sm-dmrg-cc-nat-occ"
python -u 13-sm-dmrg-cc-nat-occ.py 25
python -u 13-sm-dmrg-cc-nat-occ.py 26
python -u 13-sm-dmrg-cc-nat-occ.py 27
echo "### 14-sm-act"
python -u 14-sm-act.py
echo "### 15-sm-act-energy"
python -u 15-sm-act-energy-a-uhf.py
python -u 15-sm-act-energy-b-uccsd.py
python -u 15-sm-act-energy-c-uccsd-t.py
echo "### 16-sm-large-pop"
python -u 16-sm-large-pop-a-nuc-bs7.py
python -u 16-sm-large-pop-b-nuc-bs8.py
python -u 16-sm-large-pop-c-nuc-bs10.py
python -u 16-sm-large-pop-d-elec-bs7.py
python -u 16-sm-large-pop-e-elec-bs8.py
python -u 16-sm-large-pop-f-elec-bs10.py
echo "### 17-sm-large-pop-order"
python -u 17-sm-large-pop-order-a.py
python -u 17-sm-large-pop-order-b.py
python -u 17-sm-large-pop-order-c.py
python -u 17-sm-large-pop-order-d.py
python -u 17-sm-large-pop-order-e.py
python -u 17-sm-large-pop-order-f.py
python -u 17-sm-large-pop-order-g.py
python -u 17-sm-large-pop-order-h.py
echo "### 18-sm-large-cc-nat-occ"
python -u 18-sm-large-cc-nat-occ.py 25
python -u 18-sm-large-cc-nat-occ.py 26
python -u 18-sm-large-cc-nat-occ.py 27
echo "### 19-sm-fm-qmmm-ener"
python -u 19-sm-fm-qmmm-ener-a-uhf.py
python -u 19-sm-fm-qmmm-ener-b-uccsd.py
python -u 19-sm-fm-qmmm-ener-c-uccsd-t.py
echo "### 20-sm-fe4-spin-corr"
python -u 20-sm-fe4-spin-corr.py
echo "### 21-sm-fe4-extra-a"
python -u 21-sm-fe4-extra-a.py
echo "### 21-sm-fe4-extra-b"
python -u 21-sm-fe4-extra-b.py
echo "### 22-sm-fe4-ener"
python -u 22-sm-fe4-ener.py > 22-sm-fe4-ener.out
echo "### 23-sm-fe4-udmrg-fit"
python -u 23-sm-fe4-udmrg-fit.py > 23-sm-fe4-udmrg-fit.out
echo "### 29-sm-heis-ener"
python -u 29-sm-heis-ener.py
