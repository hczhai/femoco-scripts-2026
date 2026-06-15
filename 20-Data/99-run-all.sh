echo "### 00-llduc-uhf"
python -u 00-llduc-uhf.py | tee 00-llduc-uhf.out
echo "### 01-llduc-uccsd"
python -u 01-llduc-uccsd.py | tee 01-llduc-uccsd.out
echo "### 02-llduc-uccsdt"
python -u 02-llduc-uccsdt.py | tee 02-llduc-uccsdt.out
echo "### 03-llduc-udmrg-5k"
python -u 03-llduc-udmrg-5k.py | tee 03-llduc-udmrg-5k.out
echo "### 04-llduc-udmrg-18k"
python -u 04-llduc-udmrg-18k.py | tee 04-llduc-udmrg-18k.out
echo "### 05-llduc-udmrg-8k"
python -u 05-llduc-udmrg-8k.py | tee 05-llduc-udmrg-8k.out
echo "### 06-llduc-fno-ucc"
python -u 06-llduc-fno-ucc.py | tee 06-llduc-fno-ucc.out
echo "### 07-llduc-ghf"
python -u 07-llduc-ghf.py | tee 07-llduc-ghf.out
echo "### 08-llduc-gccsd"
python -u 08-llduc-gccsd.py | tee 08-llduc-gccsd.out
echo "### 09-llduc-udmrg-5k-ref"
python -u 09-llduc-udmrg-5k-ref.py | tee 09-llduc-udmrg-5k-ref.out
echo "### 10-fe2s2"
python -u 10-fe2s2.py  | tee 10-fe2s2.out
echo "### 11-llduc-fno-udmrg"
python -u 11-llduc-fno-udmrg.py | tee 11-llduc-fno-udmrg.out
echo "### 12-large-act"
python -u 12-large-act.py | tee 12-large-act.out
echo "### 13-large-uhf"
python -u 13-large-uhf.py | tee 13-large-uhf.out
echo "### 14-large-uccsd"
python -u 14-large-uccsd.py | tee 14-large-uccsd.out
echo "### 15-qmmm-uhf"
python -u 15-qmmm-uhf.py | tee 15-qmmm-uhf.out
echo "### 16-qmmm-uccsd"
python -u 16-qmmm-uccsd.py | tee 16-qmmm-uccsd.out
echo "### 17-large-ump2"
python -u 17-large-ump2.py | tee 17-large-ump2.out
echo "### 18-fe4s4-sa-dmrg"
python -u 18-fe4s4-sa-dmrg.py | tee 18-fe4s4-sa-dmrg.out
echo "### 19-fe4s4-uhf"
python -u 19-fe4s4-uhf.py | tee 19-fe4s4-uhf.out
echo "### 20-fe4s4-uccsd"
python -u 20-fe4s4-uccsd.py | tee 20-fe4s4-uccsd.out
echo "### 21-fe4s4-uccsdtq"
python -u 21-fe4s4-uccsdtq.py | tee 21-fe4s4-uccsdtq.out
echo "### 22-fe4s4-udmrg"
python -u 22-fe4s4-udmrg.py | tee 22-fe4s4-udmrg.out
echo "### 23-mofe3"
python -u 23-mofe3.py | tee 23-mofe3.out
echo "### 24-heisenberg"
python -u 24-heisenberg.py | tee 24-heisenberg.out
echo "### 25-fe-oxidation"
python -u 25-fe-oxidation.py | tee 25-fe-oxidation.out
echo "### 26-llduc-composite"
python -u 26-llduc-composite.py | tee 26-llduc-composite.out

echo "### 70-dmrg-arxiv-2506-13386"
python -u 70-dmrg-arxiv-2506-13386.py
echo "### 71-dmrg-arxiv-2208-02199"
python -u 71-dmrg-arxiv-2208-02199.py
echo "### 72-dft-10.1002-qua.25627"
python -u 72-dft-10.1002-qua.25627.py
echo "### 73-expr-jacs-2000-122-4926"
python -u 73-expr-jacs-2000-122-4926.py
echo "### 80-bs-states"
python -u 80-bs-states.py
echo "### 90-colors"
python -u 90-colors.py