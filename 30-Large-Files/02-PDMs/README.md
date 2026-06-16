
The 2PDMs can be downloaded using the following commands:

Note: do not run multiple ``python ../main.py download`` simultaneously for the same target folders.

## FeMo-co LLDUC Model FNO-UCC NFrozen=46/42/38/34/30/26/22 2PDMs (BS7-C, BS8-E, BS8-F)

```bash
python ../main.py download v0 "fm-fno-uccsd-f46-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/01-cc-02-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsd-f42-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/02-cc-02-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsd-f38-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/03-cc-02-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsd-f34-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/04-cc-02-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsd-f30-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/05-cc-02-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsd-f26-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/06-cc-02-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsd-f22-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/07-cc-02-%02d-pdm2.npy" 25 28

python ../main.py download v0 "fm-fno-uccsdt-f46-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/01-cc-03-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsdt-f42-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/02-cc-03-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsdt-f38-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/03-cc-03-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsdt-f34-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/04-cc-03-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsdt-f30-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/05-cc-03-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsdt-f26-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/06-cc-03-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-fno-uccsdt-f22-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/07-cc-03-%02d-pdm2.npy" 25 28

python ../main.py download v0 "fm-fno-uccsdtq-f46-2pdm" "/tmp/06-FEMOCO-LLDUC/14-FNO-UCC-NPDM/data/01-cc-04-%02d-pdm2.npy" 25 28
```

## FeMo-co LLDUC Model FNO-UDMRG NFrozen=46/42/38/34/30/26 D=5k/8k 2PDMs (BS7-C, BS8-E, BS8-F)

```bash
python ../main.py download v0 "fm-fno-udmrg-5k-f46" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/21-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-5k-f42" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/22-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-5k-f38" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/23-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-5k-f34" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/24-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-5k-f30" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/25-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-5k-f26" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/26-%02d-2pdm.npy" 25 28

python ../main.py download v0 "fm-fno-udmrg-8k-f46" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/31-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-8k-f42" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/32-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-8k-f38" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/33-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-8k-f34" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/34-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-8k-f30" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/35-%02d-2pdm.npy" 25 28
python ../main.py download v0 "fm-fno-udmrg-8k-f26" "/tmp/06-FEMOCO-LLDUC/15-FNO-UDMRG/data/36-%02d-2pdm.npy" 25 28
```

## FeMo-co LLDUC Model UCCSD/UCCSDT 2PDMs (BS7-C, BS8-E, BS8-F)

```bash
python ../main.py download v0 "fm-uccsd-30-2pdm" "/tmp/06-FEMOCO-LLDUC/16-UCC-NPDM/data/30-cc-02-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-uccsd-32-2pdm" "/tmp/06-FEMOCO-LLDUC/16-UCC-NPDM/data/32-cc-02-%02d-pdm2.npy" 25 28
python ../main.py download v0 "fm-uccsdt-30-2pdm" "/tmp/06-FEMOCO-LLDUC/16-UCC-NPDM/data/30-cc-03-%02d-pdm2.npy" 25 28
```

## FeMo-co LLDUC Model UDMRG D=5k/6k/7k/8k/9k/11k 2PDMs (BS8-E)

```bash
python ../main.py download v0 "fm-udmrg-5k-2pdm" "/tmp/06-FEMOCO-LLDUC/17-UDMRG-NPDM/data/01-%02d-2pdm.npy" 26 27
python ../main.py download v0 "fm-udmrg-6k-2pdm" "/tmp/06-FEMOCO-LLDUC/17-UDMRG-NPDM/data/02-%02d-2pdm.npy" 26 27
python ../main.py download v0 "fm-udmrg-7k-2pdm" "/tmp/06-FEMOCO-LLDUC/17-UDMRG-NPDM/data/03-%02d-2pdm.npy" 26 27
python ../main.py download v0 "fm-udmrg-8k-2pdm" "/tmp/06-FEMOCO-LLDUC/17-UDMRG-NPDM/data/04-%02d-2pdm.npy" 26 27
python ../main.py download v0 "fm-udmrg-9k-2pdm" "/tmp/06-FEMOCO-LLDUC/17-UDMRG-NPDM/data/05-%02d-2pdm.npy" 26 27
python ../main.py download v0 "fm-udmrg-11k-2pdm" "/tmp/06-FEMOCO-LLDUC/17-UDMRG-NPDM/data/06-%02d-2pdm.npy" 26 27
```
