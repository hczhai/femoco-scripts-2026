
The MPS files can be downloaded using the following commands:

Note: do not run multiple ``python ../main.py download`` simultaneously for the same target folders.

## FeMo-co LLDUC Model UDMRG D=12k/18k MPS (BS7-C, BS8-E)

```bash
python ../main.py download v0 "fm-udmrg-18k" "/tmp/06-FEMOCO-LLDUC/02-UDMRG-25/mps-18k-%02d" 25 26
python ../main.py download v0 "fm-udmrg-18k" "/tmp/06-FEMOCO-LLDUC/03-UDMRG-26/mps-18k-%02d" 26 27

python ../main.py download v0 "fm-udmrg-12k" "/tmp/06-FEMOCO-LLDUC/02-UDMRG-25/mps-12k-%02d" 25 26
python ../main.py download v0 "fm-udmrg-12k" "/tmp/06-FEMOCO-LLDUC/02-UDMRG-26/mps-12k-%02d" 26 27
```

## FeMo-co LLDUC Model UDMRG D=5k MPS

```bash
python ../main.py download v0 "fm-udmrg-5k" "/tmp/06-FEMOCO-LLDUC/04-UDMRG-5k/mps-5k-%02d" 0 35
```

## FeMo-co LLDUC Model UDMRG D=8k MPS

```bash
python ../main.py download v0 "fm-udmrg-8k" "/tmp/06-FEMOCO-LLDUC/09-UDMRG-8k/mps-8k-%02d" 4 5
python ../main.py download v0 "fm-udmrg-8k" "/tmp/06-FEMOCO-LLDUC/09-UDMRG-8k/mps-8k-%02d" 25 28
```

## Fe4S4 Model UDMRG D=3k/8k/12k MPS

```bash
python ../main.py download v0 "fe4-udmrg-3k" "/tmp/04-Fe4S4/01-UHF-UCC-UDMRG/mps-3k-%02d" 0 3
python ../main.py download v0 "fe4-udmrg-8k" "/tmp/04-Fe4S4/01-UHF-UCC-UDMRG/mps-8k-%02d" 0 3
python ../main.py download v0 "fe4-udmrg-12k" "/tmp/04-Fe4S4/01-UHF-UCC-UDMRG/mps-12k-%02d" 0 3
```

## Fe4S4 Model spin-adapted DMRG D=14k/20k MPS

```bash
python ../main.py download v0 "fe4-su2-dmrg-14k" "/tmp/04-Fe4S4/02-SU2-DMRG-20k/mps-14k-%02d" 0 1
python ../main.py download v0 "fe4-su2-dmrg-14k" "/tmp/04-Fe4S4/02-SU2-DMRG-20k/mps-14k-%02d" 2 3
python ../main.py download v0 "fe4-su2-dmrg-20k" "/tmp/04-Fe4S4/02-SU2-DMRG-20k/mps-20k-%02d" 0 3
```
