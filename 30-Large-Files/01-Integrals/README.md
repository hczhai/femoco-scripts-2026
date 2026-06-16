
The integrals for large active spaces can be downloaded using the following commands:

Note: do not run multiple ``python ../main.py download`` simultaneously for the same target folders.

## FeMo-co Two-electron Integrals for 117o, 180o, 285o, 404o Active Spaces

See ``07-FEMOCO-Large/00-Active-Space/99-README.md``.

```bash
python ../main.py download v0 "fm-large-act-g2e" "/tmp/07-FEMOCO-Large/00-Active-Space/data/18-fm-active/%02d-g2e.npy" 0 4
```

## FeMo-co Large Active Spaces MO coefficients

```bash
python ../main.py download v0 "fm-large-act-chk" "/tmp/07-FEMOCO-Large/00-Active-Space/00-ref-hs_bp86.chk" -1 0
python ../main.py download v0 "fm-large-act-mo" "/tmp/07-FEMOCO-Large/00-Active-Space/01-mo-coeff.npy" -1 0
python ../main.py download v0 "fm-large-act-dm" "/tmp/07-FEMOCO-Large/00-Active-Space/01-uhf-dm1.npy" -1 0
```
