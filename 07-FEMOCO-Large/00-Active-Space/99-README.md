
## LLDUC model (76o, 113e)

- ``18-fm-active/11-orig-mo-coeff.npy`` (2818, 2818) AO -> LMO (LLDUC, all) (obtained from the LLDUC work)
- ``00-ref-hs_bp86.chk`` (2, 2818, 2818) AO 1pdm for UKS (high-spin) (obtained from the LLDUC work)
- ``18-fm-active/10-full-overlap.npy`` (2818, 2818) AO -> AO overlap (generated in ``00-ovlp.py``)

## [4] Active space (994o, 423e)

- ``01-mo-coeff.npy`` (2, 2818, 2818) AO -> UKS-MO (high-spin) (generated in ``01-init-scf.py``)
- ``01-uhf-dm1.npy`` (2, 2818, 2818) UKS (high-spin) 1pdm in AO basis (generated in ``01-init-scf.py``)

- ``04-lo-coeff.npy`` (2818, 2818) AO -> LMO (Localized spin-averaged natural orbtial from UKS) (generated in ``04-uno.py``)

- ``05-sel-mo-occ.npy`` (994, ) LMO occupation in (994o, 423e) active spce (generated in ``05-fm-sel.py``)
- ``05-sel-mo-coeff.npy`` (2818, 994) AO -> LMO (model-4, active) (generated in ``05-fm-sel.py``)
- ``05-sel-core-mo-coeff.npy`` (2818, 362) AO -> LMO (model-4, core) (generated in ``05-fm-sel.py``)
- ``05-sel-virt-mo-coeff.npy`` (2818, 1462) AO -> LMO (model-4, active) (generated in ``05-fm-sel.py``)

- ``08-proj-dm1.npy`` (2, 994, 994) UHF initial guess projected from AO into LMO (model-4, active) (generated in ``08-dm1-proj.py``)

- ``12-fm-active/ecore.npy`` (model-4, core energy) (generated in ``12-other-ints-solv.py``)
- ``12-fm-active/h1e.npy`` (994, 994) (model-4, active h1e) (generated in ``12-other-ints-solv.py``)
- ``12-fm-active/g2e.npy`` (122272789870,) (model-4, active g2e, 912 GB) (generated in ``12-merge-g2e.py``)

- ``13-uhf-dm1.npy`` (2, 994, 994) UHF 1pdm in LMO (model-4, active) (generated in ``13-uhf.py``)

- ``14-mo-coeff.npy`` (2, 994, 994) LMO (model-4, active) -> UHF (generated in ``14-ump2.py``)
- ``14-ump2-dm1.npy`` (2, 994, 994) UMP2 1pdm in UHF basis (generated in ``14-ump2.py``)

- ``17-uno-occ.npy`` (994, ) UMP2-NO occupation (generated in ``17-mp2-uno-dmao.py``)
- ``17-uno-coeff.npy`` (994, 994) LMO (model-4, active) -> UMP2-NO (generated in ``17-mp2-uno-dmao.py``)
- ``17-mo-ovlp.npy`` (994, ) overlap (norm) between UMP2-NO and LLDUC active space (generated in ``17-mp2-uno-dmao.py``)
- ``17-mo-coeff.npy`` (``18-fm-active/10-full-mo-coeff.npy``) (2818, 2818) AO -> UMP2-NO (spin-restricted, all) (generated in ``17-mp2-uno-dmao.py``)

## [0] Active space (117o, 105e)

- ``18-fm-active/00-ecore.npy`` (model-0, core energy) (generated in ``18-act-ints.py``)
- ``18-fm-active/00-h1e.npy`` (117, 117) (model-0, active h1e) (generated in ``18-act-ints.py``)
- ``18-fm-active/00-g2e.npy`` (23829156,) (model-0, active g2e, 182 MB) (generated in ``18-act-ints.py``)
- ``18-fm-active/00-proja.npy`` (76, 117)  LMO (LLDUC, active) -> MO (model-0, active) (generated in ``20-proj-76.py``)
- ``18-fm-active/00-projc.npy`` (517, 117) LMO (LLDUC, core) -> MO (model-0, active) (generated in ``20-proj-76.py``)

## [1] Active space (180o, 143e)

- ``18-fm-active/01-ecore.npy`` (model-1, core energy) (generated in ``18-act-ints.py``)
- ``18-fm-active/01-h1e.npy`` (180, 180) (model-1, active h1e) (generated in ``18-act-ints.py``)
- ``18-fm-active/01-g2e.npy`` (132690195,) (model-1, active g2e, 1013 MB) (generated in ``18-act-ints.py``)
- ``18-fm-active/01-proja.npy`` (76, 180)  LMO (LLDUC, active) -> MO (model-1, active) (generated in ``20-proj-76.py``)
- ``18-fm-active/01-projc.npy`` (517, 180) LMO (LLDUC, core) -> MO (model-1, active) (generated in ``20-proj-76.py``)

## [2] Active space (285o, 217e)

- ``18-fm-active/02-ecore.npy`` (model-2, core energy) (generated in ``18-act-ints.py``)
- ``18-fm-active/02-h1e.npy`` (285, 285) (model-2, active h1e) (generated in ``18-act-ints.py``)
- ``18-fm-active/02-g2e.npy`` (830505390,) (model-2, active g2e, 6.2 GB) (generated in ``18-act-ints.py``)
- ``18-fm-active/02-proja.npy`` (76, 285)  LMO (LLDUC, active) -> MO (model-2, active) (generated in ``20-proj-76.py``)
- ``18-fm-active/02-projc.npy`` (517, 285) LMO (LLDUC, core) -> MO (model-2, active) (generated in ``20-proj-76.py``)

## [3] Active space (404o, 277e)

- ``18-fm-active/03-ecore.npy`` (model-3, core energy) (generated in ``18-act-ints.py``)
- ``18-fm-active/03-h1e.npy`` (404, 404) (model-3, active h1e) (generated in ``18-act-ints.py``)
- ``18-fm-active/03-g2e.npy`` (3346478955,) (model-3, active g2e, 25 GB) (generated in ``18-act-ints.py``)
- ``18-fm-active/03-proja.npy`` (76, 404)  LMO (LLDUC, active) -> MO (model-3, active) (generated in ``20-proj-76.py``)
- ``18-fm-active/03-projc.npy`` (517, 404) LMO (LLDUC, core) -> MO (model-3, active) (generated in ``20-proj-76.py``)
