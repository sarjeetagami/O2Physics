# rho00 MC closure for kstarpbpb.cxx: notes

Reviewer request: full closure with injected rho00 = 1/3, and also 0.25 and 0.40. The complete chain must recover the
injected value: reconstruction, efficiency correction, fitting and event-plane resolution correction.

## Current wagon setup (kstar_nene_mc_closure / kstar_oo_mc_closure)

- **Base wagon (MC true Psi, no smearing):**
  - processMCClosure, processMCClosureME, processMCClosurePhi, processMCClosureMEPhi = 1
  - all other process* = 0
- **Subwagon "closure" (reconstructed Psi_FT0C):**
  - processMCClosureRecoEP, processMCClosureMERecoEP, processMCClosureRecoEPPhi, processMCClosureMERecoEPPhi = 1
  - the four true-Psi switches = 0
  - needs the Q-vector task and epvector in the dependencies
- **Common settings:**
  - doClosureTest = 0, cfgSAFrame = 0, cfgEPNormalHarmonic = 2
  - closure.injectPolarization = 1, closure.rho00Inj = 0.333333 (then 0.25, 0.40)
  - closure.useTrueEP = 1, closure.splitSample = 1, closure.meRequireSingleReco = 1
  - closure.genAcceptanceCut = 0, closure.requirePrimaryTracks = 0
  - closure.fillRotBkg = 1, with nBkgRotations = 4, confMinRot = 2.7925, confMaxRot = 3.4907
  - cfgNoMixedEvents = 5
- **Shared settings:** rotation and mixing settings (nBkgRotations, confMinRot/confMaxRot, cfgNoMixedEvents,
  axisVertex, axisMultiplicityClass) are the same for K*0 and phi.
- **Efficiency:** processMC is not needed. The closure builds its own Acc x Eff from the even events
  (hSparseRecTrueEff / hSparseGenEff).
- **Resolution correction:** rho00 = 1/3 + (rho_obs - 1/3) * 4 / (1 + 3R), with R = <cos(2 n dPsi)>.
  - n = 1: R = R2
  - n = 2: R = R4{Psi2}, obtained from R2 through chi (Ollitrault)

Analysis macros (in this folder): MCClosureTrueEP.C (base wagon), MCClosureRecoEP.C (subwagon).

## Finding: the reconstructed-EP closure cannot work with the current MC (wagon change alone is not enough)

The reconstructed FT0C plane in this MC carries no information about the true plane:

- true resolution <cos 2(Psi_FT0C - Psi_true)> = -0.0004
- sub-event resolution = 0.084

No wagon setting changes that. With harmonic 1 or 2, the injected alignment is still diluted by the maximum factor
of 4, and the sub-event resolution (0.084) does not describe that dilution.

| Option | What changes | Works? |
|---|---|---|
| Same MC, harmonic set to 1 | Wagon only | No. True resolution is still 0, so the corrected result is biased (0.267 expected for 0.25). |
| An MC production with flow | Wagon only (different dataset) | Yes, if such a production exists for Ne-Ne. The sub-event method is then tested too. |
| Same MC, plane smeared with a known resolution | Code + wagon | Yes. Recommended. |

### Recommended code change (not yet implemented)

- **New option in the `closure` configurable group**, used by the four RecoEP process functions: instead of the
  FT0C plane, use the true plane plus a random angle. The spread of the random angle gives the target
  resolution (about 0.28, as on data).
- **Store the exact resolution** of that smearing in a new histogram. The macro uses it in the correction.
- **Keep the injection and the generated truth on the true plane**, as now.
- **Size:** one or two configurables, the smeared angle in the four RecoEP functions, and one histogram.
  MCClosureRecoEP.C needs a matching update to read the stored resolution.

### Wagon settings after the code change

- Switch the smearing option on, with the target resolution (about 0.28).
- Set cfgEPNormalHarmonic to the value the data wagons use, so the closure tests the same axis as the data.
- Keep closure.fillRotBkg on.
- Add a subwagon with rho00 = 0.40 (the reviewer names it), next to 1/3 and 0.25.

### Expected precision

- With R = 0.28 the correction factor drops from 4 to about 2.2, so the corrected errors shrink by almost half.
- From the present sample (pT-averaged, rotational background):
  - phi: about +-0.03, which separates 0.25 from 1/3 at about 2.5 sigma
  - K*0: about +-0.05, which separates them at about 1.5 sigma
- More MC statistics would still help, especially for K*0.

### Limitation

The smearing only tests the correction formula at a known resolution. It does not test whether the sub-event method
gives the right resolution on data; only an MC with flow can do that.
