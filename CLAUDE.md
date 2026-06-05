# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ROOT-based analysis to measure the weak mixing angle (sin²θ_W) from muon-neutrino / electron elastic scattering (ν_μ + e⁻ → ν_μ + e⁻) at the SBND detector using the BNB beam at Fermilab. All analysis code runs as ROOT macros (interpreted C++ with `.C` extension) or via PyROOT.

## Running the Code

**Main fit (with real input histograms):**
```bash
root -l -b -q 'WeinbergAngleFit_SBND.C("inputs.root")'
# Also available with different rebinnings:
root -l -b -q 'WeinbergAngleFit_SBND.C("inputs10.root")'
root -l -b -q 'WeinbergAngleFit_SBND.C("inputs20.root")'
root -l -b -q 'WeinbergAngleFit_SBND.C("inputs40.root")'
```

**Main fit (placeholder/demo mode, no input file needed):**
```bash
root -l -b -q WeinbergAngleFit_SBND.C
```

**Event selection (produces fit inputs):**
```bash
root -l -b -q NuESelection.C        # uses default path in /data/sbnd/NuEElastic/
root -l -b -q 'NuESelection.C("input.root","output.root")'
# Output: NuESelection_output.root, NuESelection.{pdf,png}
# Then use output as input to the fit:
root -l -b -q 'WeinbergAngleFit_SBND.C("NuESelection_output.root")'
```

**Fast cut tuning — skim once, optimise many times** (`NuESkim.C` + `NuEOptimize.C`):
Reading the 143 GB merged tree takes ~25 min, so cut optimisation is split into a
single slow skim pass followed by fast in-memory tuning.
```bash
# Step 1 — ONE pass over the full file → small candidate skim (~25 min):
root -l -b -q 'NuESkim.C+'                       # → NuESkim.root
#   loose preselection: ν-slice + reco-ν-vertex FV + primary e-shower (best-PDG==11, E>10 MeV)
#   stores a flat TTree of all discriminating vars + truth flags, plus the
#   efficiency denominator (h_sig_total), pre-selection cut-flow, and total POT.

# Step 2 — tune cuts on the skim (seconds) maximising POT-scaled S/√(S+B):
root -l -b -q 'NuEOptimize.C+("NuESkim.root","NuESelection_output.root")'
#   coordinate-ascent over every threshold; prints the optimal working point and
#   writes fit-ready histograms (no second pass over the 143 GB file needed).
```
The optimiser FOM is **POT-scaled** S/√(S+B): each MC event is weighted by its
per-sample factor (signal ×0.0106, ν-bkg ×1.656) so it maximises the significance
of the *real expected yields*, not the misleading raw-MC purity. It enforces a
≥50 raw-bkg-event floor so the background prediction stays statistically
meaningful — the optimum is background-MC-statistics-limited.

**Neutrino magnetic moment (NMM) analysis** (`magnetic_moment/`):
```bash
# Box-cut selection optimised for low-T_e (10–600 MeV, 120 × 5 MeV bins):
root -l -b -q magnetic_moment/NuMMSelection.C
root -l -b -q 'magnetic_moment/NuMMSelection.C("input.root","out.root")'
# Output: NuMMSelection_output.root, NuMMSelection_output.{pdf,png}

# MVA-based selection (alternative to box cuts — train first, then apply):
# Step 1 — train TMVA classifiers (BDTG, BDT, MLP); ~25 min on 9M events:
root -l -b -q 'magnetic_moment/NuMMClassifier.C+("input.root")'
# Output: magnetic_moment/TMVA_output.root
#         magnetic_moment/tmva_dataset/weights/NuMMClassifier_BDTG.weights.xml
# Inspect ROC curves and choose MVA_CUT working point:
#   root -l 'magnetic_moment/tmva_dataset/weights/TMVAGui.C'
# Step 2 — apply trained BDTG (edit MVA_CUT in NuMMSelection_MVA.C first):
root -l -b -q 'magnetic_moment/NuMMSelection_MVA.C+("input.root")'
# Output: magnetic_moment/NuMMSelection_MVA_output.root + plots

# NMM fit (placeholder/demo mode — no input needed):
root -l -b -q magnetic_moment/NuMMFit_SBND.C
# NMM fit with real selection output (merge flux first — see below):
root -l -b -q 'magnetic_moment/NuMMFit_SBND.C("NuMMFit_inputs.root")'
# Output: NuMMFit_SBND_results.root, NuMMFit_SBND.{pdf,png},
#         NuMMFit_SBND_scan.{pdf,png}
```

**Cross-section and sensitivity plots:**
```bash
root -l -b -q xsec/NuMuE_Xsec_SBND.C
root -l -b -q xsec/NuMuE_FluxAvg_sin2tW_SBND.C
```

**Utility macros:**
```bash
root -l -b -q get_sbnd_flux.C      # digitise BNB flux → sbnd_flux_converted.root
root -l -b -q rebin_histograms.C   # rebin inputs.root → outputs_reduced.root
python nu_mu_flux.py                # alternative flux digitisation via PyROOT (covers 0–2.5 GeV)
```

**Compile and run with ACLiC (faster repeated execution):**
```bash
root -l -b -q 'WeinbergAngleFit_SBND.C+("inputs.root")'
```

## Architecture

### common/NuE_common.h — shared selection infrastructure

Every macro that reads the merged MC tree (`ana/NuE`) needs the same boilerplate
to do it. That boilerplate lives **once** in `common/NuE_common.h` and is
`#include`d by all of them, rather than being copy-pasted per file:

| Symbol | Purpose |
|--------|---------|
| `Sentinel::` | missing-value sentinels (`-999000`, `-990`) + `valid()`/`validDedx()` |
| `FV::` | SBND fiducial volume bounds + `inFV()` |
| `ScaleFactors::` | per-sample POT normalisation to 1×10²¹ POT |
| `Preselection::SLICE_SCORE_MIN` | slice-score floor used by `FindNuSliceIndex` |
| `BranchVars` | superset of every `ana/NuE` branch any macro reads |
| `SetBranchAddresses` | disables all branches, re-enables only those in `BranchVars` (the read-speed optimisation the 143 GB skim needs; a no-op for correctness since every macro reads the tree only through `BranchVars`) |
| `GetPCAAngle`, `GetRecoNuVertex`, `FindNuSliceIndex` | shared lookups/helpers |

Included by `NuESelection.C`, `NuESkim.C` (`"common/NuE_common.h"`),
`magnetic_moment/NuMMSelection.C` and `magnetic_moment/NuMM_common.h`
(`"../common/NuE_common.h"`; GCC resolves `""` includes relative to the
including file, so this works regardless of CWD).

**What deliberately stays per-macro** (it genuinely differs between analyses, so
do *not* move it into the shared header): the tuned working-point thresholds
(`Cuts`/`Pre`/`PreCuts`), the reco-energy `Binning` (120 bins 0–1.2 GeV for the
Weinberg fit vs 120 bins 0–0.6 GeV for NMM), the cut-flow `enum`/labels, and the
two flavours of `ElecCandidate`/`FindElecCandidate` (box-cut `Double_t` vs MVA
`Float_t`). `magnetic_moment/NuMM_common.h` keeps only these NMM-MVA-specific
pieces and includes `NuE_common.h` for the rest. `NuEOptimize.C` reads the flat
skim tree (`NuESkim.root`), not `ana/NuE`, so it shares none of this.

### Input histogram conventions

The fit macro expects these named histograms inside the input ROOT file:

| Name      | Type  | Description |
|-----------|-------|-------------|
| `h_data`  | TH1D  | Observed reco-energy electron spectrum |
| `h_eff`   | TH1D  | Selection efficiency vs true recoil energy |
| `h_smear` | TH2D  | Smearing matrix M[reco][true], column-normalised (ΣM[:,j]=1) |
| `h_flux`  | TH1D  | BNB ν_μ flux dΦ/dEν [/cm²/POT/GeV] |
| `h_bkg`   | TH1D  | Background prediction (optional) |

If any required histogram is missing or the file cannot be opened, the macro falls back to internally generated placeholder histograms (Gaussian flux, flat efficiency, Gaussian smearing, Poisson-fluctuated pseudo-data at SM).

### Prediction pipeline

`BuildTrueSpectrum(sin2thW)` — integrates dσ/dT × flux × efficiency over true-energy bins using a trapezoidal rule with `Setup::XSEC_NSUB = 8` sub-steps per bin.

`SmearToReco(hTrue)` — applies the smearing matrix to map true → reco space.

`BuildRecoPrediction(sin2thW)` — calls both above, then adds background.

### Fit configuration (in `Setup` namespace)

```cpp
Setup::FIT_BIN_MIN = 10   // first reco bin included in χ²
Setup::FIT_BIN_MAX = 40   // last reco bin included in χ²
Setup::XSEC_NSUB   = 8    // integration substeps per T bin
Setup::N_ELECTRONS = 2.437e31  // fiducial electron targets
Setup::POT         = 1.0e21   // protons on target (3-yr run)
```

Adjust `FIT_BIN_MIN`/`FIT_BIN_MAX` to change the fitted energy range (e.g., to avoid threshold noise or low-statistics tails).

### Minimisation

TMinuit with MIGRAD → HESSE → MINOS sequence. The FCN implements Neyman χ² with data errors (sqrt(obs) for Poisson bins). HESSE gives the parabolic uncertainty; MINOS gives asymmetric profile-likelihood errors — both are saved to `WeinbergAngleFit_SBND_results.root` as a `TVectorD fit_results[6]` = [sin2thW, err_hesse, chi2, ndf, eplus, eminus].

### Output files

| File | Contents |
|------|----------|
| `WeinbergAngleFit_SBND_results.root` | Fit results vector + best-fit reco histogram |
| `WeinbergAngleFit_SBND.pdf/png` | Data vs best-fit with residuals panel |
| `WeinbergAngleFit_SBND_scan.pdf/png` | Δχ² profile scan ±4σ around best-fit |
| `xsec/NuMuE_Xsec_SBND_results.root` | Cross-section TGraphs |
| `xsec/NuMuE_FluxAvg_sin2tW_SBND.root` | Flux-averaged dσ/dT TGraphs |

### Physics constants (SM values used throughout)

```
GF     = 1.16638e-5  GeV⁻²
me     = 0.511e-3    GeV
hbarc2 = 3.8938e-28  cm²·GeV²
sin²θ_W(SM) = 0.2312
gL = −½ + sin²θ_W
gR =      sin²θ_W
```

### NuESelection.C

Selects ν_μ + e⁻ elastic scatter events from the SBND merged MC tree (`ana/NuE`). Runs over 9.1M events in ~25 min on a single core.

**Input tree branches** (key ones):
| Branch | Type | Description |
|--------|------|-------------|
| `nuEScatter` | `Int_t` | 1 = true ν-e elastic event |
| `truth_recoilElectronEnergy` | `vector<double>` | True recoil electron energy in **MeV** |
| `reco_sliceCategory` | `vector<double>` | 1 = pandora neutrino-candidate slice |
| `reco_sliceScore` | `vector<double>` | Pandora NuScore (−1 to 1; −999999 = invalid) |
| `reco_sliceInteraction` | `vector<double>` | Truth interaction code; 1098 = ν-e elastic |
| `reco_particleIsPrimary` | `vector<double>` | 1 = primary particle in neutrino vertex |
| `reco_particleRazzledPDG11` | `vector<double>` | Razzled P(electron) |
| `reco_particleRazzledPDG22` | `vector<double>` | Razzled P(photon) — used for e/γ rejection |
| `reco_particleRazzledBestPDG` | `vector<double>` | Best Razzled PDG (11/13/22/211/2212) |
| `reco_particleShowerBestPlaneEnergy` | `vector<double>` | Shower energy in **MeV** |
| `reco_particleShowerLength` | `vector<double>` | Shower length [cm] (captured in skim) |
| `reco_particleShowerOpenAngle` | `vector<double>` | Shower opening angle [rad] |
| `reco_particleTheta` | `vector<double>` | Polar angle from beam axis [rad] |
| `reco_particleBestPlanedEdx` | `vector<double>` | dE/dx at shower start [MeV/cm]; −999 = invalid |
| `reco_particlePlane{0,1,2}dEdx` | `vector<double>` | Per-plane dE/dx [MeV/cm] |
| `reco_neutrinoVX/VY/VZ` | `vector<double>` | **Reco** neutrino vertex [cm] — used for the FV cut |
| `reco_neutrinoSliceID` | `vector<double>` | Slice ID each reco neutrino belongs to |
| `angleRecalculationPCASlice_angle` | `vector<double>` | PCA-recalculated angle per slice [rad] |

Sentinel values: −999999 for most missing fields; −999 for dEdx.

**Cut sequence and achieved performance** (on merged MC, 9.508×10²² POT):

| Cut | Signal eff. |
|-----|------------|
| Has pandora ν-slice (cat=1) | 81.6% |
| Vertex in FV | 63.1% |
| ≥1 primary e⁻ candidate | 49.6% |
| Razzled P(e) > 0.50 | 46.1% |
| TrackScore < 0.50 | 44.2% |
| Shower energy 30–1500 MeV | 41.7% |
| 0 extra primaries | 38.4% |
| dE/dx 0.5–3.5 MeV/cm | 35.4% |
| θ < 0.55 rad | 34.5% |
| E×θ² < 10 MeV | **31.1%** |

Final: **19,356 selected events, 96.5% purity, 30.6% signal efficiency**.

**Output histograms** (fit-macro compatible, energy axis in GeV):
- `h_data` — selected events (signal+bkg), 40 bins, 0–1.2 GeV
- `h_bkg` — background subset
- `h_eff` — efficiency vs true T_e
- `h_smear` — column-normalised smearing matrix [reco][true]

Cuts are in the `Cuts::` namespace at the top of the file; FV in `FV::`.

### magnetic_moment/

Contains two macros for the NMM measurement:

**`NuMMSelection.C`** — same structure as `NuESelection.C` but optimised for NMM:
- `SHOWER_E_MIN = 0.010 GeV` (10 MeV, down from 30 MeV) to access the 1/T NMM excess
- `SHOWER_E_MAX = 0.600 GeV` and `N_RECO = 120` bins (0–0.6 GeV, 5 MeV/bin)
- Tighter `RAZZLE_ELEC_MIN = 0.70` (vs 0.50 in NuESelection)
- Tighter `DEDX_MIN = 0.8` MeV/cm (vs 0.5), `DEDX_MAX = 3.0` MeV/cm (vs 3.5)
- Tighter `ETHETA2_MAX = 0.003` GeV ≈ 3 MeV (vs 0.010 GeV); exploits ν-e forward collimation
- `TrackScore < 0.50` kept same as NuESelection — tightening to 0.25 catastrophically cuts low-E electrons (10–100 MeV showers are compact and score more track-like)
- Extra diagnostic histogram `h_lowT_sig/bkg` (100 bins × 1 MeV, 0–100 MeV)

**Achieved performance** (merged MC, 9.508×10²² POT):

| Cut | Signal eff. |
|-----|------------|
| Has ν-slice | 81.65% |
| Vertex in FV | 63.08% |
| ≥1 prim. e⁻ cand. | 49.57% |
| Razzled P(e) > 0.70 | 35.38% |
| TrackScore < 0.50 | 34.32% |
| E_shw ∈ [10,600] MeV | 24.01% |
| 0 extra primaries | 22.12% |
| dE/dx ∈ [0.8,3.0] MeV/cm | 20.09% |
| θ < 0.55 rad | 19.51% |
| E×θ² < 3 MeV | **11.88%** |

Final: **7,360 selected events, 98.5% raw purity, 11.9% signal efficiency**.  
Scaled to 1×10²¹ POT: 76.7 signal + 187.1 ν-bkg = 263.8 total, **29.1% purity**.  
The scaled purity is lower than the raw because the signal-enhanced MC (9.508×10²² POT) scales down by ×0.0106 while ν-bkg MC scales up by ×1.656.

**`NuMMFit_SBND.C`** — fits μ_ν (neutrino magnetic moment) using:

NMM cross-section formula (Vogel & Engel 1989):
```
dσ_NMM/dT = (πα²ℏ²c²/m_e²) × (μ_ν/μ_B)² × (1/T − 1/E_ν)   [cm²/GeV]
           = NMM_SIGMA0 × μ_par² × (1/T − 1/E_ν)
```
where `NMM_SIGMA0 = πα²ℏ²c²/m_e² × (10⁻¹⁰)² ≈ 2.495e-45 cm²` and `μ_par` is in units of 10⁻¹⁰ μ_B.

Total prediction = P_SM (fixed sin²θ_W = 0.2312) + μ_par² × P_NMM_shape + P_bkg.

P_SM and P_NMM_shape are **precomputed once** (`g_hPredSM`, `g_hNMMshape`) before TMinuit runs, making the FCN fast (only scales `g_hNMMshape` by μ_par²).

Single free parameter: `mu_nu` in units of 10⁻¹⁰ μ_B, bounded [0, 20].
90% CL upper limit found by scanning Δχ² until Δχ² = 2.71 (one-sided).

Expected SBND sensitivity (placeholder mode, 1×10²¹ POT): μ_ν < ~2 × 10⁻¹⁰ μ_B.

Fit results saved to `NuMMFit_SBND_results.root` as `TVectorD fit_results[5]` = [μ_bf, μ_err, μ_ul90, χ², ndf].

### rebin_histograms.C

The `rebin` variable controls all four histograms (data, eff, bkg, smear). The smear must be rebinned by the same factor in both X and Y — it uses `Rebin2D(rebin, rebin)`. After rebinning, efficiency is rescaled by `1/rebin` to keep it a per-bin fraction rather than a sum.

### archive/

Contains `WeinbergAngleFit_SBNDold.C` — the pre-fix version of the main macro, kept for reference. Do not run it; it contains the bugs documented in the current macro's `FIX` comments.

### Global state

`WeinbergAngleFit_SBND.C` uses five global `TH*` pointers (`g_hData`, `g_hEff`, `g_hBkg`, `g_hSmear`, `g_hFlux`) required by the TMinuit FCN callback — this is a ROOT pattern, not a design choice.
