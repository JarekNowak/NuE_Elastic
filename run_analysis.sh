#!/bin/bash
# run_analysis.sh
# Full analysis chain for ν_μ-e elastic scattering at SBND
#
# Usage (run from NuE_Elastic/):
#   bash run_analysis.sh [--weinberg | --nmm | --all] [OPTIONS] [input.root]
#
# Modes:
#   --weinberg   Run NuE selection → Weinberg angle fit         [default]
#   --nmm        Run NMM selection → NMM fit (delegates to magnetic_moment/run_nmm_analysis.sh)
#   --all        Run both analyses sequentially
#   --demo       Skip selection; run fit in placeholder/demo mode (no input MC needed)
#
# Weinberg options (only with --weinberg or --all):
#   --rebin=N    Rebin fit inputs by factor N (default: 1 = no rebin)
#
# NMM options (only with --nmm or --all):
#   --box        Box-cut selection for NMM [default]
#   --mva        MVA-based selection for NMM
#   --train      Retrain BDTG classifier before MVA selection
#
# Positional:
#   input.root   Merged MC file [default: /data/sbnd/NuEElastic/merged_nu+eIntimeBNB_DLNuE_22April.root]
#
# Examples:
#   bash run_analysis.sh --demo                          # placeholder fit, no MC needed
#   bash run_analysis.sh --weinberg input.root           # full Weinberg chain
#   bash run_analysis.sh --weinberg --rebin=2 input.root # Weinberg with 2× rebinning
#   bash run_analysis.sh --nmm --mva input.root          # NMM with MVA selection
#   bash run_analysis.sh --all input.root                # both analyses

set -e

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT="/data/sbnd/NuEElastic/merged_nu+eIntimeBNB_DLNuE_22April.root"
MODE_W=0        # run Weinberg
MODE_N=0        # run NMM
DEMO=0
REBIN=1
NMM_MODE="box"
NMM_TRAIN=0

# ── Parse arguments ───────────────────────────────────────────────────────────
for arg in "$@"; do
  case "$arg" in
    --weinberg) MODE_W=1 ;;
    --nmm)      MODE_N=1 ;;
    --all)      MODE_W=1; MODE_N=1 ;;
    --demo)     DEMO=1 ;;
    --box)      NMM_MODE="box" ;;
    --mva)      NMM_MODE="mva" ;;
    --train)    NMM_TRAIN=1 ;;
    --rebin=*)  REBIN="${arg#--rebin=}" ;;
    *.root)     INPUT="$arg" ;;
    *)
      echo "Unknown option: $arg"
      echo "Usage: $0 [--weinberg|--nmm|--all] [--demo] [--rebin=N] [--box|--mva] [--train] [input.root]"
      exit 1
      ;;
  esac
done

# Default: --weinberg if neither mode flag given
if [ "$MODE_W" = "0" ] && [ "$MODE_N" = "0" ]; then
  MODE_W=1
fi

# ── Helpers ───────────────────────────────────────────────────────────────────
FLUX_FILE="sbnd_flux_converted.root"
W_SEL_OUT="NuESelection_output.root"
W_FIT_IN="weinberg/WeinbergFit_inputs.root"

banner() {
  echo ""
  echo "══════════════════════════════════════════════════════════════════"
  echo "  $1"
  echo "══════════════════════════════════════════════════════════════════"
}

# ── Flux histogram (shared by both analyses) ──────────────────────────────────
ensure_flux() {
  if [ ! -f "$FLUX_FILE" ]; then
    banner "Generating BNB flux histogram"
    root -l -b -q 'weinberg/get_sbnd_flux.C'
    echo "[OK] $FLUX_FILE created"
  else
    echo "[SKIP] Flux: $FLUX_FILE already exists"
  fi
}

# ── Weinberg angle analysis ───────────────────────────────────────────────────
run_weinberg() {
  banner "Weinberg Angle Analysis  (sin²θ_W fit)"

  if [ "$DEMO" = "1" ]; then
    echo "  Mode: demo (placeholder histograms — no MC input)"
    root -l -b -q 'weinberg/WeinbergAngleFit_SBND.C'
    return
  fi

  echo "  Input file : $INPUT"

  if [ ! -f "$INPUT" ]; then
    echo "[ERROR] Input file not found: $INPUT"
    exit 1
  fi

  ensure_flux

  # Step 1: Event selection
  banner "Step 1: NuE event selection (~25 min)"
  root -l -b -q "NuESelection.C+(\"${INPUT}\",\"${W_SEL_OUT}\")"
  echo "[OK] Selection output: $W_SEL_OUT"

  # Step 2: Rebin (optional)
  if [ "$REBIN" -gt 1 ]; then
    banner "Step 2: Rebinning by factor $REBIN"
    root -l -b -q "weinberg/rebin_histograms.C(\"${W_SEL_OUT}\",\"${W_SEL_OUT%.root}_rebin${REBIN}.root\")"
    W_SEL_OUT="${W_SEL_OUT%.root}_rebin${REBIN}.root"
    echo "[OK] Rebinned: $W_SEL_OUT"
  fi

  # Step 3: Merge flux into fit input file
  banner "Step 3: Merging flux → fit input"
  _TMP_MERGE=$(mktemp /tmp/w_merge_XXXXXX.C)
  cat > "$_TMP_MERGE" << EOF
void $(basename "$_TMP_MERGE" .C)() {
TFile* fSel  = TFile::Open("${W_SEL_OUT}", "READ");
if (!fSel || fSel->IsZombie()) { printf("[ERROR] Cannot open %s\n","${W_SEL_OUT}"); exit(1); }
TFile* fFlux = TFile::Open("${FLUX_FILE}", "READ");
if (!fFlux || fFlux->IsZombie()) { printf("[ERROR] Cannot open %s\n","${FLUX_FILE}"); exit(1); }
TH1* hFlux = (TH1*)fFlux->Get("h_flux");
if (!hFlux) { printf("[ERROR] h_flux not found in %s\n","${FLUX_FILE}"); exit(1); }
TFile* fOut = TFile::Open("${W_FIT_IN}", "RECREATE");
fSel->Get("h_data") ->Write();
fSel->Get("h_bkg")  ->Write();
fSel->Get("h_eff")  ->Write();
fSel->Get("h_smear")->Write();
hFlux->Write();
fOut->Close(); fSel->Close(); fFlux->Close();
printf("[OK] Written: ${W_FIT_IN}\n");
}
EOF
  root -l -b -q "$_TMP_MERGE"
  rm -f "$_TMP_MERGE"

  # Step 4: Weinberg fit
  banner "Step 4: Weinberg angle fit"
  root -l -b -q "weinberg/WeinbergAngleFit_SBND.C(\"${W_FIT_IN}\")"

  # Step 5: Print result
  banner "Weinberg Result"
  _TMP_PRINT=$(mktemp /tmp/w_print_XXXXXX.C)
  _TMP_FUNC=$(basename "$_TMP_PRINT" .C)
  cat > "$_TMP_PRINT" << EOF
void ${_TMP_FUNC}() {
TFile* f = TFile::Open("WeinbergAngleFit_SBND_results.root", "READ");
if (!f || f->IsZombie()) { printf("[ERROR] Cannot open results file\n"); exit(1); }
TVectorD* v = (TVectorD*)f->Get("fit_results");
if (!v) { printf("[ERROR] fit_results not found\n"); exit(1); }
printf("\n");
printf("  sin2thW best-fit  : %.5f\n",           (*v)[0]);
printf("  HESSE error       : %.5f\n",           (*v)[1]);
printf("  chi2 / ndf        : %.1f / %.0f\n",    (*v)[2], (*v)[3]);
printf("  MINOS eplus       : %.5f\n",           (*v)[4]);
printf("  MINOS eminus      : %.5f\n",           (*v)[5]);
printf("\n");
f->Close();
}
EOF
  root -l -b -q "$_TMP_PRINT"
  rm -f "$_TMP_PRINT"

  echo "  Fit results : WeinbergAngleFit_SBND_results.root"
  echo "  Plots       : WeinbergAngleFit_SBND.{pdf,png}"
  echo "                WeinbergAngleFit_SBND_scan.{pdf,png}"
}

# ── NMM analysis ──────────────────────────────────────────────────────────────
run_nmm() {
  banner "NMM Analysis  (mode: ${NMM_MODE})"

  if [ "$DEMO" = "1" ]; then
    echo "  Mode: demo (placeholder histograms — no MC input)"
    root -l -b -q 'magnetic_moment/NuMMFit_SBND.C'
    return
  fi

  NMM_ARGS="--${NMM_MODE}"
  [ "$NMM_TRAIN" = "1" ] && NMM_ARGS="$NMM_ARGS --train"
  NMM_ARGS="$NMM_ARGS ${INPUT}"
  bash magnetic_moment/run_nmm_analysis.sh $NMM_ARGS
}

# ── Main ──────────────────────────────────────────────────────────────────────
[ "$MODE_W" = "1" ] && run_weinberg
[ "$MODE_N" = "1" ] && run_nmm

banner "Done"
