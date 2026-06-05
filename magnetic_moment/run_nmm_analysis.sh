#!/bin/bash
# run_nmm_analysis.sh
# Full NMM analysis chain: selection → flux merge → fit → 90% CL upper limit on μ_ν
#
# Usage:
#   ./magnetic_moment/run_nmm_analysis.sh [--box | --mva] [--train] [input.root]
#
# Options:
#   --box      Use box-cut selection (NuMMSelection.C)          [default]
#   --mva      Use MVA selection (NuMMSelection_MVA.C)
#   --train    Retrain the BDTG classifier before applying MVA
#   input.root Path to merged MC ROOT file
#              [default: /data/sbnd/NuEElastic/merged_nu+eIntimeBNB_DLNuE_22April.root]
#
# Run from the NuE_Elastic/ directory:
#   cd /home/t2k/nowak/NuE_Elastic
#   bash magnetic_moment/run_nmm_analysis.sh --mva

set -e  # exit on first error

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT="/data/sbnd/NuEElastic/merged_nu+eIntimeBNB_DLNuE_22April.root"
MODE="box"
TRAIN=0

# ── Parse arguments ───────────────────────────────────────────────────────────
for arg in "$@"; do
  case "$arg" in
    --box)   MODE="box"  ;;
    --mva)   MODE="mva"  ;;
    --train) TRAIN=1     ;;
    *.root)  INPUT="$arg" ;;
    *)
      echo "Unknown option: $arg"
      echo "Usage: $0 [--box | --mva] [--train] [input.root]"
      exit 1
      ;;
  esac
done

# ── Derived paths ─────────────────────────────────────────────────────────────
OUTDIR="magnetic_moment"
FLUX_FILE="sbnd_flux_converted.root"
WEIGHTS_FILE="${OUTDIR}/tmva_dataset/weights/NuMMClassifier_BDTG.weights.xml"

if [ "$MODE" = "box" ]; then
  SEL_OUTPUT="${OUTDIR}/NuMMSelection_output.root"
  FIT_INPUT="${OUTDIR}/NuMMFit_inputs.root"
else
  SEL_OUTPUT="${OUTDIR}/NuMMSelection_MVA_output.root"
  FIT_INPUT="${OUTDIR}/NuMMFit_inputs_MVA.root"
fi

FIT_RESULTS="${OUTDIR}/NuMMFit_SBND_results.root"

# ── Helper: print a section header ───────────────────────────────────────────
banner() {
  echo ""
  echo "══════════════════════════════════════════════════════════════════"
  echo "  $1"
  echo "══════════════════════════════════════════════════════════════════"
}

# ── Check prerequisites ───────────────────────────────────────────────────────
banner "NuMM Analysis Chain  (mode: ${MODE})"
echo "  Input file : $INPUT"
echo "  Flux file  : $FLUX_FILE"
echo "  Sel output : $SEL_OUTPUT"
echo "  Fit input  : $FIT_INPUT"
echo ""

if [ ! -f "$INPUT" ]; then
  echo "[ERROR] Input file not found: $INPUT"
  exit 1
fi

# ── Step 0: Flux (generate if missing) ───────────────────────────────────────
if [ ! -f "$FLUX_FILE" ]; then
  banner "Step 0: Generating BNB flux histogram"
  root -l -b -q 'weinberg/get_sbnd_flux.C'
  echo "[OK] $FLUX_FILE created"
else
  echo "[SKIP] Step 0: $FLUX_FILE already exists"
fi

# ── Step 1: Train classifier (MVA mode only, if requested or weights missing) ─
if [ "$MODE" = "mva" ]; then
  if [ "$TRAIN" = "1" ] || [ ! -f "$WEIGHTS_FILE" ]; then
    banner "Step 1a: Training BDTG classifier"
    root -l -b -q "magnetic_moment/NuMMClassifier.C+(\"${INPUT}\")"
    echo "[OK] Weights saved to $WEIGHTS_FILE"
  else
    echo "[SKIP] Step 1a: Using existing weights: $WEIGHTS_FILE"
  fi
fi

# ── Step 2: Event selection ───────────────────────────────────────────────────
if [ "$MODE" = "box" ]; then
  banner "Step 1: Box-cut selection (~25 min)"
  root -l -b -q "magnetic_moment/NuMMSelection.C+(\"${INPUT}\",\"${SEL_OUTPUT}\")"
else
  banner "Step 1b: MVA selection (~25 min)"
  root -l -b -q "magnetic_moment/NuMMSelection_MVA.C+(\"${INPUT}\",\"${SEL_OUTPUT}\")"
fi
echo "[OK] Selection output: $SEL_OUTPUT"

# ── Step 3: Merge flux into fit input file ────────────────────────────────────
banner "Step 2: Merging flux into fit input"

_TMP_MERGE=$(mktemp /tmp/nmm_merge_XXXXXX.C)
cat > "$_TMP_MERGE" << EOF
void $(basename "$_TMP_MERGE" .C)() {
TFile* fSel  = TFile::Open("${SEL_OUTPUT}", "READ");
if (!fSel || fSel->IsZombie()) { printf("[ERROR] Cannot open %s\n", "${SEL_OUTPUT}"); exit(1); }
TFile* fFlux = TFile::Open("${FLUX_FILE}", "READ");
if (!fFlux || fFlux->IsZombie()) { printf("[ERROR] Cannot open %s\n", "${FLUX_FILE}"); exit(1); }
TH1* hFlux = (TH1*)fFlux->Get("h_flux");
if (!hFlux) { printf("[ERROR] h_flux not found in %s\n", "${FLUX_FILE}"); exit(1); }
TFile* fOut = TFile::Open("${FIT_INPUT}", "RECREATE");
fSel->Get("h_data") ->Write();
fSel->Get("h_bkg")  ->Write();
fSel->Get("h_eff")  ->Write();
fSel->Get("h_smear")->Write();
hFlux->Write();
fOut->Close(); fSel->Close(); fFlux->Close();
printf("[OK] Written: ${FIT_INPUT}\n");
}
EOF
root -l -b -q "$_TMP_MERGE"
rm -f "$_TMP_MERGE"

# ── Step 4: NMM fit ───────────────────────────────────────────────────────────
banner "Step 3: NMM fit"
root -l -b -q "magnetic_moment/NuMMFit_SBND.C(\"${FIT_INPUT}\")"
echo "[OK] Fit results: $FIT_RESULTS"

# ── Step 5: Print result ──────────────────────────────────────────────────────
banner "Result"
_TMP_PRINT=$(mktemp /tmp/nmm_print_XXXXXX.C)
cat > "$_TMP_PRINT" << 'ROOTEOF'
void FUNCNAME() {
TFile* f = TFile::Open("magnetic_moment/NuMMFit_SBND_results.root", "READ");
if (!f || f->IsZombie()) { printf("[ERROR] Cannot open results file\n"); exit(1); }
TVectorD* v = (TVectorD*)f->Get("fit_results");
if (!v) { printf("[ERROR] fit_results not found\n"); exit(1); }
printf("\n");
printf("  mu_nu best-fit     : %.3f x 10^-10 mu_B\n",  (*v)[0]);
printf("  mu_nu uncertainty  : %.3f x 10^-10 mu_B\n",  (*v)[1]);
printf("  mu_nu 90%% CL UL   : %.3f x 10^-10 mu_B\n",  (*v)[2]);
printf("  chi2 / ndf         : %.1f / %.0f\n",           (*v)[3], (*v)[4]);
printf("\n");
f->Close();
}
ROOTEOF
# Patch function name to match the filename ROOT requires
_TMP_FUNC=$(basename "$_TMP_PRINT" .C)
sed -i "s/FUNCNAME/${_TMP_FUNC}/" "$_TMP_PRINT"
root -l -b -q "$_TMP_PRINT"
rm -f "$_TMP_PRINT"

banner "Done"
echo "  Output plots : ${OUTDIR}/NuMMFit_SBND.{pdf,png}"
echo "                 ${OUTDIR}/NuMMFit_SBND_scan.{pdf,png}"
echo "  Results file : $FIT_RESULTS"
echo ""
