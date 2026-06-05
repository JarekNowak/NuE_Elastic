/*=============================================================================
  NuESelection.C
  ─────────────────────────────────────────────────────────────────────────────
  Event selection for ν_μ + e⁻ → ν_μ + e⁻ elastic scattering at SBND.

  Input : merged_nu+eIntimeBNB_DLNuE_22April.root  (ana/NuE + ana/SubRun)

  Selection strategy:
    1. Find the pandora neutrino-candidate slice (category == 1).
    2. Identify a primary shower in that slice with Razzled PDG == 11.
    3. Require no extra significant primaries (single-shower topology).
    4. Apply fiducial volume, track-score, dE/dx, angle, and E×θ² cuts.

  Output ROOT file contains histograms directly usable as input to
  WeinbergAngleFit_SBND.C:
    h_data  – selected events vs reco electron kinetic energy [GeV]
    h_bkg   – background events vs reco electron kinetic energy [GeV]
    h_eff   – efficiency vs true electron kinetic energy [GeV]
    h_smear – smearing matrix M[reco][true] [GeV × GeV]

  Usage:
    root -l -b -q NuESelection.C
    root -l -b -q 'NuESelection.C("path/to/input.root","output.root")'

  Author : Jarek Nowak
  Date   : 2026-05-19
=============================================================================*/

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THStack.h"
#include "TMath.h"
#include "common/NuE_common.h"

// ──────────────────────── SELECTION CUTS ─────────────────────────────────────

// Thresholds tuned by NuEOptimize.C to maximise the POT-scaled figure of merit
// S/√(S+B) at 1×10²¹ POT (see NuESkim.C / NuEOptimize.C). Because the signal MC
// scales DOWN (×0.0106) and ν-background MC scales UP (×1.656) to reach target
// POT, each raw background event is worth ≈156× a raw signal event in the scaled
// sample — so the optimum favours much higher purity than raw-efficiency tuning.
// The working point is background-MC-statistics-limited (≥50 raw bkg events kept
// for a <15% background prediction error); see comments in NuEOptimize.C.
namespace Cuts {
  // Slice level
  const Double_t SLICE_SCORE_MIN    = -1.5;  // exclude truly invalid slices
  // Shower PID
  const Double_t RAZZLE_ELEC_MIN   =  0.86;  // P(electron) from Razzled (was 0.50)
  const Double_t RAZZLE_GAMMA_MAX  =  0.15;  // P(photon) from Razzled — e/γ rejection (new)
  const Int_t    RAZZLE_BEST_PDG   =  11;    // best-PDG must be electron
  const Double_t TRACK_SCORE_MAX   =  0.425; // shower-like (0=shower, 1=track) (was 0.50)
  // Shower kinematics (energies in GeV, angles in rad, dEdx in MeV/cm)
  const Double_t SHOWER_E_MIN      =  0.115; // GeV  (was 0.030; low-T region is bkg-dominated)
  const Double_t SHOWER_E_MAX      =  1.300; // GeV  (was 1.500)
  const Double_t THETA_MAX         =  0.15;  // rad  (was 0.55; exploits ν-e forward collimation)
  const Double_t DEDX_MIN          =  0.30;  // MeV/cm (was 0.5)
  const Double_t DEDX_MAX          =  3.75;  // MeV/cm (was 3.5)
  // E×θ² ≤ 2mₑ (= 0.001022 GeV) kinematically; allow ~10× for resolution
  const Double_t ETHETA2_MAX       =  0.010; // GeV
  // Multiplicity: at most this many extra significant primaries
  const Double_t EXTRA_PRIM_E_MIN  =  0.020; // GeV — threshold for "significant"
  const Int_t    MAX_EXTRA_PRIM    =  0;     // require single-shower topology
}

// ───────────────────── OUTPUT HISTOGRAM BINNING ───────────────────────────────

namespace Binning {
  const Int_t    N_RECO  = 120;
  const Double_t E_LO    = 0.0, E_HI = 1.2;   // GeV, reco/true T_e range
}

// ──────────────────────── CUT-FLOW LABELS ────────────────────────────────────

enum CutStep {
  kAll = 0,
  kHasNuSlice,
  kFiducialVtx,
  kOneElecCandidate,
  kRazzleScore,
  kRazzleGamma,
  kTrackScore,
  kShowerEnergy,
  kNoExtraPrim,
  kDedx,
  kTheta,
  kETheta2,
  kNCuts
};

const char* CutLabel[kNCuts] = {
  "All events",
  "Has #nu slice (cat=1)",
  "Reco #nu vertex in FV",
  "#geq1 prim. e^{-} cand.",
  "Razzled P(e)>0.86",
  "Razzled P(#gamma)<0.15",
  "TrackScore<0.425",
  "E_{shw} #in [115,1300] MeV",
  "0 extra primaries",
  "dE/dx #in [0.30,3.75] MeV/cm",
  "#theta < 0.15 rad",
  "E#times#theta^{2} < 10 MeV"
};

// Result of the electron-candidate search within one slice.
struct ElecCandidate {
  Int_t    pIdx      = -1;    // index in particle vectors
  Double_t energy    = 0.0;   // GeV
  Double_t theta     = -1.0;  // rad (PCA angle preferred, falls back to particle theta)
  Double_t dedx      = -999.0;// MeV/cm
  Double_t razzle11  = -1.0;  // P(electron) from Razzled
  Double_t razzle22  = -1.0;  // P(photon) from Razzled (e/γ separation)
  Double_t showerLen = -1.0;  // shower length [cm]
  Double_t showerOA  = -1.0;  // shower opening angle [rad]
  Double_t VX        = 0.0;   // shower start vertex [cm]
  Double_t VY        = 0.0;
  Double_t VZ        = 0.0;
  Int_t    nExtraPrim = 0;    // number of extra significant primaries in the slice
  bool     found     = false;
};

// Finds the primary electron shower candidate in a given slice.
// Also counts extra significant primaries (tracks or showers) in the same slice.
// Returns the candidate with the highest Razzled P(electron).
ElecCandidate FindElecCandidate(const BranchVars& b, Double_t sliceID,
                                Double_t pcaAngle) {
  ElecCandidate best;

  for (Int_t p = 0; p < (Int_t)b.reco_particlePDG->size(); ++p) {
    if ((*b.reco_particleSliceID)[p] != sliceID) continue;
    if ((*b.reco_particleIsPrimary)[p]  < 0.5)  continue;
    if ((*b.reco_particleClearCosmic)[p] > 0.5)  continue;

    Double_t E_MeV  = (*b.reco_particleShowerBestPlaneEnergy)[p];
    Double_t E_GeV  = E_MeV / 1000.0;
    Double_t rz11   = (*b.reco_particleRazzledPDG11)[p];
    Double_t rzBest = (*b.reco_particleRazzledBestPDG)[p];
    Double_t trkScr = (*b.reco_particleTrackScore)[p];

    // Count any significant primary as potential extra-primary
    if (Sentinel::valid(E_MeV) && E_GeV > Cuts::EXTRA_PRIM_E_MIN) {
      // Only count it as "extra" if it's not the electron candidate
      // (counted later once best candidate is known)
      ;
    }

    // Electron candidate: Razzled best PDG == 11, valid razzle score
    if (!Sentinel::valid(rz11))    continue;
    if ((Int_t)(rzBest + 0.5) != Cuts::RAZZLE_BEST_PDG) continue;
    if (!Sentinel::valid(E_MeV)) continue;

    if (!best.found || rz11 > best.razzle11) {
      best.found    = true;
      best.pIdx     = p;
      best.energy   = E_GeV;
      best.razzle11 = rz11;
      best.razzle22 = (*b.reco_particleRazzledPDG22)[p];
      best.showerLen= (*b.reco_particleShowerLength)[p];
      best.showerOA = (*b.reco_particleShowerOpenAngle)[p];
      best.dedx     = (*b.reco_particleBestPlanedEdx)[p];

      // Use PCA slice angle if valid, else particle theta
      Double_t pca = pcaAngle;
      best.theta = (Sentinel::valid(pca)) ? pca : (*b.reco_particleTheta)[p];

      best.VX = (*b.reco_particleVX)[p];
      best.VY = (*b.reco_particleVY)[p];
      best.VZ = (*b.reco_particleVZ)[p];
    }
  }

  if (!best.found) return best;

  // Count extra significant primaries (not the selected electron)
  best.nExtraPrim = 0;
  for (Int_t p = 0; p < (Int_t)b.reco_particlePDG->size(); ++p) {
    if (p == best.pIdx) continue;
    if ((*b.reco_particleSliceID)[p] != sliceID) continue;
    if ((*b.reco_particleIsPrimary)[p]  < 0.5)  continue;
    if ((*b.reco_particleClearCosmic)[p] > 0.5)  continue;

    Double_t E_MeV = (*b.reco_particleShowerBestPlaneEnergy)[p];
    if (!Sentinel::valid(E_MeV)) continue;
    if (E_MeV / 1000.0 > Cuts::EXTRA_PRIM_E_MIN)
      ++best.nExtraPrim;
  }

  return best;
}

// ─────────────────────────── MAIN FUNCTION ───────────────────────────────────

void NuESelection(
    const char* inputFile  = "/data/sbnd/NuEElastic/merged_nu+eIntimeBNB_DLNuE_22April.root",
    const char* outputFile = "NuESelection_output.root")
{
  std::cout << "\n"
            << "╔══════════════════════════════════════════════════════════════╗\n"
            << "║  ν_μ + e⁻  Elastic Scattering  —  Event Selection  —  SBND  ║\n"
            << "╚══════════════════════════════════════════════════════════════╝\n\n";

  // ── Open input file ──────────────────────────────────────────────────────
  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERROR] Cannot open input file: " << inputFile << "\n";
    return;
  }
  TTree* tNuE    = (TTree*)fIn->Get("ana/NuE");
  TTree* tSubRun = (TTree*)fIn->Get("ana/SubRun");
  if (!tNuE || !tSubRun) {
    std::cerr << "[ERROR] Trees ana/NuE or ana/SubRun not found.\n";
    fIn->Close(); return;
  }

  // ── Compute total POT from SubRun tree ──────────────────────────────────
  Double_t pot_sr = 0.0;
  tSubRun->SetBranchAddress("pot", &pot_sr);
  Double_t totalPOT = 0.0;
  for (Long64_t sr = 0; sr < tSubRun->GetEntries(); ++sr) {
    tSubRun->GetEntry(sr);
    totalPOT += pot_sr;
  }
  printf("[INFO] Total POT: %.4e\n\n", totalPOT);

  // ── Set up branch addresses ──────────────────────────────────────────────
  BranchVars b;
  SetBranchAddresses(tNuE, b);

  // ── Declare output histograms ────────────────────────────────────────────
  const Int_t   N  = Binning::N_RECO;
  const Double_t lo = Binning::E_LO, hi = Binning::E_HI;

  // Fit-macro inputs
  TH1D* h_data  = new TH1D("h_data",  "Selected events;T_{e}^{reco} [GeV];Events/bin", N,lo,hi);
  TH1D* h_bkg   = new TH1D("h_bkg",   "Background;T_{e}^{reco} [GeV];Events/bin",      N,lo,hi);
  TH1D* h_sig   = new TH1D("h_sig",   "Signal (selected);T_{e}^{reco} [GeV];Events/bin",N,lo,hi);
  TH2D* h_smear = new TH2D("h_smear",
    "Smearing matrix;T_{e}^{reco} [GeV];T_{e}^{true} [GeV]",
    N,lo,hi, N,lo,hi);
  // Denominator for efficiency: all true signal events vs true T_e
  TH1D* h_sig_total = new TH1D("h_sig_total",
    "All true signal;T_{e}^{true} [GeV];Events/bin", N,lo,hi);
  // Numerator for efficiency: SELECTED signal vs TRUE T_e (correct efficiency;
  // previously the reco-binned h_sig was divided by the true-binned h_sig_total).
  TH1D* h_sig_true = new TH1D("h_sig_true",
    "Selected signal;T_{e}^{true} [GeV];Events/bin", N,lo,hi);
  // Efficiency (filled after loop)
  TH1D* h_eff = new TH1D("h_eff",
    "Selection efficiency;T_{e}^{true} [GeV];#epsilon", N,lo,hi);
  h_eff->Sumw2();

  // Diagnostic histograms (signal vs background before and after cuts)
  TH1D* h_theta_sig  = new TH1D("h_theta_sig", "#theta (signal);#theta [rad];",       50,0,0.6);
  TH1D* h_theta_bkg  = new TH1D("h_theta_bkg", "#theta (bkg);#theta [rad];",          50,0,0.6);
  TH1D* h_etheta2_sig= new TH1D("h_etheta2_sig","E#theta^{2} (signal);E#theta^{2} [GeV];",50,0,0.015);
  TH1D* h_etheta2_bkg= new TH1D("h_etheta2_bkg","E#theta^{2} (bkg);E#theta^{2} [GeV];",  50,0,0.015);
  TH1D* h_dedx_sig   = new TH1D("h_dedx_sig",  "dE/dx (signal);dE/dx [MeV/cm];",     50,0,6);
  TH1D* h_dedx_bkg   = new TH1D("h_dedx_bkg",  "dE/dx (bkg);dE/dx [MeV/cm];",        50,0,6);
  TH1D* h_rz11_sig   = new TH1D("h_rz11_sig",  "Razzled P(e) (signal);P(e);",         50,0,1);
  TH1D* h_rz11_bkg   = new TH1D("h_rz11_bkg",  "Razzled P(e) (bkg);P(e);",            50,0,1);
  TH1D* h_rz22_sig   = new TH1D("h_rz22_sig",  "Razzled P(#gamma) (signal);P(#gamma);",50,0,1);
  TH1D* h_rz22_bkg   = new TH1D("h_rz22_bkg",  "Razzled P(#gamma) (bkg);P(#gamma);",   50,0,1);
  TH1D* h_shwOA_sig  = new TH1D("h_shwOA_sig", "Shower open angle (signal);OA [rad];", 50,0,0.6);
  TH1D* h_shwOA_bkg  = new TH1D("h_shwOA_bkg", "Shower open angle (bkg);OA [rad];",    50,0,0.6);
  TH1D* h_sliceScore_sig = new TH1D("h_sliceScore_sig","Slice score (signal);Score;",  50,-1,1);
  TH1D* h_sliceScore_bkg = new TH1D("h_sliceScore_bkg","Slice score (bkg);Score;",    50,-1,1);

  // Cut-flow counters: [cut][0=all, 1=signal]
  Double_t cutflow[kNCuts][2] = {};

  // ── Event loop ───────────────────────────────────────────────────────────
  const Long64_t nEntries = tNuE->GetEntries();
  printf("[INFO] Processing %lld events…\n\n", nEntries);

  for (Long64_t ev = 0; ev < nEntries; ++ev) {
    if (ev % 500000 == 0)
      printf("  … event %lld / %lld (%.0f%%)\n", ev, nEntries, 100.0*ev/nEntries);

    tNuE->GetEntry(ev);

    bool isTrueSignal = (b.nuEScatter == 1);

    // Fill signal denominator for efficiency (all true signal events)
    if (isTrueSignal && b.truth_recoilElectronEnergy &&
        !b.truth_recoilElectronEnergy->empty()) {
      Double_t trueE_GeV = (*b.truth_recoilElectronEnergy)[0] / 1000.0;
      h_sig_total->Fill(trueE_GeV);
    }

    // ── Cut 0: All events ─────────────────────────────────────────────────
    cutflow[kAll][0]++;
    if (isTrueSignal) cutflow[kAll][1]++;

    // ── Cut 1: Find pandora neutrino slice ────────────────────────────────
    Int_t sliceIdx = FindNuSliceIndex(b);
    if (sliceIdx < 0) continue;

    cutflow[kHasNuSlice][0]++;
    if (isTrueSignal) cutflow[kHasNuSlice][1]++;

    Double_t sliceID    = (*b.reco_sliceID)[sliceIdx];
    Double_t sliceScore = (*b.reco_sliceScore)[sliceIdx];
    Double_t sliceInter = (*b.reco_sliceInteraction)[sliceIdx];
    Double_t pcaAngle   = GetPCAAngle(b, sliceID);

    // ── Cut 2: Fiducial vertex ────────────────────────────────────────────
    // Use the GENUINE reconstructed neutrino vertex (reco_neutrinoVX/VY/VZ),
    // not reco_sliceTrueVX which is a truth-matched quantity (unavailable in
    // real data and filled for only ~38% of signal events).
    Double_t vtxX, vtxY, vtxZ;
    if (!GetRecoNuVertex(b, sliceID, vtxX, vtxY, vtxZ)) continue;
    if (!FV::inFV(vtxX, vtxY, vtxZ)) continue;

    cutflow[kFiducialVtx][0]++;
    if (isTrueSignal) cutflow[kFiducialVtx][1]++;

    // ── Cut 3: Find primary electron shower candidate ─────────────────────
    ElecCandidate cand = FindElecCandidate(b, sliceID, pcaAngle);
    if (!cand.found) continue;

    cutflow[kOneElecCandidate][0]++;
    if (isTrueSignal) cutflow[kOneElecCandidate][1]++;

    // ── Cut 4: Razzled P(electron) ────────────────────────────────────────
    if (cand.razzle11 < Cuts::RAZZLE_ELEC_MIN) continue;

    cutflow[kRazzleScore][0]++;
    if (isTrueSignal) cutflow[kRazzleScore][1]++;

    // ── Cut 5: Razzled P(photon) — e/γ separation ─────────────────────────
    if (!Sentinel::valid(cand.razzle22) || cand.razzle22 > Cuts::RAZZLE_GAMMA_MAX) continue;

    cutflow[kRazzleGamma][0]++;
    if (isTrueSignal) cutflow[kRazzleGamma][1]++;

    // ── Cut 6: Track score (shower-like) ──────────────────────────────────
    Double_t trkScore = (*b.reco_particleTrackScore)[cand.pIdx];
    if (!Sentinel::valid(trkScore) || trkScore > Cuts::TRACK_SCORE_MAX) continue;

    cutflow[kTrackScore][0]++;
    if (isTrueSignal) cutflow[kTrackScore][1]++;

    // ── Cut 7: Shower energy ──────────────────────────────────────────────
    if (cand.energy < Cuts::SHOWER_E_MIN || cand.energy > Cuts::SHOWER_E_MAX) continue;

    cutflow[kShowerEnergy][0]++;
    if (isTrueSignal) cutflow[kShowerEnergy][1]++;

    // ── Cut 8: Single-shower topology (no extra significant primaries) ─────
    if (cand.nExtraPrim > Cuts::MAX_EXTRA_PRIM) continue;

    cutflow[kNoExtraPrim][0]++;
    if (isTrueSignal) cutflow[kNoExtraPrim][1]++;

    // Fill pre-dEdx diagnostics (at this stage for all passing events)
    if (Sentinel::validDedx(cand.dedx)) {
      if (isTrueSignal) h_dedx_sig->Fill(cand.dedx);
      else              h_dedx_bkg->Fill(cand.dedx);
    }
    if (isTrueSignal) {
      h_theta_sig->Fill(cand.theta);
      h_rz11_sig->Fill(cand.razzle11);
      h_rz22_sig->Fill(cand.razzle22);
      h_shwOA_sig->Fill(cand.showerOA);
      h_sliceScore_sig->Fill(sliceScore);
    } else {
      h_theta_bkg->Fill(cand.theta);
      h_rz11_bkg->Fill(cand.razzle11);
      h_rz22_bkg->Fill(cand.razzle22);
      h_shwOA_bkg->Fill(cand.showerOA);
      h_sliceScore_bkg->Fill(sliceScore);
    }

    // ── Cut 9: Electron-like dE/dx ────────────────────────────────────────
    if (!Sentinel::validDedx(cand.dedx) ||
        cand.dedx < Cuts::DEDX_MIN || cand.dedx > Cuts::DEDX_MAX) continue;

    cutflow[kDedx][0]++;
    if (isTrueSignal) cutflow[kDedx][1]++;

    // ── Cut 10: Forward angle ─────────────────────────────────────────────
    if (cand.theta < 0 || cand.theta > Cuts::THETA_MAX) continue;

    cutflow[kTheta][0]++;
    if (isTrueSignal) cutflow[kTheta][1]++;

    if (isTrueSignal) h_etheta2_sig->Fill(cand.energy * cand.theta * cand.theta);
    else              h_etheta2_bkg->Fill(cand.energy * cand.theta * cand.theta);

    // ── Cut 11: Kinematic E×θ² constraint ─────────────────────────────────
    Double_t eTheta2 = cand.energy * cand.theta * cand.theta;
    if (eTheta2 > Cuts::ETHETA2_MAX) continue;

    cutflow[kETheta2][0]++;
    if (isTrueSignal) cutflow[kETheta2][1]++;

    // ── Event passes full selection ────────────────────────────────────────
    bool isSelectedSignal = isTrueSignal && (sliceInter == 1098.0);

    h_data->Fill(cand.energy);
    if (isSelectedSignal) {
      h_sig->Fill(cand.energy);
      // Smearing matrix M[reco][true] and true-energy numerator for efficiency
      if (b.truth_recoilElectronEnergy && !b.truth_recoilElectronEnergy->empty()) {
        Double_t trueE_GeV = (*b.truth_recoilElectronEnergy)[0] / 1000.0;
        h_smear->Fill(cand.energy, trueE_GeV);
        h_sig_true->Fill(trueE_GeV);
      }
    } else {
      h_bkg->Fill(cand.energy);
    }
  } // end event loop

  // ── Compute efficiency ───────────────────────────────────────────────────
  // Normalise smearing matrix per true-energy column
  for (Int_t iT = 1; iT <= N; ++iT) {
    Double_t colSum = 0.0;
    for (Int_t iR = 1; iR <= N; ++iR) colSum += h_smear->GetBinContent(iR, iT);
    if (colSum > 0.0)
      for (Int_t iR = 1; iR <= N; ++iR)
        h_smear->SetBinContent(iR, iT, h_smear->GetBinContent(iR, iT) / colSum);
  }

  // Efficiency = selected signal / total signal (both binned in TRUE T_e)
  h_eff->Divide(h_sig_true, h_sig_total, 1.0, 1.0, "B");  // binomial errors

  // ── Print cut-flow table ─────────────────────────────────────────────────
  printf("\n");
  printf("══════════════════════════════════════════════════════════════════\n");
  printf("  CUT-FLOW TABLE\n");
  printf("──────────────────────────────────────────────────────────────────\n");
  printf("  %-34s %10s %10s  %7s\n", "Cut", "All evts", "Signal", "Sig eff");
  printf("──────────────────────────────────────────────────────────────────\n");
  Double_t sig0 = cutflow[kAll][1];
  for (Int_t c = 0; c < kNCuts; ++c) {
    Double_t eff = (sig0 > 0) ? 100.0 * cutflow[c][1] / sig0 : 0.0;
    printf("  %-34s %10.0f %10.0f  %6.2f%%\n",
           CutLabel[c], cutflow[c][0], cutflow[c][1], eff);
  }
  printf("══════════════════════════════════════════════════════════════════\n\n");

  Double_t nSel    = h_data->Integral();
  Double_t nSig    = h_sig->Integral();
  Double_t nBkg    = h_bkg->Integral();
  Double_t purity  = (nSel > 0) ? 100.0 * nSig / nSel : 0.0;
  Double_t sigEff  = (sig0 > 0) ? 100.0 * nSig / sig0 : 0.0;
  printf("  Selected events : %.0f  (signal: %.0f  bkg: %.0f)\n", nSel, nSig, nBkg);
  printf("  Purity          : %.2f%%\n", purity);
  printf("  Signal efficiency: %.2f%%\n", sigEff);
  printf("  Total POT       : %.3e\n\n", totalPOT);

  // ── Expected yields at target exposure ─────────────────────────────────
  // Scale signal and background to 1e21 POT using per-sample factors.
  // Signal and neutrino-background MC have different effective POTs.
  Double_t nSig_sc  = nSig  * ScaleFactors::SIGNAL;
  Double_t nBkg_sc  = nBkg  * ScaleFactors::NU_BKG;
  Double_t nTot_sc  = nSig_sc + nBkg_sc;
  Double_t purity_sc = (nTot_sc > 0) ? 100.0 * nSig_sc / nTot_sc : 0.0;
  printf("  ── Scaled to %.0e POT ──────────────────────────────────────\n",
         ScaleFactors::TARGET_POT);
  printf("  Signal events   : %7.1f  (SF = %.7f)\n", nSig_sc,  ScaleFactors::SIGNAL);
  printf("  Nu bkg events   : %7.1f  (SF = %.5f)\n", nBkg_sc,  ScaleFactors::NU_BKG);
  printf("  Cosmic bkg      :    —    (separate sample, SF = %.5f)\n",
         ScaleFactors::COSMIC_BKG);
  printf("  Total (nu only) : %7.1f\n", nTot_sc);
  printf("  Purity (nu bkg) : %.2f%%\n", purity_sc);
  // POT-scaled figure of merit (the quantity NuEOptimize.C maximises).
  Double_t fom_sc = (nTot_sc > 0) ? nSig_sc / std::sqrt(nTot_sc) : 0.0;
  printf("  FOM S/#sqrt(S+B): %.3f   (POT-scaled; rawS=%.0f rawB=%.0f)\n\n",
         fom_sc, nSig, nBkg);

  // ── Save histograms ──────────────────────────────────────────────────────
  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "[ERROR] Cannot create output file: " << outputFile << "\n";
    fIn->Close(); return;
  }

  // Fit-macro inputs
  h_data ->Write();
  h_bkg  ->Write();
  h_eff  ->Write();
  h_smear->Write();
  h_sig  ->Write();
  h_sig_true->Write();
  h_sig_total->Write();

  // Diagnostic histograms
  h_theta_sig->Write();   h_theta_bkg->Write();
  h_etheta2_sig->Write(); h_etheta2_bkg->Write();
  h_dedx_sig->Write();    h_dedx_bkg->Write();
  h_rz11_sig->Write();    h_rz11_bkg->Write();
  h_rz22_sig->Write();    h_rz22_bkg->Write();
  h_shwOA_sig->Write();   h_shwOA_bkg->Write();
  h_sliceScore_sig->Write(); h_sliceScore_bkg->Write();

  fOut->Close();
  printf("[INFO] Histograms saved to %s\n\n", outputFile);

  // ── Quick summary plot ───────────────────────────────────────────────────
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

  TCanvas* c1 = new TCanvas("c_sel","NuE Selection",1600,800);
  c1->Divide(4,2);

  auto drawSigBkg = [&](Int_t pad, TH1D* hs, TH1D* hb, const char* title) {
    c1->cd(pad);
    hb->SetFillColorAlpha(kAzure-4, 0.5);
    hb->SetLineColor(kAzure+1);
    hs->SetFillColorAlpha(kRed-4, 0.5);
    hs->SetLineColor(kRed+1);
    THStack* stk = new THStack(Form("stk%d",pad), title);
    stk->Add(hb); stk->Add(hs);
    stk->Draw("HIST");
    stk->GetXaxis()->SetTitle(hs->GetXaxis()->GetTitle());
    stk->GetYaxis()->SetTitle("Events");
    TLegend* lg = new TLegend(0.55,0.72,0.90,0.90);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
    lg->AddEntry(hs,"Signal","f"); lg->AddEntry(hb,"Background","f");
    lg->Draw();
  };

  drawSigBkg(1, h_sig,      h_bkg,      "Selected events");
  drawSigBkg(2, h_rz11_sig, h_rz11_bkg, "Razzled P(e)");
  drawSigBkg(3, h_rz22_sig, h_rz22_bkg, "Razzled P(#gamma)");
  drawSigBkg(4, h_theta_sig,h_theta_bkg,"#theta [rad]");
  drawSigBkg(5, h_etheta2_sig,h_etheta2_bkg,"E#times#theta^{2} [GeV]");
  drawSigBkg(6, h_dedx_sig, h_dedx_bkg, "dE/dx [MeV/cm]");
  drawSigBkg(7, h_shwOA_sig, h_shwOA_bkg,"Shower open angle [rad]");

  c1->cd(8);
  h_eff->SetLineColor(kBlue+1); h_eff->SetLineWidth(2);
  h_eff->SetMarkerStyle(20);
  h_eff->Draw("E1");
  TLatex lt; lt.SetNDC(); lt.SetTextSize(0.045);
  lt.DrawLatex(0.15,0.92,"Selection efficiency vs T_{e}^{true} [GeV]");

  c1->SaveAs("NuESelection.pdf");
  c1->SaveAs("NuESelection.png");
  printf("[INFO] Summary plots saved: NuESelection.{pdf,png}\n\n");

  fIn->Close();
}
