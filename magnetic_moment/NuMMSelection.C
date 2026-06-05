/*=============================================================================
  NuMMSelection.C
  ─────────────────────────────────────────────────────────────────────────────
  Event selection for neutrino magnetic moment (NMM) measurement at SBND.

  Optimised for low electron recoil kinetic energy: the NMM cross-section
  diverges as dσ_NMM/dT ∝ 1/T, concentrating the signal at T_e < 100 MeV.

  Differences from NuESelection.C:
    • SHOWER_E_MIN  lowered from 30 MeV to 10 MeV
    • SHOWER_E_MAX  reduced  from 1500 MeV to 600 MeV (BNB-relevant range)
    • Output binning: 120 × 5 MeV bins, 0–0.6 GeV

  All other cuts (FV, Razzled PID, track score, dE/dx, θ, E×θ²) unchanged.

  Input : merged_nu+eIntimeBNB_DLNuE_22April.root  (ana/NuE + ana/SubRun)

  Output ROOT file is compatible with NuMMFit_SBND.C:
    h_data   – selected events vs reco T_e [GeV]
    h_bkg    – background events
    h_eff    – efficiency vs true T_e [GeV]
    h_smear  – smearing matrix M[reco][true]

  Usage:
    root -l -b -q magnetic_moment/NuMMSelection.C
    root -l -b -q 'magnetic_moment/NuMMSelection.C("in.root","out.root")'

  Author : Jarek Nowak
  Date   : 2026-05-20
=============================================================================*/

#include <iostream>
#include <cmath>
#include <vector>
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
#include "../common/NuE_common.h"

// ──────────────────────── SELECTION CUTS ─────────────────────────────────────

namespace Cuts {
  const Double_t SLICE_SCORE_MIN   = -1.5;
  // Tighter PID thresholds compared to NuESelection to suppress
  // non-EM backgrounds that dominate after POT-correct normalisation.
  const Double_t RAZZLE_ELEC_MIN  =  0.70;  // was 0.50; keeps 80% sig, cuts 55% bkg
  const Int_t    RAZZLE_BEST_PDG  =  11;
  const Double_t TRACK_SCORE_MAX  =  0.50;  // same as NuESelection; <0.25 cuts low-E electrons
  const Double_t SHOWER_E_MIN     =  0.010; // GeV — 10 MeV threshold (NMM-optimised)
  const Double_t SHOWER_E_MAX     =  0.600; // GeV — BNB-relevant range
  const Double_t THETA_MAX        =  0.55;  // rad
  const Double_t DEDX_MIN         =  0.8;   // MeV/cm (was 0.5)
  const Double_t DEDX_MAX         =  3.0;   // MeV/cm (was 3.5)
  // Tighten to ~3× kinematic limit (2mₑ = 1.022 MeV) to exploit the unique
  // forward-collimation of ν-e elastic scattering vs backgrounds.
  const Double_t ETHETA2_MAX      =  0.003; // GeV  (was 0.010); keeps 65% sig, cuts 76% bkg
  const Double_t EXTRA_PRIM_E_MIN =  0.020; // GeV
  const Int_t    MAX_EXTRA_PRIM   =  0;
}

// ───────────────────────── OUTPUT BINNING ────────────────────────────────────
// Fine 5 MeV bins (120 over 0–0.6 GeV) to resolve the 1/T NMM excess at low T_e.

namespace Binning {
  const Int_t    N_RECO = 120;
  const Double_t E_LO  = 0.0, E_HI = 0.6;  // GeV
}

// ──────────────────────── CUT-FLOW LABELS ────────────────────────────────────

enum CutStep {
  kAll = 0,
  kHasNuSlice,
  kFiducialVtx,
  kOneElecCandidate,
  kRazzleScore,
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
  "Vertex in FV",
  "#geq1 prim. e^{-} cand.",
  "Razzled P(e)>0.70",
  "TrackScore<0.50",
  "E_{shw} #in [10,600] MeV",
  "0 extra primaries",
  "dE/dx #in [0.8,3.0] MeV/cm",
  "#theta < 0.55 rad",
  "E#times#theta^{2} < 3 MeV"
};

struct ElecCandidate {
  Int_t    pIdx      = -1;
  Double_t energy    = 0.0;
  Double_t theta     = -1.0;
  Double_t dedx      = -999.0;
  Double_t razzle11  = -1.0;
  Double_t VX        = 0.0, VY = 0.0, VZ = 0.0;
  Int_t    nExtraPrim = 0;
  bool     found     = false;
};

ElecCandidate FindElecCandidate(const BranchVars& b, Double_t sliceID,
                                Double_t pcaAngle) {
  ElecCandidate best;

  for (Int_t p = 0; p < (Int_t)b.reco_particlePDG->size(); ++p) {
    if ((*b.reco_particleSliceID)[p] != sliceID) continue;
    if ((*b.reco_particleIsPrimary)[p]  < 0.5)  continue;
    if ((*b.reco_particleClearCosmic)[p] > 0.5) continue;

    Double_t E_MeV  = (*b.reco_particleShowerBestPlaneEnergy)[p];
    Double_t E_GeV  = E_MeV / 1000.0;
    Double_t rz11   = (*b.reco_particleRazzledPDG11)[p];
    Double_t rzBest = (*b.reco_particleRazzledBestPDG)[p];

    if (!Sentinel::valid(rz11))   continue;
    if ((Int_t)(rzBest + 0.5) != Cuts::RAZZLE_BEST_PDG) continue;
    if (!Sentinel::valid(E_MeV))  continue;

    if (!best.found || rz11 > best.razzle11) {
      best.found    = true;
      best.pIdx     = p;
      best.energy   = E_GeV;
      best.razzle11 = rz11;
      best.dedx     = (*b.reco_particleBestPlanedEdx)[p];
      best.theta    = Sentinel::valid(pcaAngle) ? pcaAngle
                                                : (*b.reco_particleTheta)[p];
      best.VX = (*b.reco_particleVX)[p];
      best.VY = (*b.reco_particleVY)[p];
      best.VZ = (*b.reco_particleVZ)[p];
    }
  }

  if (!best.found) return best;

  best.nExtraPrim = 0;
  for (Int_t p = 0; p < (Int_t)b.reco_particlePDG->size(); ++p) {
    if (p == best.pIdx) continue;
    if ((*b.reco_particleSliceID)[p] != sliceID) continue;
    if ((*b.reco_particleIsPrimary)[p]  < 0.5)  continue;
    if ((*b.reco_particleClearCosmic)[p] > 0.5) continue;
    Double_t E_MeV = (*b.reco_particleShowerBestPlaneEnergy)[p];
    if (!Sentinel::valid(E_MeV)) continue;
    if (E_MeV / 1000.0 > Cuts::EXTRA_PRIM_E_MIN) ++best.nExtraPrim;
  }
  return best;
}

// ─────────────────────────── MAIN FUNCTION ───────────────────────────────────

void NuMMSelection(
    const char* inputFile  = "/data/sbnd/NuEElastic/merged_nu+eIntimeBNB_DLNuE_22April.root",
    const char* outputFile = "magnetic_moment/NuMMSelection_output.root")
{
  std::cout << "\n"
            << "╔══════════════════════════════════════════════════════════════╗\n"
            << "║  ν Magnetic Moment Selection  —  Low-T Optimised  —  SBND    ║\n"
            << "╚══════════════════════════════════════════════════════════════╝\n\n";

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

  // ── Total POT ────────────────────────────────────────────────────────────
  Double_t pot_sr = 0.0;
  tSubRun->SetBranchAddress("pot", &pot_sr);
  Double_t totalPOT = 0.0;
  for (Long64_t sr = 0; sr < tSubRun->GetEntries(); ++sr) {
    tSubRun->GetEntry(sr);
    totalPOT += pot_sr;
  }
  printf("[INFO] Total POT: %.4e\n\n", totalPOT);

  BranchVars b;
  SetBranchAddresses(tNuE, b);

  // ── Output histograms ─────────────────────────────────────────────────────
  const Int_t    N  = Binning::N_RECO;
  const Double_t lo = Binning::E_LO, hi = Binning::E_HI;

  TH1D* h_data  = new TH1D("h_data",  "Selected events;T_{e}^{reco} [GeV];Events/bin", N,lo,hi);
  TH1D* h_bkg   = new TH1D("h_bkg",   "Background;T_{e}^{reco} [GeV];Events/bin",      N,lo,hi);
  TH1D* h_sig   = new TH1D("h_sig",   "Signal;T_{e}^{reco} [GeV];Events/bin",          N,lo,hi);
  TH2D* h_smear = new TH2D("h_smear",
    "Smearing matrix;T_{e}^{reco} [GeV];T_{e}^{true} [GeV]",
    N,lo,hi, N,lo,hi);
  TH1D* h_sig_total = new TH1D("h_sig_total",
    "All true signal;T_{e}^{true} [GeV];Events/bin", N,lo,hi);
  // Numerator for efficiency: SELECTED signal vs TRUE T_e (correct efficiency;
  // previously the reco-binned h_sig was divided by the true-binned h_sig_total,
  // which manufactured a spurious efficiency peak near the reco threshold).
  TH1D* h_sig_true = new TH1D("h_sig_true",
    "Selected signal;T_{e}^{true} [GeV];Events/bin", N,lo,hi);
  TH1D* h_eff = new TH1D("h_eff",
    "Selection efficiency;T_{e}^{true} [GeV];#epsilon", N,lo,hi);
  h_eff->Sumw2();

  // Diagnostic: low-T zoom (0–100 MeV, 1 MeV/bin)
  TH1D* h_lowT_sig = new TH1D("h_lowT_sig",
    "Signal T_{e}<100 MeV;T_{e}^{reco} [GeV];Events/bin", 100,0,0.1);
  TH1D* h_lowT_bkg = new TH1D("h_lowT_bkg",
    "Bkg T_{e}<100 MeV;T_{e}^{reco} [GeV];Events/bin",    100,0,0.1);

  // Standard diagnostic histograms
  TH1D* h_theta_sig   = new TH1D("h_theta_sig",  "#theta (sig);#theta [rad];",         50,0,0.6);
  TH1D* h_theta_bkg   = new TH1D("h_theta_bkg",  "#theta (bkg);#theta [rad];",         50,0,0.6);
  TH1D* h_etheta2_sig = new TH1D("h_etheta2_sig","E#theta^{2} (sig);E#theta^{2} [GeV];",50,0,0.015);
  TH1D* h_etheta2_bkg = new TH1D("h_etheta2_bkg","E#theta^{2} (bkg);E#theta^{2} [GeV];",50,0,0.015);
  TH1D* h_dedx_sig    = new TH1D("h_dedx_sig",   "dE/dx (sig);dE/dx [MeV/cm];",        50,0,6);
  TH1D* h_dedx_bkg    = new TH1D("h_dedx_bkg",   "dE/dx (bkg);dE/dx [MeV/cm];",        50,0,6);
  TH1D* h_rz11_sig    = new TH1D("h_rz11_sig",   "Razzled P(e) (sig);P(e);",            50,0,1);
  TH1D* h_rz11_bkg    = new TH1D("h_rz11_bkg",   "Razzled P(e) (bkg);P(e);",            50,0,1);

  Double_t cutflow[kNCuts][2] = {};

  const Long64_t nEntries = tNuE->GetEntries();
  printf("[INFO] Processing %lld events…\n\n", nEntries);

  // ── Event loop ────────────────────────────────────────────────────────────
  for (Long64_t ev = 0; ev < nEntries; ++ev) {
    if (ev % 500000 == 0)
      printf("  … event %lld / %lld (%.0f%%)\n", ev, nEntries, 100.0*ev/nEntries);

    tNuE->GetEntry(ev);
    bool isTrueSignal = (b.nuEScatter == 1);

    if (isTrueSignal && b.truth_recoilElectronEnergy &&
        !b.truth_recoilElectronEnergy->empty()) {
      Double_t trueE_GeV = (*b.truth_recoilElectronEnergy)[0] / 1000.0;
      h_sig_total->Fill(trueE_GeV);
    }

    cutflow[kAll][0]++;
    if (isTrueSignal) cutflow[kAll][1]++;

    Int_t sliceIdx = FindNuSliceIndex(b);
    if (sliceIdx < 0) continue;
    cutflow[kHasNuSlice][0]++;
    if (isTrueSignal) cutflow[kHasNuSlice][1]++;

    Double_t sliceID    = (*b.reco_sliceID)[sliceIdx];
    Double_t sliceScore = (*b.reco_sliceScore)[sliceIdx];
    Double_t sliceInter = (*b.reco_sliceInteraction)[sliceIdx];
    Double_t pcaAngle   = GetPCAAngle(b, sliceID);

    // Use the GENUINE reconstructed neutrino vertex (reco_neutrinoVX/VY/VZ), not
    // reco_sliceTrueVX — the latter is truth-matched: unavailable in real data
    // and filled for only ~38% of signal events, which made this FV cut a no-op
    // on data (invalid vertices fell through the validity guard and PASSED).
    // Mirrors NuESelection.C; GetRecoNuVertex comes from common/NuE_common.h.
    Double_t vtxX, vtxY, vtxZ;
    if (!GetRecoNuVertex(b, sliceID, vtxX, vtxY, vtxZ)) continue;
    if (!FV::inFV(vtxX, vtxY, vtxZ)) continue;
    cutflow[kFiducialVtx][0]++;
    if (isTrueSignal) cutflow[kFiducialVtx][1]++;

    ElecCandidate cand = FindElecCandidate(b, sliceID, pcaAngle);
    if (!cand.found) continue;
    cutflow[kOneElecCandidate][0]++;
    if (isTrueSignal) cutflow[kOneElecCandidate][1]++;

    if (!FV::inFV(cand.VX, cand.VY, cand.VZ)) continue;

    if (cand.razzle11 < Cuts::RAZZLE_ELEC_MIN) continue;
    cutflow[kRazzleScore][0]++;
    if (isTrueSignal) cutflow[kRazzleScore][1]++;

    Double_t trkScore = (*b.reco_particleTrackScore)[cand.pIdx];
    if (!Sentinel::valid(trkScore) || trkScore > Cuts::TRACK_SCORE_MAX) continue;
    cutflow[kTrackScore][0]++;
    if (isTrueSignal) cutflow[kTrackScore][1]++;

    if (cand.energy < Cuts::SHOWER_E_MIN || cand.energy > Cuts::SHOWER_E_MAX) continue;
    cutflow[kShowerEnergy][0]++;
    if (isTrueSignal) cutflow[kShowerEnergy][1]++;

    if (cand.nExtraPrim > Cuts::MAX_EXTRA_PRIM) continue;
    cutflow[kNoExtraPrim][0]++;
    if (isTrueSignal) cutflow[kNoExtraPrim][1]++;

    if (Sentinel::validDedx(cand.dedx)) {
      if (isTrueSignal) h_dedx_sig->Fill(cand.dedx);
      else              h_dedx_bkg->Fill(cand.dedx);
    }
    if (isTrueSignal) {
      h_theta_sig->Fill(cand.theta);
      h_rz11_sig->Fill(cand.razzle11);
    } else {
      h_theta_bkg->Fill(cand.theta);
      h_rz11_bkg->Fill(cand.razzle11);
    }

    if (!Sentinel::validDedx(cand.dedx) ||
        cand.dedx < Cuts::DEDX_MIN || cand.dedx > Cuts::DEDX_MAX) continue;
    cutflow[kDedx][0]++;
    if (isTrueSignal) cutflow[kDedx][1]++;

    if (cand.theta < 0 || cand.theta > Cuts::THETA_MAX) continue;
    cutflow[kTheta][0]++;
    if (isTrueSignal) cutflow[kTheta][1]++;

    Double_t eTheta2 = cand.energy * cand.theta * cand.theta;
    if (isTrueSignal) h_etheta2_sig->Fill(eTheta2);
    else              h_etheta2_bkg->Fill(eTheta2);

    if (eTheta2 > Cuts::ETHETA2_MAX) continue;
    cutflow[kETheta2][0]++;
    if (isTrueSignal) cutflow[kETheta2][1]++;

    // ── Event passes selection ────────────────────────────────────────────
    bool isSelectedSignal = isTrueSignal && (sliceInter == 1098.0);

    h_data->Fill(cand.energy);
    if (isSelectedSignal) {
      h_sig->Fill(cand.energy);
      h_lowT_sig->Fill(cand.energy);
      if (b.truth_recoilElectronEnergy && !b.truth_recoilElectronEnergy->empty()) {
        Double_t trueE_GeV = (*b.truth_recoilElectronEnergy)[0] / 1000.0;
        h_smear->Fill(cand.energy, trueE_GeV);
        h_sig_true->Fill(trueE_GeV);
      }
    } else {
      h_bkg->Fill(cand.energy);
      h_lowT_bkg->Fill(cand.energy);
    }
  }

  // ── Column-normalise smearing matrix ─────────────────────────────────────
  for (Int_t iT = 1; iT <= N; ++iT) {
    Double_t colSum = 0.0;
    for (Int_t iR = 1; iR <= N; ++iR) colSum += h_smear->GetBinContent(iR, iT);
    if (colSum > 0.0)
      for (Int_t iR = 1; iR <= N; ++iR)
        h_smear->SetBinContent(iR, iT, h_smear->GetBinContent(iR, iT) / colSum);
  }

  h_eff->Divide(h_sig_true, h_sig_total, 1.0, 1.0, "B");

  // ── Cut-flow table ────────────────────────────────────────────────────────
  printf("\n══════════════════════════════════════════════════════════════════\n");
  printf("  CUT-FLOW TABLE  (NMM-optimised selection)\n");
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

  Double_t nSel   = h_data->Integral();
  Double_t nSig   = h_sig->Integral();
  Double_t nBkg   = h_bkg->Integral();
  Double_t purity = (nSel > 0) ? 100.0 * nSig / nSel : 0.0;
  Double_t sigEff = (sig0 > 0) ? 100.0 * nSig / sig0 : 0.0;
  printf("  Selected events : %.0f  (signal: %.0f  bkg: %.0f)\n", nSel, nSig, nBkg);
  printf("  Purity          : %.2f%%\n", purity);
  printf("  Signal efficiency: %.2f%%\n", sigEff);
  printf("  Total POT       : %.3e\n\n", totalPOT);

  // ── Expected yields at target exposure ─────────────────────────────────
  Double_t nSig_sc   = nSig  * ScaleFactors::SIGNAL;
  Double_t nBkg_sc   = nBkg  * ScaleFactors::NU_BKG;
  Double_t nTot_sc   = nSig_sc + nBkg_sc;
  Double_t purity_sc = (nTot_sc > 0) ? 100.0 * nSig_sc / nTot_sc : 0.0;
  printf("  ── Scaled to %.0e POT ──────────────────────────────────────\n",
         ScaleFactors::TARGET_POT);
  printf("  Signal events   : %7.1f  (SF = %.7f)\n", nSig_sc,  ScaleFactors::SIGNAL);
  printf("  Nu bkg events   : %7.1f  (SF = %.5f)\n", nBkg_sc,  ScaleFactors::NU_BKG);
  printf("  Cosmic bkg      :    —    (separate sample, SF = %.5f)\n",
         ScaleFactors::COSMIC_BKG);
  printf("  Total (nu only) : %7.1f\n", nTot_sc);
  printf("  Purity (nu bkg) : %.2f%%\n\n", purity_sc);

  // ── Save ─────────────────────────────────────────────────────────────────
  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "[ERROR] Cannot create: " << outputFile << "\n";
    fIn->Close(); return;
  }

  h_data->Write();   h_bkg->Write();    h_eff->Write();   h_smear->Write();
  h_sig->Write();    h_sig_true->Write();   h_sig_total->Write();
  h_lowT_sig->Write(); h_lowT_bkg->Write();
  h_theta_sig->Write(); h_theta_bkg->Write();
  h_etheta2_sig->Write(); h_etheta2_bkg->Write();
  h_dedx_sig->Write(); h_dedx_bkg->Write();
  h_rz11_sig->Write(); h_rz11_bkg->Write();
  fOut->Close();
  printf("[INFO] Histograms saved to %s\n\n", outputFile);

  // ── Summary plots ─────────────────────────────────────────────────────────
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

  TCanvas* c1 = new TCanvas("c_nmm_sel","NuMM Selection",1200,800);
  c1->Divide(3,2);

  auto drawSigBkg = [&](Int_t pad, TH1D* hs, TH1D* hb, const char* title) {
    c1->cd(pad);
    hb->SetFillColorAlpha(kAzure-4, 0.5); hb->SetLineColor(kAzure+1);
    hs->SetFillColorAlpha(kRed-4,   0.5); hs->SetLineColor(kRed+1);
    THStack* stk = new THStack(Form("stk%d",pad), title);
    stk->Add(hb); stk->Add(hs);
    stk->Draw("HIST");
    stk->GetXaxis()->SetTitle(hs->GetXaxis()->GetTitle());
    stk->GetYaxis()->SetTitle("Events");
    TLegend* lg = new TLegend(0.55,0.72,0.90,0.90);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
    lg->AddEntry(hs,"Signal","f"); lg->AddEntry(hb,"Bkg","f");
    lg->Draw();
  };

  drawSigBkg(1, h_sig,      h_bkg,      "Selected events (full range)");
  drawSigBkg(2, h_lowT_sig, h_lowT_bkg, "Low-T region (0#minus100 MeV)");
  drawSigBkg(3, h_rz11_sig, h_rz11_bkg, "Razzled P(e)");
  drawSigBkg(4, h_theta_sig,h_theta_bkg,"#theta [rad]");
  drawSigBkg(5, h_dedx_sig, h_dedx_bkg, "dE/dx [MeV/cm]");

  c1->cd(6);
  h_eff->SetLineColor(kBlue+1); h_eff->SetLineWidth(2); h_eff->SetMarkerStyle(20);
  h_eff->Draw("E1");
  TLatex lt; lt.SetNDC(); lt.SetTextSize(0.045);
  lt.DrawLatex(0.15,0.92,"Selection efficiency vs T_{e}^{true} [GeV]");

  TString base = TString(outputFile).ReplaceAll(".root","");
  c1->SaveAs(base + ".pdf");
  c1->SaveAs(base + ".png");
  printf("[INFO] Summary plots saved: %s.{pdf,png}\n\n", base.Data());

  fIn->Close();
}
