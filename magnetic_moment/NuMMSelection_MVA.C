/*=============================================================================
  NuMMSelection_MVA.C
  ─────────────────────────────────────────────────────────────────────────────
  Applies the TMVA BDTG classifier trained by NuMMClassifier.C to select
  ν-e elastic scatter candidates for the neutrino magnetic moment analysis.

  Replaces the box-cut sequence in NuMMSelection.C with a single MVA score
  cut. The output histograms (h_data, h_bkg, h_eff, h_smear) are identical
  in format to NuMMSelection.C and are compatible with NuMMFit_SBND.C.

  Pre-selection (applied before MVA — physics and detector requirements):
    • Has Pandora ν-slice with valid NuScore
    • Vertex reconstructed within fiducial volume
    • ≥1 primary e⁻ candidate (razzledBestPDG=11) with valid track score
    • Shower energy ∈ [10, 600] MeV  (NMM-relevant range, hard analysis cut)

  MVA selection:
    • BDTG score > MVA_CUT  (tune based on ROC curve from NuMMClassifier.C)

  Usage:
    root -l -b -q 'magnetic_moment/NuMMSelection_MVA.C("input.root")'
    root -l -b -q 'magnetic_moment/NuMMSelection_MVA.C+("input.root")'

  To change the MVA working point, edit MVA_CUT below (BDTG range: −1 to +1).
  Tighter cuts improve purity; looser cuts improve signal efficiency.

  Author : Jarek Nowak / Claude
  Date   : 2026-05-21
=============================================================================*/

#include <iostream>
#include <cstdio>
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
#include "TString.h"
#include "TMath.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "NuMM_common.h"

// ── Working point ─────────────────────────────────────────────────────────────
// BDTG output runs from −1 (background-like) to +1 (signal-like).
// Default 0.1 is a starting point; open the TMVA GUI after training to choose
// the operating point that meets your signal efficiency or purity target.
static const Float_t MVA_CUT = 0.8f;

// ─────────────────────────── MAIN FUNCTION ───────────────────────────────────

void NuMMSelection_MVA(
    const char* inputFile  = "/data/sbnd/NuEElastic/merged_nu+eIntimeBNB_DLNuE_22April.root",
    const char* outputFile = "magnetic_moment/NuMMSelection_MVA_output.root",
    const char* weightsFile = "magnetic_moment/tmva_dataset/weights/NuMMClassifier_BDTG.weights.xml")
{
  std::cout << "\n"
            << "╔══════════════════════════════════════════════════════════════╗\n"
            << "║  NuMM MVA Selection  —  SBND                                 ║\n"
            << "╚══════════════════════════════════════════════════════════════╝\n\n";

  printf("[INFO] BDTG weights : %s\n", weightsFile);
  printf("[INFO] MVA cut      : %.3f\n\n", (Double_t)MVA_CUT);

  // ── Set up TMVA reader ───────────────────────────────────────────────────────
  // Variables MUST be added in the same order as in NuMMClassifier.C.
  TMVA::Tools::Instance();
  TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

  Float_t v_razzle_p_elec, v_track_score, v_shower_energy;
  Float_t v_dedx, v_has_dedx, v_theta, v_etheta2;
  Float_t v_n_extra_prim, v_slice_score;
  Float_t v_true_energy;

  reader->AddVariable("razzle_p_elec", &v_razzle_p_elec);
  reader->AddVariable("track_score",   &v_track_score);
  reader->AddVariable("shower_energy", &v_shower_energy);
  reader->AddVariable("dedx",          &v_dedx);
  reader->AddVariable("has_dedx",      &v_has_dedx);
  reader->AddVariable("theta",         &v_theta);
  reader->AddVariable("etheta2",       &v_etheta2);
  reader->AddVariable("n_extra_prim",  &v_n_extra_prim);
  reader->AddVariable("slice_score",   &v_slice_score);
  reader->AddSpectator("true_energy",  &v_true_energy);

  reader->BookMVA("BDTG", weightsFile);

  // ── Open input ───────────────────────────────────────────────────────────────
  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERROR] Cannot open: " << inputFile << "\n";
    delete reader; return;
  }
  TTree* tNuE    = (TTree*)fIn->Get("ana/NuE");
  TTree* tSubRun = (TTree*)fIn->Get("ana/SubRun");
  if (!tNuE || !tSubRun) {
    std::cerr << "[ERROR] Trees not found.\n"; fIn->Close(); delete reader; return;
  }

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

  // ── Output histograms ─────────────────────────────────────────────────────────
  const Int_t    N  = Binning::N_RECO;
  const Double_t lo = Binning::E_LO, hi = Binning::E_HI;

  TH1D* h_data      = new TH1D("h_data",      "Selected;T_{e}^{reco} [GeV];Events/bin",    N,lo,hi);
  TH1D* h_bkg       = new TH1D("h_bkg",       "Background;T_{e}^{reco} [GeV];Events/bin",  N,lo,hi);
  TH1D* h_sig       = new TH1D("h_sig",       "Signal;T_{e}^{reco} [GeV];Events/bin",      N,lo,hi);
  TH1D* h_sig_total = new TH1D("h_sig_total", "All sig;T_{e}^{true} [GeV];Events/bin",     N,lo,hi);
  TH1D* h_eff       = new TH1D("h_eff",       "Efficiency;T_{e}^{true} [GeV];#varepsilon", N,lo,hi);
  TH2D* h_smear     = new TH2D("h_smear",
      "Smearing;T_{e}^{reco} [GeV];T_{e}^{true} [GeV]", N,lo,hi, N,lo,hi);

  // Diagnostic: MVA score distributions and low-T zoom
  TH1D* h_mva_sig  = new TH1D("h_mva_sig",  "BDTG score (sig);BDTG score;Events", 100,-1,1);
  TH1D* h_mva_bkg  = new TH1D("h_mva_bkg",  "BDTG score (bkg);BDTG score;Events", 100,-1,1);
  TH1D* h_lowT_sig = new TH1D("h_lowT_sig", "Signal T_{e}<100 MeV;T_{e}^{reco} [GeV];Events", 100,0,0.1);
  TH1D* h_lowT_bkg = new TH1D("h_lowT_bkg", "Bkg T_{e}<100 MeV;T_{e}^{reco} [GeV];Events",   100,0,0.1);

  h_eff->Sumw2();

  // ── Event loop ────────────────────────────────────────────────────────────────
  const Long64_t nEntries = tNuE->GetEntries();
  printf("[INFO] Processing %lld events...\n\n", nEntries);

  Double_t sig0 = 0.0;

  for (Long64_t ev = 0; ev < nEntries; ++ev) {
    if (ev % 500000 == 0)
      printf("  ... event %lld / %lld (%.0f%%)\n", ev, nEntries, 100.0 * ev / nEntries);
    tNuE->GetEntry(ev);

    bool isTrueSignal = (b.nuEScatter == 1);

    if (isTrueSignal && b.truth_recoilElectronEnergy &&
        !b.truth_recoilElectronEnergy->empty()) {
      Double_t trueE = (*b.truth_recoilElectronEnergy)[0] / 1000.0;
      h_sig_total->Fill(trueE);
      sig0 += 1.0;
    }

    // ── Pre-selection ─────────────────────────────────────────────────────────
    Int_t sliceIdx = FindNuSliceIndex(b);
    if (sliceIdx < 0) continue;

    Double_t sliceID    = (*b.reco_sliceID)[sliceIdx];
    Double_t sliceScore = (*b.reco_sliceScore)[sliceIdx];
    Double_t sliceInter = (*b.reco_sliceInteraction)[sliceIdx];
    Double_t pcaAngle   = GetPCAAngle(b, sliceID);

    Double_t vtxX = (*b.reco_sliceTrueVX)[sliceIdx];
    Double_t vtxY = (*b.reco_sliceTrueVY)[sliceIdx];
    Double_t vtxZ = (*b.reco_sliceTrueVZ)[sliceIdx];
    if (Sentinel::valid(vtxX) && Sentinel::valid(vtxY) && Sentinel::valid(vtxZ))
      if (!FV::inFV(vtxX, vtxY, vtxZ)) continue;

    ElecCandidate cand = FindElecCandidate(b, sliceID, pcaAngle);
    if (!cand.found)                        continue;
    if (!FV::inFV(cand.VX, cand.VY, cand.VZ)) continue;
    if (cand.trackScore < 0.0f)             continue;  // same requirement as training
    // Energy range [10, 600] MeV enforced in FindElecCandidate via PreCuts

    // ── Evaluate BDTG ─────────────────────────────────────────────────────────
    v_razzle_p_elec = cand.razzle11;
    v_track_score   = cand.trackScore;
    v_shower_energy = cand.energy;
    v_dedx          = (cand.dedx > 0.0f) ? TMath::Min(cand.dedx, 20.0f) : 0.0f;
    v_has_dedx      = cand.hasDedx;
    v_theta         = cand.theta;
    v_etheta2       = cand.energy * cand.theta * cand.theta;
    v_n_extra_prim  = cand.nExtraPrim;
    v_slice_score   = (Float_t)sliceScore;
    v_true_energy   = 0.0f;

    Float_t mvaScore = reader->EvaluateMVA("BDTG");

    bool isSelectedSignal = isTrueSignal && (TMath::Abs(sliceInter - 1098.0) < 0.5);

    // Fill MVA score distributions before cut (for working-point inspection)
    if (isSelectedSignal)    h_mva_sig->Fill(mvaScore);
    else if (!isTrueSignal)  h_mva_bkg->Fill(mvaScore);

    if (mvaScore < MVA_CUT) continue;

    // ── Event passes selection ────────────────────────────────────────────────
    h_data->Fill(cand.energy);
    if (isSelectedSignal) {
      h_sig->Fill(cand.energy);
      h_lowT_sig->Fill(cand.energy);
      if (b.truth_recoilElectronEnergy && !b.truth_recoilElectronEnergy->empty()) {
        Double_t trueE = (*b.truth_recoilElectronEnergy)[0] / 1000.0;
        h_smear->Fill(cand.energy, trueE);
      }
    } else if (!isTrueSignal) {
      h_bkg->Fill(cand.energy);
      h_lowT_bkg->Fill(cand.energy);
    }
  }

  // ── Column-normalise smearing matrix ─────────────────────────────────────────
  for (Int_t iT = 1; iT <= N; ++iT) {
    Double_t colSum = 0.0;
    for (Int_t iR = 1; iR <= N; ++iR) colSum += h_smear->GetBinContent(iR, iT);
    if (colSum > 0.0)
      for (Int_t iR = 1; iR <= N; ++iR)
        h_smear->SetBinContent(iR, iT, h_smear->GetBinContent(iR, iT) / colSum);
  }

  h_eff->Divide(h_sig, h_sig_total, 1.0, 1.0, "B");

  // ── Summary ───────────────────────────────────────────────────────────────────
  Double_t nSel   = h_data->Integral();
  Double_t nSig   = h_sig->Integral();
  Double_t nBkg   = h_bkg->Integral();
  Double_t purity = (nSel > 0) ? 100.0 * nSig / nSel : 0.0;
  Double_t sigEff = (sig0 > 0) ? 100.0 * nSig / sig0  : 0.0;

  printf("\n══════════════════════════════════════════════════════════════════\n");
  printf("  MVA SELECTION RESULTS  (BDTG > %.2f)\n", (Double_t)MVA_CUT);
  printf("──────────────────────────────────────────────────────────────────\n");
  printf("  Selected events  : %.0f  (signal: %.0f  bkg: %.0f)\n", nSel, nSig, nBkg);
  printf("  Raw purity       : %.2f%%\n", purity);
  printf("  Signal efficiency: %.2f%%\n", sigEff);
  printf("  Total POT        : %.3e\n", totalPOT);

  Double_t nSig_sc   = nSig * ScaleFactors::SIGNAL;
  Double_t nBkg_sc   = nBkg * ScaleFactors::NU_BKG;
  Double_t nTot_sc   = nSig_sc + nBkg_sc;
  Double_t purity_sc = (nTot_sc > 0) ? 100.0 * nSig_sc / nTot_sc : 0.0;

  printf("\n  ── Scaled to %.0e POT ──────────────────────────────────────\n",
         ScaleFactors::TARGET_POT);
  printf("  Signal events   : %7.1f  (SF = %.7f)\n", nSig_sc,  ScaleFactors::SIGNAL);
  printf("  Nu bkg events   : %7.1f  (SF = %.5f)\n", nBkg_sc,  ScaleFactors::NU_BKG);
  printf("  Cosmic bkg      :    —    (separate sample, SF = %.5f)\n",
         ScaleFactors::COSMIC_BKG);
  printf("  Total (nu only) : %7.1f\n", nTot_sc);
  printf("  Purity (nu bkg) : %.2f%%\n", purity_sc);
  printf("══════════════════════════════════════════════════════════════════\n\n");

  // ── Save ─────────────────────────────────────────────────────────────────────
  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "[ERROR] Cannot create: " << outputFile << "\n";
    fIn->Close(); delete reader; return;
  }
  h_data->Write();     h_bkg->Write();      h_sig->Write();
  h_eff->Write();      h_smear->Write();    h_sig_total->Write();
  h_mva_sig->Write();  h_mva_bkg->Write();
  h_lowT_sig->Write(); h_lowT_bkg->Write();
  fOut->Close();
  printf("[INFO] Histograms saved to %s\n\n", outputFile);

  // ── Summary plot ─────────────────────────────────────────────────────────────
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

  TCanvas* c1 = new TCanvas("c_mva_sel", "NuMM MVA Selection", 1200, 800);
  c1->Divide(3, 2);

  auto drawSB = [&](Int_t pad, TH1D* hs, TH1D* hb, const char* title) {
    c1->cd(pad);
    hb->SetFillColorAlpha(kAzure-4, 0.5); hb->SetLineColor(kAzure+1);
    hs->SetFillColorAlpha(kRed-4,   0.5); hs->SetLineColor(kRed+1);
    THStack* stk = new THStack(Form("stk%d",pad), title);
    stk->Add(hb); stk->Add(hs);
    stk->Draw("HIST");
    stk->SetMinimum(0.5);
    stk->GetXaxis()->SetTitle(hs->GetXaxis()->GetTitle());
    stk->GetYaxis()->SetTitle("Events");
    gPad->SetLogy(1);
    TLegend* lg = new TLegend(0.55, 0.72, 0.90, 0.90);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
    lg->AddEntry(hs, "Signal", "f"); lg->AddEntry(hb, "Bkg", "f");
    lg->Draw();
  };

  drawSB(1, h_sig,      h_bkg,      "Selected events (0#minus600 MeV)");
  drawSB(2, h_lowT_sig, h_lowT_bkg, "Low-T_{e} region (0#minus100 MeV)");
  drawSB(3, h_mva_sig,  h_mva_bkg,  "BDTG score (all pre-selected)");

  c1->cd(4);
  h_eff->SetLineColor(kBlue+1); h_eff->SetLineWidth(2); h_eff->SetMarkerStyle(20);
  h_eff->Draw("E1");
  TLatex lt; lt.SetNDC(); lt.SetTextSize(0.045);
  lt.DrawLatex(0.12, 0.92, "Selection efficiency vs T_{e}^{true} [GeV]");

  c1->cd(5);
  // Signal/background ratio vs reco energy
  TH1D* h_ratio = (TH1D*)h_sig->Clone("h_ratio");
  h_ratio->SetTitle("S / (S+B) per bin;T_{e}^{reco} [GeV];S/(S+B)");
  TH1D* h_tot   = (TH1D*)h_sig->Clone("h_tot");
  h_tot->Add(h_bkg);
  h_ratio->Divide(h_tot);
  h_ratio->SetLineColor(kBlack); h_ratio->SetLineWidth(2);
  h_ratio->Draw("HIST");

  c1->cd(6);
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.038);
  tx.DrawLatex(0.05, 0.90, Form("BDTG cut: %.2f", (Double_t)MVA_CUT));
  tx.DrawLatex(0.05, 0.80, Form("Signal eff: %.1f%%", sigEff));
  tx.DrawLatex(0.05, 0.70, Form("Raw purity: %.1f%%", purity));
  tx.DrawLatex(0.05, 0.60, Form("Scaled purity (1e21): %.1f%%", purity_sc));
  tx.DrawLatex(0.05, 0.50, Form("Sig (1e21): %.1f events", nSig_sc));
  tx.DrawLatex(0.05, 0.40, Form("Bkg (1e21): %.1f events", nBkg_sc));

  TString base = TString(outputFile).ReplaceAll(".root", "");
  c1->SaveAs(base + ".pdf");
  c1->SaveAs(base + ".png");
  printf("[INFO] Summary plots saved: %s.{pdf,png}\n\n", base.Data());

  delete reader;
  fIn->Close();
}
