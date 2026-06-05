/*=============================================================================
  NuESkim.C
  ─────────────────────────────────────────────────────────────────────────────
  ONE-PASS skim of the SBND merged MC for the ν_μ + e⁻ elastic selection.

  Purpose:
    Reading the 143 GB merged tree takes ~25 min. To allow fast cut tuning we
    make a single pass that applies only a LOOSE preselection and writes a small
    flat TTree ("cand") with one row per electron-shower candidate, holding every
    discriminating reconstruction variable plus truth flags and per-sample POT
    weights. We also store everything needed to rebuild the final analysis
    histograms WITHOUT a second full pass:
      • h_sig_total  – efficiency denominator (all true signal vs true T_e)
      • v_preflow    – pre-selection cut-flow counts (all / signal)
      • v_meta       – total POT and bookkeeping

  Loose preselection (looser than any tunable cut, so the optimiser can explore):
    1. pandora ν-candidate slice (category==1, valid score)
    2. reco neutrino vertex in fiducial volume   (reco_neutrinoVX/VY/VZ)
    3. ≥1 primary, non-clear-cosmic shower with Razzled best-PDG==11,
       valid shower energy > 10 MeV

  Downstream:
    NuEOptimize.C  – tune cuts on the skim (POT-scaled S/√(S+B)) and emit
                     fit-ready output histograms.

  Usage:
    root -l -b -q NuESkim.C
    root -l -b -q 'NuESkim.C("in.root","NuESkim.root")'

  Author : Jarek Nowak  (skim infrastructure added 2026-05-29)
=============================================================================*/

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TVectorD.h"
#include "common/NuE_common.h"

// ──────────────── LOOSE PRESELECTION (skim level only) ───────────────────────
namespace Pre {
  const Double_t SLICE_SCORE_MIN  = -1.5;   // exclude truly invalid slices
  const Int_t    RAZZLE_BEST_PDG  =  11;    // candidate is an electron shower
  const Double_t SHOWER_E_MIN_GEV =  0.010; // 10 MeV loose floor (tunable upward)
  const Double_t EXTRA_PRIM_E_MIN =  0.020; // GeV — "significant" extra primary
}

// ─────────────────────────── MAIN ────────────────────────────────────────────
void NuESkim(
    const char* inputFile  = "/data/sbnd/NuEElastic/merged_nu+eIntimeBNB_DLNuE_22April.root",
    const char* outputFile = "NuESkim.root")
{
  std::cout << "\n=== NuESkim : one-pass loose skim for fast cut tuning ===\n\n";

  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) { std::cerr << "[ERROR] cannot open " << inputFile << "\n"; return; }
  TTree* tNuE    = (TTree*)fIn->Get("ana/NuE");
  TTree* tSubRun = (TTree*)fIn->Get("ana/SubRun");
  if (!tNuE || !tSubRun) { std::cerr << "[ERROR] trees not found\n"; fIn->Close(); return; }

  Double_t pot_sr = 0.0, totalPOT = 0.0;
  tSubRun->SetBranchAddress("pot", &pot_sr);
  for (Long64_t sr = 0; sr < tSubRun->GetEntries(); ++sr) { tSubRun->GetEntry(sr); totalPOT += pot_sr; }
  printf("[INFO] Total POT: %.4e\n", totalPOT);

  BranchVars b;
  SetBranchAddresses(tNuE, b);

  // ── Output skim file & tree ───────────────────────────────────────────────
  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  TTree* tc = new TTree("cand", "Loose-preselected e-shower candidates");

  Int_t    s_isSig, s_nuE, s_nExtraPrim;
  Double_t s_eReco, s_eTrue, s_sliceInter, s_sliceScore;
  Double_t s_rzE, s_rzG, s_rzMu, s_rzPi, s_rzP, s_bestPDG;
  Double_t s_trk, s_theta, s_thetaPart, s_etheta2;
  Double_t s_dedx, s_dedx0, s_dedx1, s_dedx2;
  Double_t s_shwLen, s_shwOA;
  Double_t s_vtxX, s_vtxY, s_vtxZ;

  tc->Branch("isSig",      &s_isSig);       // nuE==1 && sliceInter==1098
  tc->Branch("nuE",        &s_nuE);         // raw nuEScatter flag
  tc->Branch("eReco",      &s_eReco);       // GeV
  tc->Branch("eTrue",      &s_eTrue);       // GeV (-1 if N/A)
  tc->Branch("sliceInter", &s_sliceInter);
  tc->Branch("sliceScore", &s_sliceScore);
  tc->Branch("rzE",        &s_rzE);
  tc->Branch("rzG",        &s_rzG);
  tc->Branch("rzMu",       &s_rzMu);
  tc->Branch("rzPi",       &s_rzPi);
  tc->Branch("rzP",        &s_rzP);
  tc->Branch("bestPDG",    &s_bestPDG);
  tc->Branch("trk",        &s_trk);
  tc->Branch("theta",      &s_theta);       // PCA-preferred
  tc->Branch("thetaPart",  &s_thetaPart);   // particle theta
  tc->Branch("etheta2",    &s_etheta2);     // GeV
  tc->Branch("dedx",       &s_dedx);
  tc->Branch("dedx0",      &s_dedx0);
  tc->Branch("dedx1",      &s_dedx1);
  tc->Branch("dedx2",      &s_dedx2);
  tc->Branch("shwLen",     &s_shwLen);
  tc->Branch("shwOA",      &s_shwOA);
  tc->Branch("nExtraPrim", &s_nExtraPrim);
  tc->Branch("vtxX",       &s_vtxX);
  tc->Branch("vtxY",       &s_vtxY);
  tc->Branch("vtxZ",       &s_vtxZ);

  // Efficiency denominator: all true signal vs true T_e (120 bins, 0–1.2 GeV)
  TH1D* h_sig_total = new TH1D("h_sig_total",
    "All true signal;T_{e}^{true} [GeV];Events/bin", 120, 0.0, 1.2);
  h_sig_total->SetDirectory(fOut);

  // Pre-selection cut-flow: [stage][0=all,1=signal]; stages: All,Slice,FV,Cand
  Double_t pre[4][2] = {};

  const Long64_t nEntries = tNuE->GetEntries();
  printf("[INFO] Processing %lld events…\n", nEntries);

  for (Long64_t ev = 0; ev < nEntries; ++ev) {
    if (ev % 1000000 == 0) { printf("  … %lld / %lld\n", ev, nEntries); fflush(stdout); }
    tNuE->GetEntry(ev);

    const bool isTrueSignal = (b.nuEScatter == 1);

    if (isTrueSignal && b.truth_recoilElectronEnergy && !b.truth_recoilElectronEnergy->empty())
      h_sig_total->Fill((*b.truth_recoilElectronEnergy)[0] / 1000.0);

    pre[0][0]++; if (isTrueSignal) pre[0][1]++;

    Int_t sIdx = FindNuSliceIndex(b);
    if (sIdx < 0) continue;
    pre[1][0]++; if (isTrueSignal) pre[1][1]++;

    Double_t sliceID    = (*b.reco_sliceID)[sIdx];
    Double_t sliceScore = (*b.reco_sliceScore)[sIdx];
    Double_t sliceInter = (*b.reco_sliceInteraction)[sIdx];
    Double_t pcaAngle   = GetPCAAngle(b, sliceID);

    // Genuine reco neutrino vertex (not truth-matched). Require it valid & in FV.
    Double_t vx, vy, vz;
    if (!GetRecoNuVertex(b, sliceID, vx, vy, vz)) continue;
    if (!FV::inFV(vx, vy, vz)) continue;
    pre[2][0]++; if (isTrueSignal) pre[2][1]++;

    // Find best primary electron-shower candidate (highest Razzled P(e))
    Int_t bestP = -1; Double_t bestRzE = -1.0;
    for (size_t p = 0; p < b.reco_particleSliceID->size(); ++p) {
      if ((*b.reco_particleSliceID)[p] != sliceID) continue;
      if ((*b.reco_particleIsPrimary)[p] < 0.5) continue;
      if ((*b.reco_particleClearCosmic)[p] > 0.5) continue;
      Double_t E   = (*b.reco_particleShowerBestPlaneEnergy)[p];
      Double_t rzE = (*b.reco_particleRazzledPDG11)[p];
      Double_t rzB = (*b.reco_particleRazzledBestPDG)[p];
      if (!Sentinel::valid(E) || E/1000.0 < Pre::SHOWER_E_MIN_GEV) continue;
      if (!Sentinel::valid(rzE)) continue;
      if ((Int_t)(rzB + 0.5) != Pre::RAZZLE_BEST_PDG) continue;
      if (rzE > bestRzE) { bestRzE = rzE; bestP = p; }
    }
    if (bestP < 0) continue;
    pre[3][0]++; if (isTrueSignal) pre[3][1]++;

    // Count extra significant primaries (not the candidate)
    Int_t nExtra = 0;
    for (size_t p = 0; p < b.reco_particleSliceID->size(); ++p) {
      if ((Int_t)p == bestP) continue;
      if ((*b.reco_particleSliceID)[p] != sliceID) continue;
      if ((*b.reco_particleIsPrimary)[p] < 0.5) continue;
      if ((*b.reco_particleClearCosmic)[p] > 0.5) continue;
      Double_t E = (*b.reco_particleShowerBestPlaneEnergy)[p];
      if (Sentinel::valid(E) && E/1000.0 > Pre::EXTRA_PRIM_E_MIN) ++nExtra;
    }

    // Fill candidate row
    Double_t E_GeV = (*b.reco_particleShowerBestPlaneEnergy)[bestP] / 1000.0;
    Double_t pca   = pcaAngle;
    Double_t theta = Sentinel::valid(pca) ? pca : (*b.reco_particleTheta)[bestP];

    s_nuE        = b.nuEScatter;
    s_isSig      = (isTrueSignal && sliceInter == 1098.0) ? 1 : 0;
    s_eReco      = E_GeV;
    s_eTrue      = (isTrueSignal && b.truth_recoilElectronEnergy && !b.truth_recoilElectronEnergy->empty())
                   ? (*b.truth_recoilElectronEnergy)[0] / 1000.0 : -1.0;
    s_sliceInter = sliceInter;
    s_sliceScore = sliceScore;
    s_rzE        = (*b.reco_particleRazzledPDG11)[bestP];
    s_rzG        = (*b.reco_particleRazzledPDG22)[bestP];
    s_rzMu       = (*b.reco_particleRazzledPDG13)[bestP];
    s_rzPi       = (*b.reco_particleRazzledPDG211)[bestP];
    s_rzP        = (*b.reco_particleRazzledPDG2212)[bestP];
    s_bestPDG    = (*b.reco_particleRazzledBestPDG)[bestP];
    s_trk        = (*b.reco_particleTrackScore)[bestP];
    s_theta      = theta;
    s_thetaPart  = (*b.reco_particleTheta)[bestP];
    s_etheta2    = E_GeV * theta * theta;
    s_dedx       = (*b.reco_particleBestPlanedEdx)[bestP];
    s_dedx0      = (*b.reco_particlePlane0dEdx)[bestP];
    s_dedx1      = (*b.reco_particlePlane1dEdx)[bestP];
    s_dedx2      = (*b.reco_particlePlane2dEdx)[bestP];
    s_shwLen     = (*b.reco_particleShowerLength)[bestP];
    s_shwOA      = (*b.reco_particleShowerOpenAngle)[bestP];
    s_nExtraPrim = nExtra;
    s_vtxX = vx; s_vtxY = vy; s_vtxZ = vz;

    tc->Fill();
  } // event loop

  // ── Metadata vectors ──────────────────────────────────────────────────────
  TVectorD v_meta(4);
  v_meta[0] = totalPOT;
  v_meta[1] = ScaleFactors::SIGNAL;
  v_meta[2] = ScaleFactors::NU_BKG;
  v_meta[3] = (Double_t)nEntries;

  TVectorD v_preflow(8); // All,Slice,FV,Cand  ×  (all,sig)
  for (int i = 0; i < 4; ++i) { v_preflow[2*i] = pre[i][0]; v_preflow[2*i+1] = pre[i][1]; }

  fOut->cd();
  tc->Write();
  h_sig_total->Write();
  v_meta.Write("v_meta");
  v_preflow.Write("v_preflow");

  printf("\n[INFO] Skim candidates written: %lld\n", tc->GetEntries());
  printf("[INFO] Pre-flow (all / sig): All=%.0f/%.0f  Slice=%.0f/%.0f  FV=%.0f/%.0f  Cand=%.0f/%.0f\n",
         pre[0][0],pre[0][1], pre[1][0],pre[1][1], pre[2][0],pre[2][1], pre[3][0],pre[3][1]);
  printf("[INFO] Saved skim → %s\n\n", outputFile);

  fOut->Close();
  fIn->Close();
}
