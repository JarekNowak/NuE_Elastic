/*=============================================================================
  NuMMClassifier.C
  ─────────────────────────────────────────────────────────────────────────────
  Trains TMVA multivariate classifiers (BDTG, BDT, MLP) on ν-e elastic scatter
  events to replace the box-cut selection in NuMMSelection.C.

  Phase 1: Loop over the merged MC file, apply loose pre-selection (ν-slice,
           FV, ≥1 primary e⁻ candidate), and write a flat TTree with one row
           per event containing eight discriminating variables + signal label.

  Phase 2: Feed the flat trees to TMVA. Trains BDTG (primary), BDT, and MLP.
           Weights saved to magnetic_moment/tmva_dataset/weights/.
           TMVA output (ROC curves, over-training checks) saved to
           magnetic_moment/TMVA_output.root.
           Open the TMVA GUI to inspect results:
             root -l 'magnetic_moment/tmva_dataset/weights/NuMMClassifier_BDTG.weights.xml'

  MVA input variables (per best primary e⁻ candidate, no box cuts applied):
    razzle_p_elec   Razzled P(electron)              [0, 1]
    track_score     Pandora shower/track score        [0, 1]
    shower_energy   Reco shower energy                [GeV]
    dedx            dE/dx at shower start             [MeV/cm]; 0 if unavailable
    has_dedx        dE/dx validity flag               {0, 1}
    theta           PCA angle from beam axis          [rad]
    etheta2         E × θ²                            [GeV]
    n_extra_prim    Extra primary-like particles      [count as float]
    slice_score     Pandora NuScore                   [-1, 1]

  Signal     : nuEScatter == 1 && sliceInteraction == 1098
  Background : nuEScatter == 0

  Usage:
    root -l -b -q 'magnetic_moment/NuMMClassifier.C("input.root")'
    root -l -b -q 'magnetic_moment/NuMMClassifier.C+("input.root")'

  Author : Jarek Nowak / Claude
  Date   : 2026-05-21
=============================================================================*/

#include <iostream>
#include <cstdio>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"

#include "NuMM_common.h"

// ─────────────────────────── MAIN FUNCTION ───────────────────────────────────

void NuMMClassifier(
    const char* inputFile = "/data/sbnd/NuEElastic/merged_nu+eIntimeBNB_DLNuE_22April.root",
    const char* outDir    = "magnetic_moment")
{
  std::cout << "\n"
            << "╔══════════════════════════════════════════════════════════════╗\n"
            << "║  NuMM MVA Classifier Training  —  SBND                       ║\n"
            << "╚══════════════════════════════════════════════════════════════╝\n\n";

  // ── Open input ───────────────────────────────────────────────────────────────
  TFile* fIn = TFile::Open(inputFile, "READ");
  if (!fIn || fIn->IsZombie()) {
    std::cerr << "[ERROR] Cannot open: " << inputFile << "\n"; return;
  }
  TTree* tNuE = (TTree*)fIn->Get("ana/NuE");
  if (!tNuE) {
    std::cerr << "[ERROR] Tree ana/NuE not found.\n"; fIn->Close(); return;
  }

  BranchVars b;
  SetBranchAddresses(tNuE, b);

  // ── Phase 1: Build flat training trees ───────────────────────────────────────
  printf("Phase 1: Building flat training tree from %lld events...\n\n",
         tNuE->GetEntries());

  TString flatPath = TString(outDir) + "/tmva_flat.root";
  TFile* fFlat = TFile::Open(flatPath, "RECREATE");
  if (!fFlat || fFlat->IsZombie()) {
    std::cerr << "[ERROR] Cannot create: " << flatPath << "\n";
    fIn->Close(); return;
  }

  // Variables written to the flat tree — exactly what NuMMSelection_MVA.C reads
  Float_t v_razzle_p_elec, v_track_score, v_shower_energy;
  Float_t v_dedx, v_has_dedx, v_theta, v_etheta2;
  Float_t v_n_extra_prim, v_slice_score;
  Float_t v_true_energy;   // spectator — not used in classifier, for diagnostics

  TTree* sigTree = new TTree("sig_tree", "Signal (nuEScatter==1, sliceInter==1098)");
  TTree* bkgTree = new TTree("bkg_tree", "Background (nuEScatter==0)");

  // Set up branches on one tree, then replicate to the other
  TTree* trees[2] = { sigTree, bkgTree };
  for (Int_t it = 0; it < 2; ++it) {
    TTree* t = trees[it];
    t->Branch("razzle_p_elec",  &v_razzle_p_elec, "razzle_p_elec/F");
    t->Branch("track_score",    &v_track_score,   "track_score/F");
    t->Branch("shower_energy",  &v_shower_energy, "shower_energy/F");
    t->Branch("dedx",           &v_dedx,          "dedx/F");
    t->Branch("has_dedx",       &v_has_dedx,      "has_dedx/F");
    t->Branch("theta",          &v_theta,         "theta/F");
    t->Branch("etheta2",        &v_etheta2,       "etheta2/F");
    t->Branch("n_extra_prim",   &v_n_extra_prim,  "n_extra_prim/F");
    t->Branch("slice_score",    &v_slice_score,   "slice_score/F");
    t->Branch("true_energy",    &v_true_energy,   "true_energy/F");
  }

  const Long64_t nEntries = tNuE->GetEntries();
  Long64_t nSigFlat = 0, nBkgFlat = 0;

  for (Long64_t ev = 0; ev < nEntries; ++ev) {
    if (ev % 500000 == 0)
      printf("  %lld / %lld (%.0f%%)\n", ev, nEntries, 100.0 * ev / nEntries);
    tNuE->GetEntry(ev);

    bool isTrueSignal = (b.nuEScatter == 1);

    // Loose pre-selection — identical to NuMMSelection_MVA.C
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
    if (!cand.found)                   continue;
    if (!FV::inFV(cand.VX, cand.VY, cand.VZ)) continue;
    if (cand.trackScore < 0.0f)        continue;  // require valid track score

    // Fill MVA input variables
    v_razzle_p_elec = cand.razzle11;
    v_track_score   = cand.trackScore;
    v_shower_energy = cand.energy;
    // Cap at 20 MeV/cm to suppress extreme outliers that give dedx zero importance.
    // Physical EM showers are [0.5, 5] MeV/cm; values above ~20 are reconstruction artefacts.
    v_dedx          = (cand.dedx > 0.0f) ? TMath::Min(cand.dedx, 20.0f) : 0.0f;
    v_has_dedx      = cand.hasDedx;
    v_theta         = cand.theta;
    v_etheta2       = cand.energy * cand.theta * cand.theta;
    v_n_extra_prim  = cand.nExtraPrim;
    v_slice_score   = (Float_t)sliceScore;
    v_true_energy   = 0.0f;

    bool isSignal = isTrueSignal && (TMath::Abs(sliceInter - 1098.0) < 0.5);

    if (isSignal) {
      if (b.truth_recoilElectronEnergy && !b.truth_recoilElectronEnergy->empty())
        v_true_energy = (Float_t)((*b.truth_recoilElectronEnergy)[0] / 1000.0);
      sigTree->Fill();
      ++nSigFlat;
    } else if (!isTrueSignal) {
      bkgTree->Fill();
      ++nBkgFlat;
    }
    // nuEScatter==1 but sliceInter!=1098: genuinely ambiguous, excluded from training
  }

  printf("\n[INFO] Flat tree: %lld signal, %lld background\n", nSigFlat, nBkgFlat);
  fFlat->Write();

  // ── Phase 2: Train TMVA classifiers ──────────────────────────────────────────
  printf("\nPhase 2: Training TMVA classifiers...\n\n");

  TMVA::Tools::Instance();

  TString tmvaOutPath = TString(outDir) + "/TMVA_output.root";
  TFile* fTMVA = TFile::Open(tmvaOutPath, "RECREATE");
  if (!fTMVA || fTMVA->IsZombie()) {
    std::cerr << "[ERROR] Cannot create: " << tmvaOutPath << "\n";
    fFlat->Close(); fIn->Close(); return;
  }

  auto* factory = new TMVA::Factory("NuMMClassifier", fTMVA,
      "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

  TString datasetDir = TString(outDir) + "/tmva_dataset";
  auto* loader = new TMVA::DataLoader(datasetDir.Data());

  // Variables — ORDER MUST MATCH NuMMSelection_MVA.C exactly
  loader->AddVariable("razzle_p_elec", 'F');
  loader->AddVariable("track_score",   'F');
  loader->AddVariable("shower_energy", 'F');
  loader->AddVariable("dedx",          'F');
  loader->AddVariable("has_dedx",      'F');
  loader->AddVariable("theta",         'F');
  loader->AddVariable("etheta2",       'F');
  loader->AddVariable("n_extra_prim",  'F');
  loader->AddVariable("slice_score",   'F');

  loader->AddSpectator("true_energy", 'F');

  loader->AddSignalTree(sigTree,     1.0);
  loader->AddBackgroundTree(bkgTree, 1.0);

  // Normalise class weights so signal and background contribute equally.
  // Background outnumbers signal ~13:1; NormMode=NumEvents corrects for this.
  loader->PrepareTrainingAndTestTree("", "",
      "nTrain_Signal=0:nTrain_Background=0:"
      "SplitMode=Random:SplitSeed=42:NormMode=NumEvents:!V");

  // BDTG — gradient-boosted trees, primary method
  factory->BookMethod(loader, TMVA::Types::kBDT, "BDTG",
      "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=Grad:"
      "Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:"
      "nCuts=20:MaxDepth=3");

  // BDT — AdaBoost, faster training, useful baseline for comparison
  factory->BookMethod(loader, TMVA::Types::kBDT, "BDT",
      "!H:!V:NTrees=500:MinNodeSize=5%:MaxDepth=3:"
      "BoostType=AdaBoost:AdaBoostBeta=0.5:nCuts=20");

  // MLP — neural network, can capture non-linear correlations missed by trees
  factory->BookMethod(loader, TMVA::Types::kMLP, "MLP",
      "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=800:"
      "HiddenLayers=N+5,N:TestRate=5:!UseRegulator");

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  fTMVA->Close();
  fFlat->Close();
  fIn->Close();

  printf("\n══════════════════════════════════════════════════════════════════\n");
  printf("  Training complete.\n");
  printf("  Weights : %s/weights/NuMMClassifier_BDTG.weights.xml\n", datasetDir.Data());
  printf("  TMVA out: %s\n", tmvaOutPath.Data());
  printf("──────────────────────────────────────────────────────────────────\n");
  printf("  Inspect ROC curves and working-point table:\n");
  printf("    root -l '%s/weights/TMVAGui.C'\n", datasetDir.Data());
  printf("  Then choose MVA_CUT in NuMMSelection_MVA.C and run:\n");
  printf("    root -l -b -q 'magnetic_moment/NuMMSelection_MVA.C(\"%s\")'\n", inputFile);
  printf("══════════════════════════════════════════════════════════════════\n\n");
}
