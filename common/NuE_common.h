// ============================================================================
//  NuE_common.h
//  ---------------------------------------------------------------------------
//  Single source of truth for the selection infrastructure shared by every
//  macro that reads the SBND merged MC tree (ana/NuE):
//
//     NuESelection.C, NuESkim.C            (top-level Weinberg-angle pipeline)
//     magnetic_moment/NuMMSelection.C      (NMM box-cut selection)
//     magnetic_moment/NuMM_common.h        (NMM MVA pre-selection: Classifier + MVA)
//
//  Before this header existed, each of those carried its own byte-for-byte copy
//  of the blocks below. They are now defined here exactly once.
//
//  What lives here (identical for all analyses):
//     Sentinel            - missing-value sentinels / validity checks
//     FV                  - SBND fiducial volume + inFV()
//     ScaleFactors        - per-sample POT normalisation to 1e21 POT
//     Preselection        - slice/candidate thresholds used by FindNuSliceIndex
//     BranchVars          - superset of every ana/NuE branch any macro reads
//     SetBranchAddresses  - wires BranchVars (only the listed branches are read)
//     GetPCAAngle         - PCA-recalculated slice angle lookup
//     GetRecoNuVertex     - genuine reco neutrino vertex lookup (for the FV cut)
//     FindNuSliceIndex    - best pandora neutrino-candidate slice (category==1)
//
//  What stays in each macro (genuinely analysis-specific, do NOT move here):
//     Cuts / Pre / PreCuts   - the tuned working-point thresholds
//     Binning                - reco-energy axis (120/0-1.2 vs 60/0-0.6 GeV)
//     CutStep enum + labels  - cut-flow bookkeeping
//     ElecCandidate + FindElecCandidate - box-cut (Double_t) vs MVA (Float_t)
//
//  Include path (GCC resolves "" includes relative to the including file):
//     top level         : #include "common/NuE_common.h"
//     magnetic_moment/  : #include "../common/NuE_common.h"
// ============================================================================

#ifndef NUE_COMMON_H
#define NUE_COMMON_H

#include <vector>
#include "TTree.h"

// ─────────────────────────── SENTINEL VALUES ─────────────────────────────────
// Branches are set to -999999 when missing; dEdx uses a -999 sentinel.
namespace Sentinel {
  const Double_t INVALID  = -999000.0;
  const Double_t DEDX_INV = -990.0;
  inline bool valid(Double_t x)     { return x > INVALID; }
  inline bool validDedx(Double_t x) { return x > DEDX_INV; }
}

// ─────────────────────── SBND FIDUCIAL VOLUME ────────────────────────────────
// Active volume: X ∈ [−200,200], Y ∈ [−200,200], Z ∈ [0,500] cm.
// FV margins: X ±24 cm, Y ±17 cm, Z front 5 cm / back 10 cm.
namespace FV {
  const Double_t X_MIN = -176.0, X_MAX =  176.0;  // cm
  const Double_t Y_MIN = -183.0, Y_MAX =  183.0;  // cm
  const Double_t Z_MIN =    5.0, Z_MAX =  490.0;  // cm

  inline bool inFV(Double_t x, Double_t y, Double_t z) {
    return x > X_MIN && x < X_MAX &&
           y > Y_MIN && y < Y_MAX &&
           z > Z_MIN && z < Z_MAX;
  }
}

// ──────────────── POT NORMALISATION FACTORS (to 1×10²¹ POT) ─────────────────
// Signal and neutrino-background MC come from separate productions with
// different effective POT, so each is scaled by its own factor.
namespace ScaleFactors {
  const Double_t TARGET_POT  = 1.0e21;
  const Double_t SIGNAL      = 0.0105845;  // signal MC → 1e21 POT
  const Double_t NU_BKG      = 1.65585;    // neutrino background MC → 1e21 POT
  const Double_t COSMIC_BKG  = 0.85935;    // in-time cosmic MC → 1e21 POT
}

// ─────────────── SHARED PRESELECTION THRESHOLDS ──────────────────────────────
// Used by FindNuSliceIndex below. The analysis-specific working-point cuts live
// in each macro's own Cuts/Pre/PreCuts namespace.
namespace Preselection {
  const Double_t SLICE_SCORE_MIN = -1.5;  // exclude truly invalid slices
}

// ─────────────────── BRANCH VARIABLE BLOCK ───────────────────────────────────
// Superset of every ana/NuE branch any selection macro reads. SetBranchAddresses
// disables all branches first and re-enables only these, so carrying unused
// fields here costs nothing at read time.
struct BranchVars {
  // Event ID
  UInt_t    eventID, runID, subRunID;

  // Truth
  Int_t     nuEScatter;
  Double_t  nuEScatterTrueVX, nuEScatterTrueVY, nuEScatterTrueVZ;
  std::vector<Double_t>* truth_recoilElectronEnergy = nullptr;

  // Slice-level
  std::vector<Double_t>* reco_sliceID          = nullptr;
  std::vector<Double_t>* reco_sliceCategory    = nullptr;
  std::vector<Double_t>* reco_sliceScore       = nullptr;
  std::vector<Double_t>* reco_sliceInteraction = nullptr;
  std::vector<Double_t>* reco_sliceTrueVX      = nullptr;
  std::vector<Double_t>* reco_sliceTrueVY      = nullptr;
  std::vector<Double_t>* reco_sliceTrueVZ      = nullptr;

  // Particle-level
  std::vector<Double_t>* reco_particleSliceID               = nullptr;
  std::vector<Double_t>* reco_particlePDG                   = nullptr;
  std::vector<Double_t>* reco_particleIsPrimary             = nullptr;
  std::vector<Double_t>* reco_particleClearCosmic           = nullptr;
  std::vector<Double_t>* reco_particleVX                    = nullptr;
  std::vector<Double_t>* reco_particleVY                    = nullptr;
  std::vector<Double_t>* reco_particleVZ                    = nullptr;
  std::vector<Double_t>* reco_particleTrackScore            = nullptr;
  std::vector<Double_t>* reco_particleShowerBestPlaneEnergy = nullptr;
  std::vector<Double_t>* reco_particleTheta                 = nullptr;
  std::vector<Double_t>* reco_particleBestPlanedEdx         = nullptr;
  std::vector<Double_t>* reco_particlePlane0dEdx            = nullptr;
  std::vector<Double_t>* reco_particlePlane1dEdx            = nullptr;
  std::vector<Double_t>* reco_particlePlane2dEdx            = nullptr;
  std::vector<Double_t>* reco_particleRazzledPDG11          = nullptr;
  std::vector<Double_t>* reco_particleRazzledPDG13          = nullptr;
  std::vector<Double_t>* reco_particleRazzledPDG22          = nullptr;
  std::vector<Double_t>* reco_particleRazzledPDG211         = nullptr;
  std::vector<Double_t>* reco_particleRazzledPDG2212        = nullptr;
  std::vector<Double_t>* reco_particleRazzledBestPDG        = nullptr;
  std::vector<Double_t>* reco_particleShowerLength          = nullptr;
  std::vector<Double_t>* reco_particleShowerOpenAngle       = nullptr;

  // Reconstructed neutrino vertex (genuine reco, not truth-matched)
  std::vector<Double_t>* reco_neutrinoVX      = nullptr;
  std::vector<Double_t>* reco_neutrinoVY      = nullptr;
  std::vector<Double_t>* reco_neutrinoVZ      = nullptr;
  std::vector<Double_t>* reco_neutrinoSliceID = nullptr;

  // PCA-recalculated angles per slice
  std::vector<Double_t>* pcaSlice_angle   = nullptr;
  std::vector<Double_t>* pcaSlice_sliceID = nullptr;
};

// Disables every branch, then enables and wires only the ones above. This is
// the read-speed optimisation that the 143 GB skim pass needs; it is a no-op
// for correctness because every macro reads the tree solely through BranchVars.
inline void SetBranchAddresses(TTree* t, BranchVars& b) {
  t->SetBranchStatus("*", 0);
  auto on = [&](const char* n, void* p) {
    t->SetBranchStatus(n, 1);
    t->SetBranchAddress(n, p);
  };
  on("eventID",                            &b.eventID);
  on("runID",                              &b.runID);
  on("subRunID",                           &b.subRunID);
  on("nuEScatter",                         &b.nuEScatter);
  on("nuEScatterTrueVX",                   &b.nuEScatterTrueVX);
  on("nuEScatterTrueVY",                   &b.nuEScatterTrueVY);
  on("nuEScatterTrueVZ",                   &b.nuEScatterTrueVZ);
  on("truth_recoilElectronEnergy",         &b.truth_recoilElectronEnergy);
  on("reco_sliceID",                       &b.reco_sliceID);
  on("reco_sliceCategory",                 &b.reco_sliceCategory);
  on("reco_sliceScore",                    &b.reco_sliceScore);
  on("reco_sliceInteraction",              &b.reco_sliceInteraction);
  on("reco_sliceTrueVX",                   &b.reco_sliceTrueVX);
  on("reco_sliceTrueVY",                   &b.reco_sliceTrueVY);
  on("reco_sliceTrueVZ",                   &b.reco_sliceTrueVZ);
  on("reco_particleSliceID",               &b.reco_particleSliceID);
  on("reco_particlePDG",                   &b.reco_particlePDG);
  on("reco_particleIsPrimary",             &b.reco_particleIsPrimary);
  on("reco_particleClearCosmic",           &b.reco_particleClearCosmic);
  on("reco_particleVX",                    &b.reco_particleVX);
  on("reco_particleVY",                    &b.reco_particleVY);
  on("reco_particleVZ",                    &b.reco_particleVZ);
  on("reco_particleTrackScore",            &b.reco_particleTrackScore);
  on("reco_particleShowerBestPlaneEnergy", &b.reco_particleShowerBestPlaneEnergy);
  on("reco_particleTheta",                 &b.reco_particleTheta);
  on("reco_particleBestPlanedEdx",         &b.reco_particleBestPlanedEdx);
  on("reco_particlePlane0dEdx",            &b.reco_particlePlane0dEdx);
  on("reco_particlePlane1dEdx",            &b.reco_particlePlane1dEdx);
  on("reco_particlePlane2dEdx",            &b.reco_particlePlane2dEdx);
  on("reco_particleRazzledPDG11",          &b.reco_particleRazzledPDG11);
  on("reco_particleRazzledPDG13",          &b.reco_particleRazzledPDG13);
  on("reco_particleRazzledPDG22",          &b.reco_particleRazzledPDG22);
  on("reco_particleRazzledPDG211",         &b.reco_particleRazzledPDG211);
  on("reco_particleRazzledPDG2212",        &b.reco_particleRazzledPDG2212);
  on("reco_particleRazzledBestPDG",        &b.reco_particleRazzledBestPDG);
  on("reco_particleShowerLength",          &b.reco_particleShowerLength);
  on("reco_particleShowerOpenAngle",       &b.reco_particleShowerOpenAngle);
  on("reco_neutrinoVX",                    &b.reco_neutrinoVX);
  on("reco_neutrinoVY",                    &b.reco_neutrinoVY);
  on("reco_neutrinoVZ",                    &b.reco_neutrinoVZ);
  on("reco_neutrinoSliceID",               &b.reco_neutrinoSliceID);
  on("angleRecalculationPCASlice_angle",   &b.pcaSlice_angle);
  on("angleRecalculationPCASlice_sliceID", &b.pcaSlice_sliceID);
}

// ──────────────────────── HELPER FUNCTIONS ───────────────────────────────────

// PCA-recalculated angle for a given slice ID, or Sentinel::INVALID if absent.
inline Double_t GetPCAAngle(const BranchVars& b, Double_t sliceID) {
  for (Int_t i = 0; i < (Int_t)b.pcaSlice_sliceID->size(); ++i)
    if ((*b.pcaSlice_sliceID)[i] == sliceID)
      return (*b.pcaSlice_angle)[i];
  return Sentinel::INVALID;
}

// Genuine reconstructed neutrino vertex for the given slice (NOT truth-matched).
// Returns true and fills x,y,z if a valid reco neutrino is associated with the
// slice; false otherwise. Replaces use of reco_sliceTrueVX, a truth-matched
// quantity unavailable in real data and filled for only ~38% of signal events.
inline bool GetRecoNuVertex(const BranchVars& b, Double_t sliceID,
                            Double_t& x, Double_t& y, Double_t& z) {
  if (!b.reco_neutrinoSliceID) return false;
  for (Int_t i = 0; i < (Int_t)b.reco_neutrinoSliceID->size(); ++i) {
    if ((*b.reco_neutrinoSliceID)[i] != sliceID) continue;
    x = (*b.reco_neutrinoVX)[i];
    y = (*b.reco_neutrinoVY)[i];
    z = (*b.reco_neutrinoVZ)[i];
    return Sentinel::valid(x) && Sentinel::valid(y) && Sentinel::valid(z);
  }
  return false;
}

// Index of the pandora neutrino-candidate slice (category==1) with the highest
// reco_sliceScore. Returns -1 if none qualify.
inline Int_t FindNuSliceIndex(const BranchVars& b) {
  Int_t    best      = -1;
  Double_t bestScore = Sentinel::INVALID;
  for (Int_t s = 0; s < (Int_t)b.reco_sliceCategory->size(); ++s) {
    if ((*b.reco_sliceCategory)[s] < 0.5)                        continue; // not neutrino-like
    if (!Sentinel::valid((*b.reco_sliceScore)[s]))               continue; // no valid score
    if ((*b.reco_sliceScore)[s] < Preselection::SLICE_SCORE_MIN) continue; // too cosmic-like
    if ((*b.reco_sliceScore)[s] > bestScore) {
      bestScore = (*b.reco_sliceScore)[s];
      best = s;
    }
  }
  return best;
}

#endif // NUE_COMMON_H
