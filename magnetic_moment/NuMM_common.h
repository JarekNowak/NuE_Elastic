// NuMM_common.h
// Shared pre-selection infrastructure for NuMMClassifier.C and NuMMSelection_MVA.C.
// Include this header in both macros to keep the pre-selection logic identical.

#ifndef NUMM_COMMON_H
#define NUMM_COMMON_H

#include <vector>
#include "TTree.h"
#include "../common/NuE_common.h"

// Loose thresholds used only to identify the electron candidate — no box cuts.
// SHOWER_E_MIN/MAX match the analysis range so training and application see the
// same background population; training on out-of-range events dilutes the MVA.
namespace PreCuts {
  const Double_t SLICE_SCORE_MIN  = -1.5;
  const Int_t    RAZZLE_BEST_PDG  =  11;
  const Double_t SHOWER_E_MIN     =  0.010;  // GeV — 10 MeV (NMM analysis floor)
  const Double_t SHOWER_E_MAX     =  0.600;  // GeV — 600 MeV (NMM analysis ceiling)
  const Double_t EXTRA_PRIM_E_MIN =  0.020;  // GeV — threshold for counting extra primaries
}

namespace Binning {
  const Int_t    N_RECO = 60;
  const Double_t E_LO   = 0.0, E_HI = 0.6;  // GeV
}

// Per-event electron candidate with all raw variables for MVA input.
// No box cuts are applied here; variables are passed directly to the classifier.
struct ElecCandidate {
  Int_t    pIdx       = -1;
  Float_t  energy     = 0.0f;   // GeV
  Float_t  theta      = -1.0f;  // rad
  Float_t  dedx       = -1.0f;  // MeV/cm; -1 = invalid
  Float_t  hasDedx    = 0.0f;   // 1 = valid, 0 = invalid (flag for MVA)
  Float_t  razzle11   = -1.0f;
  Float_t  trackScore = -1.0f;
  Float_t  nExtraPrim = 0.0f;
  Double_t VX = 0.0, VY = 0.0, VZ = 0.0;
  bool     found      = false;
};

// Picks the primary particle with razzledBestPDG=11 that has the highest
// razzle11 score. No box cuts on individual variables — the MVA does that.
inline ElecCandidate FindElecCandidate(const BranchVars& b,
                                       Double_t sliceID,
                                       Double_t pcaAngle) {
  ElecCandidate best;

  for (Int_t p = 0; p < (Int_t)b.reco_particlePDG->size(); ++p) {
    if ((*b.reco_particleSliceID)[p]  != sliceID) continue;
    if ((*b.reco_particleIsPrimary)[p]  < 0.5)    continue;
    if ((*b.reco_particleClearCosmic)[p] > 0.5)   continue;

    Double_t rz11   = (*b.reco_particleRazzledPDG11)[p];
    Double_t rzBest = (*b.reco_particleRazzledBestPDG)[p];
    Double_t E_MeV  = (*b.reco_particleShowerBestPlaneEnergy)[p];

    if (!Sentinel::valid(rz11))  continue;
    if ((Int_t)(rzBest + 0.5) != PreCuts::RAZZLE_BEST_PDG) continue;
    if (!Sentinel::valid(E_MeV)) continue;
    Double_t E_GeV = E_MeV / 1000.0;
    if (E_GeV < PreCuts::SHOWER_E_MIN || E_GeV > PreCuts::SHOWER_E_MAX) continue;

    if (!best.found || rz11 > best.razzle11) {
      Double_t trkSc  = (*b.reco_particleTrackScore)[p];
      Double_t dx     = (*b.reco_particleBestPlanedEdx)[p];
      Double_t ang    = Sentinel::valid(pcaAngle) ? pcaAngle
                                                  : (*b.reco_particleTheta)[p];
      best.found      = true;
      best.pIdx       = p;
      best.energy     = (Float_t)E_GeV;
      best.razzle11   = (Float_t)rz11;
      best.trackScore = Sentinel::valid(trkSc) ? (Float_t)trkSc : -1.0f;
      best.theta      = (Float_t)ang;
      best.dedx       = Sentinel::validDedx(dx) ? (Float_t)dx : -1.0f;
      best.hasDedx    = Sentinel::validDedx(dx) ? 1.0f : 0.0f;
      best.VX = (*b.reco_particleVX)[p];
      best.VY = (*b.reco_particleVY)[p];
      best.VZ = (*b.reco_particleVZ)[p];
    }
  }

  if (!best.found) return best;

  best.nExtraPrim = 0.0f;
  for (Int_t p = 0; p < (Int_t)b.reco_particlePDG->size(); ++p) {
    if (p == best.pIdx)                           continue;
    if ((*b.reco_particleSliceID)[p] != sliceID) continue;
    if ((*b.reco_particleIsPrimary)[p]  < 0.5)   continue;
    if ((*b.reco_particleClearCosmic)[p] > 0.5)  continue;
    Double_t E_MeV = (*b.reco_particleShowerBestPlaneEnergy)[p];
    if (!Sentinel::valid(E_MeV)) continue;
    if (E_MeV / 1000.0 > PreCuts::EXTRA_PRIM_E_MIN) best.nExtraPrim += 1.0f;
  }
  return best;
}

#endif // NUMM_COMMON_H
