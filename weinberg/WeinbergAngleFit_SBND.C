/*=============================================================================
  WeinbergAngleFit_SBND.C
  ─────────────────────────────────────────────────────────────────────────────
  ROOT macro to fit the weak mixing angle (sin²θ_W) from
  muon-neutrino / electron elastic scattering:

        ν_μ  +  e⁻  →  ν_μ  +  e⁻

  at the SBND detector (BNB neutrino beam, Fermilab).

  Physics:
    dσ/dT = (G_F² m_e)/(2π) · [ gL²  +  gR²·(1 − T/Eν)²  −  gL·gR·(m_e T / Eν²) ]

    gL = −½ + sin²θ_W     (left-handed coupling)
    gR =      sin²θ_W     (right-handed coupling)
    T  = recoil electron kinetic energy  [GeV]
    Eν = neutrino energy                 [GeV]

  Inputs  (histogram names defined in the constants block below):
    • h_data         – observed recoil-electron energy spectrum (reco bins)
    • h_eff          – selection efficiency vs true recoil energy
    • h_smear        – 2-D smearing matrix  M[reco][true]  (normalised per true-energy column)
    • h_flux         – BNB ν_μ flux  dΦ/dEν  [/cm²/POT/GeV]
    • h_bkg          – (optional) total background prediction

  Output:
    • Best-fit sin²θ_W with statistical uncertainty
    • χ² / ndf
    • Comparison plot: data vs best-fit prediction

  Usage:
    root -l -b -q 'WeinbergAngleFit_SBND.C("inputs.root")'

  Author : Jarek Nowak
  Date   : 2026-03-28
  Fixed  : 2026-04-26
=============================================================================*/

#include <iostream>
#include <cmath>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TPad.h"
#include "TRandom3.h"

// ─────────────────────────────── CONSTANTS ───────────────────────────────────

namespace Physics {
  const Double_t GF     = 1.16638e-5;   // Fermi constant  [GeV⁻²]
  const Double_t me     = 0.511e-3;     // electron mass   [GeV]
  const Double_t hbarc2 = 3.8938e-28;  // ħc² in cm²·GeV²  (conv. GeV⁻² → cm²)
  // Standard Model prediction
  const Double_t sin2tW_SM = 0.2312;
}

namespace Setup {
  // Histogram names inside the input ROOT file
  const char* DATA_HIST   = "h_data";       // TH1D, reco-energy bins
  const char* EFF_HIST    = "h_eff";        // TH1D, true-energy bins
  const char* SMEAR_HIST  = "h_smear";      // TH2D  [reco][true]
  const char* FLUX_HIST   = "h_flux";       // TH1D  [Eν bins, GeV]
  const char* BKG_HIST    = "h_bkg";        // TH1D, reco-energy bins (optional)

  // SBND detector parameters
  const Double_t N_ELECTRONS = 2.437e31; // e⁻ targets in fiducial vol.
  const Double_t POT         = 1.0e21;   // Protons on target (3-yr run)

  // Integration sub-steps per T bin (increase for higher accuracy)
  // FIX [IMPROVEMENT]: was hardcoded magic number 4 inside BuildTrueSpectrum
  const Int_t    XSEC_NSUB  = 8;
    
  const Int_t FIT_BIN_MIN = 10;  // Skip bins 1–9 (threshold noise / very low-T region)
  const Int_t FIT_BIN_MAX = 40;  // Skip bins above 40 (low statistics at high T)
}

// ─────────────────── GLOBAL POINTERS (for TMinuit callback) ──────────────────

TH1D* g_hData   = nullptr;
TH1D* g_hEff    = nullptr;
TH1D* g_hBkg    = nullptr;
TH2D* g_hSmear  = nullptr;
TH1D* g_hFlux   = nullptr;

// ─────────────────────── FORWARD DECLARATIONS ────────────────────────────────
// FIX [BUG-5]: Missing forward declarations caused undefined-function errors
// when compiling with ACLiC (.C+) or strict parsers.

Double_t dsigma_dT(Double_t T, Double_t Enu, Double_t sin2thW);
Double_t T_max(Double_t Enu);
TH1D*    BuildTrueSpectrum(Double_t sin2thW);
TH1D*    SmearToReco(TH1D* hTrue);
TH1D*    BuildRecoPrediction(Double_t sin2thW);
void     FCN(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t iflag);
void     GeneratePlaceholderHistograms();
void     DrawResult(Double_t sin2thW_fit, Double_t sin2thW_err,Double_t chi2, Int_t ndf);
void     DrawDeltaChi2Scan(Double_t sin2thW_fit, Double_t sin2thW_err);

// ─────────────────────────── CROSS-SECTION ───────────────────────────────────

/// Differential cross-section dσ/dT [cm²/GeV] for νμ-e elastic scattering.
/// @param T       recoil electron kinetic energy [GeV]
/// @param Enu     neutrino energy                [GeV]
/// @param sin2thW sin²θ_W
Double_t dsigma_dT(Double_t T, Double_t Enu, Double_t sin2thW) {
  if (T < 0 || T > Enu || Enu <= 0) return 0.0;

  const Double_t gL = -0.5 + sin2thW;
  const Double_t gR =        sin2thW;
  const Double_t y  = T / Enu;

  // Kinematic maximum of y: y_max = 1 / (1 + me/(2Eν))
  if (y < 0 || y > 1.0 / (1.0 + 0.5 * Physics::me / Enu)) return 0.0;

  Double_t term1 = gL * gL;
  Double_t term2 = gR * gR * (1.0 - y) * (1.0 - y);
  Double_t term3 = gL * gR * Physics::me * T / (Enu * Enu);

  Double_t dsig = (Physics::GF * Physics::GF * Physics::me) / (2.0 * TMath::Pi())
                  * (term1 + term2 - term3)
                  * Physics::hbarc2;  // → cm²/GeV
  return std::max(dsig, 0.0);
}

/// Maximum recoil energy kinematically allowed:  T_max = 2Eν² / (me + 2Eν)
Double_t T_max(Double_t Enu) {
  return Enu / (1.0 + Physics::me / (2.0 * Enu));
}

// ──────────────────────── PREDICTION BUILDER ─────────────────────────────────

/// Build predicted true-T spectrum for a given sin²θ_W.
/// Returns a TH1D* (owned by caller) in the true-energy binning of g_hEff.
TH1D* BuildTrueSpectrum(Double_t sin2thW) {
  TH1D* hTrue = (TH1D*)g_hEff->Clone("h_true_pred");
  hTrue->Reset();
  hTrue->SetDirectory(nullptr);

  const Int_t nT   = hTrue->GetNbinsX();
  const Int_t nEnu = g_hFlux->GetNbinsX();

  for (Int_t iT = 1; iT <= nT; ++iT) {
    Double_t T_lo = hTrue->GetBinLowEdge(iT);
    Double_t T_hi = T_lo + hTrue->GetBinWidth(iT);
    Double_t T_c  = 0.5 * (T_lo + T_hi);

    Double_t rate = 0.0;

    for (Int_t iE = 1; iE <= nEnu; ++iE) {
      Double_t Enu  = g_hFlux->GetBinCenter(iE);
      Double_t dEnu = g_hFlux->GetBinWidth(iE);
      Double_t flux = g_hFlux->GetBinContent(iE); // /cm²/POT/GeV

      if (T_c > T_max(Enu)) continue;

      // Integrate dσ/dT over the T bin using trapezoidal rule
      // FIX [IMPROVEMENT]: nSub promoted from magic literal 4 → Setup::XSEC_NSUB
      Double_t integral = 0.0;
      for (Int_t k = 0; k < Setup::XSEC_NSUB; ++k) {
        Double_t T1 = T_lo + (T_hi - T_lo) * k         / Setup::XSEC_NSUB;
        Double_t T2 = T_lo + (T_hi - T_lo) * (k + 1.0) / Setup::XSEC_NSUB;
        integral += 0.5 * (dsigma_dT(T1, Enu, sin2thW)
                         + dsigma_dT(T2, Enu, sin2thW))
                    * (T2 - T1);
      }
      // flux   [/cm²/POT/GeV] × dEnu [GeV]  →  [/cm²/POT]
      // xsec   [cm²/GeV]      × dT   [GeV]  →  [cm²]
      // product                             →  dimensionless / POT / target
      rate += flux * dEnu * integral;
    }

    Double_t nEvents = rate * Setup::N_ELECTRONS * Setup::POT;
    Double_t eff     = g_hEff->GetBinContent(iT);
    hTrue->SetBinContent(iT, nEvents * eff);
  }

  return hTrue;
}

/// Smear a true-energy spectrum to reco-energy spectrum using the smearing matrix.
/// M[reco][true] is normalised such that ΣM[:,j] = 1 for each true bin j.
TH1D* SmearToReco(TH1D* hTrue) {
  const Int_t nReco = g_hSmear->GetNbinsX();
  const Int_t nTrue = g_hSmear->GetNbinsY();

  TH1D* hReco = (TH1D*)g_hData->Clone("h_reco_pred");
  hReco->Reset();
  hReco->SetDirectory(nullptr);

  for (Int_t iR = 1; iR <= nReco; ++iR) {
    Double_t val = 0.0;
    for (Int_t iT = 1; iT <= nTrue; ++iT) {
      val += g_hSmear->GetBinContent(iR, iT) * hTrue->GetBinContent(iT);
    }
    hReco->SetBinContent(iR, val);
  }
  return hReco;
}

/// Full prediction in reco space (signal + background).
TH1D* BuildRecoPrediction(Double_t sin2thW) {
  TH1D* hTrue = BuildTrueSpectrum(sin2thW);
  TH1D* hReco = SmearToReco(hTrue);
  delete hTrue;

  if (g_hBkg) {
    hReco->Add(g_hBkg);
  }
  return hReco;
}

// ─────────────────── CHI-SQUARED FUNCTION FOR MINUIT ─────────────────────────

/// χ² function passed to TMinuit (Neyman χ² with data errors).
void FCN(Int_t& npar, Double_t* grad, Double_t& fval,
         Double_t* par, Int_t iflag) {
  Double_t sin2thW = par[0];

  if (sin2thW < 0.0 || sin2thW > 0.5) {
    fval = 1e30;
    return;
  }

  TH1D* hPred = BuildRecoPrediction(sin2thW);

  Double_t chi2  = 0.0;
  Int_t    nBins = g_hData->GetNbinsX();

  for (Int_t i = 1; i <= nBins; ++i) {
    
    if (Setup::FIT_BIN_MIN > 0 && i < Setup::FIT_BIN_MIN) continue;
    if (Setup::FIT_BIN_MAX > 0 && i > Setup::FIT_BIN_MAX) continue;
      
    Double_t obs  = g_hData->GetBinContent(i);
    Double_t pred = hPred->GetBinContent(i);
    Double_t err  = g_hData->GetBinError(i);

    // Skip truly empty bins (no data, no usable uncertainty)
    if (err <= 0 && obs <= 0) continue;

    // FIX [BUG-2]: was "err = 1.0" — should be sqrt(obs) for Poisson statistics.
    // The comment said "use Poisson" but the original code set err=1 unconditionally.
    if (err <= 0) err = TMath::Sqrt(obs);
    if (err <= 0) err = 1.0;  // last-resort guard for truly zero-obs bins

    chi2 += (obs - pred) * (obs - pred) / (err * err);
  }

  fval = chi2;
  delete hPred;
}

// ──────────────────── PLACEHOLDER HISTOGRAM GENERATOR ────────────────────────

void GeneratePlaceholderHistograms() {
  std::cout << "[INFO] Generating placeholder histograms for demonstration.\n"
            << "       Replace with real SBND data when available.\n";

  const Int_t    nT  = 20;
  const Double_t Tlo = 0.0, Thi = 1.0;   // [GeV]
  const Int_t    nE  = 40;
  const Double_t Elo = 0.0, Ehi = 3.0;   // [GeV]

  // ── Flux (BNB ν_μ, simplified Gaussian peak) ──────────────────────────────
  // FIX [BUG-3]: original code allocated temporaries before calling here and
  // then leaked them when this function overwrote the global pointers.
  // Now this function is solely responsible for all global histogram creation.
  g_hFlux = new TH1D("h_flux",
    "BNB #nu_{#mu} flux;E_{#nu} [GeV];#Phi [cm^{-2}POT^{-1}GeV^{-1}]",
    nE, Elo, Ehi);
  for (Int_t i = 1; i <= nE; ++i) {
    Double_t E   = g_hFlux->GetBinCenter(i);
    Double_t phi = 1.2e-11 * TMath::Gaus(E, 0.7, 0.25);
    g_hFlux->SetBinContent(i, phi);
  }

  // ── Efficiency (flat-ish, falling at high T) ──────────────────────────────
  g_hEff = new TH1D("h_eff", "Efficiency;T_{e} [GeV];#epsilon", nT, Tlo, Thi);
  for (Int_t i = 1; i <= nT; ++i) {
    Double_t T   = g_hEff->GetBinCenter(i);
    Double_t eff = 0.70 * TMath::Exp(-0.5 * T);
    g_hEff->SetBinContent(i, eff);
  }

  // ── Smearing matrix (Gaussian blur σ ≈ 5% of true + 1% floor) ────────────
  g_hSmear = new TH2D("h_smear",
    "Smearing matrix;T_{e}^{reco} [GeV];T_{e}^{true} [GeV]",
    nT, Tlo, Thi, nT, Tlo, Thi);
  for (Int_t iT = 1; iT <= nT; ++iT) {
    Double_t Ttrue  = g_hSmear->GetYaxis()->GetBinCenter(iT);
    Double_t sigmaT = 0.05 * Ttrue + 0.01;
    Double_t colSum = 0.0;
    for (Int_t iR = 1; iR <= nT; ++iR) {
      Double_t Treco = g_hSmear->GetXaxis()->GetBinCenter(iR);
      Double_t val   = TMath::Gaus(Treco, Ttrue, sigmaT);
      g_hSmear->SetBinContent(iR, iT, val);
      colSum += val;
    }
    if (colSum > 0) {
      for (Int_t iR = 1; iR <= nT; ++iR)
        g_hSmear->SetBinContent(iR, iT,
          g_hSmear->GetBinContent(iR, iT) / colSum);
    }
  }

  // ── Signal prediction at SM value ────────────────────────────────────────
  Double_t sin2thW_true = Physics::sin2tW_SM;
  TH1D* hSignal  = BuildTrueSpectrum(sin2thW_true);
  TH1D* hRecoSig = SmearToReco(hSignal);
  delete hSignal;

  // ── Background (10% of signal, flat) ─────────────────────────────────────
  g_hBkg = new TH1D("h_bkg", "Background;T_{e}^{reco} [GeV];Events",
                    nT, Tlo, Thi);
  for (Int_t i = 1; i <= nT; ++i)
    g_hBkg->SetBinContent(i, 0.10 * hRecoSig->GetBinContent(i));

  // ── Data = signal + bkg + Poisson fluctuations ────────────────────────────
  g_hData = new TH1D("h_data", "SBND data;T_{e}^{reco} [GeV];Events / bin",
                     nT, Tlo, Thi);
  TRandom3 rng(42);
  for (Int_t i = 1; i <= nT; ++i) {
    Double_t mu  = hRecoSig->GetBinContent(i) + g_hBkg->GetBinContent(i);
    Double_t obs = (mu > 0) ? rng.Poisson(mu) : 0.0;
    g_hData->SetBinContent(i, obs);
    // FIX [IMPROVEMENT]: use sqrt(obs) for error; protect against obs==0
    g_hData->SetBinError(i, (obs > 0) ? TMath::Sqrt(obs) : 1.0);
  }

  delete hRecoSig;
}

// ─────────────────────────── PLOTTING HELPERS ────────────────────────────────

void DrawResult(Double_t sin2thW_fit, Double_t sin2thW_err,
                Double_t chi2, Int_t ndf) {

  TH1D* hBestFit = BuildRecoPrediction(sin2thW_fit);
  hBestFit->SetLineColor(kRed + 1);
  hBestFit->SetLineWidth(2);
  hBestFit->SetFillColorAlpha(kRed + 1, 0.15);
  hBestFit->SetFillStyle(1001);

  TH1D* hUp   = BuildRecoPrediction(sin2thW_fit + sin2thW_err);
  TH1D* hDown = BuildRecoPrediction(sin2thW_fit - sin2thW_err);
  hUp->SetLineColor(kRed - 7);   hUp->SetLineStyle(2);
  hDown->SetLineColor(kRed - 7); hDown->SetLineStyle(2);

  // FIX [BUG-1]: original code dereferenced g_hBkg without null check.
  // Background is documented as optional — must guard here.
  TH1D* hBkgPlot = nullptr;
  if (g_hBkg) {
    hBkgPlot = (TH1D*)g_hBkg->Clone("h_bkg_plot");
    hBkgPlot->SetFillColor(kAzure - 4);
    hBkgPlot->SetFillStyle(3354);
    hBkgPlot->SetLineColor(kAzure + 1);
  }

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TCanvas* c1 = new TCanvas("c_fit", "Weinberg Angle Fit — SBND", 900, 750);
  c1->SetLeftMargin(0.12);
  c1->Divide(1, 2, 0, 0);

  // ── Upper pad: spectrum comparison ───────────────────────────────────────
  TPad* pTop = (TPad*)c1->cd(1);
  pTop->SetPad(0, 0.30, 1, 1.00);
  pTop->SetLeftMargin(0.12);
  pTop->SetBottomMargin(0.02);
  pTop->SetTopMargin(0.08);
  pTop->Draw();
  pTop->cd();

  g_hData->SetMarkerStyle(20);
  g_hData->SetMarkerSize(0.9);
  g_hData->SetLineColor(kBlack);
  g_hData->GetYaxis()->SetTitle("Events / bin");
  g_hData->GetYaxis()->SetTitleSize(0.06);
  g_hData->GetYaxis()->SetLabelSize(0.055);
  g_hData->GetXaxis()->SetLabelSize(0);
  g_hData->SetTitle("");
  g_hData->Draw("E1");

  hBestFit->Draw("HIST SAME");
  if (hBkgPlot) hBkgPlot->Draw("HIST SAME");
  hUp->Draw("HIST SAME");
  hDown->Draw("HIST SAME");
  g_hData->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.052);
  leg->AddEntry(g_hData,  "SBND data",          "lpe");
  leg->AddEntry(hBestFit, "Best-fit signal + bkg", "lf");
  if (hBkgPlot) leg->AddEntry(hBkgPlot, "Background", "f");
  leg->AddEntry(hUp,      "Best-fit #pm1#sigma", "l");
  leg->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.055);
  latex.DrawLatex(0.14, 0.88, "#bf{SBND}  #nu_{#mu}-e^{-} elastic scattering");
  latex.SetTextSize(0.050);
  latex.DrawLatex(0.14, 0.80,
    Form("sin^{2}#theta_{W} = %.4f #pm %.4f", sin2thW_fit, sin2thW_err));
  latex.DrawLatex(0.14, 0.73,
    Form("#chi^{2}/ndf = %.1f / %d  (p = %.3f)",
         chi2, ndf, TMath::Prob(chi2, ndf)));
  latex.SetTextColor(kGray + 1);
  latex.DrawLatex(0.14, 0.66,
    Form("SM: sin^{2}#theta_{W} = %.4f", Physics::sin2tW_SM));

  // ── Lower pad: residuals ──────────────────────────────────────────────────
  TPad* pBot = (TPad*)c1->cd(2);
  pBot->SetPad(0, 0.00, 1, 0.30);
  pBot->SetLeftMargin(0.12);
  pBot->SetTopMargin(0.02);
  pBot->SetBottomMargin(0.28);
  pBot->Draw();
  pBot->cd();

  Int_t  nBins  = g_hData->GetNbinsX();
  TH1D*  hResid = (TH1D*)g_hData->Clone("h_resid");
  hResid->Reset();
  hResid->SetTitle(";T_{e}^{reco} [GeV];(Data-Pred)/#sigma");
  hResid->GetXaxis()->SetTitleSize(0.13);
  hResid->GetXaxis()->SetLabelSize(0.12);
  hResid->GetYaxis()->SetTitleSize(0.10);
  hResid->GetYaxis()->SetLabelSize(0.10);
  hResid->GetYaxis()->SetTitleOffset(0.45);
  hResid->GetYaxis()->SetNdivisions(505);

  for (Int_t i = 1; i <= nBins; ++i) {
    Double_t obs  = g_hData->GetBinContent(i);
    Double_t pred = hBestFit->GetBinContent(i);
    Double_t err  = g_hData->GetBinError(i);
    if (err > 0)
      hResid->SetBinContent(i, (obs - pred) / err);
  }
  hResid->SetMarkerStyle(20);
  hResid->SetMarkerSize(0.8);
  hResid->SetLineColor(kBlack);
  hResid->GetYaxis()->SetRangeUser(-3.5, 3.5);
  hResid->Draw("E1");

  TLine* zeroLine = new TLine(g_hData->GetXaxis()->GetXmin(), 0,
                               g_hData->GetXaxis()->GetXmax(), 0);
  zeroLine->SetLineColor(kRed + 1); zeroLine->SetLineWidth(2); zeroLine->Draw();

  TLine* p1sig = new TLine(g_hData->GetXaxis()->GetXmin(),  1,
                            g_hData->GetXaxis()->GetXmax(),  1);
  TLine* m1sig = new TLine(g_hData->GetXaxis()->GetXmin(), -1,
                            g_hData->GetXaxis()->GetXmax(), -1);
  p1sig->SetLineStyle(2); p1sig->SetLineColor(kGray + 1); p1sig->Draw();
  m1sig->SetLineStyle(2); m1sig->SetLineColor(kGray + 1); m1sig->Draw();

  c1->SaveAs("WeinbergAngleFit_SBND.pdf");
  c1->SaveAs("WeinbergAngleFit_SBND.png");
  std::cout << "\n[INFO] Plots saved: WeinbergAngleFit_SBND.{pdf,png}\n";

  // FIX [BUG-4]: histograms drawn to the canvas are owned/referenced by it.
  // Deleting them here leaves dangling pointers inside the pad if the canvas
  // is later redrawn or accessed interactively.  Ownership is transferred to
  // the canvas by not deleting — they will be cleaned up when c1 is deleted or
  // goes out of scope.  hResid is the only one safe to release since it is not
  // reused after SaveAs.
  // (hBestFit, hUp, hDown, hBkgPlot are left to the canvas.)
}

// ─────────────────── PROFILE-LIKELIHOOD SCAN (OPTIONAL) ──────────────────────

void DrawDeltaChi2Scan(Double_t sin2thW_fit, Double_t sin2thW_err) {
  const Int_t    nScan = 100;
  const Double_t slo   = std::max(0.0, sin2thW_fit - 4 * sin2thW_err);
  const Double_t shi   = std::min(0.5, sin2thW_fit + 4 * sin2thW_err);

  // Reference chi2 at best-fit
  Int_t    npar_dummy = 1;
  Double_t grad_dummy = 0, chi2_min = 0;
  Double_t par_tmp[1] = {sin2thW_fit};
  Int_t    flag_dummy = 4;
  FCN(npar_dummy, &grad_dummy, chi2_min, par_tmp, flag_dummy);

  TGraph* gScan = new TGraph(nScan);
  gScan->SetTitle("#Delta#chi^{2} scan;sin^{2}#theta_{W};#Delta#chi^{2}");
  gScan->SetLineColor(kBlue + 1);
  gScan->SetLineWidth(2);

  for (Int_t i = 0; i < nScan; ++i) {
    Double_t sw  = slo + (shi - slo) * i / (nScan - 1);
    Double_t chi2 = 0, grad = 0;
    Double_t par[1] = {sw};
    FCN(npar_dummy, &grad, chi2, par, flag_dummy);
    gScan->SetPoint(i, sw, chi2 - chi2_min);
  }

  TCanvas* c2 = new TCanvas("c_scan", "#Delta#chi^{2} scan", 700, 500);
  c2->SetLeftMargin(0.13);
  gScan->Draw("AL");
  gScan->GetYaxis()->SetRangeUser(0, 10);

  Double_t xlo = gScan->GetXaxis()->GetXmin();
  Double_t xhi = gScan->GetXaxis()->GetXmax();
  TLine* l1 = new TLine(xlo, 1, xhi, 1);
  TLine* l2 = new TLine(xlo, 4, xhi, 4);
  l1->SetLineStyle(2); l1->SetLineColor(kRed);    l1->Draw();
  l2->SetLineStyle(2); l2->SetLineColor(kOrange); l2->Draw();

  TLatex lt; lt.SetNDC(); lt.SetTextSize(0.04);
  lt.DrawLatex(0.65, 0.28, "1#sigma (#Delta#chi^{2}=1)");
  lt.DrawLatex(0.65, 0.52, "2#sigma (#Delta#chi^{2}=4)");
  lt.SetTextColor(kGray + 2);
  lt.DrawLatex(0.14, 0.86,
    Form("sin^{2}#theta_{W} = %.4f #pm %.4f", sin2thW_fit, sin2thW_err));

  TLine* lSM = new TLine(Physics::sin2tW_SM, 0, Physics::sin2tW_SM, 10);
  lSM->SetLineColor(kGreen + 2); lSM->SetLineWidth(2); lSM->Draw();
  TLatex lt2; lt2.SetTextSize(0.035); lt2.SetTextColor(kGreen + 2);
  lt2.DrawLatex(Physics::sin2tW_SM + 0.001, 8.0, "SM");

  c2->SaveAs("WeinbergAngleFit_SBND_scan.pdf");
  c2->SaveAs("WeinbergAngleFit_SBND_scan.png");
  std::cout << "[INFO] Scan plots saved: WeinbergAngleFit_SBND_scan.{pdf,png}\n";
}

// ─────────────────────────── MAIN FUNCTION ───────────────────────────────────

void WeinbergAngleFit_SBND(const char* inputFile = "") {

  std::cout << "\n"
            << "╔══════════════════════════════════════════════════════════════╗\n"
            << "║  Weinberg Angle Fit  —  SBND  ν_μ-e elastic scattering       ║\n"
            << "╚══════════════════════════════════════════════════════════════╝\n\n";

  // ── Load or generate histograms ───────────────────────────────────────────
  TFile* fIn = nullptr;
  bool   usePlaceholder = true;

  if (inputFile && strlen(inputFile) > 0) {
    fIn = TFile::Open(inputFile, "READ");
    if (!fIn || fIn->IsZombie()) {
      std::cerr << "[WARN] Cannot open '" << inputFile
                << "'. Falling back to placeholder histograms.\n";
    } else {
      g_hData  = dynamic_cast<TH1D*>(fIn->Get(Setup::DATA_HIST));
      g_hEff   = dynamic_cast<TH1D*>(fIn->Get(Setup::EFF_HIST));
      g_hSmear = dynamic_cast<TH2D*>(fIn->Get(Setup::SMEAR_HIST));
      g_hFlux  = dynamic_cast<TH1D*>(fIn->Get(Setup::FLUX_HIST));
      g_hBkg   = dynamic_cast<TH1D*>(fIn->Get(Setup::BKG_HIST));

      if (!g_hData || !g_hEff || !g_hSmear || !g_hFlux) {
        std::cerr << "[WARN] Required histograms not found in file. "
                  << "Falling back to placeholder.\n";
      } else {
        g_hData->SetDirectory(nullptr);
        g_hEff->SetDirectory(nullptr);
        g_hSmear->SetDirectory(nullptr);
        g_hFlux->SetDirectory(nullptr);
        if (g_hBkg) g_hBkg->SetDirectory(nullptr);
        usePlaceholder = false;

        // Normalise smearing matrix columns if not already done
        Int_t nReco = g_hSmear->GetNbinsX();
        Int_t nTrue = g_hSmear->GetNbinsY();
        for (Int_t iT = 1; iT <= nTrue; ++iT) {
          Double_t colSum = 0.0;
          for (Int_t iR = 1; iR <= nReco; ++iR)
            colSum += g_hSmear->GetBinContent(iR, iT);
          if (colSum > 0.0)
            for (Int_t iR = 1; iR <= nReco; ++iR)
              g_hSmear->SetBinContent(iR, iT,
                g_hSmear->GetBinContent(iR, iT) / colSum);
        }
      }
    }
  }

  if (usePlaceholder) {
    // FIX [BUG-3]: original code allocated temporaries for flux/eff/smear,
    // then called GeneratePlaceholderHistograms() which overwrote the global
    // pointers — leaking the temporaries.  Now we call the function directly
    // with no prior allocation; it creates all globals itself.
    GeneratePlaceholderHistograms();
  }

  // ── Sanity checks ─────────────────────────────────────────────────────────
  if (!g_hData || !g_hEff || !g_hSmear || !g_hFlux) {
    std::cerr << "[ERROR] Histograms not properly initialised. Aborting.\n";
    return;
  }

  Int_t nDataBins = g_hData->GetNbinsX();
  Int_t nDF = 0;
  for (Int_t i = 1; i <= nDataBins; ++i) {
    bool inRange = true;
    if (Setup::FIT_BIN_MIN > 0 && i < Setup::FIT_BIN_MIN) inRange = false;
    if (Setup::FIT_BIN_MAX > 0 && i > Setup::FIT_BIN_MAX) inRange = false;
    if (inRange && g_hData->GetBinContent(i) > 0) ++nDF;
  }
  --nDF; // 1 free parameter

  std::cout << "[INFO] Data bins with content: " << nDF + 1 << "\n"
            << "[INFO] Free parameters:         1  (sin²θ_W)\n"
            << "[INFO] Effective ndf:           " << nDF << "\n\n";

  // ── Minimisation with TMinuit ─────────────────────────────────────────────
  TMinuit minuit(1);
  minuit.SetFCN(FCN);
  minuit.SetPrintLevel(1);
  minuit.SetErrorDef(1.0); // Δχ² = 1 → 1σ

  Double_t arglist[10];
  Int_t    ierflg = 0;

  minuit.mnparm(0, "sin2thW",
                Physics::sin2tW_SM, // start value
                0.005,              // step
                0.10,               // lower limit
                0.40,               // upper limit
                ierflg);

  // Minimise
  arglist[0] = 5000;  // max calls
  arglist[1] = 0.01;  // tolerance
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

  // Improve with HESSE for the parabolic uncertainty
  arglist[0] = 2000;
  minuit.mnexcm("HESSE", arglist, 1, ierflg);

  // FIX [IMPROVEMENT]: run MINOS for asymmetric (profile-likelihood) errors
  arglist[0] = 2000;
  minuit.mnexcm("MINOS", arglist, 1, ierflg);

  // Extract result
  Double_t sin2thW_fit, sin2thW_err;
  minuit.GetParameter(0, sin2thW_fit, sin2thW_err);

  // Get chi2 at best-fit
  Double_t chi2_min = 0, grad = 0;
  Double_t par[1] = {sin2thW_fit};
  Int_t    npar   = 1;
  FCN(npar, &grad, chi2_min, par, 4);

  // Also retrieve MINOS asymmetric errors
  Double_t eplus = 0, eminus = 0, eparab = 0, gcc = 0;
  minuit.mnerrs(0, eplus, eminus, eparab, gcc);

  // ── Print results ─────────────────────────────────────────────────────────
  std::cout << "\n"
            << "═══════════════════════════════════════════════════════════════\n"
            << "  RESULTS\n"
            << "───────────────────────────────────────────────────────────────\n"
            << Form("  sin²θ_W (fit)  = %.5f ± %.5f  (HESSE)\n",
                    sin2thW_fit, sin2thW_err)
            << Form("               = %.5f + %.5f / %.5f  (MINOS)\n",
                    sin2thW_fit, eplus, eminus)
            << Form("  sin²θ_W (SM)   = %.5f\n",         Physics::sin2tW_SM)
            << Form("  Pull           = %+.2f σ\n",
                    (sin2thW_fit - Physics::sin2tW_SM) / sin2thW_err)
            << Form("  χ²             = %.2f\n",          chi2_min)
            << Form("  ndf            = %d\n",             nDF)
            << Form("  p-value        = %.4f\n",           TMath::Prob(chi2_min, nDF))
            << "═══════════════════════════════════════════════════════════════\n\n";

  // ── Plots ─────────────────────────────────────────────────────────────────
  DrawResult(sin2thW_fit, sin2thW_err, chi2_min, nDF);
  DrawDeltaChi2Scan(sin2thW_fit, sin2thW_err);

  // ── Save results to ROOT file ─────────────────────────────────────────────
  TFile* fOut = TFile::Open("WeinbergAngleFit_SBND_results.root", "RECREATE");
  if (fOut && !fOut->IsZombie()) {
    TVectorD results(6);
    results[0] = sin2thW_fit;
    results[1] = sin2thW_err;
    results[2] = chi2_min;
    results[3] = (Double_t)nDF;
    results[4] = eplus;   // FIX [IMPROVEMENT]: save MINOS errors too
    results[5] = eminus;
    results.Write("fit_results"); // [sin2thW, err_hesse, chi2, ndf, eplus, eminus]

    TH1D* hBestFit = BuildRecoPrediction(sin2thW_fit);
    hBestFit->SetName("h_best_fit_reco");
    hBestFit->Write();
    delete hBestFit;  // FIX [BUG-6]: was never deleted in original

    g_hData->Write("h_data_out");
    fOut->Close();
    delete fOut;
    std::cout << "[INFO] Results saved to WeinbergAngleFit_SBND_results.root\n\n";
  }

  if (fIn) { fIn->Close(); delete fIn; }
}
