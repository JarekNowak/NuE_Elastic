/*=============================================================================
  NuMMFit_SBND.C
  ─────────────────────────────────────────────────────────────────────────────
  Fit for the neutrino magnetic moment (NMM) using ν_μ + e⁻ → ν_μ + e⁻
  elastic scattering at SBND.

  Physics:
    The electromagnetic contribution from a non-zero NMM adds to the SM
    neutral-current cross-section:

      dσ/dT = dσ_SM/dT  +  dσ_NMM/dT

    where (Vogel & Engel 1989):

      dσ_NMM/dT = (πα²ℏ²c²/m_e²) × (μ_ν/μ_B)² × (1/T − 1/E_ν)

    The NMM term diverges as 1/T → signal concentrated at low recoil energy.

  Fit strategy:
    • sin²θ_W fixed at the SM value (0.2312).
    • Single free parameter: μ_ν in units of 10⁻¹⁰ μ_B.
    • Prediction = P_SM + μ_par² × P_NMM_shape + P_bkg.
    • TMinuit minimises Neyman χ² over μ_par (physical range μ_par ≥ 0).
    • 90% CL upper limit from Δχ² profile scan (Δχ² = 2.71 threshold).

  Input file must contain (produced by NuMMSelection.C):
    h_data   TH1D  reco T_e spectrum
    h_eff    TH1D  selection efficiency vs true T_e
    h_smear  TH2D  smearing matrix M[reco][true], column-normalised
    h_flux   TH1D  BNB ν_μ flux dΦ/dE_ν [/cm²/POT/GeV]
    h_bkg    TH1D  background prediction (optional)

  Usage:
    root -l -b -q 'magnetic_moment/NuMMFit_SBND.C("inputs.root")'
    root -l -b -q magnetic_moment/NuMMFit_SBND.C   # placeholder mode

  Output:
    NuMMFit_SBND_results.root   — fit result vector + histograms
    NuMMFit_SBND.{pdf,png}      — data vs best-fit overlay
    NuMMFit_SBND_scan.{pdf,png} — Δχ² profile with 90% CL mark

  Author : Jarek Nowak
  Date   : 2026-05-20
=============================================================================*/

#include <iostream>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TRandom3.h"

// ─────────────────────────── CONFIGURATION ───────────────────────────────────

namespace Setup {
  // Normalisation
  const Double_t N_ELECTRONS = 2.437e31;  // fiducial e⁻ targets in SBND
  const Double_t POT         = 1.0e21;    // target exposure [protons on target]
  const Int_t    XSEC_NSUB  = 8;         // integration substeps per true-T bin

  // Fit parameter (μ_ν in units of 10⁻¹⁰ μ_B)
  const Double_t MU_INIT  = 1.0;   // starting value
  const Double_t MU_STEP  = 0.1;
  const Double_t MU_MIN   = 0.0;   // physical lower bound
  const Double_t MU_MAX   = 20.0;  // upper bound of parameter range
  const Double_t MU_SCAN_MAX = 15.0;  // upper end of Δχ² scan
  const Int_t    MU_SCAN_STEPS = 300; // scan resolution

  // Fit range in reco bins (0 = use all)
  const Int_t FIT_BIN_MIN = 0;
  const Int_t FIT_BIN_MAX = 0;

  // Input histogram names
  const char* DATA_HIST  = "h_data";
  const char* EFF_HIST   = "h_eff";
  const char* SMEAR_HIST = "h_smear";
  const char* FLUX_HIST  = "h_flux";
  const char* BKG_HIST   = "h_bkg";
}

// ─────────────────────────── PHYSICS CONSTANTS ───────────────────────────────

namespace Physics {
  const Double_t GF       = 1.16638e-5;  // GeV⁻²
  const Double_t me       = 0.511e-3;    // GeV
  const Double_t hbarc2   = 3.8938e-28;  // cm²·GeV²
  const Double_t sin2tW   = 0.2312;      // SM value (fixed)
  const Double_t alpha    = 1.0/137.036; // fine-structure constant

  // SM coupling constants
  const Double_t gL = -0.5 + sin2tW;
  const Double_t gR =        sin2tW;

  // NMM cross-section prefactor: πα²ℏ²c²/m_e²
  // With μ_ν in units of 10⁻¹⁰ μ_B:
  //   dσ_NMM/dT = NMM_SIGMA0 × μ_par² × (1/T − 1/E_ν)   [cm²/GeV]
  const Double_t NMM_SIGMA0 = TMath::Pi() * alpha * alpha * hbarc2
                               / (me * me)          // cm² (base)
                               * 1e-20;             // × (10⁻¹⁰)² unit conversion
}

// ──────────────────────── GLOBAL HISTOGRAM POINTERS ──────────────────────────
// Required by TMinuit FCN (must be global).

TH1D* g_hData       = nullptr;
TH1D* g_hEff        = nullptr;
TH1D* g_hBkg        = nullptr;
TH2D* g_hSmear      = nullptr;
TH1D* g_hFlux       = nullptr;
TH1D* g_hPredSM     = nullptr;  // precomputed SM signal reco spectrum
TH1D* g_hNMMshape   = nullptr;  // precomputed NMM shape per μ_par² unit

// ─────────────────── CROSS-SECTION FUNCTIONS ─────────────────────────────────

/// Kinematic maximum recoil for neutrino energy Enu [GeV].
inline Double_t T_max(Double_t Enu) {
  return 2.0 * Enu * Enu / (Physics::me + 2.0 * Enu);
}

/// SM differential cross-section dσ/dT [cm²/GeV] at given T and Enu [GeV].
Double_t dsigma_SM(Double_t T, Double_t Enu) {
  if (T <= 0 || T >= T_max(Enu)) return 0.0;
  const Double_t& GF  = Physics::GF;
  const Double_t& me  = Physics::me;
  const Double_t& gL  = Physics::gL;
  const Double_t& gR  = Physics::gR;
  const Double_t& h2  = Physics::hbarc2;

  Double_t ratio  = T / Enu;
  Double_t termL  = gL * gL;
  Double_t termR  = gR * gR * (1.0 - ratio) * (1.0 - ratio);
  Double_t termX  = gL * gR * me * T / (Enu * Enu);
  Double_t xsec   = GF * GF * me / (2.0 * TMath::Pi())
                    * (termL + termR - termX);
  return xsec * h2;  // [cm²/GeV]
}

/// NMM differential cross-section dσ_NMM/dT per unit μ_par² [cm²/GeV].
/// μ_par is in units of 10⁻¹⁰ μ_B; result is the shape factor (μ_par² ×
/// return_value gives the physical cross-section contribution).
Double_t dsigma_NMM(Double_t T, Double_t Enu) {
  if (T <= 0 || T >= T_max(Enu)) return 0.0;
  return Physics::NMM_SIGMA0 * (1.0/T - 1.0/Enu);  // [cm²/GeV] per μ_par²
}

// ─────────────────── SPECTRUM BUILDER FUNCTIONS ──────────────────────────────

/// Build true-T spectrum from a given differential cross-section kernel.
/// kernel(T, Enu) returns dσ/dT in cm²/GeV.
TH1D* BuildTrueSpectrum(const char* name,
                         Double_t (*kernel)(Double_t, Double_t)) {
  const Int_t nT   = g_hEff->GetNbinsX();
  const Int_t nEnu = g_hFlux->GetNbinsX();

  TH1D* hTrue = (TH1D*)g_hEff->Clone(name);
  hTrue->Reset();
  hTrue->SetDirectory(nullptr);

  for (Int_t iT = 1; iT <= nT; ++iT) {
    Double_t T_lo = hTrue->GetBinLowEdge(iT);
    Double_t T_hi = T_lo + hTrue->GetBinWidth(iT);

    Double_t rate = 0.0;
    for (Int_t iE = 1; iE <= nEnu; ++iE) {
      Double_t Enu  = g_hFlux->GetBinCenter(iE);
      Double_t dEnu = g_hFlux->GetBinWidth(iE);
      Double_t flux = g_hFlux->GetBinContent(iE);
      if (flux <= 0) continue;
      if ((0.5*(T_lo+T_hi)) > T_max(Enu)) continue;

      Double_t integral = 0.0;
      for (Int_t k = 0; k < Setup::XSEC_NSUB; ++k) {
        Double_t T1 = T_lo + (T_hi-T_lo) * k         / Setup::XSEC_NSUB;
        Double_t T2 = T_lo + (T_hi-T_lo) * (k + 1.0) / Setup::XSEC_NSUB;
        integral += 0.5 * (kernel(T1, Enu) + kernel(T2, Enu)) * (T2 - T1);
      }
      rate += flux * dEnu * integral;
    }

    Double_t nEvents = rate * Setup::N_ELECTRONS * Setup::POT;
    Double_t eff     = g_hEff->GetBinContent(iT);
    hTrue->SetBinContent(iT, nEvents * eff);
  }
  return hTrue;
}

/// Apply smearing matrix to map true-T spectrum → reco-T spectrum.
/// Creates histogram from smear-matrix X-axis (does not require g_hData).
TH1D* SmearToReco(TH1D* hTrue, const char* name = "h_reco_tmp") {
  const Int_t nReco = g_hSmear->GetNbinsX();
  const Int_t nTrue = g_hSmear->GetNbinsY();
  TH1D* hReco = new TH1D(name, ";T_{e}^{reco} [GeV];Events/bin",
                          nReco,
                          g_hSmear->GetXaxis()->GetXmin(),
                          g_hSmear->GetXaxis()->GetXmax());
  hReco->SetDirectory(nullptr);
  for (Int_t iR = 1; iR <= nReco; ++iR) {
    Double_t val = 0.0;
    for (Int_t iT = 1; iT <= nTrue; ++iT)
      val += g_hSmear->GetBinContent(iR, iT) * hTrue->GetBinContent(iT);
    hReco->SetBinContent(iR, val);
  }
  return hReco;
}

/// Precompute SM and NMM shape spectra (call once after loading histograms).
void PrecomputeSpectra() {
  TH1D* hSMtrue  = BuildTrueSpectrum("h_sm_true",  dsigma_SM);
  TH1D* hNMMtrue = BuildTrueSpectrum("h_nmm_true", dsigma_NMM);

  if (g_hPredSM)   { delete g_hPredSM;   g_hPredSM   = nullptr; }
  if (g_hNMMshape) { delete g_hNMMshape; g_hNMMshape = nullptr; }

  g_hPredSM   = SmearToReco(hSMtrue,  "h_pred_sm");
  g_hNMMshape = SmearToReco(hNMMtrue, "h_nmm_shape");

  delete hSMtrue; delete hNMMtrue;
}

/// Build total reco prediction: SM + μ_par² × NMM_shape + background.
TH1D* BuildRecoPrediction(Double_t mu_par) {
  TH1D* hPred = (TH1D*)g_hPredSM->Clone("h_pred_total");
  hPred->SetDirectory(nullptr);

  Double_t mu2 = mu_par * mu_par;
  for (Int_t i = 1; i <= hPred->GetNbinsX(); ++i) {
    Double_t val = g_hPredSM->GetBinContent(i)
                   + mu2 * g_hNMMshape->GetBinContent(i);
    if (g_hBkg) val += g_hBkg->GetBinContent(i);
    hPred->SetBinContent(i, val);
  }
  return hPred;
}

// ─────────────────── CHI-SQUARED (Baker–Cousins Poisson LR) ──────────────────

/// Baker–Cousins (Poisson likelihood-ratio) χ² contribution for one bin:
///     2·[ pred − obs + obs·ln(obs/pred) ]   (→ 2·pred as obs → 0).
/// Unlike Neyman χ² (err² = obs) this is valid for low/zero counts and keeps
/// empty bins (obs=0 → 2·pred) — essential here, because the NMM upper limit is
/// driven by the low-T bins where the 1/T excess would appear but the SM
/// expectation is small. Asymptotically χ²-distributed, so the Δχ² = 2.71
/// (one-sided 90% CL) threshold stays valid.
inline Double_t PoissonChi2Bin(Double_t obs, Double_t pred) {
  const Double_t tiny = 1e-9;
  if (pred < tiny) pred = tiny;            // guard ln()/division; pred is ≥ 0
  if (obs <= 0.0)  return 2.0 * pred;      // obs·ln(obs/pred) → 0 as obs → 0
  return 2.0 * (pred - obs + obs * std::log(obs / pred));
}

/// Baker–Cousins χ² of the data vs the SM + μ²·NMM prediction at mu_par.
/// Single implementation shared by FCN, the best-fit χ², the 90% CL upper-limit
/// scan and the Δχ² scan plot (previously four separate copies of this loop).
Double_t Chi2AtMu(Double_t mu_par) {
  TH1D* hPred = BuildRecoPrediction(mu_par);
  Double_t chi2  = 0.0;
  Int_t    nBins = g_hData->GetNbinsX();
  for (Int_t i = 1; i <= nBins; ++i) {
    if (Setup::FIT_BIN_MIN > 0 && i < Setup::FIT_BIN_MIN) continue;
    if (Setup::FIT_BIN_MAX > 0 && i > Setup::FIT_BIN_MAX) continue;
    chi2 += PoissonChi2Bin(g_hData->GetBinContent(i), hPred->GetBinContent(i));
  }
  delete hPred;
  return chi2;
}

void FCN(Int_t& npar, Double_t* grad, Double_t& fval,
         Double_t* par, Int_t iflag) {
  fval = Chi2AtMu(par[0]);
}

// ──────────────────────── PLOTTING HELPERS ───────────────────────────────────

void DrawResult(Double_t mu_bf, Double_t mu_ul90, Double_t chi2, Int_t ndf) {
  TH1D* hBest = BuildRecoPrediction(mu_bf);
  TH1D* hSM   = BuildRecoPrediction(0.0);

  hBest->SetLineColor(kRed+1);  hBest->SetLineWidth(2);
  hBest->SetFillColorAlpha(kRed+1, 0.15); hBest->SetFillStyle(1001);
  hSM->SetLineColor(kBlue+1);   hSM->SetLineStyle(2); hSM->SetLineWidth(2);

  TH1D* hBkgPlot = nullptr;
  if (g_hBkg) {
    hBkgPlot = (TH1D*)g_hBkg->Clone("h_bkg_plot");
    hBkgPlot->SetFillColorAlpha(kGray, 0.5);
    hBkgPlot->SetLineColor(kGray+2);
  }

  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  TCanvas* c1 = new TCanvas("c_nmm_fit","NuMM Fit — SBND", 900, 750);
  c1->Divide(1, 2, 0, 0);

  // Top pad: data + prediction
  TPad* pTop = (TPad*)c1->cd(1);
  pTop->SetPad(0, 0.30, 1, 1.00);
  pTop->SetBottomMargin(0.02); pTop->SetLeftMargin(0.12);

  Double_t ymax = 1.3 * std::max(g_hData->GetMaximum(), hBest->GetMaximum());
  hBest->GetYaxis()->SetRangeUser(0, ymax);
  hBest->GetXaxis()->SetLabelSize(0);
  hBest->GetYaxis()->SetTitle("Events / bin");
  hBest->GetYaxis()->SetTitleSize(0.06);
  hBest->Draw("HIST");
  hSM->Draw("HIST SAME");
  if (hBkgPlot) { hBkgPlot->Draw("HIST SAME"); }
  g_hData->SetMarkerStyle(20); g_hData->SetMarkerSize(0.8);
  g_hData->Draw("E1 SAME");

  TLegend* lg = new TLegend(0.55, 0.62, 0.90, 0.90);
  lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.048);
  lg->AddEntry(g_hData, "Pseudo-data",        "pe");
  lg->AddEntry(hBest, Form("Best fit #mu_{#nu}=%.2f#times10^{-10}#mu_{B}",mu_bf),"lf");
  lg->AddEntry(hSM,   "SM only (#mu_{#nu}=0)","l");
  if (hBkgPlot) lg->AddEntry(hBkgPlot, "Background","f");
  lg->Draw();

  TLatex tl; tl.SetNDC(); tl.SetTextSize(0.050);
  tl.DrawLatex(0.14, 0.93,
    "#nu_{#mu}+e^{#minus}  Magnetic Moment Fit  —  SBND  1#times10^{21} POT");
  tl.SetTextSize(0.043);
  tl.DrawLatex(0.14, 0.20,
    Form("#chi^{2}/ndf = %.1f / %d   90%% CL UL: #mu_{#nu} < %.2f#times10^{-10}#mu_{B}",
         chi2, ndf, mu_ul90));

  // Bottom pad: residuals (data − SM) / error, showing NMM excess
  TPad* pBot = (TPad*)c1->cd(2);
  pBot->SetPad(0, 0.00, 1, 0.30);
  pBot->SetTopMargin(0.02); pBot->SetBottomMargin(0.30); pBot->SetLeftMargin(0.12);

  TH1D* hResid = (TH1D*)g_hData->Clone("h_resid");
  hResid->Reset();
  hResid->SetTitle(";T_{e}^{reco} [GeV];(data#minusSM)/#sigma");
  hResid->GetXaxis()->SetTitleSize(0.12); hResid->GetXaxis()->SetLabelSize(0.10);
  hResid->GetYaxis()->SetTitleSize(0.10); hResid->GetYaxis()->SetLabelSize(0.09);
  hResid->GetYaxis()->SetTitleOffset(0.5);

  for (Int_t i = 1; i <= hResid->GetNbinsX(); ++i) {
    Double_t obs  = g_hData->GetBinContent(i);
    Double_t sm   = hSM->GetBinContent(i);
    Double_t err  = g_hData->GetBinError(i);
    if (err > 0) hResid->SetBinContent(i, (obs - sm) / err);
    hResid->SetBinError(i, 1.0);
  }
  Double_t absMax = 0.0;
  for (Int_t i = 1; i <= hResid->GetNbinsX(); ++i)
    absMax = std::max(absMax, std::abs(hResid->GetBinContent(i)));
  hResid->GetYaxis()->SetRangeUser(-absMax*1.4, absMax*1.4);
  hResid->SetMarkerStyle(20); hResid->SetMarkerSize(0.7);
  hResid->Draw("E1");
  TLine* zero = new TLine(hResid->GetXaxis()->GetXmin(), 0,
                           hResid->GetXaxis()->GetXmax(), 0);
  zero->SetLineStyle(2); zero->Draw();

  c1->SaveAs("magnetic_moment/NuMMFit_SBND.pdf");
  c1->SaveAs("magnetic_moment/NuMMFit_SBND.png");
  printf("[INFO] Fit plot saved: NuMMFit_SBND.{pdf,png}\n");

  delete hBest; delete hSM;
  if (hBkgPlot) delete hBkgPlot;
  delete hResid;
}

void DrawChiScan(Double_t mu_bf, Double_t chi2_min, Double_t mu_ul90) {
  const Int_t    N    = Setup::MU_SCAN_STEPS;
  const Double_t muhi = Setup::MU_SCAN_MAX;

  TGraph* gScan = new TGraph(N);
  for (Int_t i = 0; i < N; ++i) {
    Double_t mu = muhi * i / (N - 1.0);
    gScan->SetPoint(i, mu, Chi2AtMu(mu) - chi2_min);
  }

  gStyle->SetOptStat(0);
  TCanvas* c2 = new TCanvas("c_nmm_scan","NuMM Δχ² Scan — SBND", 800, 550);
  c2->SetLeftMargin(0.12); c2->SetBottomMargin(0.13);

  gScan->SetTitle(";#mu_{#nu} [10^{-10} #mu_{B}];#Delta#chi^{2}");
  gScan->SetLineColor(kBlue+1); gScan->SetLineWidth(2);
  gScan->GetYaxis()->SetRangeUser(0, 12);
  gScan->Draw("AL");

  // 90% CL line at Δχ² = 2.71
  Double_t xmin = gScan->GetXaxis()->GetXmin();
  Double_t xmax = gScan->GetXaxis()->GetXmax();
  TLine* lCL = new TLine(xmin, 2.71, xmax, 2.71);
  lCL->SetLineColor(kRed+1); lCL->SetLineStyle(2); lCL->SetLineWidth(2);
  lCL->Draw();

  // Vertical mark at 90% CL upper limit
  TLine* lUL = new TLine(mu_ul90, 0, mu_ul90, 2.71);
  lUL->SetLineColor(kRed+1); lUL->SetLineStyle(3); lUL->SetLineWidth(2);
  lUL->Draw();

  TLatex tl; tl.SetNDC(); tl.SetTextSize(0.046);
  tl.DrawLatex(0.14, 0.93,
    "SBND  #nu_{#mu}+e^{#minus}  NMM Sensitivity  —  1#times10^{21} POT");
  tl.SetTextSize(0.040);
  tl.DrawLatex(0.55, 0.78,
    Form("90%% CL UL: #mu_{#nu} < %.2f#times10^{-10} #mu_{B}", mu_ul90));

  c2->SaveAs("magnetic_moment/NuMMFit_SBND_scan.pdf");
  c2->SaveAs("magnetic_moment/NuMMFit_SBND_scan.png");
  printf("[INFO] Scan plot saved: NuMMFit_SBND_scan.{pdf,png}\n");

  delete gScan;
}

// ──────────────────── PLACEHOLDER HISTOGRAM GENERATOR ────────────────────────

void GeneratePlaceholders() {
  std::cout << "[INFO] Generating placeholder histograms (demo mode).\n"
            << "       Run NuMMSelection.C on real data for an actual result.\n";

  const Int_t nT = 60; const Double_t Tlo = 0.0, Thi = 0.6;
  const Int_t nE = 60; const Double_t Elo = 0.0, Ehi = 3.0;

  // BNB ν_μ flux at SBND (~110 m): peak ≈ 3×10⁻⁸ /cm²/POT/GeV
  g_hFlux = new TH1D("h_flux","Flux;E_{#nu} [GeV];#Phi",nE,Elo,Ehi);
  for (Int_t i = 1; i <= nE; ++i) {
    Double_t E = g_hFlux->GetBinCenter(i);
    g_hFlux->SetBinContent(i, 3.0e-8 * TMath::Gaus(E, 0.7, 0.25));
  }

  g_hEff = new TH1D("h_eff","Efficiency;T_{e} [GeV];#epsilon",nT,Tlo,Thi);
  for (Int_t i = 1; i <= nT; ++i) {
    Double_t T = g_hEff->GetBinCenter(i);
    // Efficiency rises from threshold (10 MeV) and falls at high T
    Double_t eff = (T < 0.010) ? 0.0 :
                   0.40 * (1.0 - TMath::Exp(-(T-0.010)/0.020)) * TMath::Exp(-0.3*T);
    g_hEff->SetBinContent(i, eff);
  }

  g_hSmear = new TH2D("h_smear","Smearing;T_{e}^{reco} [GeV];T_{e}^{true} [GeV]",
                       nT,Tlo,Thi, nT,Tlo,Thi);
  for (Int_t iT = 1; iT <= nT; ++iT) {
    Double_t Ttrue = g_hSmear->GetYaxis()->GetBinCenter(iT);
    Double_t sig   = 0.08*Ttrue + 0.005;
    Double_t colSum = 0.0;
    for (Int_t iR = 1; iR <= nT; ++iR) {
      Double_t Treco = g_hSmear->GetXaxis()->GetBinCenter(iR);
      Double_t val   = TMath::Gaus(Treco, Ttrue, sig);
      g_hSmear->SetBinContent(iR, iT, val);
      colSum += val;
    }
    if (colSum > 0)
      for (Int_t iR = 1; iR <= nT; ++iR)
        g_hSmear->SetBinContent(iR, iT, g_hSmear->GetBinContent(iR,iT)/colSum);
  }

  // Generate SM pseudo-data (μ_ν = 0)
  PrecomputeSpectra();

  g_hBkg = new TH1D("h_bkg","Background;T_{e}^{reco} [GeV];Events",nT,Tlo,Thi);
  for (Int_t i = 1; i <= nT; ++i)
    g_hBkg->SetBinContent(i, 0.05 * g_hPredSM->GetBinContent(i));

  g_hData = new TH1D("h_data","Pseudo-data;T_{e}^{reco} [GeV];Events/bin",nT,Tlo,Thi);
  TRandom3 rng(42);
  for (Int_t i = 1; i <= nT; ++i) {
    Double_t mu = g_hPredSM->GetBinContent(i) + g_hBkg->GetBinContent(i);
    Double_t obs = (mu > 0) ? rng.Poisson(mu) : 0.0;
    g_hData->SetBinContent(i, obs);
    g_hData->SetBinError(i, (obs > 0) ? TMath::Sqrt(obs) : 1.0);
  }
}

// ───────────────────────────── MAIN FUNCTION ─────────────────────────────────

void NuMMFit_SBND(const char* inputFile = "",
                  const char* outputFile = "magnetic_moment/NuMMFit_SBND_results.root")
{
  std::cout << "\n"
            << "╔══════════════════════════════════════════════════════════════╗\n"
            << "║  Neutrino Magnetic Moment Fit  —  SBND  ν_μ-e elastic        ║\n"
            << "╚══════════════════════════════════════════════════════════════╝\n\n";

  bool usePlaceholders = (inputFile == nullptr || strlen(inputFile) == 0);

  if (!usePlaceholders) {
    TFile* fIn = TFile::Open(inputFile, "READ");
    if (!fIn || fIn->IsZombie()) {
      std::cerr << "[ERROR] Cannot open: " << inputFile << "\n";
      return;
    }
    g_hData  = dynamic_cast<TH1D*>(fIn->Get(Setup::DATA_HIST));
    g_hEff   = dynamic_cast<TH1D*>(fIn->Get(Setup::EFF_HIST));
    g_hSmear = dynamic_cast<TH2D*>(fIn->Get(Setup::SMEAR_HIST));
    g_hFlux  = dynamic_cast<TH1D*>(fIn->Get(Setup::FLUX_HIST));
    g_hBkg   = dynamic_cast<TH1D*>(fIn->Get(Setup::BKG_HIST));

    if (!g_hData || !g_hEff || !g_hSmear || !g_hFlux) {
      std::cerr << "[ERROR] Required histogram(s) missing from " << inputFile << "\n";
      fIn->Close(); return;
    }
    // Detach from file so they survive fIn close
    if (g_hData)  g_hData ->SetDirectory(nullptr);
    if (g_hEff)   g_hEff  ->SetDirectory(nullptr);
    if (g_hBkg)   g_hBkg  ->SetDirectory(nullptr);
    if (g_hSmear) g_hSmear->SetDirectory(nullptr);
    if (g_hFlux)  g_hFlux ->SetDirectory(nullptr);
    fIn->Close();

    printf("[INFO] Loaded histograms from %s\n", inputFile);
    PrecomputeSpectra();
  } else {
    GeneratePlaceholders();  // also calls PrecomputeSpectra internally
  }

  // Effective degrees of freedom: the Baker–Cousins sum runs over every in-range
  // bin, so count all of them (not just non-empty ones), minus the free parameter.
  Int_t nDF = 0;
  Int_t nDataBins = g_hData->GetNbinsX();
  for (Int_t i = 1; i <= nDataBins; ++i) {
    if (Setup::FIT_BIN_MIN > 0 && i < Setup::FIT_BIN_MIN) continue;
    if (Setup::FIT_BIN_MAX > 0 && i > Setup::FIT_BIN_MAX) continue;
    ++nDF;
  }
  --nDF;  // 1 free parameter (μ_ν)
  printf("[INFO] Fitted bins (in range): %d\n", nDF + 1);
  printf("[INFO] Effective ndf: %d\n\n", nDF);

  // ── TMinuit fit ───────────────────────────────────────────────────────────
  TMinuit minuit(1);
  minuit.SetFCN(FCN);
  minuit.SetPrintLevel(1);

  Int_t ierflg = 0;
  Double_t arglist[10];

  arglist[0] = 1.0;
  minuit.mnexcm("SET ERRDEF", arglist, 1, ierflg);  // Δχ² = 1 for 1σ

  minuit.mnparm(0, "mu_nu", Setup::MU_INIT, Setup::MU_STEP,
                Setup::MU_MIN, Setup::MU_MAX, ierflg);

  arglist[0] = 5000; arglist[1] = 0.01;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
  minuit.mnexcm("HESSE",  arglist, 1, ierflg);
  minuit.mnexcm("MINOS",  arglist, 1, ierflg);

  Double_t mu_bf, mu_err_hesse;
  Double_t mu_eplus, mu_eminus, mu_gcc;
  Int_t    mu_ivar;
  minuit.GetParameter(0, mu_bf, mu_err_hesse);
  minuit.mnerrs(0, mu_eplus, mu_eminus, mu_err_hesse, mu_gcc);

  // Physical best-fit: enforce non-negative
  Double_t mu_phys = std::max(0.0, mu_bf);

  // χ² at best-fit (Baker–Cousins)
  Double_t chi2_min = Chi2AtMu(mu_phys);

  // ── 90% CL upper limit via Δχ² scan ──────────────────────────────────────
  // Scan mu from mu_phys upward until Δχ² ≥ 2.71 (one-sided 90% CL).
  Double_t mu_ul90 = Setup::MU_SCAN_MAX;
  {
    const Int_t    Nscan = Setup::MU_SCAN_STEPS;
    const Double_t muhi  = Setup::MU_SCAN_MAX;
    for (Int_t i = 0; i < Nscan; ++i) {
      Double_t mu = mu_phys + (muhi - mu_phys) * i / (Nscan - 1.0);
      if (Chi2AtMu(mu) - chi2_min >= 2.71) { mu_ul90 = mu; break; }
    }
  }

  // ── Print results ─────────────────────────────────────────────────────────
  printf("\n");
  printf("═══════════════════════════════════════════════════════════════\n");
  printf("  NEUTRINO MAGNETIC MOMENT FIT RESULTS\n");
  printf("───────────────────────────────────────────────────────────────\n");
  printf("  μ_ν best-fit  = %.3f × 10⁻¹⁰ μ_B  ± %.3f (HESSE)\n",
         mu_phys, mu_err_hesse);
  printf("  90%% CL UL     < %.3f × 10⁻¹⁰ μ_B  (Δχ² = 2.71)\n", mu_ul90);
  printf("  χ²  at best   = %.2f\n", chi2_min);
  printf("  ndf           = %d\n", nDF);
  if (nDF > 0)
    printf("  p-value       = %.4f\n", TMath::Prob(chi2_min, nDF));
  printf("═══════════════════════════════════════════════════════════════\n\n");

  // ── Save results ──────────────────────────────────────────────────────────
  TFile* fOut = TFile::Open(outputFile, "RECREATE");
  if (fOut && !fOut->IsZombie()) {
    TVectorD res(5);
    res[0] = mu_phys;
    res[1] = mu_err_hesse;
    res[2] = mu_ul90;
    res[3] = chi2_min;
    res[4] = (Double_t)nDF;
    res.Write("fit_results");

    TH1D* hBest = BuildRecoPrediction(mu_phys);
    TH1D* hSM   = BuildRecoPrediction(0.0);
    hBest->SetName("h_best_fit");  hBest->Write();
    hSM->SetName("h_sm_pred");     hSM->Write();
    g_hNMMshape->SetName("h_nmm_shape"); g_hNMMshape->Write();
    g_hData->Write("h_data_out");
    delete hBest; delete hSM;
    fOut->Close();
    printf("[INFO] Results saved to %s\n", outputFile);
  }

  // ── Plots ─────────────────────────────────────────────────────────────────
  DrawResult(mu_phys, mu_ul90, chi2_min, nDF);
  DrawChiScan(mu_phys, chi2_min, mu_ul90);
}
