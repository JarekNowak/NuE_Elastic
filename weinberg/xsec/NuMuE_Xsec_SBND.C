/*=============================================================================
  NuMuE_Xsec_SBND.C
  ─────────────────────────────────────────────────────────────────────────────
  Plots the differential and total cross-section for

        ν_μ  +  e⁻  →  ν_μ  +  e⁻   (elastic, NC)

  as a function of the recoil electron kinetic energy T, for neutrino
  energies representative of the BNB beam at SBND.

  Physics:
    dσ/dT = (G_F² mₑ / 2π) [ gL² + gR²(1 − T/Eν)² − gL·gR·mₑT/Eν² ]
    gL = −½ + sin²θ_W
    gR =      sin²θ_W
    T ∈ [0,  Eν / (1 + mₑ/2Eν)]   (kinematic limit)

  Outputs:
    NuMuE_Xsec_SBND_diff.pdf/png   – dσ/dT vs T  for several Eν + flux avg
    NuMuE_Xsec_SBND_total.pdf/png  – total σ(Eν)
    NuMuE_Xsec_SBND_ratio.pdf/png  – (σ_SM ± Δsin²θ_W) / σ_SM  ratio band
    NuMuE_Xsec_SBND_results.root   – all histograms / TF1s saved

  Usage:
    root -l -b -q NuMuE_Xsec_SBND.C

  Author : <your name>
  Date   : 2026-03-28
=============================================================================*/

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TAxis.h"
#include "TFile.h"
#include "TPad.h"
#include "TColor.h"
#include "TArrow.h"
#include "TBox.h"
#include "TFrame.h"

// ─────────────────────────── CONSTANTS ───────────────────────────────────────
namespace Phys {
  const Double_t GF     = 1.16638e-5;    // Fermi constant     [GeV⁻²]
  const Double_t me     = 0.511e-3;      // electron mass      [GeV]
  const Double_t hbarc2 = 3.8938e-28;   // ħc²               [cm²·GeV²]
  const Double_t sin2tW = 0.2312;        // SM weak mixing angle
}

// BNB energy range relevant to SBND
namespace BNB {
  const Double_t Enu_min  = 0.10;   // [GeV]
  const Double_t Enu_max  = 3.00;   // [GeV]
  const Double_t Enu_peak = 0.70;   // approximate flux peak [GeV]
}

// ─────────────────────────── PHYSICS FUNCTIONS ───────────────────────────────

/// Maximum recoil kinetic energy (kinematic limit) [GeV]
inline Double_t T_max(Double_t Enu) {
  return Enu / (1.0 + Phys::me / (2.0 * Enu));
}

/// dσ/dT [cm²/GeV]
///   @param T       recoil electron kinetic energy [GeV]
///   @param Enu     incident neutrino energy        [GeV]
///   @param sin2thW sin²θ_W
Double_t dsigma_dT(Double_t T, Double_t Enu, Double_t sin2thW) {
  if (T < 0.0 || Enu <= 0.0) return 0.0;
  if (T > T_max(Enu))        return 0.0;

  Double_t gL = -0.5 + sin2thW;
  Double_t gR =        sin2thW;
  Double_t y  = T / Enu;

  Double_t xs = (Phys::GF * Phys::GF * Phys::me / (2.0 * TMath::Pi())) * Phys::hbarc2
                * (gL*gL + gR*gR*(1-y)*(1-y) - gL*gR*Phys::me*T/(Enu*Enu));
  return std::max(xs, 0.0);
}

/// Total cross-section σ(Eν) by numerical integration [cm²]
Double_t sigma_total(Double_t Enu, Double_t sin2thW, Int_t nSteps = 500) {
  Double_t Tm  = T_max(Enu);
  Double_t dT  = Tm / nSteps;
  Double_t sum = 0.0;
  for (Int_t i = 0; i < nSteps; ++i) {
    Double_t T1 = i       * dT;
    Double_t T2 = (i + 1) * dT;
    sum += 0.5 * (dsigma_dT(T1, Enu, sin2thW) + dsigma_dT(T2, Enu, sin2thW)) * dT;
  }
  return sum;
}

/// Approximate BNB ν_μ flux shape (for flux-averaged plots) [arb. units/GeV]
Double_t bnb_flux(Double_t Enu) {
  if (Enu < 0) return 0.0;
  // Gaussian primary peak + high-energy tail
  Double_t phi = TMath::Gaus(Enu, 0.70, 0.22)
               + 0.20 * TMath::Gaus(Enu, 1.30, 0.35)
               + 0.05 * TMath::Gaus(Enu, 2.00, 0.50);
  return std::max(phi, 0.0);
}

/// Flux-averaged dσ/dT [cm²/GeV]  〈dσ/dT〉_Φ = ∫ Φ(Eν) dσ/dT dEν / ∫ Φ dEν
Double_t flux_avg_dsigma_dT(Double_t T, Double_t sin2thW,
                             Int_t nE = 200) {
  Double_t dE    = (BNB::Enu_max - BNB::Enu_min) / nE;
  Double_t num   = 0.0;
  Double_t denom = 0.0;
  for (Int_t k = 0; k < nE; ++k) {
    Double_t Enu = BNB::Enu_min + (k + 0.5) * dE;
    Double_t phi = bnb_flux(Enu) * dE;
    num   += phi * dsigma_dT(T, Enu, sin2thW);
    denom += phi;
  }
  return (denom > 0) ? num / denom : 0.0;
}

// ─────────────────────────── STYLE SETUP ─────────────────────────────────────
void SetSBNDStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFillStyle(0);
  gStyle->SetLegendFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.052, "XYZ");
  gStyle->SetLabelSize(0.046, "XYZ");
  gStyle->SetTitleOffset(1.05, "X");
  gStyle->SetTitleOffset(1.10, "Y");
}

// ─────────────── PLOT 1: DIFFERENTIAL CROSS SECTION vs T ─────────────────────

void Plot_Differential(TFile* fOut) {

  // Neutrino energies to show (GeV) — span the BNB spectrum
  std::vector<Double_t> ENU = {0.20, 0.40, 0.70, 1.00, 1.50, 2.50};

  // Colours (palette-friendly, colour-blind safe)
  std::vector<Int_t> colours = {
    TColor::GetColor("#E8673A"),   // 0.20 GeV  warm red
    TColor::GetColor("#E8A700"),   // 0.40 GeV  amber
    TColor::GetColor("#3DAA6E"),   // 0.70 GeV  green  (BNB peak)
    TColor::GetColor("#185FA5"),   // 1.00 GeV  blue
    TColor::GetColor("#8E44AD"),   // 1.50 GeV  purple
    TColor::GetColor("#2C3E50"),   // 2.50 GeV  dark slate
  };
  std::vector<Int_t> styles = {1, 1, 1, 1, 1, 1};

  const Int_t nT   = 800;
  const Double_t Tlo_GeV = 1e-4;
  const Double_t Thi_GeV = 2.80;   // up to BNB max ~3 GeV

  // Build TGraphs for each Eν
  std::vector<TGraph*> graphs;
  for (Size_t ie = 0; ie < ENU.size(); ++ie) {
    Double_t Enu = ENU[ie];
    Double_t Tm  = T_max(Enu);

    TGraph* g = new TGraph();
    g->SetLineColor(colours[ie]);
    g->SetLineWidth(2);
    g->SetLineStyle(styles[ie]);

    Int_t np = 0;
    for (Int_t it = 0; it < nT; ++it) {
      // Logarithmic spacing in T (better for log-x plot)
      Double_t logT = TMath::Log10(Tlo_GeV) +
                      (TMath::Log10(Thi_GeV) - TMath::Log10(Tlo_GeV)) * it / (nT - 1);
      Double_t T = TMath::Power(10.0, logT);
      if (T > Tm) break;
      Double_t xs = dsigma_dT(T, Enu, Phys::sin2tW);
      if (xs > 0) { g->SetPoint(np++, T * 1000.0, xs); }  // T in MeV
    }
    graphs.push_back(g);
  }

  // Flux-averaged curve (dashed black)
  TGraph* gFlux = new TGraph();
  gFlux->SetLineColor(kBlack);
  gFlux->SetLineWidth(3);
  gFlux->SetLineStyle(2);   // dashed
  {
    Int_t np = 0;
    for (Int_t it = 0; it < nT; ++it) {
      Double_t logT = TMath::Log10(Tlo_GeV) +
                      (TMath::Log10(Thi_GeV) - TMath::Log10(Tlo_GeV)) * it / (nT-1);
      Double_t T  = TMath::Power(10.0, logT);
      Double_t xs = flux_avg_dsigma_dT(T, Phys::sin2tW);
      if (xs > 0) gFlux->SetPoint(np++, T * 1000.0, xs); // MeV
    }
  }

  // ── Draw ──
  TCanvas* c1 = new TCanvas("c_diff", "dSigma/dT", 900, 680);

  // Invisible frame to set axes
  TH1D* hFrame = new TH1D("hFrame_diff",
                           ";T_{e}^{kin}  [MeV];d#sigma/dT  [cm^{2}/GeV]",
                           100, 0.1, 2500.0);
  hFrame->GetXaxis()->SetMoreLogLabels();
  hFrame->GetXaxis()->SetNoExponent();
  hFrame->GetYaxis()->SetMaxDigits(2);
  hFrame->GetYaxis()->SetRangeUser(5e-47, 1.5e-42);
  hFrame->Draw("AXIS");
  c1->SetLogx();
  c1->SetLogy();

  // Draw individual Eν curves
  for (auto g : graphs) g->Draw("L SAME");

  // Draw flux-average
  gFlux->Draw("L SAME");

  // Legend
  TLegend* leg = new TLegend(0.55, 0.40, 0.94, 0.91);
  leg->SetTextSize(0.040);
  leg->SetHeader("E_{#nu}  [GeV]  (BNB #nu_{#mu})", "C");
  for (Size_t ie = 0; ie < ENU.size(); ++ie) {
    leg->AddEntry(graphs[ie], Form("  %.2f GeV", ENU[ie]), "l");
  }
  leg->AddEntry(gFlux, "  BNB flux average", "l");
  leg->Draw();

  // SBND label
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(62);
  ltx.SetTextSize(0.055);
  ltx.DrawLatex(0.16, 0.88, "SBND");
  ltx.SetTextFont(42);
  ltx.SetTextSize(0.043);
  ltx.DrawLatex(0.16, 0.82, "#nu_{#mu} + e^{-} #rightarrow #nu_{#mu} + e^{-}");
  ltx.DrawLatex(0.16, 0.76, Form("sin^{2}#theta_{W} = %.4f (SM)", Phys::sin2tW));

  // Kinematic edge arrows for BNB peak
  Double_t Tm_peak_MeV = T_max(BNB::Enu_peak) * 1000.0;
  TArrow* arr = new TArrow(Tm_peak_MeV, 3e-46, Tm_peak_MeV, 8e-45, 0.018, ">");
  arr->SetLineColor(TColor::GetColor("#3DAA6E"));
  arr->SetLineWidth(2);
  arr->Draw();
  TLatex ltx2;
  ltx2.SetTextFont(42);
  ltx2.SetTextSize(0.036);
  ltx2.SetTextColor(TColor::GetColor("#3DAA6E"));
  ltx2.DrawLatex(Tm_peak_MeV * 0.42, 1.5e-45, "T_{max}(0.70 GeV)");

  c1->RedrawAxis();
  c1->SaveAs("NuMuE_Xsec_SBND_diff.pdf");
  c1->SaveAs("NuMuE_Xsec_SBND_diff.png");
  std::cout << "[INFO] Saved: NuMuE_Xsec_SBND_diff.{pdf,png}\n";

  // Save to ROOT file
  fOut->cd();
  for (Size_t ie = 0; ie < ENU.size(); ++ie) {
    graphs[ie]->SetName(Form("g_dsig_Enu%.2f", ENU[ie]));
    graphs[ie]->Write();
  }
  gFlux->SetName("g_dsig_fluxavg");
  gFlux->Write();
}

// ─────────────── PLOT 2: TOTAL CROSS SECTION vs Eν ───────────────────────────

void Plot_Total(TFile* fOut) {

  const Int_t nE   = 300;
  const Double_t Elo = 0.05, Ehi = 3.5;   // [GeV]

  // SM prediction
  TGraph* gSM = new TGraph(nE);
  gSM->SetLineColor(TColor::GetColor("#185FA5"));
  gSM->SetLineWidth(3);

  // Analytical approximation at high energy: σ ≈ GF² mₑ Eν (gL² + gR²/3) / π
  // We compute numerically for full accuracy
  for (Int_t i = 0; i < nE; ++i) {
    Double_t Enu = Elo + (Ehi - Elo) * i / (nE - 1);
    Double_t sig = sigma_total(Enu, Phys::sin2tW);
    gSM->SetPoint(i, Enu, sig);
  }

  // ± 1σ uncertainty band on sin²θ_W (world avg δ ≈ 0.0007)
  Double_t dSW = 0.0020;   // conservative SBND uncertainty goal
  TGraph* gUp   = new TGraph(nE);
  TGraph* gDown = new TGraph(nE);
  gUp->SetLineColor(TColor::GetColor("#185FA5"));
  gUp->SetLineStyle(2); gUp->SetLineWidth(2);
  gDown->SetLineColor(TColor::GetColor("#185FA5"));
  gDown->SetLineStyle(2); gDown->SetLineWidth(2);

  for (Int_t i = 0; i < nE; ++i) {
    Double_t Enu = Elo + (Ehi - Elo) * i / (nE - 1);
    gUp->SetPoint(i,   Enu, sigma_total(Enu, Phys::sin2tW + dSW));
    gDown->SetPoint(i, Enu, sigma_total(Enu, Phys::sin2tW - dSW));
  }

  // Shade the band between Up and Down using a TGraph closed polygon
  const Int_t nBand = 2 * nE;
  TGraph* gBand = new TGraph(nBand);
  for (Int_t i = 0; i < nE; ++i) {
    Double_t x, y;
    gUp->GetPoint(i, x, y);
    gBand->SetPoint(i, x, y);
  }
  for (Int_t i = 0; i < nE; ++i) {
    Double_t x, y;
    gDown->GetPoint(nE - 1 - i, x, y);
    gBand->SetPoint(nE + i, x, y);
  }
  gBand->SetFillColorAlpha(TColor::GetColor("#185FA5"), 0.18);
  gBand->SetLineColor(0);
  gBand->SetFillStyle(1001);

  // Approximate linear scaling σ ∝ Eν (shown as dashed reference)
  TGraph* gLin = new TGraph(2);
  gLin->SetPoint(0, Elo, sigma_total(Elo, Phys::sin2tW));
  gLin->SetPoint(1, Ehi, sigma_total(Elo, Phys::sin2tW) * (Ehi / Elo));
  gLin->SetLineColor(kGray + 1);
  gLin->SetLineStyle(3);
  gLin->SetLineWidth(1);

  // ── Draw ──
  TCanvas* c2 = new TCanvas("c_total", "Total cross-section", 900, 680);

  TH1D* hF2 = new TH1D("hFrame_tot", ";E_{#nu}  [GeV];#sigma  [cm^{2}]",
                        100, Elo, Ehi);
  Double_t sigLo = sigma_total(Elo, Phys::sin2tW - dSW) * 0.4;
  Double_t sigHi = sigma_total(Ehi, Phys::sin2tW + dSW) * 1.5;
  hF2->GetYaxis()->SetRangeUser(sigLo, sigHi);
  hF2->Draw("AXIS");
  c2->SetLogy();

  gLin->Draw("L SAME");
  gBand->Draw("F SAME");
  gUp->Draw("L SAME");
  gDown->Draw("L SAME");
  gSM->Draw("L SAME");

  // BNB flux distribution (right-aligned secondary axis, scaled to fit)
  {
    Double_t fluxPeak = 0.0;
    for (Int_t i = 0; i < 200; ++i)
      fluxPeak = std::max(fluxPeak, bnb_flux(Elo + (Ehi-Elo)*i/199));

    TGraph* gFluxBkg = new TGraph();
    Double_t scale   = sigHi * 0.35 / fluxPeak;   // visual scale
    Int_t np = 0;
    for (Int_t i = 0; i < 200; ++i) {
      Double_t E = Elo + (Ehi - Elo) * i / 199.0;
      Double_t f = bnb_flux(E) * scale + sigLo * 1.5;
      gFluxBkg->SetPoint(np++, E, f);
    }
    gFluxBkg->SetFillColorAlpha(TColor::GetColor("#E8A700"), 0.20);
    gFluxBkg->SetLineColor(TColor::GetColor("#BA7517"));
    gFluxBkg->SetLineWidth(1);
    gFluxBkg->SetLineStyle(1);
    gFluxBkg->Draw("LF SAME");

    TLatex ltxF;
    ltxF.SetTextSize(0.036);
    ltxF.SetTextColor(TColor::GetColor("#BA7517"));
    ltxF.DrawLatex(1.4, sigLo * 2.5, "BNB #nu_{#mu} flux (scaled)");
  }

  // Linear-E reference label
  TLatex ltxLin;
  ltxLin.SetTextSize(0.032);
  ltxLin.SetTextColor(kGray + 1);
  ltxLin.DrawLatex(0.20, sigma_total(0.15, Phys::sin2tW)*1.3, "#sigma #propto E_{#nu}");

  TLegend* leg2 = new TLegend(0.16, 0.68, 0.56, 0.91);
  leg2->SetTextSize(0.041);
  leg2->AddEntry(gSM,   "#nu_{#mu} e^{-} elastic (SM, NC)", "l");
  leg2->AddEntry(gBand, Form("#pm#delta sin^{2}#theta_{W} = %.4f", dSW), "f");
  leg2->Draw();

  TLatex ltx;
  ltx.SetNDC(); ltx.SetTextFont(62); ltx.SetTextSize(0.055);
  ltx.DrawLatex(0.55, 0.88, "SBND");
  ltx.SetTextFont(42); ltx.SetTextSize(0.043);
  ltx.DrawLatex(0.55, 0.82, "#nu_{#mu} + e^{-} #rightarrow #nu_{#mu} + e^{-}");
  ltx.DrawLatex(0.55, 0.76, Form("sin^{2}#theta_{W} = %.4f", Phys::sin2tW));

  c2->RedrawAxis();
  c2->SaveAs("NuMuE_Xsec_SBND_total.pdf");
  c2->SaveAs("NuMuE_Xsec_SBND_total.png");
  std::cout << "[INFO] Saved: NuMuE_Xsec_SBND_total.{pdf,png}\n";

  fOut->cd();
  gSM->SetName("g_sigma_total_SM"); gSM->Write();
  gBand->SetName("g_sigma_band");   gBand->Write();
}

// ─────────────── PLOT 3: SENSITIVITY — σ RATIO BAND ──────────────────────────

void Plot_Ratio(TFile* fOut) {

  // Show how the total cross-section changes with sin²θ_W around SM
  const Int_t    nSW = 200;
  const Double_t swLo = 0.10, swHi = 0.40;
  const Double_t Enu  = BNB::Enu_peak;   // evaluate at BNB flux peak

  Double_t sig_SM = sigma_total(Enu, Phys::sin2tW);

  TGraph* gRatio = new TGraph(nSW);
  gRatio->SetLineColor(TColor::GetColor("#E8673A"));
  gRatio->SetLineWidth(3);

  for (Int_t i = 0; i < nSW; ++i) {
    Double_t sw = swLo + (swHi - swLo) * i / (nSW - 1);
    Double_t r  = sigma_total(Enu, sw) / sig_SM;
    gRatio->SetPoint(i, sw, r);
  }

  // Mark SM value
  TGraph* gSMpt = new TGraph(1);
  gSMpt->SetPoint(0, Phys::sin2tW, 1.0);
  gSMpt->SetMarkerStyle(20);
  gSMpt->SetMarkerSize(1.3);
  gSMpt->SetMarkerColor(TColor::GetColor("#185FA5"));

  // ── dσ/dT vs T comparison for low vs high sin²θ_W ──────────────────────
  // (two-panel canvas: left = ratio band, right = dσ/dT for 3 sin²θ_W vals)

  TCanvas* c3 = new TCanvas("c_ratio", "Sensitivity", 1200, 580);
  c3->Divide(2, 1);

  // ── Left pad: σ/σ_SM vs sin²θ_W ──
  c3->cd(1);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);

  TH1D* hFR = new TH1D("hFR", ";sin^{2}#theta_{W};#sigma / #sigma_{SM}", 100, swLo, swHi);
  hFR->GetYaxis()->SetRangeUser(0.0, 4.0);
  hFR->Draw("AXIS");

  // Unity line
  TLine* lOne = new TLine(swLo, 1.0, swHi, 1.0);
  lOne->SetLineStyle(2); lOne->SetLineColor(kGray + 1); lOne->Draw();

  gRatio->Draw("L SAME");
  gSMpt->Draw("P SAME");

  // SM marker text
  TLatex ltxSM;
  ltxSM.SetTextSize(0.040);
  ltxSM.SetTextColor(TColor::GetColor("#185FA5"));
  ltxSM.DrawLatex(Phys::sin2tW + 0.010, 1.08, "SM");

  // Shade SBND sensitivity goal (δsin²θ_W ≈ ±0.005)
  Double_t swGoal = 0.005;
  TBox* bGoal = new TBox(Phys::sin2tW - swGoal, 0.0,
                          Phys::sin2tW + swGoal, 4.0);
  bGoal->SetFillColorAlpha(TColor::GetColor("#185FA5"), 0.10);
  bGoal->SetLineColorAlpha(TColor::GetColor("#185FA5"), 0.50);
  bGoal->SetLineStyle(2);
  bGoal->Draw();

  TLegend* leg3 = new TLegend(0.16, 0.76, 0.92, 0.91);
  leg3->SetTextSize(0.037);
  leg3->AddEntry(gRatio, Form("#sigma(sin^{2}#theta_{W}) / #sigma_{SM}  @ E_{#nu}=%.2f GeV", Enu), "l");
  leg3->AddEntry(bGoal,  "SBND sensitivity goal (#pm0.005)", "f");
  leg3->Draw();

  TLatex ltxL;
  ltxL.SetNDC(); ltxL.SetTextFont(62); ltxL.SetTextSize(0.052);
  ltxL.DrawLatex(0.16, 0.68, "SBND");
  ltxL.SetTextFont(42); ltxL.SetTextSize(0.038);
  ltxL.DrawLatex(0.16, 0.62, "#nu_{#mu} + e^{-} #rightarrow #nu_{#mu} + e^{-}");

  gPad->RedrawAxis();

  // ── Right pad: dσ/dT for 3 sin²θ_W values ──
  c3->cd(2);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04);
  gPad->SetLogx();
  gPad->SetLogy();

  TH1D* hFD = new TH1D("hFD_ratio",
                        ";T_{e}^{kin}  [MeV];d#sigma/dT  [cm^{2}/GeV]",
                        100, 0.1, 700.0);
  hFD->GetXaxis()->SetMoreLogLabels();
  hFD->GetXaxis()->SetNoExponent();
  hFD->GetYaxis()->SetRangeUser(1e-46, 2e-42);
  hFD->Draw("AXIS");

  std::vector<Double_t> swVals = {0.18, Phys::sin2tW, 0.28};
  std::vector<Int_t>    swCols = {
    TColor::GetColor("#E8673A"),
    TColor::GetColor("#185FA5"),
    TColor::GetColor("#3DAA6E"),
  };
  std::vector<TString> swLabels = {"0.18", "0.2312 (SM)", "0.28"};

  std::vector<TGraph*> gSW;
  const Int_t nT2 = 500;

  for (Size_t is = 0; is < swVals.size(); ++is) {
    TGraph* g = new TGraph();
    g->SetLineColor(swCols[is]);
    g->SetLineWidth(2);
    Double_t Tm = T_max(Enu);
    Int_t np = 0;
    for (Int_t it = 0; it < nT2; ++it) {
      Double_t logT = TMath::Log10(0.0001) +
                      (TMath::Log10(0.700) - TMath::Log10(0.0001)) * it / (nT2-1);
      Double_t T  = TMath::Power(10.0, logT);
      if (T > Tm) break;
      Double_t xs = dsigma_dT(T, Enu, swVals[is]);
      if (xs > 0) g->SetPoint(np++, T * 1000.0, xs);
    }
    g->Draw("L SAME");
    gSW.push_back(g);
  }

  TLegend* leg4 = new TLegend(0.16, 0.62, 0.84, 0.91);
  leg4->SetTextSize(0.037);
  leg4->SetHeader(Form("E_{#nu} = %.2f GeV (BNB peak)", Enu), "C");
  for (Size_t is = 0; is < swVals.size(); ++is)
    leg4->AddEntry(gSW[is], Form("sin^{2}#theta_{W} = %s", swLabels[is].Data()), "l");
  leg4->Draw();

  gPad->RedrawAxis();
  c3->SaveAs("NuMuE_Xsec_SBND_ratio.pdf");
  c3->SaveAs("NuMuE_Xsec_SBND_ratio.png");
  std::cout << "[INFO] Saved: NuMuE_Xsec_SBND_ratio.{pdf,png}\n";

  fOut->cd();
  gRatio->SetName("g_sigma_ratio_vs_sw"); gRatio->Write();
}

// ─────────────── PLOT 4: T_max KINEMATIC CURVE ───────────────────────────────

void Plot_TMax(TFile* fOut) {

  const Int_t    nE  = 300;
  const Double_t Elo = 0.05, Ehi = 3.5;

  TGraph* gTmax = new TGraph(nE);
  gTmax->SetLineColor(TColor::GetColor("#8E44AD"));
  gTmax->SetLineWidth(3);

  for (Int_t i = 0; i < nE; ++i) {
    Double_t Enu = Elo + (Ehi - Elo) * i / (nE-1);
    gTmax->SetPoint(i, Enu, T_max(Enu) * 1000.0); // MeV
  }

  // High-energy asymptote T_max → Eν  (m_e/2 correction negligible)
  TGraph* gAsym = new TGraph(2);
  gAsym->SetPoint(0, Elo, Elo * 1000.0);
  gAsym->SetPoint(1, Ehi, Ehi * 1000.0);
  gAsym->SetLineColor(kGray + 1);
  gAsym->SetLineStyle(2);
  gAsym->SetLineWidth(2);

  TCanvas* c4 = new TCanvas("c_tmax", "T_max", 700, 580);

  TH1D* hFT = new TH1D("hFT",
                        ";E_{#nu}  [GeV];T_{max}  [MeV]", 100, Elo, Ehi);
  hFT->GetYaxis()->SetRangeUser(0, Ehi * 1010.0);
  hFT->Draw("AXIS");

  gAsym->Draw("L SAME");
  gTmax->Draw("L SAME");

  // Mark BNB peak
  Double_t Ep   = BNB::Enu_peak;
  Double_t Tmax_p = T_max(Ep) * 1000.0;
  TLine* lEp = new TLine(Ep, 0, Ep, Tmax_p);
  lEp->SetLineColor(TColor::GetColor("#3DAA6E"));
  lEp->SetLineStyle(2); lEp->SetLineWidth(2); lEp->Draw();
  TLine* lTp = new TLine(Elo, Tmax_p, Ep, Tmax_p);
  lTp->SetLineColor(TColor::GetColor("#3DAA6E"));
  lTp->SetLineStyle(2); lTp->SetLineWidth(2); lTp->Draw();

  TLatex ltx;
  ltx.SetTextSize(0.038);
  ltx.SetTextColor(TColor::GetColor("#3DAA6E"));
  ltx.DrawLatex(Ep + 0.05, Tmax_p + 30,
    Form("(%.2f GeV, %.0f MeV)", Ep, Tmax_p));

  TLatex ltxA;
  ltxA.SetTextSize(0.036);
  ltxA.SetTextColor(kGray + 1);
  ltxA.DrawLatex(2.0, 2100, "T_{max} #rightarrow E_{#nu}  (m_{e}#rightarrow0)");

  TLegend* legT = new TLegend(0.15, 0.78, 0.65, 0.91);
  legT->SetTextSize(0.040);
  legT->AddEntry(gTmax, "T_{max} = E_{#nu} / (1 + m_{e}/2E_{#nu})", "l");
  legT->Draw();

  TLatex ltxH;
  ltxH.SetNDC(); ltxH.SetTextFont(62); ltxH.SetTextSize(0.055);
  ltxH.DrawLatex(0.55, 0.22, "SBND");
  ltxH.SetTextFont(42); ltxH.SetTextSize(0.040);
  ltxH.DrawLatex(0.55, 0.16, "#nu_{#mu} e^{-} kinematics");

  c4->RedrawAxis();
  c4->SaveAs("NuMuE_Xsec_SBND_tmax.pdf");
  c4->SaveAs("NuMuE_Xsec_SBND_tmax.png");
  std::cout << "[INFO] Saved: NuMuE_Xsec_SBND_tmax.{pdf,png}\n";

  fOut->cd();
  gTmax->SetName("g_tmax_vs_enu"); gTmax->Write();
}

// ─────────────────────────── MAIN ────────────────────────────────────────────

void NuMuE_Xsec_SBND() {

  std::cout << "\n"
    << "╔══════════════════════════════════════════════════════════════════╗\n"
    << "║  ν_μ + e⁻ → ν_μ + e⁻  Cross-Section  — SBND / BNB                ║\n"
    << "╚══════════════════════════════════════════════════════════════════╝\n\n";

  SetSBNDStyle();

  // Open output ROOT file
  TFile* fOut = TFile::Open("NuMuE_Xsec_SBND_results.root", "RECREATE");
  if (!fOut || fOut->IsZombie()) {
    std::cerr << "[ERROR] Cannot create output ROOT file.\n";
    return;
  }

  // Print table of key values
  std::cout << "  Eν [GeV]  |  T_max [MeV]  |  σ [cm²]      |  dσ/dT|_{T=0} [cm²/GeV]\n"
            << "  ─────────────────────────────────────────────────────────────────────\n";
  for (Double_t Enu : {0.20, 0.40, 0.70, 1.00, 1.50, 2.00, 3.00}) {
    printf("  %7.2f    |  %9.1f    |  %10.3e   |  %10.3e\n",
           Enu,
           T_max(Enu) * 1000.0,
           sigma_total(Enu, Phys::sin2tW),
           dsigma_dT(1e-6, Enu, Phys::sin2tW));
  }
  std::cout << "\n";

  // Print gL, gR
  Double_t gL = -0.5 + Phys::sin2tW;
  Double_t gR =        Phys::sin2tW;
  printf("  sin²θ_W = %.4f   gL = % .4f   gR = %.4f\n\n",
         Phys::sin2tW, gL, gR);

  // Generate all plots
  Plot_Differential(fOut);
  Plot_Total(fOut);
  Plot_Ratio(fOut);
  Plot_TMax(fOut);

  fOut->Close();
  std::cout << "\n[INFO] All results saved to NuMuE_Xsec_SBND_results.root\n\n";
}
