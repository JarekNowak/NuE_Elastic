/*=============================================================================
  NuMuE_FluxAvg_sin2tW_SBND.C
  ─────────────────────────────────────────────────────────────────────────────
  Plots the BNB ν_μ flux-averaged differential cross-section

        〈dσ/dT〉_Φ  =  ∫ Φ(Eν) · dσ/dT · dEν  /  ∫ Φ(Eν) dEν

  for ν_μ + e⁻ → ν_μ + e⁻  as a function of the recoil electron kinetic
  energy T, for several values of sin²θ_W spanning the SBND sensitivity range.

  Physics:
    dσ/dT = (G_F² mₑ / 2π) [ gL² + gR²(1 − T/Eν)² − gL·gR·mₑT/Eν² ]
    gL = −½ + sin²θ_W ,   gR = sin²θ_W
    T_max(Eν) = Eν / (1 + mₑ/2Eν)

  Outputs:
    NuMuE_FluxAvg_sin2tW_SBND.pdf / .png   – main plot
    NuMuE_FluxAvg_sin2tW_SBND_ratio.pdf / .png  – ratio to SM
    NuMuE_FluxAvg_sin2tW_SBND.root          – TGraphs saved

  Usage:
    root -l -b -q NuMuE_FluxAvg_sin2tW_SBND.C

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
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TFile.h"
#include "TPad.h"
#include "TColor.h"
#include "TAxis.h"
#include "TFrame.h"

// ─────────────────────────── CONSTANTS ───────────────────────────────────────
namespace Phys {
  const Double_t GF     = 1.16638e-5;    // Fermi constant  [GeV⁻²]
  const Double_t me     = 0.511e-3;      // electron mass   [GeV]
  const Double_t hbarc2 = 3.8938e-28;   // ħc²             [cm²·GeV²]
  const Double_t sin2tW_SM = 0.2312;     // SM value
}
namespace BNB {
  const Double_t Enu_min = 0.10;         // [GeV]
  const Double_t Enu_max = 3.00;         // [GeV]
}

// ─────────────────────────── PHYSICS ─────────────────────────────────────────

inline Double_t T_max(Double_t Enu) {
  return Enu / (1.0 + Phys::me / (2.0 * Enu));
}

Double_t dsigma_dT(Double_t T, Double_t Enu, Double_t sw2) {
  if (T < 0.0 || Enu <= 0.0 || T > T_max(Enu)) return 0.0;
  Double_t gL = -0.5 + sw2, gR = sw2;
  Double_t y  = T / Enu;
  Double_t xs = (Phys::GF * Phys::GF * Phys::me / (2.0 * TMath::Pi())) * Phys::hbarc2
                * (gL*gL + gR*gR*(1-y)*(1-y) - gL*gR*Phys::me*T/(Enu*Enu));
  return std::max(xs, 0.0);
}

/// Approximate BNB ν_μ flux shape [arb. units / GeV]
Double_t bnb_flux(Double_t Enu) {
  if (Enu < 0) return 0.0;
  Double_t phi = TMath::Gaus(Enu, 0.70, 0.22)
               + 0.20 * TMath::Gaus(Enu, 1.30, 0.35)
               + 0.05 * TMath::Gaus(Enu, 2.00, 0.50);
  return std::max(phi, 0.0);
}

/// 〈dσ/dT〉_Φ  [cm²/GeV]  — flux-averaged differential cross-section
Double_t flux_avg_dsig(Double_t T, Double_t sw2, Int_t nE = 200) {
  Double_t dE    = (BNB::Enu_max - BNB::Enu_min) / nE;
  Double_t num   = 0.0;
  Double_t denom = 0.0;
  for (Int_t k = 0; k < nE; ++k) {
    Double_t Enu = BNB::Enu_min + (k + 0.5) * dE;
    Double_t phi = bnb_flux(Enu) * dE;
    num   += phi * dsigma_dT(T, Enu, sw2);
    denom += phi;
  }
  return (denom > 0.0) ? num / denom : 0.0;
}

/// 〈σ〉_Φ  [cm²] — flux-averaged total cross-section (for table)
Double_t flux_avg_sigma(Double_t sw2, Int_t nE = 200, Int_t nT = 400) {
  Double_t dE    = (BNB::Enu_max - BNB::Enu_min) / nE;
  Double_t num   = 0.0;
  Double_t denom = 0.0;
  for (Int_t k = 0; k < nE; ++k) {
    Double_t Enu = BNB::Enu_min + (k + 0.5) * dE;
    Double_t phi = bnb_flux(Enu) * dE;
    Double_t Tm  = T_max(Enu);
    Double_t dT  = Tm / nT;
    Double_t sig = 0.0;
    for (Int_t j = 0; j < nT; ++j) {
      Double_t T1 = j * dT, T2 = (j+1) * dT;
      sig += 0.5 * (dsigma_dT(T1, Enu, sw2) + dsigma_dT(T2, Enu, sw2)) * dT;
    }
    num   += phi * sig;
    denom += phi;
  }
  return (denom > 0.0) ? num / denom : 0.0;
}

// ─────────────────────────── STYLE ───────────────────────────────────────────
void SetStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.135);
  gStyle->SetPadRightMargin(0.040);
  gStyle->SetPadTopMargin(0.060);
  gStyle->SetPadBottomMargin(0.130);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillStyle(0);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.052, "XYZ");
  gStyle->SetLabelSize(0.046, "XYZ");
  gStyle->SetTitleOffset(1.10, "X");
  gStyle->SetTitleOffset(1.15, "Y");
}

// ─────────────────────────── MAIN ────────────────────────────────────────────
void NuMuE_FluxAvg_sin2tW_SBND() {

  std::cout << "\n"
    << "╔═══════════════════════════════════════════════════════════════════╗\n"
    << "║  BNB flux-averaged dσ/dT  —  various sin²θ_W  —  SBND          ║\n"
    << "╚═══════════════════════════════════════════════════════════════════╝\n\n";

  SetStyle();

  // ── sin²θ_W values to show ────────────────────────────────────────────────
  // SM = 0.2312; span ±0.08 to bracket SBND sensitivity and BSM scenarios
  const std::vector<Double_t> SW = {
    0.15, 0.18, 0.21, 0.2312, 0.25, 0.28, 0.31, 0.35
  };

  // Colour palette (colour-blind friendly, warm→cool ordered)
  const std::vector<Int_t> COLS = {
    TColor::GetColor("#C0392B"),   // 0.15  deep red
    TColor::GetColor("#E8673A"),   // 0.18  coral
    TColor::GetColor("#E8A700"),   // 0.21  amber
    TColor::GetColor("#185FA5"),   // 0.2312  blue (SM) — bold
    TColor::GetColor("#3DAA6E"),   // 0.25  green
    TColor::GetColor("#2980B9"),   // 0.28  mid-blue
    TColor::GetColor("#8E44AD"),   // 0.31  purple
    TColor::GetColor("#2C3E50"),   // 0.35  dark slate
  };

  // ── Build TGraphs ─────────────────────────────────────────────────────────
  const Int_t    NP   = 600;           // points along T axis
  const Double_t TloG = 3e-5;          // [GeV]  lower edge (~30 eV — very forward)
  const Double_t ThiG = 2.50;          // [GeV]  upper edge

  std::vector<TGraph*> graphs(SW.size(), nullptr);

  for (Size_t is = 0; is < SW.size(); ++is) {
    TGraph* g = new TGraph();
    Int_t   np = 0;

    for (Int_t it = 0; it < NP; ++it) {
      // Logarithmic T spacing for better resolution near threshold
      Double_t logT = TMath::Log10(TloG) +
                      (TMath::Log10(ThiG) - TMath::Log10(TloG)) * it / (NP - 1);
      Double_t T  = TMath::Power(10.0, logT);
      Double_t xs = flux_avg_dsig(T, SW[is]);
      if (xs > 0.0) g->SetPoint(np++, T * 1000.0, xs);   // T → MeV on x-axis
    }

    Bool_t isSM = (SW[is] == Phys::sin2tW_SM);
    g->SetLineColor(COLS[is]);
    g->SetLineWidth(isSM ? 4 : 2);
    g->SetLineStyle(isSM ? 1 : 1);
    graphs[is] = g;
  }

  // ── Print summary table ───────────────────────────────────────────────────
  Double_t sig_SM = flux_avg_sigma(Phys::sin2tW_SM);
  std::cout << "  sin²θ_W  |  〈σ〉_Φ [cm²]    |  ratio to SM\n"
            << "  ──────────────────────────────────────────────\n";
  for (Size_t is = 0; is < SW.size(); ++is) {
    Double_t sig = flux_avg_sigma(SW[is]);
    Bool_t   isSM = (SW[is] == Phys::sin2tW_SM);
    printf("  %.4f   |  %10.3e      |  %.4f%s\n",
           SW[is], sig, sig / sig_SM, isSM ? "  ← SM" : "");
  }
  std::cout << "\n";

  // ─────────────────── CANVAS 1: dσ/dT vs T ────────────────────────────────
  TCanvas* c1 = new TCanvas("c_fluxavg", "Flux-avg dSigma/dT vs sin2tW", 960, 720);

  // Invisible frame to define axes
  TH1D* hF = new TH1D("hFrm",
                       ";T_{e}^{kin}  [MeV];"
                       "#langle d#sigma/dT #rangle_{#Phi}  [cm^{2}/GeV]",
                       100, 100, 2000.0);
  hF->GetXaxis()->SetMoreLogLabels();
  hF->GetXaxis()->SetNoExponent();
  hF->GetYaxis()->SetRangeUser(2e-47, 1e-42);
  hF->Draw("AXIS");
//  c1->SetLogx();
//  c1->SetLogy();

  // Draw all curves
  for (auto g : graphs) g->Draw("L SAME");

  // ── Legend ────────────────────────────────────────────────────────────────
  TLegend* leg = new TLegend(0.655, 0.34, 0.820, 0.62);
  leg->SetTextSize(0.038);
  leg->SetTextFont(42);
  leg->SetHeader("sin^{2}#theta_{W}", "C");
  for (Size_t is = 0; is < SW.size(); ++is) {
    Bool_t isSM = (SW[is] == Phys::sin2tW_SM);
    TString label = Form("%.4f%s", SW[is], isSM ? "  (SM)" : "");
    leg->AddEntry(graphs[is], label, "l");
  }
  leg->Draw();

  // ── Labels ────────────────────────────────────────────────────────────────
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(62);
  ltx.SetTextSize(0.058);
  ltx.DrawLatex(0.60, 0.88, "SBND");

  ltx.SetTextFont(42);
  ltx.SetTextSize(0.042);
  ltx.DrawLatex(0.60, 0.82, "#nu_{#mu} + e^{#minus} #rightarrow #nu_{#mu} + e^{#minus}");
  ltx.DrawLatex(0.60, 0.76, "BNB #nu_{#mu}  flux average");

  // Mark SM curve with a floating label on the canvas
  {
    // find y value of SM curve at T ~ 5 MeV
    Double_t T_label_MeV = 5.0;
    Double_t y_label = flux_avg_dsig(T_label_MeV * 1e-3, Phys::sin2tW_SM);
    TLatex ltxSM;
    ltxSM.SetTextSize(0.034);
    ltxSM.SetTextColor(COLS[3]);
    ltxSM.DrawLatex(T_label_MeV * 1.3, y_label * 0.40, "SM");
  }

  // ── Arrow indicating direction of increasing sin²θ_W ─────────────────────
  {
    // draw a small curved annotation bracket on the plot
    // At T = 200 MeV, read off min and max curve y-values
    Double_t T_anno = 200e-3;   // 200 MeV in GeV
    Double_t y_lo   = flux_avg_dsig(T_anno, SW.front());
    Double_t y_hi   = flux_avg_dsig(T_anno, SW.back());
    TLatex ltxA;
    ltxA.SetTextSize(0.036);
    ltxA.SetTextColor(kGray + 1);
    ltxA.DrawLatex(250.0, y_lo * 0.60,
                   "sin^{2}#theta_{W} #uparrow");
  }

  c1->RedrawAxis();
  c1->SaveAs("NuMuE_FluxAvg_sin2tW_SBND.pdf");
  c1->SaveAs("NuMuE_FluxAvg_sin2tW_SBND.png");
  std::cout << "[INFO] Saved: NuMuE_FluxAvg_sin2tW_SBND.{pdf,png}\n";

  // ─────────────────── CANVAS 2: ratio to SM ───────────────────────────────
  TCanvas* c2 = new TCanvas("c_ratio", "Ratio to SM", 960, 720);

  // Build ratio graphs (skip the SM curve itself — it's unity)
  std::vector<TGraph*> gRatios;

  // SM reference
  std::vector<Double_t> Tpts_SM, Ypts_SM;
  {
    Int_t np0 = graphs[3]->GetN();  // SM index is 3 (sin2tW = 0.2312)
    for (Int_t i = 0; i < np0; ++i) {
      Double_t x, y;
      graphs[3]->GetPoint(i, x, y);
      if (y > 0.0) { Tpts_SM.push_back(x); Ypts_SM.push_back(y); }
    }
  }

  for (Size_t is = 0; is < SW.size(); ++is) {
    if (SW[is] == Phys::sin2tW_SM) { gRatios.push_back(nullptr); continue; }

    TGraph* gr = new TGraph();
    gr->SetLineColor(COLS[is]);
    gr->SetLineWidth(2);

    Int_t npR = graphs[is]->GetN();
    Int_t nrp = 0;
    for (Int_t i = 0; i < npR; ++i) {
      Double_t x, y;
      graphs[is]->GetPoint(i, x, y);
      if (y <= 0.0) continue;

      // Interpolate SM value at same x
      Double_t ySM = 0.0;
      for (Size_t j = 1; j < Tpts_SM.size(); ++j) {
        if (Tpts_SM[j] >= x) {
          Double_t frac = (x - Tpts_SM[j-1]) / (Tpts_SM[j] - Tpts_SM[j-1]);
          ySM = Ypts_SM[j-1] + frac * (Ypts_SM[j] - Ypts_SM[j-1]);
          break;
        }
      }
      if (ySM > 0.0) gr->SetPoint(nrp++, x, y / ySM);
    }
    gRatios.push_back(gr);
  }

  TH1D* hFR = new TH1D("hFrmR",
                        ";T_{e}^{kin}  [MeV];"
                        "#langle d#sigma/dT #rangle_{#Phi} / #langle d#sigma/dT #rangle_{#Phi}^{SM}",
                        100, 100, 2000.0);
  hFR->GetXaxis()->SetMoreLogLabels();
  hFR->GetXaxis()->SetNoExponent();
  hFR->GetYaxis()->SetRangeUser(0.0, 3.6);
  hFR->Draw("AXIS");
  c2->SetLogx();

  // Unity (SM) line
  TLine* lOne = new TLine(100, 1.0, 2000.0, 1.0);
  lOne->SetLineColor(COLS[3]);
  lOne->SetLineWidth(3);
  lOne->SetLineStyle(1);
  lOne->Draw();

  for (Size_t is = 0; is < SW.size(); ++is) {
    if (gRatios[is]) gRatios[is]->Draw("L SAME");
  }

  // Legend
  TLegend* leg2 = new TLegend(0.155, 0.58, 0.560, 0.91);
  leg2->SetTextSize(0.038);
  leg2->SetHeader("sin^{2}#theta_{W}", "C");
  for (Size_t is = 0; is < SW.size(); ++is) {
    Bool_t isSM = (SW[is] == Phys::sin2tW_SM);
    TGraph* gEntry = isSM ? (TGraph*)lOne : gRatios[is];
    TString label  = Form("%.4f%s", SW[is], isSM ? "  (SM = 1)" : "");
    if (gEntry) leg2->AddEntry(gEntry, label, "l");
  }
  leg2->Draw();

  TLatex ltx2;
  ltx2.SetNDC();
  ltx2.SetTextFont(62); ltx2.SetTextSize(0.058);
  ltx2.DrawLatex(0.60, 0.88, "SBND");
  ltx2.SetTextFont(42); ltx2.SetTextSize(0.042);
  ltx2.DrawLatex(0.60, 0.82, "#nu_{#mu} + e^{#minus} #rightarrow #nu_{#mu} + e^{#minus}");
  ltx2.DrawLatex(0.60, 0.76, "BNB #nu_{#mu}  flux average");
  ltx2.DrawLatex(0.60, 0.70, "Ratio to SM prediction");

  c2->RedrawAxis();
  c2->SaveAs("NuMuE_FluxAvg_sin2tW_SBND_ratio.pdf");
  c2->SaveAs("NuMuE_FluxAvg_sin2tW_SBND_ratio.png");
  std::cout << "[INFO] Saved: NuMuE_FluxAvg_sin2tW_SBND_ratio.{pdf,png}\n";

  // ── Save all to ROOT file ─────────────────────────────────────────────────
  TFile* fOut = TFile::Open("NuMuE_FluxAvg_sin2tW_SBND.root", "RECREATE");
  if (fOut && !fOut->IsZombie()) {
    for (Size_t is = 0; is < SW.size(); ++is) {
      graphs[is]->SetName(Form("g_fluxavg_sw%.4f", SW[is]));
      graphs[is]->Write();
      if (gRatios[is]) {
        gRatios[is]->SetName(Form("g_ratio_sw%.4f", SW[is]));
        gRatios[is]->Write();
      }
    }
    fOut->Close();
    std::cout << "[INFO] Saved: NuMuE_FluxAvg_sin2tW_SBND.root\n\n";
  }
}
