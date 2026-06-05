/*=============================================================================
  NuEOptimize.C
  ─────────────────────────────────────────────────────────────────────────────
  Fast cut tuner for the ν_μ + e⁻ elastic selection, run on the candidate skim
  produced by NuESkim.C.

  Figure of merit:
      FOM = S / √(S + B)      (POT-scaled to 1×10²¹ POT)
  with S, B the POT-scaled selected signal / background event counts. Each MC
  event is weighted by its per-sample factor (signal ×SIGNAL, ν-bkg ×NU_BKG),
  so the optimiser maximises the statistical significance of the *real* expected
  yields — not the misleading raw-MC purity.

  Method: coordinate-ascent. Start from the current NuESelection cuts, scan each
  threshold over a grid holding the others fixed, adopt the best, and iterate
  until stable. A minimum raw-MC-statistics guard prevents tuning into MC noise.

  It then writes fit-ready output histograms (h_data, h_bkg, h_sig, h_eff,
  h_smear) for the chosen working point — no second pass over the 143 GB file.

  Usage:
    root -l -b -q NuEOptimize.C
    root -l -b -q 'NuEOptimize.C("NuESkim.root","NuESelection_output.root")'
=============================================================================*/

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVectorD.h"

// ─────────────────────── candidate record (in memory) ───────────────────────
struct Cand {
  int    isSig, nuE, nExtra;
  float  eReco, eTrue, sliceInter, sliceScore;
  float  rzE, rzG, rzMu, rzPi, rzP, bestPDG;
  float  trk, theta, thetaPart, etheta2;
  float  dedx, dedx0, dedx1, dedx2;
  float  shwLen, shwOA;
};

// ─────────────────────────── tunable cut set ────────────────────────────────
struct CutSet {
  double rzE_min   = 0.50;   // P(e) ≥
  double rzG_max   = 1.01;   // P(γ) ≤   (1.01 = no cut)
  double trk_max   = 0.50;   // track score ≤
  double eMin      = 0.030;  // shower E ≥ [GeV]
  double eMax      = 1.500;  // shower E ≤ [GeV]
  double dedx_min  = 0.50;   // MeV/cm
  double dedx_max  = 3.50;   // MeV/cm
  double theta_max = 0.55;   // rad
  double et2_max   = 0.010;  // E×θ² ≤ [GeV]
  double oa_min    = 0.00;   // shower open angle ≥ [rad] (0 = no cut)
  int    maxExtra  = 0;      // extra significant primaries ≤
};

static double SF_SIG = 0.0105845, SF_BKG = 1.65585;
// Guards against tuning into MC-statistics noise: require enough surviving raw
// MC events that both the signal efficiency AND the background prediction are
// statistically meaningful (≥50 bkg events ⇒ <15% stat error on B).
static const int MIN_RAW_SIG = 80;
static const int MIN_RAW_BKG = 50;

// pass / fail for a single candidate
inline bool passAll(const Cand& c, const CutSet& x) {
  if (c.rzE < x.rzE_min) return false;
  if (c.rzG > x.rzG_max) return false;
  if (c.trk < -999000 || c.trk > x.trk_max) return false;
  if (c.eReco < x.eMin || c.eReco > x.eMax) return false;
  if (c.dedx < -990 || c.dedx < x.dedx_min || c.dedx > x.dedx_max) return false;
  if (c.theta < 0 || c.theta > x.theta_max) return false;
  if (c.etheta2 > x.et2_max) return false;
  if (c.shwOA < x.oa_min) return false;
  if (c.nExtra > x.maxExtra) return false;
  return true;
}

struct Yield { double S=0, B=0; long rawS=0, rawB=0; };

Yield evalYield(const std::vector<Cand>& v, const CutSet& x) {
  Yield y;
  for (const auto& c : v) {
    if (!passAll(c, x)) continue;
    if (c.isSig) { y.S += SF_SIG; y.rawS++; }
    else         { y.B += SF_BKG; y.rawB++; }
  }
  return y;
}

inline double fom(const Yield& y) {
  if (y.rawS < MIN_RAW_SIG) return -1.0;          // not enough MC signal stats
  if (y.rawB < MIN_RAW_BKG) return -1.0;          // background prediction too noisy
  double tot = y.S + y.B;
  return (tot > 0) ? y.S / std::sqrt(tot) : 0.0;
}

// scan one cut over a grid, return value maximising FOM
template <class Setter>
double scanCut(const std::vector<Cand>& v, CutSet& x,
               const std::vector<double>& grid, Setter set, const char* name) {
  double bestF = -1, bestV = 0; bool any=false;
  for (double g : grid) {
    CutSet t = x; set(t, g);
    double f = fom(evalYield(v, t));
    if (f > bestF) { bestF = f; bestV = g; any=true; }
  }
  if (any) set(x, bestV);
  printf("    %-12s → %8.4f   (FOM=%.3f)\n", name, bestV, bestF);
  return bestF;
}

std::vector<double> arange(double lo, double hi, double step) {
  std::vector<double> g;
  for (double v = lo; v <= hi + 1e-9; v += step) g.push_back(v);
  return g;
}

void report(const char* tag, const std::vector<Cand>& v, const CutSet& x) {
  Yield y = evalYield(v, x);
  double tot = y.S + y.B;
  double pur = tot>0 ? 100.0*y.S/tot : 0;
  printf("  %-10s  S=%8.1f  B=%8.1f  tot=%8.1f  purity=%5.1f%%  "
         "FOM=%.3f   (rawS=%ld rawB=%ld)\n",
         tag, y.S, y.B, tot, pur, (tot>0?y.S/std::sqrt(tot):0), y.rawS, y.rawB);
}

// ─────────────────────────────── MAIN ────────────────────────────────────────
void NuEOptimize(const char* skimFile = "NuESkim.root",
                 const char* outFile  = "NuESelection_output.root")
{
  TFile* f = TFile::Open(skimFile, "READ");
  if (!f || f->IsZombie()) { std::cerr << "[ERROR] cannot open " << skimFile << "\n"; return; }
  TTree* t = (TTree*)f->Get("cand");
  TH1D*  h_sig_total = (TH1D*)f->Get("h_sig_total");
  TVectorD* v_meta    = (TVectorD*)f->Get("v_meta");
  TVectorD* v_preflow = (TVectorD*)f->Get("v_preflow");
  if (!t || !h_sig_total || !v_meta || !v_preflow) { std::cerr << "[ERROR] bad skim file\n"; return; }

  SF_SIG = (*v_meta)[1];
  SF_BKG = (*v_meta)[2];
  double preAll_sig = (*v_preflow)[1];   // all true signal (efficiency denom count)

  // ── load tree into memory ─────────────────────────────────────────────────
  Cand c;
  Int_t isSig, nuE, nExtra;
  Double_t eReco,eTrue,sliceInter,sliceScore,rzE,rzG,rzMu,rzPi,rzP,bestPDG,
           trk,theta,thetaPart,etheta2,dedx,dedx0,dedx1,dedx2,shwLen,shwOA;
  t->SetBranchAddress("isSig",&isSig);  t->SetBranchAddress("nuE",&nuE);
  t->SetBranchAddress("nExtraPrim",&nExtra);
  t->SetBranchAddress("eReco",&eReco);  t->SetBranchAddress("eTrue",&eTrue);
  t->SetBranchAddress("sliceInter",&sliceInter); t->SetBranchAddress("sliceScore",&sliceScore);
  t->SetBranchAddress("rzE",&rzE); t->SetBranchAddress("rzG",&rzG);
  t->SetBranchAddress("rzMu",&rzMu); t->SetBranchAddress("rzPi",&rzPi); t->SetBranchAddress("rzP",&rzP);
  t->SetBranchAddress("bestPDG",&bestPDG); t->SetBranchAddress("trk",&trk);
  t->SetBranchAddress("theta",&theta); t->SetBranchAddress("thetaPart",&thetaPart);
  t->SetBranchAddress("etheta2",&etheta2);
  t->SetBranchAddress("dedx",&dedx); t->SetBranchAddress("dedx0",&dedx0);
  t->SetBranchAddress("dedx1",&dedx1); t->SetBranchAddress("dedx2",&dedx2);
  t->SetBranchAddress("shwLen",&shwLen); t->SetBranchAddress("shwOA",&shwOA);

  std::vector<Cand> V; V.reserve(t->GetEntries());
  for (Long64_t i=0;i<t->GetEntries();++i){
    t->GetEntry(i);
    c.isSig=isSig; c.nuE=nuE; c.nExtra=nExtra;
    c.eReco=eReco; c.eTrue=eTrue; c.sliceInter=sliceInter; c.sliceScore=sliceScore;
    c.rzE=rzE; c.rzG=rzG; c.rzMu=rzMu; c.rzPi=rzPi; c.rzP=rzP; c.bestPDG=bestPDG;
    c.trk=trk; c.theta=theta; c.thetaPart=thetaPart; c.etheta2=etheta2;
    c.dedx=dedx; c.dedx0=dedx0; c.dedx1=dedx1; c.dedx2=dedx2;
    c.shwLen=shwLen; c.shwOA=shwOA;
    V.push_back(c);
  }
  printf("[INFO] Loaded %zu candidates. SF_sig=%.6f SF_bkg=%.5f  allSig=%.0f\n\n",
         V.size(), SF_SIG, SF_BKG, preAll_sig);

  // ── current (baseline) cuts ───────────────────────────────────────────────
  CutSet cur;  // defaults above = current NuESelection.C values
  printf("══════════════ BASELINE (current NuESelection cuts) ══════════════\n");
  report("baseline", V, cur);
  printf("\n");

  // ── coordinate ascent ─────────────────────────────────────────────────────
  CutSet x = cur;
  auto gRzE  = arange(0.30,0.98,0.02);
  auto gRzG  = arange(0.05,1.00,0.05); gRzG.push_back(1.01);          // 1.01 = off
  auto gTrk  = arange(0.10,0.60,0.025);
  auto gEmin = arange(0.010,0.120,0.005);
  auto gEmax = arange(0.50,1.50,0.05);
  auto gDdMin= arange(0.00,1.80,0.10);
  auto gDdMax= arange(2.50,6.00,0.25);
  auto gTh   = arange(0.05,0.70,0.025);
  auto gEt2  = arange(0.0005,0.0150,0.0005);
  auto gOA   = arange(0.00,0.16,0.01);
  std::vector<double> gExtra = {0,1,2,9};

  double prevF = -1;
  for (int pass=1; pass<=5; ++pass) {
    printf("──────── coordinate-ascent pass %d ────────\n", pass);
    scanCut(V,x,gRzE ,[](CutSet&t,double g){t.rzE_min=g;}  ,"P(e)>=");
    scanCut(V,x,gRzG ,[](CutSet&t,double g){t.rzG_max=g;}  ,"P(g)<=");
    scanCut(V,x,gTrk ,[](CutSet&t,double g){t.trk_max=g;}  ,"trk<=");
    scanCut(V,x,gEmin,[](CutSet&t,double g){t.eMin=g;}     ,"Emin");
    scanCut(V,x,gEmax,[](CutSet&t,double g){t.eMax=g;}     ,"Emax");
    scanCut(V,x,gDdMin,[](CutSet&t,double g){t.dedx_min=g;},"dEdx_min");
    scanCut(V,x,gDdMax,[](CutSet&t,double g){t.dedx_max=g;},"dEdx_max");
    scanCut(V,x,gTh  ,[](CutSet&t,double g){t.theta_max=g;},"theta<=");
    scanCut(V,x,gEt2 ,[](CutSet&t,double g){t.et2_max=g;}  ,"Eth2<=");
    scanCut(V,x,gOA  ,[](CutSet&t,double g){t.oa_min=g;}   ,"OA>=");
    double f = scanCut(V,x,gExtra,[](CutSet&t,double g){t.maxExtra=(int)g;},"extraPrim<=");
    report("pass", V, x);
    printf("\n");
    if (std::fabs(f-prevF) < 1e-4) { printf("[INFO] converged.\n\n"); break; }
    prevF = f;
  }

  // ── final summary ─────────────────────────────────────────────────────────
  Yield yB = evalYield(V, cur), yO = evalYield(V, x);
  printf("══════════════════════ OPTIMISED WORKING POINT ══════════════════════\n");
  report("baseline", V, cur);
  report("optimised", V, x);
  double effB = preAll_sig>0 ? 100.0*yB.rawS/preAll_sig : 0;
  double effO = preAll_sig>0 ? 100.0*yO.rawS/preAll_sig : 0;
  printf("\n  signal eff:  baseline %.2f%%   optimised %.2f%%\n", effB, effO);
  printf("  FOM gain  :  %.3f → %.3f  (×%.2f)\n",
         fom(yB), fom(yO), fom(yB)>0?fom(yO)/fom(yB):0);
  printf("\n  Optimised cuts:\n");
  printf("    RAZZLE_ELEC_MIN = %.3f\n", x.rzE_min);
  printf("    RAZZLE_GAMMA_MAX= %.3f %s\n", x.rzG_max, x.rzG_max>1.0?"(off)":"");
  printf("    TRACK_SCORE_MAX = %.3f\n", x.trk_max);
  printf("    SHOWER_E_MIN    = %.3f GeV\n", x.eMin);
  printf("    SHOWER_E_MAX    = %.3f GeV\n", x.eMax);
  printf("    DEDX_MIN        = %.2f MeV/cm\n", x.dedx_min);
  printf("    DEDX_MAX        = %.2f MeV/cm\n", x.dedx_max);
  printf("    THETA_MAX       = %.3f rad\n", x.theta_max);
  printf("    ETHETA2_MAX     = %.4f GeV\n", x.et2_max);
  printf("    SHOWER_OA_MIN   = %.3f rad %s\n", x.oa_min, x.oa_min<=0?"(off)":"");
  printf("    MAX_EXTRA_PRIM  = %d\n", x.maxExtra);
  printf("════════════════════════════════════════════════════════════════════\n\n");

  // ── build fit-ready output histograms at optimised working point ───────────
  const int N=120; const double lo=0.0, hi=1.2;
  TH1D* h_data = new TH1D("h_data","Selected;T_{e}^{reco} [GeV];Events/bin",N,lo,hi);
  TH1D* h_bkg  = new TH1D("h_bkg" ,"Background;T_{e}^{reco} [GeV];Events/bin",N,lo,hi);
  TH1D* h_sig  = new TH1D("h_sig" ,"Signal (reco);T_{e}^{reco} [GeV];Events/bin",N,lo,hi);
  TH1D* h_sig_true = new TH1D("h_sig_true","Selected signal;T_{e}^{true} [GeV];",N,lo,hi);
  TH2D* h_smear= new TH2D("h_smear","Smearing;T_{e}^{reco} [GeV];T_{e}^{true} [GeV]",N,lo,hi,N,lo,hi);
  TH1D* h_eff  = new TH1D("h_eff","Efficiency;T_{e}^{true} [GeV];#epsilon",N,lo,hi);
  h_eff->Sumw2();

  for (const auto& cc : V) {
    if (!passAll(cc, x)) continue;
    if (cc.isSig) {
      h_sig->Fill(cc.eReco);
      if (cc.eTrue >= 0) { h_sig_true->Fill(cc.eTrue); h_smear->Fill(cc.eReco, cc.eTrue); }
    } else {
      h_bkg->Fill(cc.eReco);
    }
    h_data->Fill(cc.eReco);
  }
  // column-normalise smearing matrix (per true-energy column)
  for (int iT=1;iT<=N;++iT){
    double col=0; for(int iR=1;iR<=N;++iR) col+=h_smear->GetBinContent(iR,iT);
    if(col>0) for(int iR=1;iR<=N;++iR) h_smear->SetBinContent(iR,iT,h_smear->GetBinContent(iR,iT)/col);
  }
  // efficiency vs TRUE energy: selected signal(true) / all signal(true)
  h_eff->Divide(h_sig_true, h_sig_total, 1.0, 1.0, "B");

  TFile* fo = TFile::Open(outFile,"RECREATE");
  h_data->Write(); h_bkg->Write(); h_sig->Write(); h_sig_true->Write();
  h_eff->Write(); h_smear->Write(); h_sig_total->Write();
  fo->Close();
  printf("[INFO] Fit-ready histograms written → %s\n\n", outFile);

  f->Close();
}
