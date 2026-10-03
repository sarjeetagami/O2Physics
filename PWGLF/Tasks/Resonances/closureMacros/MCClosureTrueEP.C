// MC closure test of the rho00 (spin alignment) analysis chain of kstarpbpb.cxx
//
// MC TRUE Psi (no event-plane smearing): for the wagon with
//   processMCClosure, processMCClosureME, processMCClosurePhi, processMCClosureMEPhi = true
// Reads <TaskDir>/closureKstar and <TaskDir>/closurePhi (the reconstructed-Psi case is MCClosureRecoEP.C)
// and, for each centrality and pT bin, runs the full chain:
//   1. raw yield vs cos(theta*): SE - normalised ME, residual background fitted and subtracted, bin counting in a mass window
//   2. compare raw yield with the true reconstructed yield (hSparseRecTrue, same window)
//   3. Acc x Eff(cos theta*) = hSparseRecTrueEff (window) / hSparseGenEff (all masses), from the independent efficiency sample
//   4. fit the corrected dN/dcos(theta*) with (1 - rho) + (3 rho - 1) cos^2 -> rho00_obs
//   5. compare with the generated input (fit of hSparseGen) and with the injected rho00
//      (true Psi: no resolution correction, rho00_obs must equal the injected value directly)
//
// Usage (one call per injected rho00, i.e. per (sub)wagon output):
//   root -l -b -q 'MCClosureTrueEP.C+("AnalysisResults.root", 0.333333, "rho0333")'
//   root -l -b -q 'MCClosureTrueEP.C+("AnalysisResults.root", 0.25,     "rho0250")'
//   root -l -b -q 'MCClosureTrueEP.C+("AnalysisResults.root", 0.40,     "rho0400")'
//
// Output: MCClosureTrueEP_<tag>.root (all histograms and graphs) and MCClosureTrueEP_<tag>.pdf (fits, distributions, rho00 vs pT).
// Final result per centrality: graph <folder>/cent.../gRhoFinal (= gRhoMeas).

#include <TAxis.h>
#include <TCanvas.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

// ======================================================================================
// User settings: adapt to your wagon
// ======================================================================================
namespace closure_cfg
{
// top-level directory in AnalysisResults.root (check with: TFile f("AnalysisResults.root"); f.ls())
const char* TaskDir = "kstarpbpb";

const int HarmonicN = 2; // cfgEPNormalHarmonic used in the wagon (only stored in the output: no resolution correction for true Psi)

// analysis binning: edges must coincide with bin edges of the THnSparse axes
const std::vector<double> CentEdges = {10., 30., 50.};
const std::vector<double> PtEdgesKstar = {1., 2., 3., 4., 6., 10.};
const std::vector<double> PtEdgesPhi = {1., 2., 3., 4., 6.};

const bool DrawMassFits = true; // one PDF page per (cent, pT) with the invariant-mass fits in each cos(theta*) bin
} // namespace closure_cfg

struct SpeciesConfig {
  std::string name;     // "Kstar" or "Phi"
  std::string label;    // for plots
  bool isKstar;         // Breit-Wigner + pol2 (K*) or Voigtian + pol1 (phi)
  double normLo, normHi; // ME normalisation region
  double fitLo, fitHi;   // mass fit range
  double winLo, winHi;   // signal (bin counting) window
  double mass, width;    // PDG values (start / fixed)
  std::vector<double> ptEdges;
};

SpeciesConfig kstarConfig()
{
  return {"Kstar", "K*^{0}", true, 1.10, 1.30, 0.75, 1.10, 0.80, 1.00, 0.8955, 0.0473, closure_cfg::PtEdgesKstar};
}

SpeciesConfig phiConfig()
{
  return {"Phi", "#phi(1020)", false, 1.04, 1.08, 0.995, 1.06, 1.010, 1.030, 1.019461, 0.004249, closure_cfg::PtEdgesPhi};
}

// ======================================================================================
// helpers
// ======================================================================================
namespace closure_util
{

// navigate "a/b/c" through TDirectory and TList levels (HistogramRegistry folders can be either)
TObject* getObject(TDirectory* top, const std::string& path)
{
  TObject* cur = top;
  size_t start = 0;
  while (start <= path.size()) {
    size_t end = path.find('/', start);
    std::string token = path.substr(start, end == std::string::npos ? std::string::npos : end - start);
    if (!token.empty()) {
      if (auto* dir = dynamic_cast<TDirectory*>(cur)) {
        cur = dir->Get(token.c_str());
      } else if (auto* coll = dynamic_cast<TCollection*>(cur)) {
        cur = coll->FindObject(token.c_str());
      } else {
        return nullptr;
      }
      if (!cur) {
        return nullptr;
      }
    }
    if (end == std::string::npos) {
      break;
    }
    start = end + 1;
  }
  return cur;
}

void listTopLevel(TFile* f)
{
  printf("  top-level keys in %s:\n", f->GetName());
  TIter next(f->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(next())) {
    printf("    %s (%s)\n", key->GetName(), key->GetClassName());
  }
}

// [lo, hi] -> first and last axis bin, warns if the edges do not coincide with bin edges
std::pair<int, int> binRange(const TAxis* ax, double lo, double hi)
{
  const double eps = 1e-6;
  int b1 = ax->FindBin(lo + eps);
  int b2 = ax->FindBin(hi - eps);
  b1 = std::max(b1, 1);
  b2 = std::min(b2, ax->GetNbins());
  if (std::abs(ax->GetBinLowEdge(b1) - lo) > 1e-4 || std::abs(ax->GetBinUpEdge(b2) - hi) > 1e-4) {
    printf("  WARNING: range [%g, %g] on axis %s is used as [%g, %g] (bin edges)\n", lo, hi, ax->GetTitle(), ax->GetBinLowEdge(b1), ax->GetBinUpEdge(b2));
  }
  return {b1, b2};
}

// axes of the closure sparses: 0 mass, 1 pT, 2 cos(theta*), 3 centrality
TH1D* projectMass(THnSparse* hs, std::pair<int, int> ptBins, int cosBin, std::pair<int, int> centBins, const TString& name)
{
  hs->GetAxis(1)->SetRange(ptBins.first, ptBins.second);
  if (cosBin > 0) {
    hs->GetAxis(2)->SetRange(cosBin, cosBin);
  } else {
    hs->GetAxis(2)->SetRange(1, hs->GetAxis(2)->GetNbins());
  }
  hs->GetAxis(3)->SetRange(centBins.first, centBins.second);
  TH1D* h = hs->Projection(0, "E");
  h->SetName(name);
  h->SetTitle("");
  h->SetDirectory(nullptr);
  for (int i = 0; i < hs->GetNdimensions(); i++) {
    hs->GetAxis(i)->SetRange();
  }
  return h;
}

double integralWindow(const TH1* h, double lo, double hi, double& err)
{
  const int b1 = h->GetXaxis()->FindBin(lo + 1e-7);
  const int b2 = h->GetXaxis()->FindBin(hi - 1e-7);
  return h->IntegralAndError(b1, b2, err);
}

// rho00 fit of a cos(theta*) distribution (integral of the function over each bin)
TF1* fitRho(TH1D* h, const TString& name, double& rho, double& rhoErr)
{
  rho = -1;
  rhoErr = 0;
  TF1* f = new TF1(name, "[0]*((1.-[1])+(3.*[1]-1.)*x*x)", -1., 1.);
  f->SetParNames("norm", "rho00");
  double sum = 0;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    sum += h->GetBinContent(i);
  }
  if (sum <= 0) {
    return f;
  }
  f->SetParameters(sum / h->GetNbinsX(), 1. / 3.);
  TFitResultPtr r = h->Fit(f, "QRNSI0");
  if (static_cast<int>(r) == 0 || static_cast<int>(r) == 4000) {
    rho = f->GetParameter(1);
    rhoErr = f->GetParError(1);
  }
  return f;
}

void styleGraph(TGraphErrors* g, int color, int marker)
{
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.2);
}

} // namespace closure_util

// ======================================================================================
// signal extraction in one (cent, pT, cos) bin
// ======================================================================================
struct SignalResult {
  double yield = 0, yieldErr = 0;
  bool fitOk = false;
};

SignalResult extractSignal(TH1D* hSE, TH1D* hME, const SpeciesConfig& s, TH1D*& hSub, TF1*& fTot, TF1*& fBkg, const TString& tag)
{
  SignalResult res;
  double e;
  const double seNorm = closure_util::integralWindow(hSE, s.normLo, s.normHi, e);
  const double meNorm = closure_util::integralWindow(hME, s.normLo, s.normHi, e);
  hSub = static_cast<TH1D*>(hSE->Clone("hSub_" + tag));
  hSub->SetDirectory(nullptr);
  if (meNorm > 0) {
    hSub->Add(hME, -seNorm / meNorm);
  }

  const double bw = hSub->GetBinWidth(1);
  double winErr;
  const double winCounts = closure_util::integralWindow(hSub, s.winLo, s.winHi, winErr);
  if (s.isKstar) {
    fTot = new TF1("fTot_" + tag, "[0]*[2]/(2.*TMath::Pi())/((x-[1])*(x-[1])+[2]*[2]/4.)+[3]+[4]*x+[5]*x*x", s.fitLo, s.fitHi);
    fTot->SetParameters(std::max(winCounts, 1.) * bw, s.mass, s.width, 0., 0., 0.);
    fTot->SetParLimits(1, s.mass - 0.03, s.mass + 0.03);
    fTot->SetParLimits(2, 0.02, 0.10);
    fBkg = new TF1("fBkg_" + tag, "[0]+[1]*x+[2]*x*x", s.fitLo, s.fitHi);
  } else {
    fTot = new TF1("fTot_" + tag, "[0]*TMath::Voigt(x-[1],[2],[3])+[4]+[5]*x", s.fitLo, s.fitHi);
    fTot->SetParameters(std::max(winCounts, 1.) * bw, s.mass, 0.001, s.width, 0., 0.);
    fTot->SetParLimits(1, s.mass - 0.005, s.mass + 0.005);
    fTot->SetParLimits(2, 0.0002, 0.005);
    fTot->FixParameter(3, s.width);
    fBkg = new TF1("fBkg_" + tag, "[0]+[1]*x", s.fitLo, s.fitHi);
  }
  fTot->SetLineColor(kRed);
  fBkg->SetLineColor(kBlue);
  fBkg->SetLineStyle(2);

  TFitResultPtr r = hSub->Fit(fTot, "QRNS0");
  res.fitOk = (static_cast<int>(r) == 0 || static_cast<int>(r) == 4000);
  const int nSigPar = s.isKstar ? 3 : 4;
  for (int i = 0; i < fBkg->GetNpar(); i++) {
    fBkg->SetParameter(i, res.fitOk ? fTot->GetParameter(nSigPar + i) : 0.);
  }

  // bin counting in the window, residual background (fit) subtracted
  const int b1 = hSub->GetXaxis()->FindBin(s.winLo + 1e-7);
  const int b2 = hSub->GetXaxis()->FindBin(s.winHi - 1e-7);
  double y = 0, e2 = 0;
  for (int i = b1; i <= b2; i++) {
    const double lo = hSub->GetXaxis()->GetBinLowEdge(i);
    const double hi = hSub->GetXaxis()->GetBinUpEdge(i);
    y += hSub->GetBinContent(i) - fBkg->Integral(lo, hi) / (hi - lo);
    e2 += hSub->GetBinError(i) * hSub->GetBinError(i);
  }
  res.yield = y;
  res.yieldErr = std::sqrt(e2);
  return res;
}

// ======================================================================================
// one closure configuration (species x true / reconstructed Psi)
// ======================================================================================
void runConfig(TFile* fin, const std::string& taskDir, const std::string& folder, const SpeciesConfig& s, bool recoEP, double rhoInj, TFile* fout, const TString& pdf)
{
  using namespace closure_util;
  const std::string base = taskDir + "/" + folder + "/";
  auto get = [&](const char* name) { return getObject(fin, base + name); };
  auto* hsSE = dynamic_cast<THnSparse*>(get("hSparseSE"));
  auto* hsME = dynamic_cast<THnSparse*>(get("hSparseME"));
  auto* hsRecTrue = dynamic_cast<THnSparse*>(get("hSparseRecTrue"));
  auto* hsGen = dynamic_cast<THnSparse*>(get("hSparseGen"));
  auto* hsRecTrueEff = dynamic_cast<THnSparse*>(get("hSparseRecTrueEff"));
  auto* hsGenEff = dynamic_cast<THnSparse*>(get("hSparseGenEff"));
  if (!hsSE || !hsME || !hsRecTrue || !hsGen || !hsRecTrueEff || !hsGenEff) {
    printf("\n=== %s: not found (or incomplete) in %s, skipped\n", base.c_str(), fin->GetName());
    return;
  }
  printf("\n=== %s (%s, %s), injected rho00 = %.4f\n", base.c_str(), s.label.c_str(), recoEP ? "reconstructed Psi_FT0C" : "MC true Psi", rhoInj);
  if (auto* hMC = dynamic_cast<TH1*>(get("hMC"))) {
    printf("  MC events: all %.0f, no reco %.0f, >1 reco %.0f, selected eff sample %.0f, selected test sample %.0f\n",
           hMC->GetBinContent(1), hMC->GetBinContent(2), hMC->GetBinContent(3), hMC->GetBinContent(4), hMC->GetBinContent(5));
  }

  TDirectory* outDir = fout->mkdir(folder.c_str());
  const TAxis* ptAxis = hsSE->GetAxis(1);
  const TAxis* cosAxis = hsSE->GetAxis(2);
  const TAxis* centAxis = hsSE->GetAxis(3);
  const int nCos = cosAxis->GetNbins();
  const int nPt = static_cast<int>(s.ptEdges.size()) - 1;
  const double rhoUnpol = 1. / 3.;

  for (size_t ic = 0; ic + 1 < closure_cfg::CentEdges.size(); ic++) {
    const std::pair<double, double> cent{closure_cfg::CentEdges[ic], closure_cfg::CentEdges[ic + 1]};
    const auto centBins = binRange(centAxis, cent.first, cent.second);
    const TString centTag = Form("cent%.0f_%.0f", cent.first, cent.second);
    TDirectory* centDir = outDir->mkdir(centTag);

    const double r2 = std::nan(""), rCorr = std::nan(""), r2True = std::nan(""), rCorrTrue = std::nan("");

    auto* gGen = new TGraphErrors();
    auto* gIdeal = new TGraphErrors();
    auto* gMeas = new TGraphErrors();
    auto* gCorr = new TGraphErrors();
    auto* gCorrTrueR = new TGraphErrors();
    auto* gYieldRatio = new TGraphErrors();
    gGen->SetName("gRhoGen");
    gGen->SetTitle("generated (input, true axis)");
    gIdeal->SetName("gRhoIdeal");
    gIdeal->SetTitle("true reco / Acc#timesEff (no signal extraction)");
    gMeas->SetName("gRhoMeas");
    gMeas->SetTitle("(SE#minusME) / Acc#timesEff");
    gCorr->SetName("gRhoCorr");
    gCorr->SetTitle("(SE#minusME) / Acc#timesEff, EP-resolution corrected");
    gCorrTrueR->SetName("gRhoCorrTrueR");
    gCorrTrueR->SetTitle("as gRhoCorr, with R from the MC true plane");
    gYieldRatio->SetName("gRawOverRecTrue");
    gYieldRatio->SetTitle("raw yield / true reco yield (cos-integrated)");

    printf("  %-10s %-10s | %-17s %-17s %-17s %-17s | raw/true\n", centTag.Data(), "pT", "rho gen", "rho ideal", "rho meas", recoEP ? "rho corrected" : "-");

    for (int ip = 0; ip < nPt; ip++) {
      const double ptLo = s.ptEdges[ip], ptHi = s.ptEdges[ip + 1];
      const auto ptBins = binRange(ptAxis, ptLo, ptHi);
      const TString tag = Form("%s_%s_%s_pt%.1f_%.1f", s.name.c_str(), recoEP ? "reco" : "true", centTag.Data(), ptLo, ptHi);

      TH1D* hRaw = new TH1D("hRaw_" + tag, ";cos(#vartheta*);raw yield", nCos, cosAxis->GetXmin(), cosAxis->GetXmax());
      TH1D* hRecTrue = static_cast<TH1D*>(hRaw->Clone("hRecTrue_" + tag));
      TH1D* hRecTrueEff = static_cast<TH1D*>(hRaw->Clone("hRecTrueEff_" + tag));
      TH1D* hGenEff = static_cast<TH1D*>(hRaw->Clone("hGenEff_" + tag));
      TH1D* hGen = static_cast<TH1D*>(hRaw->Clone("hGen_" + tag));
      for (auto* h : {hRaw, hRecTrue, hRecTrueEff, hGenEff, hGen}) {
        h->SetDirectory(nullptr);
        h->Sumw2();
      }

      TCanvas* cFit = nullptr;
      if (closure_cfg::DrawMassFits) {
        cFit = new TCanvas("cFit_" + tag, "", 1400, 900);
        const int nx = static_cast<int>(std::ceil(std::sqrt(nCos)));
        cFit->Divide(nx, static_cast<int>(std::ceil(static_cast<double>(nCos) / nx)));
      }

      int nFailedFits = 0;
      for (int ib = 1; ib <= nCos; ib++) {
        const TString btag = tag + Form("_cos%d", ib);
        TH1D* se = projectMass(hsSE, ptBins, ib, centBins, "se_" + btag);
        TH1D* me = projectMass(hsME, ptBins, ib, centBins, "me_" + btag);
        TH1D* sub = nullptr;
        TF1 *fTot = nullptr, *fBkg = nullptr;
        const SignalResult sig = extractSignal(se, me, s, sub, fTot, fBkg, btag);
        nFailedFits += sig.fitOk ? 0 : 1;
        hRaw->SetBinContent(ib, sig.yield);
        hRaw->SetBinError(ib, sig.yieldErr);

        double e;
        TH1D* rt = projectMass(hsRecTrue, ptBins, ib, centBins, "rt_" + btag);
        hRecTrue->SetBinContent(ib, integralWindow(rt, s.winLo, s.winHi, e));
        hRecTrue->SetBinError(ib, e);
        TH1D* rte = projectMass(hsRecTrueEff, ptBins, ib, centBins, "rte_" + btag);
        hRecTrueEff->SetBinContent(ib, integralWindow(rte, s.winLo, s.winHi, e));
        hRecTrueEff->SetBinError(ib, e);
        TH1D* ge = projectMass(hsGenEff, ptBins, ib, centBins, "ge_" + btag);
        hGenEff->SetBinContent(ib, ge->IntegralAndError(1, ge->GetNbinsX(), e));
        hGenEff->SetBinError(ib, e);
        TH1D* g = projectMass(hsGen, ptBins, ib, centBins, "g_" + btag);
        hGen->SetBinContent(ib, g->IntegralAndError(1, g->GetNbinsX(), e));
        hGen->SetBinError(ib, e);

        if (cFit) {
          cFit->cd(ib);
          sub->GetXaxis()->SetRangeUser(s.fitLo, s.fitHi);
          sub->SetMarkerStyle(20);
          sub->SetMarkerSize(0.5);
          sub->SetTitle(Form("cos#vartheta* [%.2f, %.2f]%s", cosAxis->GetBinLowEdge(ib), cosAxis->GetBinUpEdge(ib), sig.fitOk ? "" : " FIT FAILED"));
          sub->DrawCopy("E");
          fTot->DrawCopy("same");
          fBkg->DrawCopy("same");
          TLine l;
          l.SetLineStyle(3);
          l.DrawLine(s.winLo, sub->GetMinimum(), s.winLo, sub->GetMaximum());
          l.DrawLine(s.winHi, sub->GetMinimum(), s.winHi, sub->GetMaximum());
        }
        for (TObject* o : std::initializer_list<TObject*>{se, me, sub, fTot, fBkg, rt, rte, ge, g}) {
          delete o;
        }
      }
      if (nFailedFits > 0) {
        printf("  WARNING: %s: %d mass fits failed (no residual background subtracted there), check the PDF\n", tag.Data(), nFailedFits);
      }
      if (cFit) {
        cFit->cd(0);
        TLatex t;
        t.SetNDC();
        t.SetTextSize(0.02);
        t.DrawLatex(0.01, 0.985, Form("%s %s, %s, %.1f < p_{T} < %.1f GeV/c: SE #minus norm. ME (red: fit, blue: residual bkg, dotted: window)",
                                      s.label.c_str(), recoEP ? "reco #Psi" : "true #Psi", centTag.Data(), ptLo, ptHi));
        cFit->Print(pdf);
        delete cFit;
      }

      // Acc x Eff and corrected distributions
      TH1D* hEff = static_cast<TH1D*>(hRecTrueEff->Clone("hAccEff_" + tag));
      hEff->SetDirectory(nullptr);
      hEff->Divide(hRecTrueEff, hGenEff, 1., 1., "B");
      TH1D* hCorr = static_cast<TH1D*>(hRaw->Clone("hCorrected_" + tag));
      hCorr->SetDirectory(nullptr);
      hCorr->Divide(hEff);
      TH1D* hIdeal = static_cast<TH1D*>(hRecTrue->Clone("hRecTrueCorrected_" + tag));
      hIdeal->SetDirectory(nullptr);
      hIdeal->Divide(hEff);

      double rhoGen, eGen, rhoIdeal, eIdeal, rhoMeas, eMeas;
      TF1* fGen = fitRho(hGen, "fRhoGen_" + tag, rhoGen, eGen);
      TF1* fIdeal = fitRho(hIdeal, "fRhoIdeal_" + tag, rhoIdeal, eIdeal);
      TF1* fMeas = fitRho(hCorr, "fRhoMeas_" + tag, rhoMeas, eMeas);

      double rhoCorr = std::nan(""), eCorr = 0, rhoCorrT = std::nan(""), eCorrT = 0;
      if (recoEP && rhoMeas >= 0) {
        if (std::isfinite(rCorr)) {
          const double k = 4. / (1. + 3. * rCorr);
          rhoCorr = rhoUnpol + (rhoMeas - rhoUnpol) * k;
          eCorr = eMeas * k;
        }
        if (std::isfinite(rCorrTrue)) {
          const double k = 4. / (1. + 3. * rCorrTrue);
          rhoCorrT = rhoUnpol + (rhoMeas - rhoUnpol) * k;
          eCorrT = eMeas * k;
        }
      }

      double eRaw, eTrue;
      const double sRaw = hRaw->IntegralAndError(1, nCos, eRaw);
      const double sTrue = hRecTrue->IntegralAndError(1, nCos, eTrue);
      const double ratio = sTrue > 0 ? sRaw / sTrue : 0;
      const double ratioErr = (sTrue > 0 && sRaw != 0) ? ratio * std::sqrt(std::pow(eRaw / sRaw, 2) + std::pow(eTrue / sTrue, 2)) : 0;

      const double ptMid = 0.5 * (ptLo + ptHi), ptErr = 0.5 * (ptHi - ptLo);
      auto addPoint = [&](TGraphErrors* gr, double v, double ev, double shift) {
        if (std::isfinite(v) && v >= -1) {
          const int n = gr->GetN();
          gr->SetPoint(n, ptMid + shift * ptErr, v);
          gr->SetPointError(n, shift == 0 ? ptErr : 0, ev);
        }
      };
      addPoint(gGen, rhoGen, eGen, -0.3);
      addPoint(gIdeal, rhoIdeal, eIdeal, -0.1);
      addPoint(gMeas, rhoMeas, eMeas, 0.1);
      if (recoEP) {
        addPoint(gCorr, rhoCorr, eCorr, 0.3);
        addPoint(gCorrTrueR, rhoCorrT, eCorrT, 0.4);
      }
      addPoint(gYieldRatio, ratio, ratioErr, 0);

      printf("  %-10s %4.1f-%-5.1f | %.4f +- %.4f   %.4f +- %.4f   %.4f +- %.4f   ", "", ptLo, ptHi, rhoGen, eGen, rhoIdeal, eIdeal, rhoMeas, eMeas);
      if (recoEP) {
        printf("%.4f +- %.4f", rhoCorr, eCorr);
      } else {
        printf("%-17s", "");
      }
      printf(" | %.3f +- %.3f\n", ratio, ratioErr);

      // cos(theta*) distributions page
      {
        TCanvas c("cCos_" + tag, "", 1400, 500);
        c.Divide(3, 1);
        c.cd(1);
        TH1D* r1 = static_cast<TH1D*>(hRaw->DrawCopy("E"));
        r1->SetTitle("raw (SE#minusME, black) vs true reco (red)");
        r1->SetMinimum(0);
        r1->SetMarkerStyle(20);
        TH1D* r2c = static_cast<TH1D*>(hRecTrue->DrawCopy("E same"));
        r2c->SetLineColor(kRed);
        r2c->SetMarkerColor(kRed);
        r2c->SetMarkerStyle(24);
        c.cd(2);
        TH1D* e1 = static_cast<TH1D*>(hEff->DrawCopy("E"));
        e1->SetTitle("Acc #times Eff (efficiency sample)");
        e1->SetMinimum(0);
        c.cd(3);
        TH1D* hcN = static_cast<TH1D*>(hCorr->Clone("hcN"));
        TH1D* hgN = static_cast<TH1D*>(hGen->Clone("hgN"));
        if (hcN->Integral() > 0) {
          hcN->Scale(1. / hcN->Integral(), "width");
        }
        if (hgN->Integral() > 0) {
          hgN->Scale(1. / hgN->Integral(), "width");
        }
        hcN->SetTitle(Form("corrected (black, #rho_{00}=%.3f) vs generated (red, #rho_{00}=%.3f);cos(#vartheta*);1/N dN/dcos", rhoMeas, rhoGen));
        hcN->SetMarkerStyle(20);
        hcN->SetMinimum(0);
        hcN->SetMaximum(1.5 * std::max(hcN->GetMaximum(), hgN->GetMaximum()));
        hcN->DrawCopy("E");
        hgN->SetLineColor(kRed);
        hgN->SetMarkerColor(kRed);
        hgN->SetMarkerStyle(24);
        hgN->DrawCopy("E same");
        delete hcN;
        delete hgN;
        c.cd(0);
        TLatex t;
        t.SetNDC();
        t.SetTextSize(0.03);
        t.DrawLatex(0.01, 0.97, Form("%s %s, %s, %.1f < p_{T} < %.1f GeV/c, injected #rho_{00} = %.3f", s.label.c_str(), recoEP ? "reco #Psi" : "true #Psi",
                                     centTag.Data(), ptLo, ptHi, rhoInj));
        c.Print(pdf);
      }

      centDir->cd();
      for (TObject* o : std::initializer_list<TObject*>{hRaw, hRecTrue, hRecTrueEff, hGenEff, hGen, hEff, hCorr, hIdeal, fGen, fIdeal, fMeas}) {
        o->Write();
        delete o;
      }
    } // pT

    // rho00 vs pT page
    {
      TCanvas c("cRho_" + centTag, "", 900, 700);
      gPad->SetLeftMargin(0.13);
      TH1* frame = c.DrawFrame(s.ptEdges.front(), std::min(0.1, rhoInj - 0.15), s.ptEdges.back(), std::max(0.55, rhoInj + 0.15));
      frame->SetTitle(Form("%s %s, %s;p_{T} (GeV/c);#rho_{00}", s.label.c_str(), recoEP ? "reconstructed #Psi_{FT0C}" : "MC true #Psi", centTag.Data()));
      TLine l;
      l.SetLineStyle(2);
      l.SetLineColor(kGray + 2);
      l.DrawLine(s.ptEdges.front(), rhoInj, s.ptEdges.back(), rhoInj);
      l.SetLineStyle(3);
      l.DrawLine(s.ptEdges.front(), rhoUnpol, s.ptEdges.back(), rhoUnpol);
      styleGraph(gGen, kBlack, 24);
      styleGraph(gIdeal, kGreen + 2, 25);
      styleGraph(gMeas, recoEP ? kOrange + 7 : kRed, 20);
      styleGraph(gCorr, kRed, 21);
      styleGraph(gCorrTrueR, kMagenta + 1, 26);
      TLegend leg(0.15, 0.70, 0.88, 0.89);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.028);
      for (TGraphErrors* gr : {gGen, gIdeal, gMeas}) {
        if (gr->GetN()) {
          gr->Draw("P same");
          leg.AddEntry(gr, gr->GetTitle(), "p");
        }
      }
      if (recoEP) {
        for (TGraphErrors* gr : {gCorr, gCorrTrueR}) {
          if (gr->GetN()) {
            gr->Draw("P same");
            leg.AddEntry(gr, gr->GetTitle(), "p");
          }
        }
      }
      leg.AddEntry((TObject*)nullptr, Form("dashed: injected #rho_{00} = %.3f, dotted: 1/3", rhoInj), "");
      leg.Draw();
      c.Print(pdf);
    }

    centDir->cd();
    // final result of this configuration: corrected for reconstructed Psi, measured for true Psi
    TGraphErrors* gFinal = static_cast<TGraphErrors*>((recoEP ? gCorr : gMeas)->Clone("gRhoFinal"));
    for (TGraphErrors* gr : {gGen, gIdeal, gMeas, gCorr, gCorrTrueR, gYieldRatio, gFinal}) {
      gr->Write();
    }
    TH1D hInfo("hInfo", "injected rho00, n, R2(sub), R(corr), R2(true), R(corr, true)", 6, 0, 6);
    hInfo.SetBinContent(1, rhoInj);
    hInfo.SetBinContent(2, closure_cfg::HarmonicN);
    hInfo.SetBinContent(3, r2);
    hInfo.SetBinContent(4, rCorr);
    hInfo.SetBinContent(5, r2True);
    hInfo.SetBinContent(6, rCorrTrue);
    hInfo.Write();
  } // cent
}

// ======================================================================================
void MCClosureTrueEP(const char* inFile = "AnalysisResults.root", double rhoInj = 1. / 3., const char* tag = "rho0333")
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TH1::AddDirectory(false);

  TFile* fin = TFile::Open(inFile);
  if (!fin || fin->IsZombie()) {
    printf("ERROR: cannot open %s\n", inFile);
    return;
  }
  if (!fin->Get(closure_cfg::TaskDir)) {
    printf("ERROR: directory '%s' not found in %s, set TaskDir at the top of the macro\n", closure_cfg::TaskDir, inFile);
    closure_util::listTopLevel(fin);
    return;
  }

  TFile* fout = TFile::Open(Form("MCClosureTrueEP_%s.root", tag), "RECREATE");
  const TString pdf = Form("MCClosureTrueEP_%s.pdf", tag);
  TCanvas cOpen("cOpen", "", 900, 700);
  cOpen.Print(pdf + "[");

  runConfig(fin, closure_cfg::TaskDir, "closureKstar", kstarConfig(), false, rhoInj, fout, pdf);
  runConfig(fin, closure_cfg::TaskDir, "closurePhi", phiConfig(), false, rhoInj, fout, pdf);

  cOpen.Print(pdf + "]");
  fout->Close();
  fin->Close();
  printf("\nOutput: MCClosureTrueEP_%s.root, %s\n", tag, pdf.Data());
}
