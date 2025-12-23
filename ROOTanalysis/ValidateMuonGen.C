#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TColor.h>
#include <TApplication.h>
#include <TSystem.h>

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace {
std::string ResolveFilename(const char* explicitArg) {
  if (explicitArg && std::strlen(explicitArg) > 0) {
    return explicitArg;
  }

  std::string candidate;
  TApplication* app = gROOT->GetApplication();
  if (app) {
    const int argc = app->Argc();
    char** argv = app->Argv();
    bool captureNext = false;
    for (int i = 1; i < argc; ++i) {
      if (!argv[i]) {
        continue;
      }
      if (captureNext) {
        if (std::strlen(argv[i]) > 0) {
          candidate = argv[i];
          break;
        }
      }
      if (std::strcmp(argv[i], "--args") == 0) {
        captureNext = true;
      }
    }
  }

  if (candidate.empty()) {
    candidate = "B02ntuples.root";
  }

  return candidate;
}

std::pair<double, double> ComputeRange(TTree* tree,
                                       const char* branchName,
                                       Long64_t nEntries,
                                       double fallbackHalfWidth = 5.0) {
  double minVal = (tree && nEntries > 0) ? tree->GetMinimum(branchName) : 0.0;
  double maxVal = (tree && nEntries > 0) ? tree->GetMaximum(branchName) : 0.0;
  if (!std::isfinite(minVal) || !std::isfinite(maxVal)) {
    minVal = -fallbackHalfWidth;
    maxVal = fallbackHalfWidth;
  }

  double span = maxVal - minVal;
  if (span <= 0.0) {
    minVal -= fallbackHalfWidth;
    maxVal += fallbackHalfWidth;
    span = maxVal - minVal;
  } else {
    const double pad = 0.1 * span;
    minVal -= pad;
    maxVal += pad;
  }

  if (maxVal <= minVal) {
    maxVal = minVal + 1.0;
  }

  return {minVal, maxVal};
}
}  // namespace

void RunValidateMuonGen(const std::string& filename, bool makeImpact) {
  if (filename.empty()) {
    std::cerr << "[ValidateMuonGen] Error: empty file name resolved.\n";
    return;
  }

  gROOT->SetBatch(true);

  std::cout << "[ValidateMuonGen] Opening file: " << filename << std::endl;
  TFile* file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[ValidateMuonGen] Error: could not open '" << filename << "'.\n";
    delete file;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(file->Get("B02Evts"));
  if (!tree) {
    std::cerr << "[ValidateMuonGen] Error: TTree 'B02Evts' not found in file.\n";
    file->Close();
    delete file;
    return;
  }

  std::vector<std::string> missingBranches;
  const char* requiredBranches[] = {"thetaPri", "EevtPri",
                                    "muonX0",   "muonY0",
                                    "muonZ0"};
  const auto branchChecker = [&](const char* branchName) {
    if (!tree->GetBranch(branchName)) {
      missingBranches.emplace_back(branchName);
      return false;
    }
    return true;
  };
  for (const char* branchName : requiredBranches) {
    branchChecker(branchName);
  }

  if (!missingBranches.empty()) {
    std::cerr << "[ValidateMuonGen] Missing required branches in B02Evts: ";
    for (size_t i = 0; i < missingBranches.size(); ++i) {
      std::cerr << missingBranches[i];
      if (i + 1 != missingBranches.size()) {
        std::cerr << ", ";
      }
    }
    std::cerr << std::endl;
    file->Close();
    delete file;
    return;
  }

  const Long64_t nEntries = tree->GetEntries();
  if (nEntries == 0) {
    std::cout << "[ValidateMuonGen] Warning: B02Evts tree has no entries.\n";
  }

  gStyle->SetOptStat(0);
  TCanvas canvas("cValidate", "Muon generator validation", 900, 700);
  canvas.SetGrid();

  // 1. Downward zenith cosine: -cos(thetaPri)
  canvas.Clear();
  TH1D hCostheta("hCostheta",
                 "Muon angular distribution;cos_{zenith}^{#downarrow};Events",
                 60, 0.0, 1.0);
  double thetaPri = 0.0;
  tree->SetBranchAddress("thetaPri", &thetaPri);
  Long64_t clampLow = 0;
  Long64_t clampHigh = 0;
  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);
    double c = -std::cos(thetaPri);
    if (!std::isfinite(c)) {
      ++clampLow;
      c = 0.0;
    } else if (c < 0.0) {
      ++clampLow;
      c = 0.0;
    } else if (c > 1.0) {
      ++clampHigh;
      c = 1.0;
    }
    hCostheta.Fill(c);
  }
  hCostheta.SetLineColor(kAzure + 2);
  hCostheta.SetLineWidth(2);
  canvas.cd();
  hCostheta.Draw("hist");
  canvas.SaveAs("validate_costheta.pdf");
  std::cout << "[ValidateMuonGen] costheta clamp: low=" << clampLow
            << " high=" << clampHigh << " total=" << nEntries << "\n";

  // 2. EevtPri (energy)
  const double eMinRaw = (nEntries > 0) ? tree->GetMinimum("EevtPri") : 0.0;
  const double eMaxRaw = (nEntries > 0) ? tree->GetMaximum("EevtPri") : 0.0;

  // The tree stores EevtPri in GeV; pad the observed range slightly.
  double eMin = std::isfinite(eMinRaw) ? eMinRaw : 0.0;
  double eMax = std::isfinite(eMaxRaw) ? eMaxRaw : 0.0;
  if (eMin < 0.0) eMin = 0.0;
  if (eMax <= eMin) eMax = eMin + 1.0;
  const double eSpan = eMax - eMin;
  const double eLow = std::max(0.0, eMin - 0.1 * eSpan);
  const double eHigh = eMax + 0.1 * eSpan;

  canvas.Clear();
  TH1D hEnergy("hEnergy",
               "Primary muon energy;E_{#mu} [GeV];Events",
               80, eLow, eHigh);
  double ePri = 0.0;
  tree->SetBranchAddress("EevtPri", &ePri);
  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);
    if (!std::isfinite(ePri)) {
      continue;
    }
    hEnergy.Fill(ePri);
  }
  hEnergy.SetLineColor(kOrange + 7);
  hEnergy.SetLineWidth(2);
  canvas.cd();
  hEnergy.Draw("hist");
  canvas.SaveAs("validate_energy.pdf");
  std::cout << "[ValidateMuonGen] energy underflow=" << hEnergy.GetBinContent(0)
            << " overflow=" << hEnergy.GetBinContent(hEnergy.GetNbinsX() + 1)
            << " entries=" << hEnergy.GetEntries() << "\n";

  // 3. muonX0 vs muonY0
  auto [xLow, xHigh] = ComputeRange(tree, "muonX0", nEntries, 5.0);
  auto [yLow, yHigh] = ComputeRange(tree, "muonY0", nEntries, 5.0);
  canvas.Clear();
  canvas.SetRightMargin(0.15);
  TH2D hXY("hXY", "Muon initial position;x0 [cm];y0 [cm]",
           100, xLow, xHigh,
           100, yLow, yHigh);
  tree->Draw("muonY0:muonX0>>hXY", "", "goff");
  canvas.cd();
  hXY.Draw("COLZ");
  canvas.SaveAs("validate_xy0.pdf");
  canvas.SetRightMargin(0.08);

  // Optional: impact-plane sampling (muonXImp, muonYImp).
  if (makeImpact) {
    const bool haveImpact = tree->GetBranch("muonXImp") && tree->GetBranch("muonYImp");
    if (!haveImpact) {
      std::cout << "[ValidateMuonGen] validate_xyImpact.pdf requested but impact branches "
                   "(muonXImp/muonYImp) are missing; skipping.\n";
    } else {
      auto [xImpLow, xImpHigh] = ComputeRange(tree, "muonXImp", nEntries, 1.0);
      auto [yImpLow, yImpHigh] = ComputeRange(tree, "muonYImp", nEntries, 1.0);
      canvas.Clear();
      canvas.SetRightMargin(0.15);
      TH2D hXYImp("hXYImp", "Muon impact-plane sampling;x_{imp} [cm];y_{imp} [cm]",
                  100, xImpLow, xImpHigh,
                  100, yImpLow, yImpHigh);
      tree->Draw("muonYImp:muonXImp>>hXYImp", "", "goff");
      canvas.cd();
      hXYImp.Draw("COLZ");
      canvas.SaveAs("validate_xyImpact.pdf");
      canvas.SetRightMargin(0.08);
    }
  }

  // 4. muonZ0
  const double zMinRaw = (nEntries > 0) ? tree->GetMinimum("muonZ0") : 0.0;
  const double zMaxRaw = (nEntries > 0) ? tree->GetMaximum("muonZ0") : 0.0;
  const bool zIsConstant = (nEntries > 0) && std::isfinite(zMinRaw) &&
                           std::isfinite(zMaxRaw) &&
                           std::abs(zMaxRaw - zMinRaw) < 1e-9;
  auto [zLow, zHigh] = ComputeRange(tree, "muonZ0", nEntries, 0.5);
  canvas.Clear();
  TString zTitle = zIsConstant
                     ? Form("Muon source height (constant z0 = %.3g cm);z0 [cm];Events",
                            zMinRaw)
                     : "Muon source height;z0 [cm];Events";
  TH1D hZ("hZ", zTitle, 80, zLow, zHigh);
  double z0 = 0.0;
  tree->SetBranchAddress("muonZ0", &z0);
  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);
    if (!std::isfinite(z0)) {
      continue;
    }
    hZ.Fill(z0);
  }
  hZ.SetLineColor(kGreen + 2);
  hZ.SetLineWidth(2);
  canvas.cd();
  hZ.Draw("hist");
  canvas.SaveAs("validate_z0.pdf");
  std::cout << "[ValidateMuonGen] z0 underflow=" << hZ.GetBinContent(0)
            << " overflow=" << hZ.GetBinContent(hZ.GetNbinsX() + 1)
            << " entries=" << hZ.GetEntries() << "\n";

  file->Close();
  delete file;
  std::cout << "[ValidateMuonGen] Validation plots saved: "
               "validate_costheta.pdf, validate_energy.pdf, "
               "validate_xy0.pdf, validate_z0.pdf";
  if (makeImpact) {
    std::cout << ", validate_xyImpact.pdf";
  }
  std::cout << "\n";
}

void ValidateMuonGen(const char* filename = "", bool makeImpact = false) {
  const bool envImpact = []() {
    const char* env = gSystem ? gSystem->Getenv("SIMCCD_VALIDATE_IMPACT") : nullptr;
    if (!env || std::strlen(env) == 0) return false;
    return std::atoi(env) != 0;
  }();
  RunValidateMuonGen(ResolveFilename(filename), makeImpact || envImpact);
}
