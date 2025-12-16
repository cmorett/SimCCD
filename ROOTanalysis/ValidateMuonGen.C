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

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
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
    candidate = "OUTPUT.root";
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

void RunValidateMuonGen(const std::string& filename) {
  if (filename.empty()) {
    std::cerr << "[ValidateMuonGen] Error: empty file name resolved.\n";
    return;
  }

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

  // 1. cos(thetaPri)
  canvas.Clear();
  TH1D hCostheta("hCostheta",
                 "Muon angular distribution;cos(#theta_{pri});Events",
                 60, 0.0, 1.0);
  tree->Draw("cos(thetaPri)>>hCostheta", "", "goff");
  hCostheta.SetLineColor(kAzure + 2);
  hCostheta.SetLineWidth(2);
  canvas.cd();
  hCostheta.Draw("hist");
  canvas.SaveAs("validate_costheta.pdf");

  // 2. EevtPri (energy)
  const double eMinRaw = (nEntries > 0) ? tree->GetMinimum("EevtPri") : 0.0;
  const double eMaxRaw = (nEntries > 0) ? tree->GetMaximum("EevtPri") : 0.0;
  bool assumeMeV = std::isfinite(eMaxRaw) ? (eMaxRaw > 200.0) : false;
  const double energyScale = assumeMeV ? 1.0 / 1000.0 : 1.0;
  double eMin = std::isfinite(eMinRaw) ? eMinRaw * energyScale : 0.0;
  double eMax = std::isfinite(eMaxRaw) ? eMaxRaw * energyScale : 20.0;
  if (eMin < 0.0) {
    eMin = 0.0;
  }
  if (eMax <= eMin) {
    eMax = eMin + 1.0;
  }
  const double eSpan = eMax - eMin;
  double eLow = std::max(0.0, eMin - 0.1 * eSpan);
  double eHigh = eMax + 0.1 * eSpan;
  if (eHigh <= eLow) {
    eLow = 0.0;
    eHigh = std::max(20.0, eHigh + 1.0);
  }
  const char* energyAxisLabel = "E_{#mu} [GeV]";
  TString energyExpr = assumeMeV ? "EevtPri/1000.0" : "EevtPri";

  canvas.Clear();
  TH1D hEnergy("hEnergy",
               Form("Primary muon energy;%s;Events", energyAxisLabel),
               80, eLow, eHigh);
  tree->Draw(energyExpr + ">>hEnergy", "", "goff");
  hEnergy.SetLineColor(kOrange + 7);
  hEnergy.SetLineWidth(2);
  canvas.cd();
  hEnergy.Draw("hist");
  canvas.SaveAs("validate_energy.pdf");

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

  // 4. muonZ0
  auto [zLow, zHigh] = ComputeRange(tree, "muonZ0", nEntries, 10.0);
  canvas.Clear();
  TH1D hZ("hZ", "Muon source height;z0 [cm];Events",
          80, zLow, zHigh);
  tree->Draw("muonZ0>>hZ", "", "goff");
  hZ.SetLineColor(kGreen + 2);
  hZ.SetLineWidth(2);
  canvas.cd();
  hZ.Draw("hist");
  canvas.SaveAs("validate_z0.pdf");

  file->Close();
  delete file;
  std::cout << "[ValidateMuonGen] Validation plots saved: "
               "validate_costheta.pdf, validate_energy.pdf, "
               "validate_xy0.pdf, validate_z0.pdf\n";
}

void ValidateMuonGen(const char* filename) {
  RunValidateMuonGen(ResolveFilename(filename));
}

void ValidateMuonGen() {
  RunValidateMuonGen(ResolveFilename(nullptr));
}
