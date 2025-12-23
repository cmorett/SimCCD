#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <cstdio>
#include <cmath>

// Basic sanity checks for B02ntuples.root produced by the validation run.
void check_ntuple_sanity(const char* filename = "build/vs2022/main/Release/B02ntuples.root")
{
  TFile file(filename);
  if (file.IsZombie()) {
    printf("ERROR: could not open %s\n", filename);
    return;
  }

  auto* t = dynamic_cast<TTree*>(file.Get("B02Evts"));
  if (!t) {
    printf("ERROR: TTree B02Evts not found in %s\n", filename);
    return;
  }

  const char* required[] = {
    "thetaPri", "phiPri", "EevtPri", "muonX0", "muonY0", "muonZ0",
    "muonXImp", "muonYImp", "muonZImp",
    "EdepCCD", "nStepsCCD", "xEntryCCD", "yEntryCCD", "zEntryCCD",
    "xExitCCD", "yExitCCD", "zExitCCD", "trackLenCCD", "dirX", "dirY", "dirZ"
  };
  for (auto name : required) {
    printf("Branch %-8s : %s\n", name, t->GetBranch(name) ? "OK" : "MISSING");
  }

  const Long64_t n = t->GetEntries();
  printf("Entries: %lld\n", n);
  if (n == 0) {
    return;
  }

  auto printMinMax = [&](const char* br) {
    const double mn = t->GetMinimum(br);
    const double mx = t->GetMaximum(br);
    printf("  %-8s min=%g max=%g\n", br, mn, mx);
  };

  printMinMax("EevtPri");
  printMinMax("thetaPri");
  printMinMax("muonZ0");
  printMinMax("muonXImp");
  printMinMax("muonYImp");
  printMinMax("muonZImp");
  printMinMax("EdepCCD");
  printMinMax("trackLenCCD");

  // Clamp stats for downward zenith cosine.
  double theta = 0.0;
  t->SetBranchAddress("thetaPri", &theta);
  Long64_t clampLow = 0, clampHigh = 0;
  double z0 = 0.0;
  t->SetBranchAddress("muonZ0", &z0);
  long double zSum = 0.0L, zSum2 = 0.0L;

  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    double c = -std::cos(theta);
    if (!std::isfinite(c) || c < 0.0) {
      ++clampLow;
    } else if (c > 1.0) {
      ++clampHigh;
    }
    zSum += z0;
    zSum2 += static_cast<long double>(z0) * static_cast<long double>(z0);
  }

  const double zMean = static_cast<double>(zSum / n);
  const double zVar = static_cast<double>(zSum2 / n - zSum / n * zSum / n);
  printf("Clamp counts for -cos(thetaPri): low=%lld high=%lld total=%lld\n",
         clampLow, clampHigh, n);
  printf("muonZ0 mean=%g variance=%g (std=%g)\n", zMean, zVar, std::sqrt(std::max(0.0, zVar)));
}
