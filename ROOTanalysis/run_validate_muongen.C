#include <string>

void run_validate_muongen(const char* filename = "build/vs2022/main/Release/B02ntuples.root")
{
  gROOT->ProcessLine(".L ROOTanalysis/ValidateMuonGen.C");
  const std::string cmd = std::string("ValidateMuonGen(\"") + filename + "\")";
  gROOT->ProcessLine(cmd.c_str());
}
