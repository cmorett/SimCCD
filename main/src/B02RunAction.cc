//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4ANALYSIS_USE
//#include "B02AnalysisManager.hh"
#endif

#include "B02RunAction.hh"

#include "B02BarSD.hh"
#include "B02DetectorConstruction.hh"
#include "B02EventAction.hh"
#include "B02PrimaryGeneratorAction.hh"
#include "Randomize.hh"

#include "G4AnalysisManager.hh"
#include "G4GenericMessenger.hh"
#include "G4ProductionCuts.hh"
#include "G4RegionStore.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4VUserPhysicsList.hh"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <typeinfo>

namespace {

std::string Trim(const std::string& s) {
  const auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  const auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

}  // namespace

B02RunAction::B02RunAction(B02EventAction* b02eventAction )
  : fB02EventAction(b02eventAction)
{
  fMessenger = new G4GenericMessenger(this, "/runAction/", "Run control");
  fMessenger->DeclareProperty("useTimeSeed", fUseTimeSeed,
                              "Seed the random engine from wall-clock time at run start.");
  fMessenger->DeclareProperty("useFixedSeeds", fUseFixedSeeds,
                              "Seed the random engine from seed1/seed2 at run start.");
  fMessenger->DeclareProperty("seed1", fSeed1, "First random seed used when useFixedSeeds=true.");
  fMessenger->DeclareProperty("seed2", fSeed2, "Second random seed used when useFixedSeeds=true.");
  fMessenger->DeclareProperty("macroPath", fMacroPath,
                              "Macro path (stored in run metadata).");
  fMessenger->DeclareProperty("macroHash", fMacroHash,
                              "Macro content hash (stored in run metadata).");
  fMessenger->DeclareProperty("provenanceTag", fProvenanceTag,
                              "Optional provenance tag recorded in run metadata.");

  BuildAnalysis();
}

B02RunAction::~B02RunAction()
{
  delete fMessenger;
}

void B02RunAction::SetMacroPath(const G4String& path) { fMacroPath = path; }
void B02RunAction::SetMacroHash(const G4String& hash) { fMacroHash = hash; }
void B02RunAction::SetProvenanceTag(const G4String& tag) { fProvenanceTag = tag; }

void B02RunAction::BuildAnalysis() {
  auto analysisManager = G4AnalysisManager::Instance();

  fRunInfoNtupleId = -1;

  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(G4Threading::IsMultithreadedApplication());
  analysisManager->SetFileName("B02");

  // Creating 1D histograms
  analysisManager->CreateH1("Epri","Primary muon kinetic energy (GeV)", 100, 0., 200.);
  analysisManager->CreateH1("ThAng","cos(#theta_{pri})", 100, -1.0, 1.0);
  analysisManager->CreateH1("Ebar", "Ebar", 100,0,700);
  analysisManager->CreateH1("XBar", "Xbar", 100, -1500, 1500);
  analysisManager->CreateH1("LDis" ,"Length traced", 100, 0., 500);
  analysisManager->CreateH1("MuonDecayEdepLAr"," Energy deposited by decay muons", 100, 0, 1000);

  if ( fB02EventAction ) {
    // Ntuple 0: hits
    analysisManager->CreateNtuple("B02Hits", "B02Hits");
    analysisManager->CreateNtupleIColumn("evtId");
    analysisManager->CreateNtupleDColumn("XhitBar");
    analysisManager->CreateNtupleDColumn("YhitBar");
    analysisManager->CreateNtupleDColumn("ZhitBar");
    analysisManager->CreateNtupleDColumn("RhitBar");
    analysisManager->CreateNtupleDColumn("EhitBar");
    analysisManager->CreateNtupleDColumn("WhitBar");
    analysisManager->FinishNtuple();

    // Ntuple 1: event summary (indices must match EventAction filling)
    analysisManager->CreateNtuple("B02Evts", "B02Evts");
    analysisManager->CreateNtupleIColumn("evtID");              // 0
    analysisManager->CreateNtupleDColumn("EevtBar");            // 1
    analysisManager->CreateNtupleDColumn("WevtBar");            // 2
    analysisManager->CreateNtupleDColumn("GevtBar");            // 3
    analysisManager->CreateNtupleDColumn("ADCevtBar");          // 4
    analysisManager->CreateNtupleDColumn("EevtPri");            // 5
    analysisManager->CreateNtupleDColumn("thetaPri");           // 6
    analysisManager->CreateNtupleDColumn("phiPri");             // 7
    analysisManager->CreateNtupleIColumn("nHitBar");            // 8
    analysisManager->CreateNtupleDColumn("LengthMuLAr");        // 9
    analysisManager->CreateNtupleDColumn("MuonDecayEdepLAr");   // 10
    analysisManager->CreateNtupleDColumn("muonX0");             // 11
    analysisManager->CreateNtupleDColumn("muonY0");             // 12
    analysisManager->CreateNtupleDColumn("muonZ0");             // 13
    analysisManager->CreateNtupleDColumn("muonXImp");           // 14
    analysisManager->CreateNtupleDColumn("muonYImp");           // 15
    analysisManager->CreateNtupleDColumn("muonZImp");           // 16
    analysisManager->CreateNtupleDColumn("EdepCCD");            // 17
    analysisManager->CreateNtupleDColumn("EdepOther");          // 18
    analysisManager->CreateNtupleIColumn("nStepsCCD");          // 19
    analysisManager->CreateNtupleDColumn("xEntryCCD");          // 20
    analysisManager->CreateNtupleDColumn("yEntryCCD");          // 21
    analysisManager->CreateNtupleDColumn("zEntryCCD");          // 22
    analysisManager->CreateNtupleDColumn("xExitCCD");           // 23
    analysisManager->CreateNtupleDColumn("yExitCCD");           // 24
    analysisManager->CreateNtupleDColumn("zExitCCD");           // 25
    analysisManager->CreateNtupleDColumn("trackLenCCD");        // 26
    analysisManager->CreateNtupleDColumn("dirX");               // 27
    analysisManager->CreateNtupleDColumn("dirY");               // 28
    analysisManager->CreateNtupleDColumn("dirZ");               // 29
    analysisManager->CreateNtupleDColumn("muonCosTheta");       // 30
    analysisManager->CreateNtupleDColumn("muonWeight_s");       // 31
    analysisManager->CreateNtupleDColumn("eventLivetime_s");    // 32
    analysisManager->CreateNtupleIColumn("muonModeCode");       // 33
    analysisManager->CreateNtupleDColumn("muonEnergySampledGeV");  // 34
    analysisManager->CreateNtupleIColumn("firstHitIsCCD");      // 35
    analysisManager->CreateNtupleDColumn("geomIntersectsCCD");  // 36
    analysisManager->CreateNtupleIColumn("fluxModelCode");      // 37
    analysisManager->CreateNtupleDColumn("cfg_sourcePlaneZ_cm");// 38
    analysisManager->CreateNtupleDColumn("cfg_sourcePlaneLx_cm");// 39
    analysisManager->CreateNtupleDColumn("cfg_sourcePlaneLy_cm");// 40
    analysisManager->CreateNtupleDColumn("cfg_thetaMax_deg");   // 41
    analysisManager->CreateNtupleDColumn("cfg_EminGeV_eff");    // 42
    analysisManager->CreateNtupleDColumn("cfg_EmaxGeV_eff");    // 43
    analysisManager->CreateNtupleIColumn("muonPDG");            // 44 (new)
    analysisManager->CreateNtupleIColumn("muonChargeSign");     // 45 (new)
    analysisManager->CreateNtupleIColumn("prov_seed1");         // 46
    analysisManager->CreateNtupleIColumn("prov_seed2");         // 47
    analysisManager->CreateNtupleIColumn("prov_useTimeSeed");   // 48
    analysisManager->CreateNtupleIColumn("prov_overburdenEnabled"); // 49
    analysisManager->CreateNtupleDColumn("prov_overburdenThickness_cm"); // 50
    analysisManager->CreateNtupleDColumn("prov_overburdenZTop_cm"); // 51
    analysisManager->CreateNtupleDColumn("prov_overburdenMaterialHash"); // 52
    analysisManager->CreateNtupleDColumn("prov_ccdGammaCut_cm"); // 53
    analysisManager->CreateNtupleDColumn("prov_ccdElectronCut_cm"); // 54
    analysisManager->CreateNtupleDColumn("prov_ccdPositronCut_cm"); // 55
    analysisManager->CreateNtupleDColumn("prov_ccdMaxStep_cm"); // 56
    analysisManager->CreateNtupleDColumn("prov_ccdThickness_cm"); // 57
    analysisManager->CreateNtupleDColumn("prov_gitHashCode");   // 58
    analysisManager->CreateNtupleDColumn("prov_macroHashCode"); // 59
    analysisManager->CreateNtupleDColumn("prov_macroPathHash"); // 60
    analysisManager->CreateNtupleDColumn("prov_physicsListHash"); // 61
    analysisManager->CreateNtupleDColumn("prov_muonChargeRatio"); // 62
    analysisManager->CreateNtupleIColumn("isTargeted");        // 63
    analysisManager->FinishNtuple();

    // Ntuple 2: run-level provenance
    fRunInfoNtupleId =
        analysisManager->CreateNtuple("B02RunInfo", "Run-level provenance and configuration");
    analysisManager->CreateNtupleSColumn("gitHash");               // 0
    analysisManager->CreateNtupleIColumn("gitDirty");              // 1
    analysisManager->CreateNtupleSColumn("macroPath");             // 2
    analysisManager->CreateNtupleSColumn("macroHash");             // 3
    analysisManager->CreateNtupleSColumn("provenanceTag");         // 4
    analysisManager->CreateNtupleSColumn("physicsList");           // 5
    analysisManager->CreateNtupleIColumn("useTimeSeed");           // 6
    analysisManager->CreateNtupleIColumn("seed1");                 // 7
    analysisManager->CreateNtupleIColumn("seed2");                 // 8
    analysisManager->CreateNtupleSColumn("muonMode");              // 9
    analysisManager->CreateNtupleSColumn("fluxModel");             //10
    analysisManager->CreateNtupleSColumn("muonChargeMode");        //11
    analysisManager->CreateNtupleDColumn("muonChargeRatio");       //12
    analysisManager->CreateNtupleDColumn("cfg_sourcePlaneZ_cm");   //13
    analysisManager->CreateNtupleDColumn("cfg_sourcePlaneLx_cm");  //14
    analysisManager->CreateNtupleDColumn("cfg_sourcePlaneLy_cm");  //15
    analysisManager->CreateNtupleDColumn("cfg_thetaMax_deg");      //16
    analysisManager->CreateNtupleDColumn("cfg_EminGeV_eff");       //17
    analysisManager->CreateNtupleDColumn("cfg_EmaxGeV_eff");       //18
    analysisManager->CreateNtupleIColumn("overburdenEnabled");     //19
    analysisManager->CreateNtupleDColumn("overburdenThickness_cm");//20
    analysisManager->CreateNtupleDColumn("overburdenZTop_cm");     //21
    analysisManager->CreateNtupleSColumn("overburdenMaterial");    //22
    analysisManager->CreateNtupleDColumn("ccdGammaCut_cm");        //23
    analysisManager->CreateNtupleDColumn("ccdElectronCut_cm");     //24
    analysisManager->CreateNtupleDColumn("ccdPositronCut_cm");     //25
    analysisManager->CreateNtupleDColumn("ccdMaxStep_cm");         //26
    analysisManager->CreateNtupleDColumn("ccdThickness_cm");       //27
    analysisManager->FinishNtuple();

    analysisManager->SetNtupleFileName(0, "B02ntuples");
    analysisManager->SetNtupleFileName(1, "B02ntuples");
    if (fRunInfoNtupleId >= 0) {
      analysisManager->SetNtupleFileName(fRunInfoNtupleId, "B02ntuples");
    }
  }
}

void B02RunAction::BeginOfRunAction(const G4Run* aRun)
{
  long seeds[2] = {0L, 0L};
  bool seedsSet = false;
  if (fUseFixedSeeds) {
    seeds[0] = static_cast<long>(fSeed1);
    seeds[1] = static_cast<long>(fSeed2);
    G4Random::setTheSeeds(seeds);
    seedsSet = true;
  } else if (fUseTimeSeed) {
    time_t systime = time(nullptr);
    seeds[0] = static_cast<long>(systime);
    seeds[1] = static_cast<long>(systime * G4UniformRand());
    G4Random::setTheSeeds(seeds);
    seedsSet = true;
    fSeed1 = seeds[0];
    fSeed2 = seeds[1];
  }

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  if (seedsSet) {
    G4cout << "### Random seeds: " << seeds[0] << " " << seeds[1] << G4endl;
  } else {
    const long* current = G4Random::getTheSeeds();
    if (current) {
      G4cout << "### Random seeds (unchanged): " << current[0] << " " << current[1] << G4endl;
    }
  }
  if (const auto* detector = dynamic_cast<const B02DetectorConstruction*>(
          G4RunManager::GetRunManager()->GetUserDetectorConstruction())) {
    if (detector->GetCCDPrintInfo()) {
      detector->PrintCCDInfo();
    }
  }

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Reset histograms from previous run
  analysisManager->Reset();

  // Open an output file
  analysisManager->OpenFile();

  FillRunInfoNtuple(aRun);
}

void B02RunAction::EndOfRunAction(const G4Run* /*run*/)
{
   auto analysisManager = G4AnalysisManager::Instance();

   if(analysisManager->GetH1(2)){
     G4cout<<"Histo: EpriHis (GeV): mean: "<<analysisManager->GetH1(0)->mean()<< " rms: "<<analysisManager->GetH1(0)->rms()<<G4endl;
     G4cout<<"Histo: cos(theta): mean: "<<analysisManager->GetH1(1)->mean()<< " rms: "<<analysisManager->GetH1(1)->rms()<<G4endl;
     G4cout<<"Histo: EhitBar (MeV): mean: "<<analysisManager->GetH1(2)->mean()<< " rms: "<<analysisManager->GetH1(2)->rms()<<G4endl;
     G4cout<<"Histo: XhitBar (cm): mean: "<<analysisManager->GetH1(3)->mean()<< " rms: "<<analysisManager->GetH1(3)->rms()<<G4endl;
     G4cout<<"Histo: LDis (cm): mean: "<<analysisManager->GetH1(4)->mean()<< " rms: "<<analysisManager->GetH1(4)->rms()<<G4endl;
     G4cout<<"Histo: Muon decay energy deposited (MeV): mean: "<<analysisManager->GetH1(5)->mean()<< " rms: "<<analysisManager->GetH1(5)->rms()<<G4endl;
   }

  analysisManager->Write();
  analysisManager->CloseFile(false);
}

std::string B02RunAction::RunCommand(const std::string& cmd) const {
#if defined(_WIN32)
  FILE* pipe = _popen(cmd.c_str(), "r");
#else
  FILE* pipe = popen(cmd.c_str(), "r");
#endif
  if (!pipe) {
    return {};
  }
  std::array<char, 256> buffer{};
  std::string result;
  while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe)) {
    result += buffer.data();
  }
#if defined(_WIN32)
  _pclose(pipe);
#else
  pclose(pipe);
#endif
  return Trim(result);
}

std::string B02RunAction::DetectGitHash() {
  if (!fGitHash.empty()) {
    return fGitHash;
  }
  const std::string hash = RunCommand("git rev-parse --short HEAD");
  if (!hash.empty()) {
    fGitHash = hash;
  } else {
    fGitHash = "unknown";
  }
  return fGitHash;
}

bool B02RunAction::DetectGitDirty() {
  fGitDirty = false;
  const std::string status = RunCommand("git status --porcelain");
  if (!status.empty()) {
    fGitDirty = true;
  }
  return fGitDirty;
}

std::string B02RunAction::ComputeFileHash(const std::string& path) const {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return {};
  }
  std::ostringstream ss;
  ss << in.rdbuf();
  const std::string contents = ss.str();
  const auto h = std::hash<std::string>{}(contents);
  std::ostringstream hex;
  hex << std::hex << h;
  return hex.str();
}

void B02RunAction::FillRunInfoNtuple(const G4Run* run) {
  auto analysisManager = G4AnalysisManager::Instance();
  if (!analysisManager || fRunInfoNtupleId < 0) {
    return;
  }
  const auto ntupleId = fRunInfoNtupleId;

  DetectGitHash();
  DetectGitDirty();
  if (fMacroHash.empty() && !fMacroPath.empty()) {
    fMacroHash = ComputeFileHash(fMacroPath);
  }

  const auto* detector = dynamic_cast<const B02DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  const auto* generator = dynamic_cast<const B02PrimaryGeneratorAction*>(
      G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  const auto* physList = G4RunManager::GetRunManager()->GetUserPhysicsList();

  double cfgSourceZ = generator ? generator->GetSourcePlaneZ() / cm : 0.0;
  double cfgSourceLx = generator ? generator->GetSourcePlaneLx() / cm : 0.0;
  double cfgSourceLy = generator ? generator->GetSourcePlaneLy() / cm : 0.0;
  double cfgThetaMax = generator ? generator->GetThetaMaxDeg() : 0.0;
  double cfgEmin = generator ? generator->GetEffectiveEminGeV() : 0.0;
  double cfgEmax = generator ? generator->GetEffectiveEmaxGeV() : 0.0;
  std::string muonMode = generator ? std::string(generator->GetMuonModeString()) : "";
  std::string fluxModel = generator ? std::string(generator->GetFluxModelString()) : "";
  std::string muonChargeMode = generator ? std::string(generator->GetMuonChargeMode()) : "equal";
  double muonChargeRatio = generator ? generator->GetMuonChargeRatio() : 0.0;

  fPhysicsListName = physList ? typeid(*physList).name() : "unknown";

  // Geometry knobs from detector
  double ccdGammaCut_cm = 0.0;
  double ccdElectronCut_cm = 0.0;
  double ccdPositronCut_cm = 0.0;
  double ccdMaxStep_cm = 0.0;
  double ccdThickness_cm = 0.0;
  if (detector) {
    ccdGammaCut_cm = detector->GetCCDGammaCut() / cm;
    ccdElectronCut_cm = detector->GetCCDElectronCut() / cm;
    ccdPositronCut_cm = detector->GetCCDPositronCut() / cm;
    ccdMaxStep_cm = detector->GetCCDMaxStep() / cm;
    fCCDThickness = detector->GetCCDThickness();
    ccdThickness_cm = fCCDThickness / cm;
    fOverburdenEnabled = detector->IsOverburdenEnabled();
    fOverburdenThickness = detector->GetOverburdenThickness();
    fOverburdenZTop = detector->GetOverburdenZTop();
    fOverburdenMaterial = detector->GetOverburdenMaterialName();
  }

  analysisManager->FillNtupleSColumn(ntupleId, 0, fGitHash);
  analysisManager->FillNtupleIColumn(ntupleId, 1, fGitDirty ? 1 : 0);
  analysisManager->FillNtupleSColumn(ntupleId, 2, fMacroPath);
  analysisManager->FillNtupleSColumn(ntupleId, 3, fMacroHash);
  analysisManager->FillNtupleSColumn(ntupleId, 4, fProvenanceTag);
  analysisManager->FillNtupleSColumn(ntupleId, 5, fPhysicsListName);
  analysisManager->FillNtupleIColumn(ntupleId, 6, fUseTimeSeed ? 1 : 0);
  analysisManager->FillNtupleIColumn(ntupleId, 7, static_cast<G4int>(fSeed1));
  analysisManager->FillNtupleIColumn(ntupleId, 8, static_cast<G4int>(fSeed2));
  analysisManager->FillNtupleSColumn(ntupleId, 9, muonMode);
  analysisManager->FillNtupleSColumn(ntupleId,10, fluxModel);
  analysisManager->FillNtupleSColumn(ntupleId,11, muonChargeMode);
  analysisManager->FillNtupleDColumn(ntupleId,12, muonChargeRatio);
  analysisManager->FillNtupleDColumn(ntupleId,13, cfgSourceZ);
  analysisManager->FillNtupleDColumn(ntupleId,14, cfgSourceLx);
  analysisManager->FillNtupleDColumn(ntupleId,15, cfgSourceLy);
  analysisManager->FillNtupleDColumn(ntupleId,16, cfgThetaMax);
  analysisManager->FillNtupleDColumn(ntupleId,17, cfgEmin);
  analysisManager->FillNtupleDColumn(ntupleId,18, cfgEmax);
  analysisManager->FillNtupleIColumn(ntupleId,19, fOverburdenEnabled ? 1 : 0);
  analysisManager->FillNtupleDColumn(ntupleId,20, fOverburdenThickness / cm);
  analysisManager->FillNtupleDColumn(ntupleId,21, fOverburdenZTop / cm);
  analysisManager->FillNtupleSColumn(ntupleId,22, fOverburdenMaterial);
  analysisManager->FillNtupleDColumn(ntupleId,23, ccdGammaCut_cm);
  analysisManager->FillNtupleDColumn(ntupleId,24, ccdElectronCut_cm);
  analysisManager->FillNtupleDColumn(ntupleId,25, ccdPositronCut_cm);
  analysisManager->FillNtupleDColumn(ntupleId,26, ccdMaxStep_cm);
  analysisManager->FillNtupleDColumn(ntupleId,27, ccdThickness_cm);
  analysisManager->AddNtupleRow(ntupleId);

  G4cout << "[RunInfo] git=" << fGitHash << (fGitDirty ? " (dirty)" : "")
         << ", macro=" << fMacroPath << ", seeds=(" << fSeed1 << "," << fSeed2 << ")"
         << ", physicsList=" << fPhysicsListName << G4endl;
}
