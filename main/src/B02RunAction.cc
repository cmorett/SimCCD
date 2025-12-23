//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4ANALYSIS_USE
//#include "B02AnalysisManager.hh"
#endif

#include "B02RunAction.hh"
#include <ctime>
#include "Randomize.hh"
#include "B02EventAction.hh"
// new headers, sep3 2024

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "B02BarSD.hh"
#include "B02DetectorConstruction.hh"
#include "G4Threading.hh"
#include "G4GenericMessenger.hh"
// modification to a new implementation based on runAction.cc class example B5


B02RunAction::B02RunAction(B02EventAction* b02eventAction )
:fB02EventAction(b02eventAction)

{
  fMessenger = new G4GenericMessenger(this, "/runAction/", "Run control");
  fMessenger->DeclareProperty("useTimeSeed", fUseTimeSeed,
                              "Seed the random engine from wall-clock time at run start.");
  fMessenger->DeclareProperty("useFixedSeeds", fUseFixedSeeds,
                              "Seed the random engine from seed1/seed2 at run start.");
  fMessenger->DeclareProperty("seed1", fSeed1, "First random seed used when useFixedSeeds=true.");
  fMessenger->DeclareProperty("seed2", fSeed2, "Second random seed used when useFixedSeeds=true.");

  // Create the generic analysis manager  // new modification sep 6 
  auto analysisManager = G4AnalysisManager::Instance();
  
   
  analysisManager->SetDefaultFileType("root");
     // If the filename extension is not provided, the default file type (root)
     // will be used for all files specified without extension.
  analysisManager->SetVerboseLevel(1);

  // Default settings
  analysisManager->SetNtupleMerging(G4Threading::IsMultithreadedApplication());
     // Note: merging ntuples is available only with Root output
  analysisManager->SetFileName("B02");

  // Book histograms, ntuple

  // Creating 1D histograms
  analysisManager->CreateH1("Epri","Primary muon kinetic energy (GeV)", 100, 0., 200.);        // h1 Id = 0
  analysisManager->CreateH1("ThAng","cos(#theta_{pri})", 100, -1.0, 1.0);  // h1 Id = 1
  analysisManager->CreateH1("Ebar", "Ebar", 100,0,700);          // h1 Id = 2
  analysisManager->CreateH1("XBar", "Xbar", 100, -1500, 1500);   // h1 Id = 3
  analysisManager->CreateH1("LDis" ,"Length traced", 100, 0., 500); // h1 Id = 4 
  analysisManager->CreateH1("MuonDecayEdepLAr"," Energy deposited by decay muons", 100, 0, 1000);


  // Create 1st ntuple (id = 0)
  if ( fB02EventAction ) {
    analysisManager->CreateNtuple("B02Hits", "B02Hits");
    analysisManager->CreateNtupleIColumn("evtId");     // column Id = 0
    analysisManager->CreateNtupleDColumn("XhitBar");   // column Id = 1
    analysisManager->CreateNtupleDColumn("YhitBar");   // column Id = 2
    analysisManager->CreateNtupleDColumn("ZhitBar");   // column Id = 3
    analysisManager->CreateNtupleDColumn("RhitBar");   // column Id = 4
    analysisManager->CreateNtupleDColumn("EhitBar");   // column Id = 5
    analysisManager->CreateNtupleDColumn("WhitBar");   // column Id = 6
   //analysisManager->CreateNtupleDColumn("LDisHits");  // column Id = 7
    analysisManager->FinishNtuple();


    // Set ntuple output file
  //analysisManager->SetNtupleFileName(0, "B02ntupleHits");

  
  
  // Create 2nd ntuple (id = 1)
    //
    analysisManager->CreateNtuple("B02Evts", "B02Evts");
    // Column order: 0 evtID, 1 EevtBar, 2 WevtBar, 3 GevtBar, 4 ADCevtBar,
    // 5 EevtPri, 6 thetaPri, 7 phiPri, 8 nHitBar, 9 LengthMuLAr,
    // 10 MuonDecayEdepLAr, 11 muonX0, 12 muonY0, 13 muonZ0.
    // 14 muonXImp, 15 muonYImp, 16 muonZImp.
    // 17 EdepCCD, 18 nStepsCCD, 19-21 entry CCD (x,y,z),
    // 22-24 exit CCD (x,y,z), 25 trackLenCCD, 26-28 dirX/Y/Z.
    analysisManager->CreateNtupleIColumn("evtID");
    analysisManager->CreateNtupleDColumn("EevtBar");
    analysisManager->CreateNtupleDColumn("WevtBar");
    analysisManager->CreateNtupleDColumn("GevtBar");
    analysisManager->CreateNtupleDColumn("ADCevtBar");
    analysisManager->CreateNtupleDColumn("EevtPri");
    analysisManager->CreateNtupleDColumn("thetaPri");
    analysisManager->CreateNtupleDColumn("phiPri");
    analysisManager->CreateNtupleIColumn("nHitBar");
    analysisManager->CreateNtupleDColumn("LengthMuLAr");
    analysisManager->CreateNtupleDColumn("MuonDecayEdepLAr");
    analysisManager->CreateNtupleDColumn("muonX0");
    analysisManager->CreateNtupleDColumn("muonY0");
    analysisManager->CreateNtupleDColumn("muonZ0");
    analysisManager->CreateNtupleDColumn("muonXImp");
    analysisManager->CreateNtupleDColumn("muonYImp");
    analysisManager->CreateNtupleDColumn("muonZImp");
    analysisManager->CreateNtupleDColumn("EdepCCD");
    analysisManager->CreateNtupleIColumn("nStepsCCD");
    analysisManager->CreateNtupleDColumn("xEntryCCD");
    analysisManager->CreateNtupleDColumn("yEntryCCD");
    analysisManager->CreateNtupleDColumn("zEntryCCD");
    analysisManager->CreateNtupleDColumn("xExitCCD");
    analysisManager->CreateNtupleDColumn("yExitCCD");
    analysisManager->CreateNtupleDColumn("zExitCCD");
    analysisManager->CreateNtupleDColumn("trackLenCCD");
    analysisManager->CreateNtupleDColumn("dirX");
    analysisManager->CreateNtupleDColumn("dirY");
    analysisManager->CreateNtupleDColumn("dirZ");
    analysisManager->FinishNtuple();

 
  // Set ntuple output file
  //analysisManager->SetNtupleFileName(1, "B02ntupleEvts");
  
  analysisManager->SetNtupleFileName(0, "B02ntuples");
  analysisManager->SetNtupleFileName(1, "B02ntuples");
  
  }

 
}

B02RunAction::~B02RunAction()
{
  delete fMessenger;
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
 
  //inform the runManager to save random number seed
  
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();


  // Reset histograms from previous run
  analysisManager->Reset();

  // Open an output file
  // The default file name is set in RunAction::RunAction(),
  // it can be overwritten in a macro
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02RunAction::EndOfRunAction(const G4Run* /*run*/)
{

   auto analysisManager = G4AnalysisManager::Instance();

//if(fTree) fTree->commit();
 if(analysisManager->GetH1(2)){
 // if(fEhitBar) {
    
    G4cout<<"Histo: EpriHis (GeV): mean: "<<analysisManager->GetH1(0)->mean()<< " rms: "<<analysisManager->GetH1(0)->rms()<<G4endl;
    //G4cout<<"Histo: EpriHis (GeV): mean: "<<fEpriHis->mean()/GeV << " rms: "<<fEpriHis->rms()/GeV <<G4endl;
    G4cout<<"Histo: cos(theta): mean: "<<analysisManager->GetH1(1)->mean()<< " rms: "<<analysisManager->GetH1(1)->rms()<<G4endl;
	//	G4cout<<"Histo: ThAng: mean: "<<fThHis->mean() << " rms: "<<fThHis->rms() <<G4endl;
    G4cout<<"Histo: EhitBar (MeV): mean: "<<analysisManager->GetH1(2)->mean()<< " rms: "<<analysisManager->GetH1(2)->rms()<<G4endl;
    //G4cout<<"Histo: EhitBar (MeV): mean: "<<fEhitBar->mean()/MeV << " rms: "<<fEhitBar->rms()/MeV <<G4endl;
    G4cout<<"Histo: XhitBar (cm): mean: "<<analysisManager->GetH1(3)->mean()<< " rms: "<<analysisManager->GetH1(3)->rms()<<G4endl;
    //G4cout<<"Histo: XhitBar (cm) : mean: "<<fXhitBar->mean()/cm  << " rms: "<<fXhitBar->rms()/cm  <<G4endl;
    G4cout<<"Histo: LDis (cm): mean: "<<analysisManager->GetH1(4)->mean()<< " rms: "<<analysisManager->GetH1(4)->rms()<<G4endl;
    G4cout<<"Histo: Muon decay energy deposited (MeV): mean: "<<analysisManager->GetH1(5)->mean()<< " rms: "<<analysisManager->GetH1(5)->rms()<<G4endl;

}


  // save histograms & ntuple
  //
  //auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile(false);
    // Keep content of histos so that they are plotted.
    // The content will be reset at start of the next run.




}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
