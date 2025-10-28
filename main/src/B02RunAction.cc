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
// modification to a new implementation based on runAction.cc class example B5


B02RunAction::B02RunAction(B02EventAction* b02eventAction )
:fB02EventAction(b02eventAction)

{
  // Create the generic analysis manager  // new modification sep 6 
  auto analysisManager = G4AnalysisManager::Instance();
  
   
  analysisManager->SetDefaultFileType("root");
     // If the filename extension is not provided, the default file type (root)
     // will be used for all files specified without extension.
  analysisManager->SetVerboseLevel(1);

  // Default settings
  analysisManager->SetNtupleMerging(true);
     // Note: merging ntuples is available only with Root output
  analysisManager->SetFileName("B02");

  // Book histograms, ntuple

  // Creating 1D histograms
  analysisManager->CreateH1("Epri","Epri", 100, 0., 100);        // h1 Id = 0
  analysisManager->CreateH1("ThAng","ThAng", 100, 0., 3.14115);  // h1 Id = 1
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
    analysisManager->CreateNtupleIColumn("evtID"); // column Id = 0
    analysisManager->CreateNtupleDColumn("EevtBar"); // column Id = 1
    analysisManager->CreateNtupleDColumn("WevtBar");
    analysisManager->CreateNtupleDColumn("GevtBar");
    analysisManager->CreateNtupleDColumn("ADCevtBar");
    analysisManager->CreateNtupleDColumn("EevtPri");
    analysisManager->CreateNtupleDColumn("thetaPri");
    analysisManager->CreateNtupleDColumn("phiPri");
    analysisManager->CreateNtupleIColumn("nHitBar");
    analysisManager->CreateNtupleDColumn("LengthMuLAr");  // column Id = 7
    analysisManager->CreateNtupleDColumn("MuonDecayEdepLAr"); 
    analysisManager->FinishNtuple();

 
  // Set ntuple output file
  //analysisManager->SetNtupleFileName(1, "B02ntupleEvts");
  
  analysisManager->SetNtupleFileName(0, "B02ntuples");
  analysisManager->SetNtupleFileName(1, "B02ntuples");
  
  }

 
}

void B02RunAction::BeginOfRunAction(const G4Run* aRun)
{

long seeds[2];
time_t systime = time(NULL);
seeds[0] = (long) systime;
seeds[1] = (long) (systime*G4UniformRand());
G4Random::setTheSeeds(seeds);

G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
 
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
    
    G4cout<<"Histo: EpriHis (MeV):mean: "<<analysisManager->GetH1(0)->mean()<< " rms: "<<analysisManager->GetH1(0)->rms()<<G4endl;
    //G4cout<<"Histo: EpriHis (MeV): mean: "<<fEpriHis->mean()/MeV << " rms: "<<fEpriHis->rms()/MeV <<G4endl;
    G4cout<<"Histo: ThAng: mean: "<<analysisManager->GetH1(1)->mean()<< " rms: "<<analysisManager->GetH1(1)->rms()<<G4endl;
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
