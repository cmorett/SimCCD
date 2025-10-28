//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifdef G4ANALYSIS_USE
#include "B02AnalysisManager.hh"
#endif
 
#include "B02EventAction.hh"

// New headers sep 2024
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
// incorporate headers from B02AnalysisManager 
#include "B02BarHit.hh"
#include "G4PrimaryVertex.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"


B02EventAction::B02EventAction() 

{
  // length distribution of muons 
  fLength = 0.0;

  // energy deposited by muons that will decay
  fEnergy = 0.0;

  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1); 
}



void B02EventAction::BeginOfEventAction(const G4Event* aEvent)
{

  fLength = 0.0; // mm
  fEnergy = 0.0;
  
  if(fBarCollID==-1) {
  
  auto analysisManager = G4AnalysisManager::Instance();
    
     G4SDManager* SDman = G4SDManager::GetSDMpointer();
     fBarCollID = SDman->GetCollectionID("barCollection");
  } 

  //if(!fEpriHis)   return; // No histo booked !
  auto analysisManager = G4AnalysisManager::Instance();
  
  if(!analysisManager->GetH1(0)) return; // no histo Booked !

  G4PrimaryVertex* primaryVertex = aEvent->GetPrimaryVertex();
  G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
  fEpri = primaryParticle->GetKineticEnergy();
	fpx = primaryParticle->GetPx();
	fpy = primaryParticle->GetPy();
	fpz = primaryParticle->GetPz();
	fp = sqrt(fpx*fpx+fpy*fpy+fpz*fpz);
	fth = acos(-fpz/fp);
	fph = atan2(fpy,fpx);
	analysisManager->FillH1(0,fEpri/MeV);
        //fEpriHis->fill(fEpri/GeV);
	analysisManager->FillH1(1, cos(fth));
	//fThHis->fill(cos(fth));
	
}

void B02EventAction::EndOfEventAction(const G4Event* aEvent)

{

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
 
 // accumulated length by muons.
 
  // G4cout << "TotalLength " << fLength/cm << G4endl;
  //G4cout << "TotalEnergyMuonDecay " << fEnergy/MeV << G4endl;
  analysisManager->FillH1(4, fLength/cm);
  analysisManager->FillH1(5, fEnergy/MeV);
////    
    
  if(!analysisManager->GetH1(2)) return; // No histo booked !
  if(!analysisManager->GetNtuple(0)) return;  
  
  G4int event_id = aEvent->GetEventID();  
  //
  // get number of stored trajectories
  //
  G4TrajectoryContainer* trajectoryContainer = aEvent->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  //
  // periodic printing
  //
  if (event_id < 1000 || event_id%100 == 0) {
    G4cout << ">>> Evento " << aEvent->GetEventID() << G4endl;
    // G4cout << "TotalLength " << fLength/mm << G4endl;
    // G4cout << "    " << n_trajectories << " trajectories stored in this event." << G4endl;
  }
  
  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  B02BarHitsCollection* BHC = 
  HCE ? (B02BarHitsCollection*)(HCE->GetHC(fBarCollID)) : 0;

  if(BHC) {
    G4int n_hit = BHC->entries();
    G4double EevtBar = 0.;
    G4double WevtBar = 0.;
    G4double GevtBar = 0.;
    G4double ADCevtBar = 0.;
    for (G4int i=0;i<n_hit;i++) {
      G4double EhitBar = (*BHC)[i]->GetEdep();
      G4double WhitBar = (*BHC)[i]->GetEvis();
      G4double XhitBar = (*BHC)[i]->GetPos().getX();
      G4double YhitBar = (*BHC)[i]->GetPos().getY();
      G4double ZhitBar = (*BHC)[i]->GetPos().getZ();
      G4double RhitBar = sqrt(XhitBar*XhitBar+YhitBar*YhitBar+ZhitBar*ZhitBar);
      analysisManager->FillH1(2,EhitBar/MeV);
      //fEhitBar->fill(EhitBar/MeV);
      analysisManager->FillH1(3,XhitBar/cm);
      //fXhitBar->fill(XhitBar/cm);     
      analysisManager->FillNtupleIColumn(0,0,event_id);
      //fHitTuple->fill(0,event_id);
      analysisManager->FillNtupleDColumn(0,1,XhitBar/cm);
      //fHitTuple->fill(1,XhitBar/cm);
      analysisManager->FillNtupleDColumn(0,2,YhitBar/cm);
      //fHitTuple->fill(2,YhitBar/cm);
      analysisManager->FillNtupleDColumn(0,3,ZhitBar/cm);
      //fHitTuple->fill(3,ZhitBar/cm);
      analysisManager->FillNtupleDColumn(0,4,RhitBar/cm);
      //fHitTuple->fill(4,RhitBar/cm);
      analysisManager->FillNtupleDColumn(0,5,EhitBar/MeV);
      //fHitTuple->fill(5,EhitBar/MeV);
      analysisManager->FillNtupleDColumn(0,6,WhitBar/MeV);
      //fHitTuple->fill(6,WhitBar/MeV);
      
      //analysisManager->FillNtupleDColumn(0,7,fLength/cm);
      
      analysisManager->AddNtupleRow(0);
      //fHitTuple->addRow();
   

      EevtBar += EhitBar;
      WevtBar += WhitBar;
    } // for i hits of hit collection

    if (n_hit==0) {EevtBar = -5.; 
                   WevtBar = -5.;
    
    } // if n_hit==0

		//double res = 0.05;                      //
		//double sig = res*sqrt(550.0*WevtBar);   //
		//double xg = G4RandGauss::shoot(0, sig); //
		//GevtBar = WevtBar+xg;

		//double p1 =  0.0;      // No offset
		//double p2 = 7.3;      // 15 PE/MEV
		//double p3 = 1.0/2000; // NL turns on at ~550

		//ADCevtBar = p1+p2*GevtBar/(1+p3*GevtBar);


      analysisManager->FillNtupleIColumn(1,0,event_id);
	  //fEvtTuple->fill(0,event_id);
      analysisManager->FillNtupleDColumn(1,1,EevtBar/MeV);
	  //fEvtTuple->fill(1,EevtBar/MeV);
      analysisManager->FillNtupleDColumn(1,2,WevtBar/MeV);
	  //fEvtTuple->fill(2,WevtBar/MeV);  
      analysisManager->FillNtupleDColumn(1,3,GevtBar/MeV);
	  //fEvtTuple->fill(3,GevtBar/MeV); no
      analysisManager->FillNtupleDColumn(1,4,ADCevtBar);
	  //fEvtTuple->fill(4,ADCevtBar); no 
      analysisManager->FillNtupleDColumn(1,5, fEpri/GeV);      
          //fEvtTuple->fill(5,fEpri/GeV);
      analysisManager->FillNtupleDColumn(1,6, fth);
	  //fEvtTuple->fill(6,fth);
      analysisManager->FillNtupleDColumn(1,7, fph);
	  //fEvtTuple->fill(7,fph);
      analysisManager->FillNtupleIColumn(1,8, n_hit);
	  //fEvtTuple->fill(8, n_hit);
      analysisManager->FillNtupleDColumn(1,9,fLength/cm);
      analysisManager->FillNtupleDColumn(1,10,fEnergy/MeV);
      analysisManager->AddNtupleRow(1);

} // if BHC

// analysisManager->FillNtupleDColumn(1,9,fLength/cm);
// analysisManager->FillNtupleDColumn(1,10,fEnergy/MeV);

}  //B02EndOfEventAction



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
