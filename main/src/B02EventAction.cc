//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifdef G4ANALYSIS_USE
#include "B02AnalysisManager.hh"
#endif

#include "B02EventAction.hh"

#include <cmath>

// New headers sep 2024
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
// incorporate headers from B02AnalysisManager 
#include "B02BarHit.hh"


B02EventAction::B02EventAction()
  : fB02EventAction(nullptr),
    fBarCollID(-1),
    fEpri(0.0),
    fpx(0.0),
    fpy(0.0),
    fpz(0.0),
    fp(0.0),
    fth(0.0),
    fph(0.0),
    fLength(0.0),
    fEnergy(0.0)
{
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

  const auto primaryVertex = aEvent->GetPrimaryVertex(0);
  if (!primaryVertex) {
    fEpri = 0.0;
    fth = 0.0;
    fph = 0.0;
    return;
  }

  const auto primaryParticle = primaryVertex->GetPrimary(0);
  if (!primaryParticle) {
    fEpri = 0.0;
    fth = 0.0;
    fph = 0.0;
    return;
  }

  fEpri = primaryParticle->GetKineticEnergy();
  const auto momentum = primaryParticle->GetMomentum();
  const auto momentumMag2 = momentum.mag2();
  if (momentumMag2 > 0.) {
    fth = momentum.theta();
    fph = momentum.phi();
  } else {
    fth = 0.0;
    fph = 0.0;
  }

  analysisManager->FillH1(0,fEpri/MeV);
  analysisManager->FillH1(1,std::cos(fth));
}

void B02EventAction::EndOfEventAction(const G4Event* aEvent)

{

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Extract primary truth for ntuple filling
  G4double muonX0 = 0.0;
  G4double muonY0 = 0.0;
  G4double muonZ0 = 0.0;
  G4double thetaPri = 0.0;
  G4double phiPri = 0.0;
  G4double primaryEnergy = 0.0;

  if (const auto primaryVertex = aEvent->GetPrimaryVertex(0)) {
    muonX0 = primaryVertex->GetX0();
    muonY0 = primaryVertex->GetY0();
    muonZ0 = primaryVertex->GetZ0();

    if (const auto primaryParticle = primaryVertex->GetPrimary(0)) {
      primaryEnergy = primaryParticle->GetKineticEnergy();
      const auto momentum = primaryParticle->GetMomentum();
      if (momentum.mag2() > 0.) {
        thetaPri = momentum.theta();
        phiPri = momentum.phi();
      }
    }
  }

  fEpri = primaryEnergy;
  fth = thetaPri;
  fph = phiPri;
 
 // accumulated length by muons.
 
  // G4cout << "TotalLength " << fLength/cm << G4endl;
  //G4cout << "TotalEnergyMuonDecay " << fEnergy/MeV << G4endl;
  analysisManager->FillH1(4, fLength/cm);
  analysisManager->FillH1(5, fEnergy/MeV);
////    
    
  if(!analysisManager->GetH1(2)) return; // No histo booked !
  
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

  G4int n_hit = 0;
  G4double EevtBar = -5.;
  G4double WevtBar = -5.;
  G4double GevtBar = 0.;
  G4double ADCevtBar = 0.;

  if(BHC) {
    n_hit = BHC->entries();
    EevtBar = 0.;
    WevtBar = 0.;
    for (G4int i=0;i<n_hit;i++) {
      G4double EhitBar = (*BHC)[i]->GetEdep();
      G4double WhitBar = (*BHC)[i]->GetEvis();
      G4double XhitBar = (*BHC)[i]->GetPos().getX();
      G4double YhitBar = (*BHC)[i]->GetPos().getY();
      G4double ZhitBar = (*BHC)[i]->GetPos().getZ();
      G4double RhitBar = sqrt(XhitBar*XhitBar+YhitBar*YhitBar+ZhitBar*ZhitBar);
      analysisManager->FillH1(2,EhitBar/MeV);
      analysisManager->FillH1(3,XhitBar/cm);
      analysisManager->FillNtupleIColumn(0,0,event_id);
      analysisManager->FillNtupleDColumn(0,1,XhitBar/cm);
      analysisManager->FillNtupleDColumn(0,2,YhitBar/cm);
      analysisManager->FillNtupleDColumn(0,3,ZhitBar/cm);
      analysisManager->FillNtupleDColumn(0,4,RhitBar/cm);
      analysisManager->FillNtupleDColumn(0,5,EhitBar/MeV);
      analysisManager->FillNtupleDColumn(0,6,WhitBar/MeV);
      analysisManager->AddNtupleRow(0);

      EevtBar += EhitBar;
      WevtBar += WhitBar;
    } // for i hits of hit collection

    if (n_hit==0) {
      EevtBar = -5.;
      WevtBar = -5.;
    }
  }

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
      analysisManager->FillNtupleDColumn(1,11,muonX0/cm);
      analysisManager->FillNtupleDColumn(1,12,muonY0/cm);
      analysisManager->FillNtupleDColumn(1,13,muonZ0/cm);
      analysisManager->AddNtupleRow(1);

// analysisManager->FillNtupleDColumn(1,9,fLength/cm);
// analysisManager->FillNtupleDColumn(1,10,fEnergy/MeV);

}  //B02EndOfEventAction



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
