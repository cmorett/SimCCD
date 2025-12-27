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
#include "B02PrimaryGeneratorAction.hh"
#include "B02RunAction.hh"
#include <functional>


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
    fEnergy(0.0),
    fEdepCCDGeV(0.0),
    fNstepsCCD(0),
    fEntryCCD(),
    fExitCCD(),
    fTrackLenCCD(0.0),
    fHasEntryCCD(false),
    fDirX(0.0),
    fDirY(0.0),
    fDirZ(-1.0)
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}



void B02EventAction::BeginOfEventAction(const G4Event* aEvent)
{

  fLength = 0.0; // mm
  fEnergy = 0.0;
  fEdepCCDGeV = 0.0;
  fEdepOtherGeV = 0.0;
  fEdepGeomFlag = 0.0;
  fFirstHitIsCCD = 0;
  fNstepsCCD = 0;
  fEntryCCD = G4ThreeVector();
  fExitCCD = G4ThreeVector();
  fTrackLenCCD = 0.0;
  fHasEntryCCD = false;
  fDirX = 0.0;
  fDirY = 0.0;
  fDirZ = -1.0;
  
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

  analysisManager->FillH1(0,fEpri/GeV);
  analysisManager->FillH1(1,std::cos(fth));
}

void B02EventAction::EndOfEventAction(const G4Event* aEvent)

{
  if (fHasEntryCCD) {
    fTrackLenCCD = (fExitCCD - fEntryCCD).mag();
  } else {
    fEntryCCD = G4ThreeVector(0.,0.,0.);
    fExitCCD = G4ThreeVector(0.,0.,0.);
    fTrackLenCCD = 0.0;
  }

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Extract primary truth for ntuple filling
  G4double muonX0 = 0.0;
  G4double muonY0 = 0.0;
  G4double muonZ0 = 0.0;
  G4double muonXImp = 0.0;
  G4double muonYImp = 0.0;
  G4double muonZImp = 0.0;
  G4double thetaPri = 0.0;
  G4double phiPri = 0.0;
  G4double primaryEnergy = 0.0;
  G4double muonWeight = 0.0;
  G4double eventLivetime = 0.0;
  G4double muonCosTheta = 1.0;
  G4double muonEnergySampledGeV = 0.0;
  G4int muonModeCode = 0;
  G4int muonIsTargeted = 0;
  G4int fluxModelCode = 0;
  G4double cfgSourceZ = 0.0;
  G4double cfgSourceLx = 0.0;
  G4double cfgSourceLy = 0.0;
  G4double cfgThetaMaxDeg = 0.0;
  G4double cfgEminEff = 0.0;
  G4double cfgEmaxEff = 0.0;
  G4ThreeVector primaryDir(0.,0.,-1.);
  G4int muonPDGCode = 13;
  G4int muonChargeSign = -1;
  G4int seed1 = 0;
  G4int seed2 = 0;
  G4int useTimeSeed = 0;
  G4int overburdenEnabled = 0;
  G4double overburdenThickness = 0.0;
  G4double overburdenZTop = 0.0;
  G4double overburdenMaterialHash = 0.0;
  G4double ccdGammaCut = 0.0;
  G4double ccdElectronCut = 0.0;
  G4double ccdPositronCut = 0.0;
  G4double ccdMaxStep = 0.0;
  G4double ccdThickness = 0.0;
  G4double gitHashCode = 0.0;
  G4double macroHashCode = 0.0;
  G4double macroPathHash = 0.0;
  G4double physicsListHash = 0.0;
  G4double muonChargeRatio = 0.0;

  if (const auto primaryVertex = aEvent->GetPrimaryVertex(0)) {
    muonX0 = primaryVertex->GetX0();
    muonY0 = primaryVertex->GetY0();
    muonZ0 = primaryVertex->GetZ0();

    if (const auto primaryParticle = primaryVertex->GetPrimary(0)) {
      primaryEnergy = primaryParticle->GetKineticEnergy();
      muonPDGCode = primaryParticle->GetPDGcode();
      muonChargeSign = (muonPDGCode < 0) ? 1 : -1;
      const auto momentum = primaryParticle->GetMomentum();
      if (momentum.mag2() > 0.) {
        thetaPri = momentum.theta();
        phiPri = momentum.phi();
        primaryDir = momentum.unit();
        muonCosTheta = -primaryDir.z();
      }
    }
  }

  fEpri = primaryEnergy;
  fth = thetaPri;
  fph = phiPri;
  muonEnergySampledGeV = primaryEnergy / GeV;
  fDirX = primaryDir.x();
  fDirY = primaryDir.y();
  fDirZ = primaryDir.z();

  if (const auto* gen = dynamic_cast<const B02PrimaryGeneratorAction*>(
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())) {
    muonXImp = gen->GetMuonXImpact();
    muonYImp = gen->GetMuonYImpact();
    muonZImp = gen->GetMuonZImpact();
    muonWeight = gen->GetMuonWeight();
    eventLivetime = gen->GetEventLivetime();
    muonCosTheta = gen->GetMuonCosTheta();
    muonEnergySampledGeV = gen->GetMuonEnergySampledGeV();
    G4String mode = gen->GetMuonModeString();
    mode.toLower();
    if (mode == "tierb_plane_flux") {
      muonModeCode = 1;
    } else {
      muonModeCode = 0;
    }
    fEdepGeomFlag = gen->GetGeomIntersectsCCD() ? 1.0 : 0.0;
    fluxModelCode = gen->GetFluxModelCode();
    cfgSourceZ = gen->GetSourcePlaneZ() / cm;
    cfgSourceLx = gen->GetSourcePlaneLx() / cm;
    cfgSourceLy = gen->GetSourcePlaneLy() / cm;
    cfgThetaMaxDeg = gen->GetThetaMaxDeg();
    cfgEminEff = gen->GetEffectiveEminGeV();
    cfgEmaxEff = gen->GetEffectiveEmaxGeV();
    muonPDGCode = gen->GetMuonPDGCode();
    muonChargeSign = gen->GetMuonChargeSign();
    muonChargeRatio = gen->GetMuonChargeRatio();
    muonIsTargeted = gen->GetIsTargeted() ? 1 : 0;
  }

  if (const auto* runAction = dynamic_cast<const B02RunAction*>(
          G4RunManager::GetRunManager()->GetUserRunAction())) {
    seed1 = static_cast<G4int>(runAction->GetSeed1());
    seed2 = static_cast<G4int>(runAction->GetSeed2());
    useTimeSeed = runAction->GetUseTimeSeed() ? 1 : 0;
    overburdenEnabled = runAction->GetOverburdenEnabled() ? 1 : 0;
    overburdenThickness = runAction->GetOverburdenThickness() / cm;
    overburdenZTop = runAction->GetOverburdenZTop() / cm;
    overburdenMaterialHash =
        static_cast<G4double>(std::hash<std::string>{}(runAction->GetOverburdenMaterialName()));
    ccdGammaCut = runAction->GetCCDGammaCut() / cm;
    ccdElectronCut = runAction->GetCCDElectronCut() / cm;
    ccdPositronCut = runAction->GetCCDPositronCut() / cm;
    ccdMaxStep = runAction->GetCCDMaxStep() / cm;
    ccdThickness = runAction->GetCCDThicknessCached() / cm;
    gitHashCode = static_cast<G4double>(std::hash<std::string>{}(runAction->GetGitHash()));
    macroHashCode = static_cast<G4double>(std::hash<std::string>{}(runAction->GetMacroHash()));
    macroPathHash = static_cast<G4double>(std::hash<std::string>{}(runAction->GetMacroPath()));
    physicsListHash =
        static_cast<G4double>(std::hash<std::string>{}(runAction->GetPhysicsListName()));
  }
 
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
      analysisManager->FillNtupleDColumn(1,14,muonXImp/cm);
      analysisManager->FillNtupleDColumn(1,15,muonYImp/cm);
      analysisManager->FillNtupleDColumn(1,16,muonZImp/cm);
      analysisManager->FillNtupleDColumn(1,17,fEdepCCDGeV);
      analysisManager->FillNtupleDColumn(1,18,fEdepOtherGeV);
      analysisManager->FillNtupleIColumn(1,19,fNstepsCCD);
      analysisManager->FillNtupleDColumn(1,20,fEntryCCD.x()/cm);
      analysisManager->FillNtupleDColumn(1,21,fEntryCCD.y()/cm);
      analysisManager->FillNtupleDColumn(1,22,fEntryCCD.z()/cm);
      analysisManager->FillNtupleDColumn(1,23,fExitCCD.x()/cm);
      analysisManager->FillNtupleDColumn(1,24,fExitCCD.y()/cm);
      analysisManager->FillNtupleDColumn(1,25,fExitCCD.z()/cm);
      analysisManager->FillNtupleDColumn(1,26,fTrackLenCCD/cm);
      analysisManager->FillNtupleDColumn(1,27,fDirX);
      analysisManager->FillNtupleDColumn(1,28,fDirY);
      analysisManager->FillNtupleDColumn(1,29,fDirZ);
      analysisManager->FillNtupleDColumn(1,30,muonCosTheta);
      analysisManager->FillNtupleDColumn(1,31,muonWeight);
      analysisManager->FillNtupleDColumn(1,32,eventLivetime);
      analysisManager->FillNtupleIColumn(1,33,muonModeCode);
      analysisManager->FillNtupleDColumn(1,34,muonEnergySampledGeV);
      analysisManager->FillNtupleIColumn(1,35,fFirstHitIsCCD);
      analysisManager->FillNtupleDColumn(1,36,fEdepGeomFlag);
      analysisManager->FillNtupleIColumn(1,37,fluxModelCode);
      analysisManager->FillNtupleDColumn(1,38,cfgSourceZ);
      analysisManager->FillNtupleDColumn(1,39,cfgSourceLx);
      analysisManager->FillNtupleDColumn(1,40,cfgSourceLy);
      analysisManager->FillNtupleDColumn(1,41,cfgThetaMaxDeg);
      analysisManager->FillNtupleDColumn(1,42,cfgEminEff);
      analysisManager->FillNtupleDColumn(1,43,cfgEmaxEff);
      analysisManager->FillNtupleIColumn(1,44,muonPDGCode);
      analysisManager->FillNtupleIColumn(1,45,muonChargeSign);
      analysisManager->FillNtupleIColumn(1,46,seed1);
      analysisManager->FillNtupleIColumn(1,47,seed2);
      analysisManager->FillNtupleIColumn(1,48,useTimeSeed);
      analysisManager->FillNtupleIColumn(1,49,overburdenEnabled);
      analysisManager->FillNtupleDColumn(1,50,overburdenThickness);
      analysisManager->FillNtupleDColumn(1,51,overburdenZTop);
      analysisManager->FillNtupleDColumn(1,52,overburdenMaterialHash);
      analysisManager->FillNtupleDColumn(1,53,ccdGammaCut);
      analysisManager->FillNtupleDColumn(1,54,ccdElectronCut);
      analysisManager->FillNtupleDColumn(1,55,ccdPositronCut);
      analysisManager->FillNtupleDColumn(1,56,ccdMaxStep);
      analysisManager->FillNtupleDColumn(1,57,ccdThickness);
      analysisManager->FillNtupleDColumn(1,58,gitHashCode);
      analysisManager->FillNtupleDColumn(1,59,macroHashCode);
      analysisManager->FillNtupleDColumn(1,60,macroPathHash);
      analysisManager->FillNtupleDColumn(1,61,physicsListHash);
      analysisManager->FillNtupleDColumn(1,62,muonChargeRatio);
      analysisManager->FillNtupleIColumn(1,63,muonIsTargeted);
      analysisManager->AddNtupleRow(1);

// analysisManager->FillNtupleDColumn(1,9,fLength/cm);
// analysisManager->FillNtupleDColumn(1,10,fEnergy/MeV);

}  //B02EndOfEventAction



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void B02EventAction::AddCCDStep(double edepGeV,
                                const G4ThreeVector& prePos,
                                const G4ThreeVector& postPos,
                                G4bool isPrimary,
                                G4bool isCCDVolume)
{
  if (isCCDVolume) {
    fEdepCCDGeV += edepGeV;
    if (fFirstHitIsCCD == 0 && edepGeV > 0) {
      fFirstHitIsCCD = 1;
    }
  } else {
    fEdepOtherGeV += edepGeV;
  }

  if (isPrimary && isCCDVolume) {
    fNstepsCCD++;
    if (!fHasEntryCCD) {
      fEntryCCD = prePos;
      fHasEntryCCD = true;
    }
    fExitCCD = postPos;
  }
}
