//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B02BarSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Scintillation.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4EmSaturation.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


B02BarSD::B02BarSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="barCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02BarSD::~B02BarSD()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02BarSD::Initialize(G4HCofThisEvent* HCE)
{
  barCollection = new B02BarHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, barCollection ); 

  //femSaturation = G4LossTableManager::Instance()->EmSaturation(); //test
  
   // Birks correction in G4Scintillation
   
     emSaturation = G4LossTableManager::Instance()->EmSaturation();
     
      // Create an instance of G4Scintillation
    G4Scintillation* scintillationProcess = new G4Scintillation();

     // Add Birks correction (saturation effect) to the scintillation process
      scintillationProcess->AddSaturation(emSaturation);
     
       //G4Scintillation::AddSaturation(emSaturation);
    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B02BarSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
 
  G4double evis = emSaturation->VisibleEnergyDepositionAtAStep(aStep);
  
  //G4double evis = femSaturation->VisibleEnergyDeposition(aStep);  //giviubg trouble in G4 11.00
  //G4double evis = edep;
  
  //G4double evis = aStep->GetTotalEnergyDeposit();



  if(edep==0.) return false;

  B02BarHit* newHit = new B02BarHit();
  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetEdep     (edep);
  newHit->SetEvis     (evis);
  newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
  barCollection->insert( newHit );
  
  //newHit->Print();
  newHit->Draw();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02BarSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel>0) { 
     G4int NbHits = barCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
            << " hits in the bar: " << G4endl;
     for (G4int i=0;i<NbHits;i++) (*barCollection)[i]->Print();
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

