//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef B02SteppingAction_h
#define B02SteppingAction_h 1

#include "B02DetectorConstruction.hh"
#include "B02EventAction.hh"

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"


//class B02AnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//class B02EventAction;

class B02SteppingAction : public G4UserSteppingAction
{
  public:
   
 //B02SteppingAction(B02DetectorConstruction* b02detectorconstruction, 
  //                 B02EventAction* B02eventAction);
  B02SteppingAction(B02EventAction* B02eventAction);
  
  ~B02SteppingAction() override = default;

  void UserSteppingAction(const G4Step* aStep) override;

  private:
    // B02DetectorConstruction* fB02DetectorConstruction = nullptr;
     B02EventAction* fEventAction = nullptr;
     G4double fLength;
     G4double fEnergy;
     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
