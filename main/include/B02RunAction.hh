//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef B02RunAction_h
#define B02RunAction_h 1

#include "G4UserRunAction.hh"
//#include "B02EventAction.hh" // AAA
#include "B02BarSD.hh"

class G4Run;
class G4GenericMessenger;


class B02EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//modification 2024 

class B02RunAction : public G4UserRunAction
{
  public:
    B02RunAction(B02EventAction* b02eventAction);
   ~B02RunAction() override;

  public:
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;
  // modification 2024
  private:
    B02EventAction* fB02EventAction = nullptr;
    G4GenericMessenger* fMessenger = nullptr;
    G4bool fUseTimeSeed = true;
    G4bool fUseFixedSeeds = false;
    G4long fSeed1 = 12345;
    G4long fSeed2 = 67890;
    
};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
