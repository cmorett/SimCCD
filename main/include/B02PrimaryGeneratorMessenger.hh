//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef B02PrimaryGeneratorMessenger_h
#define B02PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class B02PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class B02PrimaryGeneratorMessenger: public G4UImessenger
{
public:
  B02PrimaryGeneratorMessenger(B02PrimaryGeneratorAction*);
  virtual ~B02PrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  B02PrimaryGeneratorAction* B02Action;
  G4UIdirectory*               gunDir; 
  G4UIcmdWithAString*          RndmCmd;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
