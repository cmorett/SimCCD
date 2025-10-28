//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef B02DetectorMessenger_h
#define B02DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class B02DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class B02DetectorMessenger: public G4UImessenger
{
  public:
   B02DetectorMessenger(B02DetectorConstruction*);
   ~B02DetectorMessenger();
   
  
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    B02DetectorConstruction* myDetector;
    
    G4UIdirectory*             B02Dir;
    G4UIdirectory*             detDir;
    G4UIcmdWithAString*        RockMatCmd;
    G4UIcmdWithAString*        Rock2MatCmd;
    G4UIcmdWithAString*        BarMatCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
