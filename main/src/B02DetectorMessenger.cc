//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B02DetectorMessenger.hh"

#include "B02DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02DetectorMessenger::B02DetectorMessenger(B02DetectorConstruction* myDet)
:myDetector(myDet)
{

  B02Dir = new G4UIdirectory("/B02/");
  B02Dir->SetGuidance("UI commands specific to this example.");
  
  detDir = new G4UIdirectory("/B02/det/");
  detDir->SetGuidance("detector control.");
  
  RockMatCmd = new G4UIcmdWithAString("/B02/det/setRockMate",this);
  RockMatCmd->SetGuidance("Select Material of the Rock.");
  RockMatCmd->SetParameterName("choice",false);
  RockMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Rock2MatCmd = new G4UIcmdWithAString("/B02/det/setRock2Mate",this);
  Rock2MatCmd->SetGuidance("Select Material of the Rock2.");
  Rock2MatCmd->SetParameterName("choice",false);
  Rock2MatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BarMatCmd = new G4UIcmdWithAString("/B02/det/setBarMate",this);
  BarMatCmd->SetGuidance("Select Material of the Bar.");
  BarMatCmd->SetParameterName("choice",false);
  BarMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02DetectorMessenger::~B02DetectorMessenger()
{
  delete RockMatCmd;
  delete detDir;
  delete B02Dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
/*  if( command == RockMatCmd )
   { myDetector->setRockMaterial(newValue);}*/

/*  if( command == Rock2MatCmd )
   { myDetector->setRock2Material(newValue);}*/

 // if( command == BarMatCmd )
//   { myDetector->setBarMaterial(newValue);}

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
