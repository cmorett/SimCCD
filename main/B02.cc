// 
// --------------------------------------------------------------
//      GEANT 4 - B02 
// --------------------------------------------------------------

#include "B02DetectorConstruction.hh"
#include "B02PhysicsList.hh"
#include "B02PrimaryGeneratorAction.hh"
#include "B02RunAction.hh"
#include "B02EventAction.hh"
#include "B02SteppingAction.hh"
#include "B02SteppingVerbose.hh"
#include "B02ActionInitialization.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4SteppingVerbose.hh"


int main(int argc,char** argv)
{
 
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }
  
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new B02SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // The Run manager
  //
  G4RunManager* runManager = new G4RunManager;
  
  // Mandatory user initialization classes
    
    
   /// Detector construction
  runManager->SetUserInitialization(new B02DetectorConstruction()); 
  // Physics list
   G4VUserPhysicsList* physics = new B02PhysicsList;
   runManager->SetUserInitialization(physics);
   //  User action initialization
  runManager->SetUserInitialization(new B02ActionInitialization());
  
   
 //////////////////////////////////////////////////////////////////////////////////////  
  /*
   auto B02DetConstruction = new B02DetectorConstruction();
   runManager->SetUserInitialization(B02DetConstruction);
   
   G4VUserPhysicsList* physics = new B02PhysicsList;
   runManager->SetUserInitialization(physics);
    
   G4VUserPrimaryGeneratorAction* gen_action = new B02PrimaryGeneratorAction(new B02DetectorConstruction);
   runManager->SetUserAction(gen_action); 
    
   G4UserEventAction* event_action = new B02EventAction();
   runManager->SetUserAction(event_action);
  
   G4UserRunAction* run_action = new B02RunAction(new B02EventAction());
   runManager->SetUserAction(run_action);
 
   //G4UserSteppingAction* stepping_action = new B02SteppingAction(new B02DetectorConstruction(), new B02EventAction());
   G4UserSteppingAction* stepping_action = new B02SteppingAction(new B02EventAction());
   runManager->SetUserAction(stepping_action);     
*/
////////////////////////////////////////////////////////////////////////////////

  // Initialize G4 kernel
  //
    runManager->Initialize();
    
  // Initialize visualization with the default graphics system
  auto visManager = new G4VisExecutive(argc, argv);
  // Constructors can also take optional arguments:
  // - a graphics system of choice, eg. "OGL"
  // - and a verbosity argument - see /vis/verbose guidance.
  // auto visManager = new G4VisExecutive(argc, argv, "OGL", "Quiet");
  // auto visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
  //
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  
  delete visManager;
  delete runManager;
 
}

