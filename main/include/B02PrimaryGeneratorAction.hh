//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifndef B02PrimaryGeneratorAction_h
#define B02PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
//#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"

class B02DetectorConstruction;
class G4ParticleGun;
class G4Event;
class B02PrimaryGeneratorMessenger;
class G4GenericMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class B02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    B02PrimaryGeneratorAction(B02DetectorConstruction*);
   ~B02PrimaryGeneratorAction() override;

  public:
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(G4String val) { rndmFlag = val;}

  private:
    void DefineCommands();
    G4GenericMessenger* fMessenger = nullptr; 
    G4ThreeVector fParticlePosition;
    G4ParticleGun* particleGun;
    G4GeneralParticleSource* fParticleSource;
    B02DetectorConstruction* myDetector;

    B02PrimaryGeneratorMessenger* gunMessenger; //messenger of this class
    G4String                         rndmFlag;	   //flag for a rndm impact point
    
    // ==== Dimensiones para caja de 10cm (ICN) ==== //
    // double R = 70.;
    // double px = 18.;
    // double py = 18.; 

    // ==== Dimensiones para CCD de 250x529 (ICN) ==== //
    double R = 8.; // cm
    double px = 1.5; // cm
    double py = 1.5; // cm

    // ==== Dimensiones para CONNIE (1022x420) ==== //
    // double R = 7.; // cm
    // double px = 1.7; // cm
    // double py = 1.7; // cm

    bool doSmith;
    //G4double energy = 10*MeV;
   
   /* 
   G4double x;
   G4double y;
   G4double z;
   G4double Px;
   G4double Py;
   G4double Pz; 
   */

};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


