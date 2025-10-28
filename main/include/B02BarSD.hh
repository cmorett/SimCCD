//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef B02BarSD_h
#define B02BarSD_h 1

#include "G4VSensitiveDetector.hh"
#include "B02BarHit.hh"
//#include "G4EMSaturation.hh"
//#include "G4LossTableManager.hh"
#include <G4EmSaturation.hh>
#include <G4LossTableManager.hh>
#include <G4Track.hh>
#include <G4SystemOfUnits.hh>


class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class B02BarSD : public G4VSensitiveDetector
{
  public:
      B02BarSD(G4String);
     ~B02BarSD();
     
     //G4double GetTotalTrackLength() const { return totalTrackLength; }

      void   Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void   EndOfEvent(G4HCofThisEvent*);

  private:
      B02BarHitsCollection* barCollection;
      G4EmSaturation* emSaturation;
      //G4EmSaturation* femSaturation; // AAA-23Ago2024

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
