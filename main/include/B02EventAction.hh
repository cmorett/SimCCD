//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef B02EventAction_h
#define B02EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

// New modification based in B5 example 


/// Event action


class B02EventAction : public G4UserEventAction
{
public:
    B02EventAction();
    ~B02EventAction() override = default;
    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;
    	
    void AddLength(double dl){fLength += dl;};
    void AddEdep(double de){fEnergy += de;};
    //G4double GetTracedLength() const;
    
private:
    B02EventAction* fB02EventAction;
    int fBarCollID;                
    double fEpri;
    double fpx, fpy, fpz, fp, fth, fph;
    double fLength;
    double fEnergy;
};

//inline void B02EventAction::AddLength(G4double dl) {
//  fLength += dl;
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
