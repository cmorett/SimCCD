//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef B02RunAction_h
#define B02RunAction_h 1

#include "G4UserRunAction.hh"
#include "B02BarSD.hh"
#include <string>

class G4Run;
class G4GenericMessenger;

class B02EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class B02RunAction : public G4UserRunAction
{
  public:
    B02RunAction(B02EventAction* b02eventAction);
   ~B02RunAction() override;

  public:
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    // Provenance accessors (useful for validation scripts)
    std::string GetGitHash() const { return fGitHash; }
    bool IsGitDirty() const { return fGitDirty; }
    std::string GetMacroPath() const { return fMacroPath; }
    std::string GetMacroHash() const { return fMacroHash; }
    std::string GetPhysicsListName() const { return fPhysicsListName; }
    G4long GetSeed1() const { return fSeed1; }
    G4long GetSeed2() const { return fSeed2; }
    G4bool GetUseTimeSeed() const { return fUseTimeSeed; }
    bool GetOverburdenEnabled() const { return fOverburdenEnabled; }
    G4double GetOverburdenThickness() const { return fOverburdenThickness; }
    G4double GetOverburdenZTop() const { return fOverburdenZTop; }
    std::string GetOverburdenMaterialName() const { return fOverburdenMaterial; }
    G4double GetCCDGammaCut() const { return fCCDGammaCut; }
    G4double GetCCDElectronCut() const { return fCCDElectronCut; }
    G4double GetCCDPositronCut() const { return fCCDPositronCut; }
    G4double GetCCDMaxStep() const { return fCCDMaxStep; }
    G4double GetCCDThicknessCached() const { return fCCDThickness; }

    void SetMacroPath(const G4String& path);
    void SetMacroHash(const G4String& hash);
    void SetProvenanceTag(const G4String& tag);

  private:
    void BuildAnalysis();
    void FillRunInfoNtuple(const G4Run* run);
    std::string RunCommand(const std::string& cmd) const;
    std::string DetectGitHash();
    bool DetectGitDirty();
    std::string ComputeFileHash(const std::string& path) const;

    B02EventAction* fB02EventAction = nullptr;
    G4GenericMessenger* fMessenger = nullptr;
    G4bool fUseTimeSeed = true;
    G4bool fUseFixedSeeds = false;
    G4long fSeed1 = 12345;
    G4long fSeed2 = 67890;

    // Provenance
    std::string fGitHash;
    bool fGitDirty = false;
    std::string fMacroPath;
    std::string fMacroHash;
    std::string fProvenanceTag;
    std::string fPhysicsListName;

    // Geometry/physics knobs captured for run metadata
    double fCCDGammaCut = 0.0;
    double fCCDElectronCut = 0.0;
    double fCCDPositronCut = 0.0;
    double fCCDMaxStep = 0.0;
    double fCCDThickness = 0.0;
    bool fOverburdenEnabled = false;
    double fOverburdenThickness = 0.0;
    double fOverburdenZTop = 0.0;
    std::string fOverburdenMaterial;
};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
