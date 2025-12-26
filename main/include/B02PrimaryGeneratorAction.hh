//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifndef B02PrimaryGeneratorAction_h
#define B02PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
//#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"
#include "MuonFluxModel.hh"

#include <utility>

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
    G4double GetMuonX0() const { return fMuonX0; }
    G4double GetMuonY0() const { return fMuonY0; }
    G4double GetMuonZ0() const { return fMuonZ0; }
    G4double GetMuonXImpact() const { return fMuonXImp; }
    G4double GetMuonYImpact() const { return fMuonYImp; }
    G4double GetMuonZImpact() const { return fMuonZImp; }
    G4double GetMuonTheta() const { return fMuonTheta; }
    G4double GetMuonPhi() const { return fMuonPhi; }
    G4double GetMuonWeight() const { return fEventWeight; }
    G4double GetEventLivetime() const { return fEventLivetime; }
    G4String GetMuonModeString() const { return fMuonModeString; }
    G4double GetMuonCosTheta() const { return fMuonCosTheta; }
    G4double GetMuonEnergySampledGeV() const { return fSampledEnergyGeV; }
    G4double GetEffectiveEminGeV() const { return fEffectiveEminGeV; }
    G4double GetEffectiveEmaxGeV() const { return fEffectiveEmaxGeV; }
    G4double GetSourcePlaneZ() const { return fSourcePlaneZ; }
    G4double GetSourcePlaneLx() const { return fSourcePlaneLx; }
    G4double GetSourcePlaneLy() const { return fSourcePlaneLy; }
    G4double GetThetaMaxDeg() const { return fThetaMaxRad / deg; }
    G4int GetFluxModelCode() const { return static_cast<G4int>(fFluxModel); }
    G4String GetFluxModelString() const { return ToString(fFluxModel).c_str(); }
    G4int GetMuonModeCode() const;
    G4bool GetGeomIntersectsCCD() const { return fGeomIntersectsCCD; }
    G4int GetMuonPDGCode() const { return fMuonPDGCode; }
    G4int GetMuonChargeSign() const { return fMuonChargeSign; }
    G4double GetMuonChargeRatio() const { return fMuonPlusToMinusRatio; }
    G4String GetMuonChargeMode() const { return fMuonChargeMode; }

    // Command handlers (invoked via the generic messenger)
    void SetMuonMode(const G4String& mode);
    void SetFluxModel(const G4String& model);
    void SetSourcePlaneZ(G4double value);
    void SetSourcePlaneLx(G4double value);
    void SetSourcePlaneLy(G4double value);
    void SetSourcePlaneAutoSize(G4bool value);
    void SetSourcePlaneMargin(G4double value);
    void SetEminGeV(G4double value);
    void SetEmaxGeV(G4double value);
    void SetMuonChargeMode(const G4String& mode);
    void SetMuonChargeRatio(G4double value);

  private:
    void DefineCommands();
    void GenerateForcedFootprint(G4Event* anEvent);
    void GenerateTierBFlux(G4Event* anEvent);
    void SelectMuonDefinition();
    std::pair<G4double, G4double> ComputeSourcePlaneSize() const;
    void EnsureSamplerReady();
    void MarkSamplerDirty();
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

    G4double fThetaMaxRad = 85.*deg;
    G4bool fUseCosmicMuons = true;
    G4String fMuonModeString = "forced_footprint";
    MuonFluxModelType fFluxModel = MuonFluxModelType::Guan2015;
    G4double fSourcePlaneZ = 8.0 * cm;
    G4double fSourcePlaneLx = 10.0 * cm;
    G4double fSourcePlaneLy = 10.0 * cm;
    G4bool fSourcePlaneAutoSize = false;
    G4double fSourcePlaneMargin = 1.0 * cm;
    G4double fRequestedEminGeV = 1.0;
    G4double fRequestedEmaxGeV = 1.0e4;
    G4double fEffectiveEminGeV = 1.0;
    G4double fEffectiveEmaxGeV = 1.0e4;
    G4double fEminGeV = 1.0;
    G4double fEmaxGeV = 1.0e4;
    G4String fMuonChargeMode = "fixedRatio";
    G4double fMuonPlusToMinusRatio = 1.25;

    G4double fZImpactPlane = 0.0; // cm reference plane for impact point
    G4double fMuonX0 = 0.0;
    G4double fMuonY0 = 0.0;
    G4double fMuonZ0 = 0.0;
    G4double fMuonXImp = 0.0;
    G4double fMuonYImp = 0.0;
    G4double fMuonZImp = 0.0;
    G4double fMuonTheta = 0.0;
    G4double fMuonPhi = 0.0;
    G4double fMuonEnergyGeV = 4.0;
    G4bool fUseFixedEnergy = true;
    G4double fMuonCosTheta = 1.0;
    G4double fSampledEnergyGeV = 0.0;
    G4int fMuonPDGCode = 13;
    G4int fMuonChargeSign = -1;
    G4double fEventWeight = 0.0;    // seconds represented by a single event
    G4double fEventLivetime = 0.0;  // identical to weight for now
    G4bool fGeomIntersectsCCD = false;

    MuonFluxSampler fFluxSampler;
    bool fSamplerDirty = true;
    G4double fCachedThetaMax = -1.0;
    G4double fCachedEminGeV = -1.0;
    G4double fCachedEmaxGeV = -1.0;
    MuonFluxModelType fCachedFluxModel = MuonFluxModelType::Guan2015;
    G4double fCachedPlaneZ = -1.0;
    G4double fCachedPlaneLx = -1.0;
    G4double fCachedPlaneLy = -1.0;
    bool fPrintedFluxInfo = false;
   
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
