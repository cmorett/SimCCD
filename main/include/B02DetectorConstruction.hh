//

#ifndef B02DetectorConstruction_h
#define B02DetectorConstruction_h 1

#include "B02DetectorMessenger.hh"
#include "B02SteppingAction.hh"
#include "B02BarSD.hh"
#include "B02DetectorConstruction.hh"

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
//Required explicitly since G4 version 10.0
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "tls.hh"
#include "G4GenericMessenger.hh"
#include "G4SDManager.hh"

class G4Box;
class G4Tubs;
class G4Paraboloid;
class G4Orb;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;
class B0DetectorMessenger;
class G4Element;
class G4VSensitiveDetector;

class B02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    B02DetectorConstruction();
   ~B02DetectorConstruction() override;

  public:

    G4VPhysicalVolume* Construct() override;
    void DefineMaterials();
  
    G4double GetWorldFullLength() const; 

    void setBarMaterial  (G4String);
    void SetCheckOverlaps(G4bool );
    
    // new 
    void SetHigh(G4double val);
    G4double GetZHalfLength() 
    { return zhigh;}
    
    G4LogicalVolume *GetScoringVolume() const { return fSiLogic; }

    

  private:
    
    // methods
    
    G4VPhysicalVolume* DefineVolumes();

    // World volume
    G4Box*             solidWorld;    // pointer to the solid envelope 
    G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
    G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope

    //active cylinder                 // pointer to the solid Cylinder
    G4Box*             Sibox;     // pointer to the solid Cylinder
    G4LogicalVolume*   SiLogic;     // pointer to the logical Cylinder
    G4VPhysicalVolume* PhysiLAr;     // pointer to the physical Cylinder 


    B02DetectorMessenger* detectorMessenger;  // pointer to the Messenger

    G4double fWorldLength;            // Full length of the world volume
    G4double fRockLength;             // Full length of the Rock volume
    G4double fRock2Length;            // Full length of the Rock2 volume
    G4double fBarDiameter;           // Full diameter of the Bar volume

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    
    
    void DefineCommands();
    
    G4GenericMessenger *fMessenger = nullptr; // implementation Messenger 
    
    //shielding
    G4Box             *xwallL1, *xwallL2, *xwallP1, *xwallP2, *yf1, *yf2, *yfP, *zroofS, *yfS, *yb1, *xwallL1N, *xwallL2N, *xwallP1N, *xwallP2N, *yf3; 
    G4LogicalVolume   *ywallLV, *xwallLV1, *xwallLV2, *xwallLVP1, *xwallLVP2, *yfLV1, *yfLV2, *yfLVP, *zroofLVS, *yfLVS, *ybLV, *xwallLV1N, *xwallLV2N, *xwallLVP1N, *xwallLVP2N, *yfLV3;
    
    //non active cylinder
    G4Tubs* cylinder;
    G4LogicalVolume* cyLV;
    
    // parameter      
    G4double zhigh = 123.0*cm;
    
    // materials and elements
    G4Element *elFe, *elCr, *elNi, *elC;
    
    //connect with stepping action clss 
    B02DetectorConstruction* fB02DetectorConstruction;
    
    // G4LogicalVolume *fLogicLAr;
    G4LogicalVolume *fSiLogic;
     
};

// inline functions
/*
inline const G4VPhysicalVolume* B02DetectorConstruction::GetActiveLAr() const {
  return physiLAr;

}
*/


#endif

