//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B02DetectorConstruction.hh"
#include "B02DetectorMessenger.hh"
#include "B02BarSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Paraboloid.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include "G4GenericMessenger.hh"
#include "G4VisAttributes.hh"
#include "G4VSensitiveDetector.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "CADMesh.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



B02DetectorConstruction::B02DetectorConstruction()

{
  DefineCommands();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02DetectorConstruction::~B02DetectorConstruction()
{
delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double B02DetectorConstruction::GetWorldFullLength() const
{return fWorldLength;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//added 


G4VPhysicalVolume* B02DetectorConstruction::Construct()
{ 
  //------------------------------------------------------ volumes
  //--------- Sizes of the principal geometrical components (solids)
  
  // new introductions to better comprehension
  DefineMaterials();
  
  // ============ Materials ======================== //
  auto Air = G4Material::GetMaterial("G4_AIR");
  auto LAr = new G4Material("LAr",   18.,  39.95*g/mole,  1.393*g/cm3);
  //LAr->GetIonisation()->SetBirksConstant(0.069*cm/MeV); //LarSoft documentation C. Zhang
  
  auto nistManager = G4NistManager::Instance();
  auto Si = nistManager->FindOrBuildMaterial("G4_Si");
  Si->GetIonisation()->SetBirksConstant(0.09*cm/MeV); // El valor de 0.1 parece ajustar bastante bien el espectro pero se necesita mas estudio
  
  
  //auto Si = new G4Material("Si",   14.,  28.0849*g/mole,  2.3396*g/cm3);
  //auto Si = G4Material::GetMaterial("G4_Si");
  //G4Element *elSi = new G4Element("Silicon", "Si", 14., 28.0849*g/mole);
  
  // === Vacuum definition == //
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  auto Vacuum = new G4Material("interGalactic", atomicNumber, massOfMole, density, kStateGas, temperature, pressure);
  // ======================= //
  
  
   //Steel (density ~8.0 g/cm³)
  auto Steel = new G4Material("StainlessSteel", 8.0*g/cm3, 4);
  // Add elements with their fractional mass
    Steel->AddElement(elFe, 0.70);  // 70% Iron
    Steel->AddElement(elCr, 0.18);  // 18% Chromium
    Steel->AddElement(elNi, 0.10);  // 10% Nickel
    Steel->AddElement(elC,  0.02);  // 2% Carbon
    
  fWorldLength = 10*cm;
  fRockLength = 0.80*fWorldLength;
  fRock2Length = 0.50*fWorldLength;
  fBarDiameter = 0.3*fRock2Length; 
     
  // World
  //------------------------------ 
  G4double HalfWorldLength = 1.*fWorldLength;
  //G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  //G4cout << "Computed tolerance = "
    //     << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
    //     << " mm" << G4endl;
  solidWorld= new G4Box("World",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  //logicWorld= new G4LogicalVolume( solidWorld, Vacuum, "World", 0, 0, 0);
  logicWorld= new G4LogicalVolume( solidWorld, Vacuum, "World");
  physiWorld = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.), logicWorld, "World", 0, false, 0, true);   
  // =========================================== //
  

  // =============== Constructor of CCD (non active volume) ===================== //

  //G4double pixel_size = 0.0015; // cm

  // G4double XLength = 1.917; // cm
  // G4double YLength = 1.587; // cm
  // G4double ZLength = 0.0725; // cm
  
  // G4double XLength = 0.6; // cm
  // G4double YLength = 0.7875; // cm
  // G4double ZLength = 0.0725*1.; // cm

  // ========= Para el ICN ======================= //
  // G4double XLength = 300 * pixel_size; // cm
  // G4double YLength = 529 * pixel_size; // cm (Debe ser la dimensión con mayor tamaño)
  // G4double ZLength = 0.0725*1.; // cm

  //G4double XLength = 250 * pixel_size; // cm
  //G4double YLength = 529 * pixel_size; // cm (Debe ser la dimensión con mayor tamaño)
  //G4double ZLength = 0.0725; // cm
  // ============================================================================ //

  // ================= Para CONNIE ================== //
  // G4double XLength = 420 * pixel_size; // cm
  // G4double YLength = 1022 * pixel_size; // cm (Debe ser la dimensión con mayor tamaño)
  // G4double ZLength = 0.068; // cm

  // G4double XLength = 420 * pixel_size; // cm
  // G4double YLength = 700 * pixel_size; // cm (Debe ser la dimensión con mayor tamaño)
  // G4double ZLength = 0.068; // cm

  // G4double XLength = 420 * pixel_size; // cm
  // G4double YLength = 600 * pixel_size; // cm (Debe ser la dimensión con mayor tamaño)
  // G4double ZLength = 0.068; // cm

  
  //Sibox = new G4Box("ccd", HalfWorldLength, HalfWorldLength, HalfWorldLength);
  //Sibox = new G4Box("CCD", 0.5*XLength*cm, 0.5*YLength*cm, 0.5*ZLength*cm);
  //SiLogic = new G4LogicalVolume(Sibox, Si, "CCD", 0, 0, 0);
  // new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), SiLogic, "CCD", vacuum_stLV, false, 0,true);
  //new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), SiLogic, "CCD", logicWorld, false, 0,true);

  //fSiLogic = SiLogic;

  // ======================================================== //

    // ======== CAD enclosure (Path A: CADMesh • STL in mm) ========
  {
    // Materials for the enclosure. Adjust to your real hardware if needed:
    // - If your lids/shell/tube are stainless, use the custom "StainlessSteel" you already create above.
    // - Otherwise, this uses Aluminum as a sensible default.
    auto nist = G4NistManager::Instance();
    G4Material* enclosureMat = nist->FindOrBuildMaterial("G4_Al"); // change to StainlessSteel if appropriate

    // Optional: world bigger, if needed for large shells (your current fWorldLength is fine; bump if you outgrow it)
    // fWorldLength = 1.0*m;  // if you decide to enlarge world, also rebuild solidWorld with the new size.

    const G4bool checkOverlaps = true;

    struct Part {
      const char* file;
      const char* name;
    } parts[] = {
      {"assets/cad/back_lid.stl",  "back_lid"},
      {"assets/cad/front_lid.stl", "front_lid"},
      {"assets/cad/left_lid.stl",  "left_lid"},
      {"assets/cad/right_lid.stl", "right_lid"},
      {"assets/cad/outer_shell.stl","outer_shell"},
      {"assets/cad/tube.stl",      "tube"},
    };

    // Because you exported from Onshape in mm **and** the parts are already in assembly position,
    // we place each one at the world origin with no extra rotation/offset.
    // If anything looks misaligned later, we can switch to a single OBJ-with-groups export, or apply per-part transforms.
    for (const auto& p : parts) {
      auto mesh  = CADMesh::TessellatedMesh::FromSTL(p.file);
      // mesh->SetScale(1.0);  // mm → mm (no scale needed). Uncomment if you ever export in different units.
      // mesh->SetOffset(0.0*mm, 0.0*mm, 0.0*mm); // not needed if assembly coords are baked-in

      auto solid = mesh->GetSolid(); // G4TessellatedSolid*
      auto lv    = new G4LogicalVolume(solid, enclosureMat, p.name);

      // Nice-to-have: make it visible and solid in the viewer
      auto vis = new G4VisAttributes(G4Colour(0.8,0.8,0.9));
      vis->SetForceSolid(true);
      lv->SetVisAttributes(vis);

      new G4PVPlacement(
        /*pRot=*/nullptr,
        /*tlate=*/G4ThreeVector(),  // (0,0,0)
        /*pCurrentLogical=*/lv,
        /*pName=*/p.name,
        /*pMotherLogical=*/logicWorld,
        /*pMany=*/false,
        /*pCopyNo=*/0,
        /*checkOverlaps=*/checkOverlaps
      );
    }
  }

 return physiWorld;
 
}


void B02DetectorConstruction::DefineMaterials()
{
 //------------------------------------------------------ materials
  
  
 //Some elements
    elFe = new G4Element("Iron", "Fe", 26., 55.85*g/mole);
    elCr = new G4Element("Chromium", "Cr", 24., 51.9961*g/mole);
    elNi = new G4Element("Nickel", "Ni", 28., 58.69*g/mole);
    elC  = new G4Element("Carbon", "C", 6., 12.01*g/mole);

//Some Nist materials
  
   auto nistManager = G4NistManager::Instance();

  // Air
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  //H2O
  nistManager->FindOrBuildMaterial("G4_WATER");
  
  //concrete 
  nistManager->FindOrBuildMaterial("G4_CONCRETE");
  
  //polyethylene
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  
  
  //G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
 
}

void B02DetectorConstruction::SetHigh(G4double val)
{
  if (!Sibox) {
      G4cerr << "Detector has not yet been constructed." << G4endl;
      return;
  }
  // zhigh = val;
  // SolidLAr-> SetZHalfLength(0.5*zhigh);
  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void B02DetectorConstruction::DefineCommands()
{
  
  // detectorMessenger = new B02DetectorMessenger(this);
  fMessenger = new G4GenericMessenger(this, "/detector/", "B02DetectorConstruction");
  
  fMessenger->DeclareMethodWithUnit("zhigh", "cm",  &B02DetectorConstruction::SetHigh, "high cylinder");
  
}


