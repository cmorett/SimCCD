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
#include "G4Exception.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <string>
#include <limits>
#include <algorithm>


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
  
  // === Vacuum definition == //
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  auto Vacuum = new G4Material("interGalactic", atomicNumber, massOfMole, density, kStateGas, temperature, pressure);
  
  
  // ===== Steel (density ~8.0 g/cm³)
  auto Steel = new G4Material("StainlessSteel", 8.0*g/cm3, 4);
  // Add elements with their fractional mass
    Steel->AddElement(elFe, 0.70);  // 70% Iron
    Steel->AddElement(elCr, 0.18);  // 18% Chromium
    Steel->AddElement(elNi, 0.10);  // 10% Nickel
    Steel->AddElement(elC,  0.02);  // 2% Carbon
    
  fWorldLength = 2.0*m;
  fRockLength = 0.80*fWorldLength;
  fRock2Length = 0.50*fWorldLength;
  fBarDiameter = 0.3*fRock2Length; 

  // Copper // 
  auto Cu = nistManager->FindOrBuildMaterial("G4_Cu");

  // Aluminum // 
  auto Al = nistManager->FindOrBuildMaterial("G4_Al");

  // AlN //
  G4Element* elAl = nistManager->FindOrBuildElement("Al");
  G4Element* elN = nistManager->FindOrBuildElement("N");
  G4Material* AlN = new G4Material("AluminumNitride", 3.26 * g/cm3, 2);
  AlN->AddElement(elAl, 1);
  AlN->AddElement(elN, 1);
     
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

  #include <fstream>
  #include <sstream>

  // ======== Load assembly.dae via Assimp and create one logical volume per mesh ========
  // ...existing code...
  {
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile("assets/cad/assembly.dae",
                               aiProcess_Triangulate |
                               aiProcess_JoinIdenticalVertices |
                               aiProcess_PreTransformVertices |
                               aiProcess_FlipWindingOrder);

    if (!scene) {
      std::string err = std::string("Failed to load assembly.dae: ") + importer.GetErrorString();
      G4Exception("B02DetectorConstruction::Construct", "AssimpImport",
                  FatalException, G4String(err));
    }

    G4cout << "Assimp: mNumMeshes = " << scene->mNumMeshes << G4endl;

    G4Material* defaultMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    const G4bool checkOverlaps = true;
    // OnShape exported units -> meters, convert to mm for Geant4
    const G4double meshScale = 1000.0 * mm;

    // debug bboxes (optional)
    for (unsigned int mdbg = 0; mdbg < scene->mNumMeshes; ++mdbg) {
      aiMesh* amdbg = scene->mMeshes[mdbg];
      if (!amdbg || amdbg->mNumVertices==0) continue;
      double minx=1e300, miny=1e300, minz=1e300;
      double maxx=-1e300, maxy=-1e300, maxz=-1e300;
      for (unsigned int v=0; v<amdbg->mNumVertices; ++v) {
        const aiVector3D &vv = amdbg->mVertices[v];
        minx = std::min(minx, (double)vv.x); maxx = std::max(maxx, (double)vv.x);
        miny = std::min(miny, (double)vv.y); maxy = std::max(maxy, (double)vv.y);
        minz = std::min(minz, (double)vv.z); maxz = std::max(maxz, (double)vv.z);
      }
      G4cout << "  mesh[" << mdbg << "] name=\"" 
             << (amdbg->mName.length?amdbg->mName.C_Str():"<no-name>")
             << "\" raw bbox x["<<minx<<","<<maxx<<"] y["<<miny<<","<<maxy<<"] z["<<minz<<","<<maxz<<"]"
             << G4endl;
      G4cout << "    scaled bbox (mm): x["<<(minx*meshScale/mm)<<" mm,"<<(maxx*meshScale/mm)<<" mm] "
             << " y["<<(miny*meshScale/mm)<<" mm,"<<(maxy*meshScale/mm)<<" mm] "
             << " z["<<(minz*meshScale/mm)<<" mm,"<<(maxz*meshScale/mm)<<" mm]"
             << G4endl;
    }

    // create a G4TessellatedSolid for each aiMesh and place it at origin (PreTransformVertices baked transforms)
    for (unsigned int m = 0; m < scene->mNumMeshes; ++m) {
      aiMesh* am = scene->mMeshes[m];
      if (!am || am->mNumFaces==0) continue;

      G4String meshName = (am->mName.length > 0) ? G4String(am->mName.C_Str()) :
                                                   G4String("mesh_") + G4String(std::to_string(m));
      G4TessellatedSolid* tess = new G4TessellatedSolid(meshName);

      for (unsigned int f = 0; f < am->mNumFaces; ++f) {
        const aiFace& face = am->mFaces[f];
        if (face.mNumIndices != 3) continue;
        const aiVector3D &v0 = am->mVertices[face.mIndices[0]];
        const aiVector3D &v1 = am->mVertices[face.mIndices[1]];
        const aiVector3D &v2 = am->mVertices[face.mIndices[2]];

        G4ThreeVector a(v0.x * meshScale, v0.y * meshScale, v0.z * meshScale);
        G4ThreeVector b(v1.x * meshScale, v1.y * meshScale, v1.z * meshScale);
        G4ThreeVector c(v2.x * meshScale, v2.y * meshScale, v2.z * meshScale);

        tess->AddFacet(new G4TriangularFacet(a, b, c, ABSOLUTE));
      }

      tess->SetSolidClosed(true);
      if (tess->GetNumberOfFacets() == 0) { delete tess; continue; }

      G4LogicalVolume* lv = new G4LogicalVolume(tess, defaultMat, meshName + "_lv");
      G4VisAttributes* vis = new G4VisAttributes(G4Colour(0.6,0.6,0.9));
      vis->SetForceSolid(true);
      lv->SetVisAttributes(vis);

      new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), lv, meshName, logicWorld, false, static_cast<G4int>(m), checkOverlaps);
    }
  }

  return physiWorld;
// ...existing code...
 
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


