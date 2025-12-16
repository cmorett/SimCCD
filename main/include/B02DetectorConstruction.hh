//

#ifndef B02DetectorConstruction_h
#define B02DetectorConstruction_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"

#include <string>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class B02BarSD;
struct aiScene;

struct GeometryOptions {
  G4String geometryMode = "cad";   // "cad" or "primitive"
  G4String cadMode = "merged";     // "merged" or "parts"
  G4String cadFile;
  G4String assetsDir;
  G4bool checkOverlaps = true;
  G4double overlapTolerance = 0.05 * mm;
  G4int overlapSamples = 1000;
  G4bool verbose = false;
};

class B02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    explicit B02DetectorConstruction(GeometryOptions opts = {});
   ~B02DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    G4double GetWorldFullLength() const { return fWorldLength; }
    G4LogicalVolume* GetScoringVolume() const { return fPrimaryScoring; }
    const GeometryOptions& GetOptions() const { return fOptions; }
    bool IsSensitiveVolume(const G4LogicalVolume* volume) const;

  private:
    void DefineMaterials();
    G4VPhysicalVolume* BuildWorld();
    G4LogicalVolume* BuildPrimitiveGeometry();
    G4LogicalVolume* BuildCadGeometry();
    G4LogicalVolume* BuildCadMergedTessellated(const aiScene& scene, G4Material* material, G4double scale);
    G4LogicalVolume* BuildCadBoundingUnion(const aiScene& scene, G4Material* material, G4double scale);
    std::vector<G4LogicalVolume*> BuildCadParts(const aiScene& scene, G4Material* material, G4double scale);
    void RegisterSensitiveVolume(G4LogicalVolume* lv);
    std::string ResolveCadPath() const;
    std::string ResolveAssetsDir() const;

    GeometryOptions fOptions;
    G4double fWorldLength = 4.0 * m;
    G4LogicalVolume* fWorldLogical = nullptr;
    G4VPhysicalVolume* fWorldPhysical = nullptr;
    std::vector<G4LogicalVolume*> fSensitiveVolumes;
    G4LogicalVolume* fPrimaryScoring = nullptr;
    B02BarSD* fBarSD = nullptr;
};

#endif
