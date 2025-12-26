//
// Detector construction with robust CAD handling and primitive fallback.
//

#include "B02DetectorConstruction.hh"

#include "B02BarSD.hh"
#include "SimCCDConfig.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ProductionCuts.hh"
#include "G4GenericMessenger.hh"
#include "G4MultiUnion.hh"
#include "G4RotationMatrix.hh"
#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4ThreeVector.hh"
#include "G4TriangularFacet.hh"
#include "G4Transform3D.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/scene.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <limits>

namespace {
namespace fs = std::filesystem;

fs::path TryCanonical(const fs::path& input) {
  std::error_code ec;
  fs::path result = fs::weakly_canonical(input, ec);
  return ec ? input : result;
}

std::size_t FixMeshOrientation(aiMesh& mesh, const aiVector3D& referencePoint) {
  if (!mesh.HasFaces() || mesh.mNumVertices == 0) {
    return 0;
  }

  std::size_t flipped = 0;
  for (unsigned int f = 0; f < mesh.mNumFaces; ++f) {
    aiFace& face = mesh.mFaces[f];
    if (face.mNumIndices != 3) {
      continue;
    }
    const aiVector3D& v0 = mesh.mVertices[face.mIndices[0]];
    const aiVector3D& v1 = mesh.mVertices[face.mIndices[1]];
    const aiVector3D& v2 = mesh.mVertices[face.mIndices[2]];

    aiVector3D normal = (v1 - v0) ^ (v2 - v0);
    const float norm2 = normal.SquareLength();
    if (norm2 < 1e-18f) {
      continue;  // Degenerate face
    }
    const aiVector3D toReference = referencePoint - v0;
    const float alignment = normal * toReference;
    if (alignment > 0.f) {
      std::swap(face.mIndices[1], face.mIndices[2]);
      ++flipped;
    }
  }
  return flipped;
}

std::size_t FixMeshOrientation(aiMesh& mesh) {
  aiVector3D centroid(0.f, 0.f, 0.f);
  if (mesh.mNumVertices > 0) {
    for (unsigned int i = 0; i < mesh.mNumVertices; ++i) {
      centroid += mesh.mVertices[i];
    }
    centroid /= static_cast<float>(mesh.mNumVertices);
  }
  return FixMeshOrientation(mesh, centroid);
}

G4VisAttributes* MakeVisAttributes(const G4Colour& color, bool solid = true) {
  auto* vis = new G4VisAttributes(color);
  vis->SetForceSolid(solid);
  vis->SetVisibility(true);
  return vis;
}
}  // namespace

B02DetectorConstruction::B02DetectorConstruction(GeometryOptions opts)
    : fOptions(std::move(opts)) {
  if (fOptions.geometryMode.empty()) {
    fOptions.geometryMode = "cad";
  }
  if (fOptions.cadMode.empty()) {
    fOptions.cadMode = "merged";
  }

  fOverburdenMessenger =
      new G4GenericMessenger(this, "/sim/overburden/", "Overburden configuration");
  fOverburdenMessenger->DeclareProperty("enable", fOverburdenEnabled,
                                        "Enable a slab above the detector to emulate roof/overburden.");
  fOverburdenMessenger->DeclarePropertyWithUnit("thickness", "cm", fOverburdenThickness,
                                                "Overburden thickness");
  fOverburdenMessenger->DeclarePropertyWithUnit("zTop", "cm", fOverburdenZTop,
                                                "Z position of the overburden top surface.");
  fOverburdenMessenger->DeclareProperty("material", fOverburdenMaterial,
                                        "G4 material name for the overburden (e.g. G4_CONCRETE).");

  fCutsMessenger = new G4GenericMessenger(this, "/sim/cuts/",
                                          "Production cuts for the CCD scoring region");
  fCutsMessenger->DeclarePropertyWithUnit("ccdGammaCut", "cm", fCCDGammaCut,
                                          "Production cut for gammas in the CCD region (0 => default).");
  fCutsMessenger->DeclarePropertyWithUnit("ccdElectronCut", "cm", fCCDElectronCut,
                                          "Production cut for electrons in the CCD region (0 => default).");
  fCutsMessenger->DeclarePropertyWithUnit("ccdPositronCut", "cm", fCCDPositronCut,
                                          "Production cut for positrons in the CCD region (0 => default).");

  fCCDMessenger = new G4GenericMessenger(this, "/sim/ccd/", "CCD-specific controls");
  fCCDMessenger->DeclarePropertyWithUnit("maxStep", "cm", fCCDMaxStep,
                                         "Optional step limit inside the CCD scoring slab (0 disables).");
}

B02DetectorConstruction::~B02DetectorConstruction() {
  delete fOverburdenMessenger;
  delete fCutsMessenger;
  delete fCCDMessenger;
}

void B02DetectorConstruction::DefineMaterials() {
  auto* nist = G4NistManager::Instance();
  nist->FindOrBuildMaterial("G4_AIR");
  nist->FindOrBuildMaterial("G4_Galactic");
  nist->FindOrBuildMaterial("G4_Si");
  nist->FindOrBuildMaterial("G4_Al");
  nist->FindOrBuildMaterial("G4_Cu");
  nist->FindOrBuildMaterial("G4_Fe");
  nist->FindOrBuildMaterial("G4_POLYETHYLENE");
}

G4VPhysicalVolume* B02DetectorConstruction::Construct() {
  DefineMaterials();
  BuildWorld();

  if (fOptions.geometryMode == "primitive") {
    fPrimaryScoring = BuildPrimitiveGeometry();
  } else {
    fPrimaryScoring = BuildCadGeometry();
  }

  BuildOverburden();
  BuildCCDOverlay();
  ConfigureCCDRegion();
  return fWorldPhysical;
}

G4VPhysicalVolume* B02DetectorConstruction::BuildWorld() {
  auto* nist = G4NistManager::Instance();
  auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");

  auto* solidWorld =
      new G4Box("World", 0.5 * fWorldLength, 0.5 * fWorldLength, 0.5 * fWorldLength);
  fWorldLogical = new G4LogicalVolume(solidWorld, vacuum, "WorldLogical");
  fWorldLogical->SetVisAttributes(nullptr);

  fWorldPhysical = new G4PVPlacement(nullptr, {}, fWorldLogical, "World", nullptr, false, 0,
                                     false);
  return fWorldPhysical;
}

G4LogicalVolume* B02DetectorConstruction::BuildPrimitiveGeometry() {
  auto* nist = G4NistManager::Instance();
  auto* si = nist->FindOrBuildMaterial("G4_Si");

  constexpr G4double halfX = 25.0 * mm;
  constexpr G4double halfY = 25.0 * mm;
  constexpr G4double halfZ = 0.5 * mm;

  auto* solid = new G4Box("PrimitiveSensor", halfX, halfY, halfZ);
  auto* logical = new G4LogicalVolume(solid, si, "PrimitiveSensorLogical");
  logical->SetVisAttributes(MakeVisAttributes(G4Colour(0.2, 0.6, 0.8)));

  new G4PVPlacement(nullptr, {}, logical, "PrimitiveSensorPV", fWorldLogical, false, 0,
                    fOptions.checkOverlaps);

  RegisterSensitiveVolume(logical);
  return logical;
}

G4LogicalVolume* B02DetectorConstruction::BuildCadGeometry() {
  const std::string cadPath = ResolveCadPath();
  Assimp::Importer importer;
  importer.SetPropertyInteger(AI_CONFIG_PP_SBP_REMOVE,
                              aiPrimitiveType_POINT | aiPrimitiveType_LINE);

  unsigned int flags = aiProcess_Triangulate | aiProcess_JoinIdenticalVertices |
                       aiProcess_PreTransformVertices | aiProcess_ImproveCacheLocality |
                       aiProcess_RemoveRedundantMaterials | aiProcess_FindDegenerates |
                       aiProcess_FindInvalidData | aiProcess_OptimizeMeshes |
                       aiProcess_GenNormals;
#ifdef aiProcess_FixInfacingNormals
  flags |= aiProcess_FixInfacingNormals;
#endif

  const aiScene* scene = importer.ReadFile(cadPath, flags);
  if (!scene || !scene->HasMeshes()) {
    std::ostringstream msg;
    msg << "Failed to load CAD file at '" << cadPath << "': "
        << importer.GetErrorString();
    G4Exception("B02DetectorConstruction::BuildCadGeometry", "AssimpImport",
                FatalException, msg.str().c_str());
  }

  aiVector3D assemblyCentroid(0.f, 0.f, 0.f);
  std::size_t totalVertices = 0;
  for (unsigned int i = 0; i < scene->mNumMeshes; ++i) {
    if (const aiMesh* mesh = scene->mMeshes[i]) {
      for (unsigned int v = 0; v < mesh->mNumVertices; ++v) {
        assemblyCentroid += mesh->mVertices[v];
      }
      totalVertices += mesh->mNumVertices;
    }
  }
  if (totalVertices > 0) {
    assemblyCentroid /= static_cast<float>(totalVertices);
  }

  std::size_t flippedFaces = 0;
  for (unsigned int i = 0; i < scene->mNumMeshes; ++i) {
    if (scene->mMeshes[i]) {
      flippedFaces += FixMeshOrientation(*scene->mMeshes[i]);
    }
  }
  if (totalVertices > 0) {
    for (unsigned int i = 0; i < scene->mNumMeshes; ++i) {
      if (scene->mMeshes[i]) {
        flippedFaces += FixMeshOrientation(*scene->mMeshes[i], assemblyCentroid);
      }
    }
  }
  G4cout << "[geometry] Loaded CAD scene with " << scene->mNumMeshes
         << " mesh(es) from: " << cadPath << G4endl;
  G4cout << "[geometry] Flipped " << flippedFaces
         << " faces to enforce outward normals." << G4endl;
  G4cout << "[geometry] CAD mode: " << fOptions.cadMode << G4endl;

  auto* nist = G4NistManager::Instance();
  auto* defaultMat = nist->FindOrBuildMaterial("G4_Al");
  const G4double meshScale = 1000.0 * mm;  // meters -> mm

  std::vector<G4LogicalVolume*> cadVolumes;
  if (fOptions.cadMode == "parts") {
    cadVolumes = BuildCadParts(*scene, defaultMat, meshScale);
  } else if (fOptions.cadMode == "tessellated") {
    auto* merged = BuildCadMergedTessellated(*scene, defaultMat, meshScale);
    if (merged) {
      cadVolumes.push_back(merged);
    }
  } else {
    auto* merged = BuildCadBoundingUnion(*scene, defaultMat, meshScale);
    if (merged) {
      cadVolumes.push_back(merged);
    }
  }

  if (cadVolumes.empty()) {
    G4Exception("B02DetectorConstruction::BuildCadGeometry", "CadEmpty", FatalException,
                "No CAD volumes were created.");
  }

  for (std::size_t i = 0; i < cadVolumes.size(); ++i) {
    auto* lv = cadVolumes[i];
    auto name = lv->GetName();
    auto* pv = new G4PVPlacement(nullptr, {}, lv, name + "_pv", fWorldLogical, false,
                                 static_cast<G4int>(i), fOptions.checkOverlaps);
    if (fOptions.checkOverlaps) {
      pv->CheckOverlaps(fOptions.overlapSamples, fOptions.overlapTolerance, true, 1);
    }
  }

  return cadVolumes.front();
}

G4LogicalVolume* B02DetectorConstruction::BuildCadMergedTessellated(const aiScene& scene,
                                                                    G4Material* material,
                                                                    G4double meshScale) {
  auto* solid = new G4TessellatedSolid("CADAssemblyMerged");

  for (unsigned int m = 0; m < scene.mNumMeshes; ++m) {
    const aiMesh* mesh = scene.mMeshes[m];
    if (!mesh) {
      continue;
    }
    for (unsigned int f = 0; f < mesh->mNumFaces; ++f) {
      const aiFace& face = mesh->mFaces[f];
      if (face.mNumIndices != 3) {
        continue;
      }
      const aiVector3D& v0 = mesh->mVertices[face.mIndices[0]];
      const aiVector3D& v1 = mesh->mVertices[face.mIndices[1]];
      const aiVector3D& v2 = mesh->mVertices[face.mIndices[2]];

      const G4ThreeVector a(v0.x * meshScale, v0.y * meshScale, v0.z * meshScale);
      const G4ThreeVector b(v1.x * meshScale, v1.y * meshScale, v1.z * meshScale);
      const G4ThreeVector c(v2.x * meshScale, v2.y * meshScale, v2.z * meshScale);

      solid->AddFacet(new G4TriangularFacet(a, b, c, ABSOLUTE));
    }
  }

  solid->SetSolidClosed(true);
  if (solid->GetNumberOfFacets() == 0) {
    delete solid;
    return nullptr;
  }

  auto* lv = new G4LogicalVolume(solid, material, "CADAssemblyMerged_lv");
  lv->SetVisAttributes(MakeVisAttributes(G4Colour(0.6, 0.6, 0.8)));
  return lv;
}

void B02DetectorConstruction::BuildOverburden() {
  if (!fOverburdenEnabled || fOverburdenThickness <= 0.0) {
    return;
  }

  auto* nist = G4NistManager::Instance();
  G4Material* material = nist->FindOrBuildMaterial(fOverburdenMaterial);
  if (!material) {
    G4Exception("B02DetectorConstruction::BuildOverburden", "OverburdenMatMissing",
                JustWarning,
                ("Could not find material '" + fOverburdenMaterial +
                 "'. Overburden disabled.")
                    .c_str());
    return;
  }

  const G4double halfX = 0.45 * fWorldLength;
  const G4double halfY = 0.45 * fWorldLength;
  const G4double halfZ = 0.5 * fOverburdenThickness;
  const G4double zCenter = fOverburdenZTop - halfZ;
  const G4double halfWorld = 0.5 * fWorldLength;
  if (std::abs(zCenter) + halfZ > halfWorld) {
    G4Exception("B02DetectorConstruction::BuildOverburden", "OverburdenOutsideWorld",
                JustWarning,
                "Overburden exceeds world bounds; skipping placement.");
    return;
  }

  auto* solid = new G4Box("OverburdenSolid", halfX, halfY, halfZ);
  auto* logical = new G4LogicalVolume(solid, material, "OverburdenLogical");
  logical->SetVisAttributes(MakeVisAttributes(G4Colour(0.5, 0.5, 0.5), false));
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zCenter), logical, "OverburdenPV",
                    fWorldLogical, false, 0, fOptions.checkOverlaps);
}

void B02DetectorConstruction::BuildCCDOverlay() {
  // Create a thin scoring box representing the CCD active silicon. This is intentionally
  // small so that misses are not counted as CCD hits even when the CAD model is large.
  // Dimensions match the legacy px/py footprint (1.5 cm x 1.5 cm) and a thin thickness.
  if (fCCDOverlayLogical || fOptions.geometryMode == "primitive") {
    return;
  }

  auto* nist = G4NistManager::Instance();
  auto* si = nist->FindOrBuildMaterial("G4_Si");

  const G4double halfX = 0.75 * cm;
  const G4double halfY = 0.75 * cm;
  const G4double halfZ = 0.025 * cm;  // 0.5 mm total thickness

  auto* solid = new G4Box("CCDOverlay", halfX, halfY, halfZ);
  fCCDOverlayLogical = new G4LogicalVolume(solid, si, "CCDOverlayLogical");
  fCCDOverlayLogical->SetVisAttributes(MakeVisAttributes(G4Colour(0.8, 0.2, 0.2), false));

  G4LogicalVolume* mother = fWorldLogical;
  if (fPrimaryScoring) {
    mother = fPrimaryScoring;
  }

  new G4PVPlacement(nullptr, {}, fCCDOverlayLogical, "CCDOverlayPV", mother, false, 0,
                    false);
  RegisterSensitiveVolume(fCCDOverlayLogical);
}

G4LogicalVolume* B02DetectorConstruction::BuildCadBoundingUnion(const aiScene& scene,
                                                                G4Material* material,
                                                                G4double meshScale) {
  auto* multi = new G4MultiUnion("CADAssemblyMerged");
  std::size_t added = 0;

  for (unsigned int m = 0; m < scene.mNumMeshes; ++m) {
    const aiMesh* mesh = scene.mMeshes[m];
    if (!mesh || mesh->mNumVertices == 0) {
      continue;
    }

    double minx = std::numeric_limits<double>::max();
    double miny = std::numeric_limits<double>::max();
    double minz = std::numeric_limits<double>::max();
    double maxx = -std::numeric_limits<double>::max();
    double maxy = -std::numeric_limits<double>::max();
    double maxz = -std::numeric_limits<double>::max();

    for (unsigned int v = 0; v < mesh->mNumVertices; ++v) {
      const aiVector3D& vv = mesh->mVertices[v];
      minx = std::min(minx, static_cast<double>(vv.x));
      miny = std::min(miny, static_cast<double>(vv.y));
      minz = std::min(minz, static_cast<double>(vv.z));
      maxx = std::max(maxx, static_cast<double>(vv.x));
      maxy = std::max(maxy, static_cast<double>(vv.y));
      maxz = std::max(maxz, static_cast<double>(vv.z));
    }

    const double dx = (maxx - minx) * meshScale;
    const double dy = (maxy - miny) * meshScale;
    const double dz = (maxz - minz) * meshScale;
    if (dx <= 0 || dy <= 0 || dz <= 0) {
      continue;
    }

    const auto center = G4ThreeVector(0.5 * (minx + maxx) * meshScale,
                                      0.5 * (miny + maxy) * meshScale,
                                      0.5 * (minz + maxz) * meshScale);
    auto* box = new G4Box("cad_bbox", 0.5 * dx, 0.5 * dy, 0.5 * dz);
    multi->AddNode(*box, G4Transform3D(G4RotationMatrix(), center));
    ++added;
  }

  if (added == 0) {
    delete multi;
    return nullptr;
  }

  multi->Voxelize();
  auto* lv = new G4LogicalVolume(multi, material, "CADAssemblyMerged_lv");
  lv->SetVisAttributes(MakeVisAttributes(G4Colour(0.6, 0.6, 0.8)));
  return lv;
}

std::vector<G4LogicalVolume*> B02DetectorConstruction::BuildCadParts(
    const aiScene& scene, G4Material* material, G4double meshScale) {
  std::vector<G4LogicalVolume*> result;
  result.reserve(scene.mNumMeshes);

  for (unsigned int m = 0; m < scene.mNumMeshes; ++m) {
    const aiMesh* mesh = scene.mMeshes[m];
    if (!mesh) {
      continue;
    }

    G4String meshName =
        (mesh->mName.length > 0) ? G4String(mesh->mName.C_Str()) : G4String("mesh_") +
                                                                    G4String(std::to_string(m));
    auto* tess = new G4TessellatedSolid(meshName);

    for (unsigned int f = 0; f < mesh->mNumFaces; ++f) {
      const aiFace& face = mesh->mFaces[f];
      if (face.mNumIndices != 3) {
        continue;
      }
      const aiVector3D& v0 = mesh->mVertices[face.mIndices[0]];
      const aiVector3D& v1 = mesh->mVertices[face.mIndices[1]];
      const aiVector3D& v2 = mesh->mVertices[face.mIndices[2]];

      const G4ThreeVector a(v0.x * meshScale, v0.y * meshScale, v0.z * meshScale);
      const G4ThreeVector b(v1.x * meshScale, v1.y * meshScale, v1.z * meshScale);
      const G4ThreeVector c(v2.x * meshScale, v2.y * meshScale, v2.z * meshScale);

      tess->AddFacet(new G4TriangularFacet(a, b, c, ABSOLUTE));
    }

    tess->SetSolidClosed(true);
    if (tess->GetNumberOfFacets() == 0) {
      delete tess;
      continue;
    }

    auto* lv = new G4LogicalVolume(tess, material, meshName + "_lv");
    lv->SetVisAttributes(MakeVisAttributes(G4Colour(0.6, 0.6, 0.8)));
    result.push_back(lv);
  }

  return result;
}

void B02DetectorConstruction::RegisterSensitiveVolume(G4LogicalVolume* lv) {
  if (!lv) {
    return;
  }
  if (std::find(fSensitiveVolumes.begin(), fSensitiveVolumes.end(), lv) ==
      fSensitiveVolumes.end()) {
    fSensitiveVolumes.push_back(lv);
  }
  if (!fPrimaryScoring) {
    fPrimaryScoring = lv;
  }
}

void B02DetectorConstruction::ConstructSDandField() {
  if (fSensitiveVolumes.empty()) {
    return;
  }

  auto* sdManager = G4SDManager::GetSDMpointer();
  if (!fBarSD) {
    fBarSD = new B02BarSD("BarSD");
    sdManager->AddNewDetector(fBarSD);
  }

  for (auto* lv : fSensitiveVolumes) {
    if (lv) {
      lv->SetSensitiveDetector(fBarSD);
    }
  }
}

bool B02DetectorConstruction::IsSensitiveVolume(const G4LogicalVolume* volume) const {
  return volume &&
         std::find(fSensitiveVolumes.begin(), fSensitiveVolumes.end(), volume) !=
             fSensitiveVolumes.end();
}

bool B02DetectorConstruction::IsCCDScoringVolume(const G4LogicalVolume* volume) const {
  return volume == fCCDOverlayLogical || volume == fPrimaryScoring;
}

G4double B02DetectorConstruction::GetCCDThickness() const {
  if (!fCCDOverlayLogical) {
    return 0.0;
  }
  if (const auto* box = dynamic_cast<const G4Box*>(fCCDOverlayLogical->GetSolid())) {
    return 2.0 * box->GetZHalfLength();
  }
  return 0.0;
}

void B02DetectorConstruction::ConfigureCCDRegion() {
  if (!fCCDOverlayLogical) {
    return;
  }

  auto* region = new G4Region("CCDRegion");
  region->AddRootLogicalVolume(fCCDOverlayLogical);

  if (fCCDGammaCut > 0.0 || fCCDElectronCut > 0.0 || fCCDPositronCut > 0.0) {
    auto* cuts = new G4ProductionCuts();
    if (fCCDGammaCut > 0.0) {
      cuts->SetProductionCut(fCCDGammaCut, "gamma");
    }
    if (fCCDElectronCut > 0.0) {
      cuts->SetProductionCut(fCCDElectronCut, "e-");
    }
    if (fCCDPositronCut > 0.0) {
      cuts->SetProductionCut(fCCDPositronCut, "e+");
    }
    region->SetProductionCuts(cuts);
  }

  if (fCCDMaxStep > 0.0) {
    auto* limits = new G4UserLimits(fCCDMaxStep);
    fCCDOverlayLogical->SetUserLimits(limits);
  }
}

std::string B02DetectorConstruction::ResolveAssetsDir() const {
  std::vector<fs::path> candidates;
  if (!fOptions.assetsDir.empty()) {
    candidates.emplace_back(std::string(fOptions.assetsDir));
  }
  if (const char* env = std::getenv("SIMCCD_ASSETS_DIR")) {
    candidates.emplace_back(env);
  }
  if (simccd::kBinaryAssetsDir && std::strlen(simccd::kBinaryAssetsDir) > 0) {
    candidates.emplace_back(simccd::kBinaryAssetsDir);
  }
  if (simccd::kSourceAssetsDir && std::strlen(simccd::kSourceAssetsDir) > 0) {
    candidates.emplace_back(simccd::kSourceAssetsDir);
  }
  candidates.emplace_back(fs::current_path() / "assets");

  for (const auto& raw : candidates) {
    const fs::path normalized = TryCanonical(raw);
    if (fs::exists(normalized / "cad")) {
      return normalized.string();
    }
  }
  return {};
}

std::string B02DetectorConstruction::ResolveCadPath() const {
  std::vector<fs::path> candidates;
  std::vector<std::string> checked;

  if (!fOptions.cadFile.empty()) {
    candidates.emplace_back(std::string(fOptions.cadFile));
  }
  if (const char* env = std::getenv("SIMCCD_CAD_FILE")) {
    candidates.emplace_back(env);
  }

  const std::string assetsDir = ResolveAssetsDir();
  if (!assetsDir.empty()) {
    candidates.emplace_back(fs::path(assetsDir) / simccd::kCadRelativePath);
  }
  if (simccd::kBinaryAssetsDir && std::strlen(simccd::kBinaryAssetsDir) > 0) {
    candidates.emplace_back(
        fs::path(simccd::kBinaryAssetsDir) / simccd::kCadRelativePath);
  }
  if (simccd::kSourceAssetsDir && std::strlen(simccd::kSourceAssetsDir) > 0) {
    candidates.emplace_back(
        fs::path(simccd::kSourceAssetsDir) / simccd::kCadRelativePath);
  }
  candidates.emplace_back(fs::current_path() / "assets" / "cad" / "assembly.dae");

  for (const auto& raw : candidates) {
    const fs::path normalized = TryCanonical(raw);
    checked.push_back(normalized.string());
    if (fs::exists(normalized)) {
      return normalized.string();
    }
  }

  std::ostringstream msg;
  msg << "Could not locate CAD file 'assembly.dae'. Checked paths:";
  for (const auto& pathStr : checked) {
    msg << "\n  - " << pathStr;
  }
  msg << "\nSet SIMCCD_CAD_FILE or SIMCCD_ASSETS_DIR, or pass --cad-file/--assets-dir.";

  G4Exception("B02DetectorConstruction::ResolveCadPath", "CadNotFound", FatalException,
              msg.str().c_str());
  return {};
}
