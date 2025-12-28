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
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"

#include <assimp/Importer.hpp>
#include <assimp/config.h>
#include <assimp/postprocess.h>
#include <assimp/scene.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

bool NodeNameMatches(const aiNode* node, const std::string& target) {
  if (!node) {
    return false;
  }
  const std::string name = node->mName.C_Str();
  return name.find(target) != std::string::npos;
}

void CollectMeshIndicesForNodeName(const aiNode* node,
                                   const std::string& target,
                                   std::vector<bool>& mask,
                                   bool inTarget);

std::unordered_set<std::string> FindMeshNamesForNode(const std::string& cadPath,
                                                     const std::string& target) {
  std::unordered_set<std::string> names;
  std::ifstream infile(cadPath);
  if (!infile) {
    return names;
  }

  struct DaeNode {
    std::string id;
    std::string name;
    std::vector<std::string> instanceGeometry;
    std::vector<std::string> instanceNode;
    std::vector<int> children;
  };

  auto parse_attr = [](const std::string& line, const std::string& key) -> std::string {
    const std::string needle = key + "=\"";
    const std::size_t start = line.find(needle);
    if (start == std::string::npos) {
      return {};
    }
    const std::size_t value_start = start + needle.size();
    const std::size_t end = line.find('"', value_start);
    if (end == std::string::npos) {
      return {};
    }
    return line.substr(value_start, end - value_start);
  };

  std::vector<DaeNode> nodes;
  std::unordered_map<std::string, int> id_to_index;
  std::vector<int> stack;

  std::string line;
  while (std::getline(infile, line)) {
    if (line.find("<node") != std::string::npos && line.find("</node") == std::string::npos) {
      DaeNode node;
      node.id = parse_attr(line, "id");
      node.name = parse_attr(line, "name");
      const int idx = static_cast<int>(nodes.size());
      nodes.push_back(std::move(node));
      if (!nodes.back().id.empty()) {
        id_to_index[nodes.back().id] = idx;
      }
      if (!stack.empty()) {
        nodes[stack.back()].children.push_back(idx);
      }
      stack.push_back(idx);
    }

    if (line.find("<instance_geometry") != std::string::npos && !stack.empty()) {
      std::string url = parse_attr(line, "url");
      if (!url.empty() && url.front() == '#') {
        url.erase(0, 1);
      }
      if (!url.empty()) {
        nodes[stack.back()].instanceGeometry.push_back(url);
      }
    }

    if (line.find("<instance_node") != std::string::npos && !stack.empty()) {
      std::string url = parse_attr(line, "url");
      if (!url.empty() && url.front() == '#') {
        url.erase(0, 1);
      }
      if (!url.empty()) {
        nodes[stack.back()].instanceNode.push_back(url);
      }
    }

    if (line.find("</node>") != std::string::npos && !stack.empty()) {
      stack.pop_back();
    }
  }

  std::unordered_set<int> visited;
  std::vector<int> targets;
  for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
    if (nodes[i].name.find(target) != std::string::npos) {
      targets.push_back(i);
    }
  }

  std::vector<int> stack_nodes;
  for (int idx : targets) {
    stack_nodes.push_back(idx);
  }

  while (!stack_nodes.empty()) {
    const int idx = stack_nodes.back();
    stack_nodes.pop_back();
    if (!visited.insert(idx).second) {
      continue;
    }
    for (const auto& geom : nodes[idx].instanceGeometry) {
      names.insert(geom);
    }
    for (const auto& ref : nodes[idx].instanceNode) {
      auto it = id_to_index.find(ref);
      if (it != id_to_index.end()) {
        stack_nodes.push_back(it->second);
      }
    }
    for (int child : nodes[idx].children) {
      stack_nodes.push_back(child);
    }
  }
  return names;
}

void CollectMeshIndicesForNodeName(const aiNode* node,
                                   const std::string& target,
                                   std::vector<bool>& mask,
                                   bool inTarget = false) {
  if (!node) {
    return;
  }

  const bool matched = inTarget || NodeNameMatches(node, target);
  if (matched) {
    for (unsigned int i = 0; i < node->mNumMeshes; ++i) {
      const unsigned int meshIndex = node->mMeshes[i];
      if (meshIndex < mask.size()) {
        mask[meshIndex] = true;
      }
    }
  }

  for (unsigned int i = 0; i < node->mNumChildren; ++i) {
    CollectMeshIndicesForNodeName(node->mChildren[i], target, mask, matched);
  }
}

G4VisAttributes* MakeVisAttributes(const G4Colour& color, bool solid = true);

G4LogicalVolume* BuildCadBoundingUnionForMeshes(const aiScene& scene,
                                                G4Material* material,
                                                G4double meshScale,
                                                const std::string& name,
                                                const std::vector<bool>* meshMask,
                                                bool includeMask,
                                                const G4Colour& color) {
  auto* multi = new G4MultiUnion(name);
  std::size_t added = 0;

  for (unsigned int m = 0; m < scene.mNumMeshes; ++m) {
    if (meshMask) {
      const bool selected = (m < meshMask->size()) ? (*meshMask)[m] : false;
      if (includeMask != selected) {
        continue;
      }
    }

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
  auto* lv = new G4LogicalVolume(multi, material, name + "_lv");
  lv->SetVisAttributes(MakeVisAttributes(color));
  return lv;
}

bool ComputeMeshBounds(const aiScene& scene,
                       const std::vector<bool>* meshMask,
                       bool includeMask,
                       G4double meshScale,
                       G4ThreeVector& outMin,
                       G4ThreeVector& outMax) {
  double minx = std::numeric_limits<double>::max();
  double miny = std::numeric_limits<double>::max();
  double minz = std::numeric_limits<double>::max();
  double maxx = -std::numeric_limits<double>::max();
  double maxy = -std::numeric_limits<double>::max();
  double maxz = -std::numeric_limits<double>::max();
  bool hasVertex = false;

  for (unsigned int m = 0; m < scene.mNumMeshes; ++m) {
    if (meshMask) {
      const bool selected = (m < meshMask->size()) ? (*meshMask)[m] : false;
      if (includeMask != selected) {
        continue;
      }
    }

    const aiMesh* mesh = scene.mMeshes[m];
    if (!mesh || mesh->mNumVertices == 0) {
      continue;
    }
    for (unsigned int v = 0; v < mesh->mNumVertices; ++v) {
      const aiVector3D& vv = mesh->mVertices[v];
      minx = std::min(minx, static_cast<double>(vv.x));
      miny = std::min(miny, static_cast<double>(vv.y));
      minz = std::min(minz, static_cast<double>(vv.z));
      maxx = std::max(maxx, static_cast<double>(vv.x));
      maxy = std::max(maxy, static_cast<double>(vv.y));
      maxz = std::max(maxz, static_cast<double>(vv.z));
      hasVertex = true;
    }
  }

  if (!hasVertex) {
    return false;
  }

  outMin = G4ThreeVector(minx * meshScale, miny * meshScale, minz * meshScale);
  outMax = G4ThreeVector(maxx * meshScale, maxy * meshScale, maxz * meshScale);
  return true;
}

G4LogicalVolume* BuildCadMergedTessellatedForMeshes(const aiScene& scene,
                                                    G4Material* material,
                                                    G4double meshScale,
                                                    const std::string& name,
                                                    const std::vector<bool>* meshMask,
                                                    bool includeMask,
                                                    const G4Colour& color) {
  auto* solid = new G4TessellatedSolid(name);

  for (unsigned int m = 0; m < scene.mNumMeshes; ++m) {
    if (meshMask) {
      const bool selected = (m < meshMask->size()) ? (*meshMask)[m] : false;
      if (includeMask != selected) {
        continue;
      }
    }

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

  auto* lv = new G4LogicalVolume(solid, material, name + "_lv");
  lv->SetVisAttributes(MakeVisAttributes(color));
  return lv;
}

bool SolidContainsBox(const G4VSolid* solid,
                      const G4ThreeVector& center,
                      const G4ThreeVector& half) {
  if (!solid) {
    return false;
  }
  const G4ThreeVector corners[8] = {
      center + G4ThreeVector( half.x(),  half.y(),  half.z()),
      center + G4ThreeVector( half.x(),  half.y(), -half.z()),
      center + G4ThreeVector( half.x(), -half.y(),  half.z()),
      center + G4ThreeVector( half.x(), -half.y(), -half.z()),
      center + G4ThreeVector(-half.x(),  half.y(),  half.z()),
      center + G4ThreeVector(-half.x(),  half.y(), -half.z()),
      center + G4ThreeVector(-half.x(), -half.y(),  half.z()),
      center + G4ThreeVector(-half.x(), -half.y(), -half.z()),
  };
  for (const auto& corner : corners) {
    const auto inside = solid->Inside(corner);
    if (inside == kOutside) {
      return false;
    }
  }
  return true;
}

G4VisAttributes* MakeVisAttributes(const G4Colour& color, bool solid) {
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
  fCCDMessenger->DeclareProperty("printInfo", fCCDPrintInfo,
                                 "Print CCD thickness and bounding box at initialization.");
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
  if (fCCDPrintInfo) {
    PrintCCDInfo();
  }
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

  constexpr G4double kCCDActiveSide = 5.461 * cm;
  constexpr G4double halfX = 0.5 * kCCDActiveSide;
  constexpr G4double halfY = 0.5 * kCCDActiveSide;
  constexpr G4double halfZ = 0.3625 * mm;

  auto* solid = new G4Box("PrimitiveSensor", halfX, halfY, halfZ);
  auto* logical = new G4LogicalVolume(solid, si, "PrimitiveSensorLogical");
  logical->SetVisAttributes(MakeVisAttributes(G4Colour(0.2, 0.6, 0.8)));

  new G4PVPlacement(nullptr, {}, logical, "PrimitiveSensorPV", fWorldLogical, false, 0,
                    fOptions.checkOverlaps);

  RegisterSensitiveVolume(logical);
  return logical;
}

G4LogicalVolume* B02DetectorConstruction::BuildCadGeometry() {
  fCadSteelLogical = nullptr;
  fCadCopperLogical = nullptr;
  if (fOptions.cadMode == "none") {
    G4cout << "[geometry] CAD mode: none (no CAD volumes)." << G4endl;
    fCCDOverlayCenter = G4ThreeVector();
    return nullptr;
  }

  const std::string cadPath = ResolveCadPath();
  Assimp::Importer importer;
  importer.SetPropertyInteger(AI_CONFIG_PP_SBP_REMOVE,
                              aiPrimitiveType_POINT | aiPrimitiveType_LINE);
  importer.SetPropertyInteger(AI_CONFIG_PP_PTV_KEEP_HIERARCHY, 1);

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
  auto* defaultMat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  if (!defaultMat) {
    defaultMat = nist->FindOrBuildMaterial("G4_Fe");
  }
  auto* copperMat = nist->FindOrBuildMaterial("G4_Cu");
  const G4double meshScale = 1000.0 * mm;  // meters -> mm

  const std::string kCopperNodeName = "Assembly 1 New";
  const auto copperMeshNames = FindMeshNamesForNode(cadPath, kCopperNodeName);

  std::vector<bool> copperMeshes(scene->mNumMeshes, false);
  if (scene->mRootNode) {
    CollectMeshIndicesForNodeName(scene->mRootNode, kCopperNodeName, copperMeshes);
  }
  if (!copperMeshNames.empty()) {
    for (unsigned int i = 0; i < scene->mNumMeshes; ++i) {
      if (scene->mMeshes[i] &&
          copperMeshNames.count(scene->mMeshes[i]->mName.C_Str()) > 0) {
        copperMeshes[i] = true;
      }
    }
  }
  const auto copperCount = static_cast<std::size_t>(
      std::count(copperMeshes.begin(), copperMeshes.end(), true));
  if (copperCount > 0) {
    G4cout << "[geometry] CAD material override: " << copperCount
           << " mesh(es) set to copper (node match: 'Assembly 1 New')." << G4endl;
  } else {
    G4cout << "[geometry] CAD material override: no meshes matched 'Assembly 1 New'; "
              "using steel for all CAD parts." << G4endl;
    if (std::getenv("SIMCCD_CAD_DUMP_NAMES")) {
      for (unsigned int i = 0; i < scene->mNumMeshes; ++i) {
        const aiMesh* mesh = scene->mMeshes[i];
        if (!mesh) {
          continue;
        }
        const std::string meshName = mesh->mName.C_Str();
        G4cout << "[geometry] mesh[" << i << "] name='" << meshName
               << "' vertices=" << mesh->mNumVertices << G4endl;
      }
    }
  }

  G4ThreeVector boundsMin;
  G4ThreeVector boundsMax;
  if (copperCount > 0 &&
      ComputeMeshBounds(*scene, &copperMeshes, true, meshScale, boundsMin, boundsMax)) {
    fCCDOverlayCenter = 0.5 * (boundsMin + boundsMax);
    G4cout << "[geometry] CCD overlay center from copper bounds: ("
           << fCCDOverlayCenter.x() / mm << ", " << fCCDOverlayCenter.y() / mm << ", "
           << fCCDOverlayCenter.z() / mm << ") mm" << G4endl;
    const G4ThreeVector span = boundsMax - boundsMin;
    const G4ThreeVector halfSpan = 0.5 * span;
    const G4ThreeVector overlayHalf(0.5 * 5.461 * cm, 0.5 * 5.461 * cm, 0.3625 * mm);
    const bool overlayFits = (overlayHalf.x() <= halfSpan.x()) &&
                             (overlayHalf.y() <= halfSpan.y()) &&
                             (overlayHalf.z() <= halfSpan.z());
    G4cout << "[geometry] Copper bounds span: (" << span.x() / mm << ", " << span.y() / mm
           << ", " << span.z() / mm << ") mm; overlay half size: (" << overlayHalf.x() / mm
           << ", " << overlayHalf.y() / mm << ", " << overlayHalf.z() / mm
           << ") mm; fits=" << (overlayFits ? "yes" : "no") << G4endl;
  } else if (ComputeMeshBounds(*scene, nullptr, true, meshScale, boundsMin, boundsMax)) {
    fCCDOverlayCenter = 0.5 * (boundsMin + boundsMax);
    G4cout << "[geometry] CCD overlay center from assembly bounds: ("
           << fCCDOverlayCenter.x() / mm << ", " << fCCDOverlayCenter.y() / mm << ", "
           << fCCDOverlayCenter.z() / mm << ") mm" << G4endl;
  }

  std::vector<G4LogicalVolume*> cadVolumes;
  if (fOptions.cadMode == "parts") {
    cadVolumes = BuildCadParts(*scene, defaultMat, copperMat, meshScale, copperMeshes);
  } else if (fOptions.cadMode == "tessellated") {
    if (copperCount > 0 && copperMat) {
      auto* steel = BuildCadMergedTessellatedForMeshes(
          *scene, defaultMat, meshScale, "CADAssemblySteel", &copperMeshes, false,
          G4Colour(0.6, 0.6, 0.7));
      auto* copper = BuildCadMergedTessellatedForMeshes(
          *scene, copperMat, meshScale, "CADAssemblyCopper", &copperMeshes, true,
          G4Colour(0.8, 0.4, 0.2));
      if (steel) {
        fCadSteelLogical = steel;
        cadVolumes.push_back(steel);
      }
      if (copper) {
        fCadCopperLogical = copper;
        cadVolumes.push_back(copper);
      }
    } else {
      auto* merged = BuildCadMergedTessellated(*scene, defaultMat, meshScale);
      if (merged) {
        fCadSteelLogical = merged;
        cadVolumes.push_back(merged);
      }
    }
  } else {
    if (copperCount > 0 && copperMat) {
      auto* steel = BuildCadBoundingUnionForMeshes(
          *scene, defaultMat, meshScale, "CADAssemblySteel", nullptr, true,
          G4Colour(0.6, 0.6, 0.7));
      auto* copper = BuildCadBoundingUnionForMeshes(
          *scene, copperMat, meshScale, "CADAssemblyCopper", &copperMeshes, true,
          G4Colour(0.8, 0.4, 0.2));
      if (steel) {
        fCadSteelLogical = steel;
        cadVolumes.push_back(steel);
      }
      if (copper) {
        fCadCopperLogical = copper;
        cadVolumes.push_back(copper);
      }
    } else {
      auto* merged = BuildCadBoundingUnion(*scene, defaultMat, meshScale);
      if (merged) {
        fCadSteelLogical = merged;
        cadVolumes.push_back(merged);
      }
    }
  }

  if (cadVolumes.empty()) {
    G4Exception("B02DetectorConstruction::BuildCadGeometry", "CadEmpty", FatalException,
                "No CAD volumes were created.");
  }
  if (!fCadSteelLogical && !cadVolumes.empty()) {
    fCadSteelLogical = cadVolumes.front();
  }

  G4LogicalVolume* steelMother = nullptr;
  if ((fOptions.cadMode == "tessellated" || fOptions.cadMode == "merged") &&
      copperCount > 0 && copperMat) {
    for (auto* lv : cadVolumes) {
      if (lv && lv->GetMaterial() != copperMat) {
        steelMother = lv;
        break;
      }
    }
  }

  G4int copyNo = 0;
  for (auto* lv : cadVolumes) {
    if (!lv) {
      continue;
    }
    auto name = lv->GetName();
    G4LogicalVolume* mother = fWorldLogical;
    if (steelMother && lv->GetMaterial() == copperMat) {
      mother = steelMother;
    }
    auto* pv = new G4PVPlacement(nullptr, {}, lv, name + "_pv", mother, false, copyNo,
                                 fOptions.checkOverlaps);
    if (fOptions.checkOverlaps) {
      pv->CheckOverlaps(fOptions.overlapSamples, fOptions.overlapTolerance, true, 1);
    }
    ++copyNo;
  }

  return cadVolumes.front();
}

G4LogicalVolume* B02DetectorConstruction::BuildCadMergedTessellated(const aiScene& scene,
                                                                    G4Material* material,
                                                                    G4double meshScale) {
  return BuildCadMergedTessellatedForMeshes(scene, material, meshScale, "CADAssemblyMerged",
                                            nullptr, true, G4Colour(0.6, 0.6, 0.8));
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
  // Dimensions match the active CCD size (5.461 cm x 5.461 cm) and a thin thickness.
  if (fCCDOverlayLogical || fOptions.geometryMode == "primitive") {
    return;
  }

  auto* nist = G4NistManager::Instance();
  auto* si = nist->FindOrBuildMaterial("G4_Si");

  const G4double halfX = 0.5 * 5.461 * cm;
  const G4double halfY = 0.5 * 5.461 * cm;
  const G4double halfZ = 0.3625 * mm;  // 0.725 mm total thickness

  auto* solid = new G4Box("CCDOverlay", halfX, halfY, halfZ);
  fCCDOverlayLogical = new G4LogicalVolume(solid, si, "CCDOverlayLogical");
  fCCDOverlayLogical->SetVisAttributes(MakeVisAttributes(G4Colour(0.8, 0.2, 0.2), false));

  const G4ThreeVector overlayHalf(halfX, halfY, halfZ);
  G4LogicalVolume* mother = fWorldLogical;
  if (fCadCopperLogical &&
      SolidContainsBox(fCadCopperLogical->GetSolid(), fCCDOverlayCenter, overlayHalf)) {
    mother = fCadCopperLogical;
  } else if (fCadSteelLogical &&
             SolidContainsBox(fCadSteelLogical->GetSolid(), fCCDOverlayCenter, overlayHalf)) {
    mother = fCadSteelLogical;
  } else if (fPrimaryScoring &&
             SolidContainsBox(fPrimaryScoring->GetSolid(), fCCDOverlayCenter, overlayHalf)) {
    mother = fPrimaryScoring;
  } else {
    G4cout << "[geometry][WARN] CCD overlay not contained in CAD solids; placing in world."
           << G4endl;
  }
  if (mother) {
    G4cout << "[geometry] CCD overlay mother volume: " << mother->GetName() << G4endl;
  }

  new G4PVPlacement(nullptr, fCCDOverlayCenter, fCCDOverlayLogical, "CCDOverlayPV", mother, false,
                    0,
                    false);
  RegisterSensitiveVolume(fCCDOverlayLogical);
}

G4LogicalVolume* B02DetectorConstruction::BuildCadBoundingUnion(const aiScene& scene,
                                                                G4Material* material,
                                                                G4double meshScale) {
  return BuildCadBoundingUnionForMeshes(scene, material, meshScale, "CADAssemblyMerged", nullptr,
                                        true, G4Colour(0.6, 0.6, 0.8));
}

std::vector<G4LogicalVolume*> B02DetectorConstruction::BuildCadParts(
    const aiScene& scene,
    G4Material* defaultMaterial,
    G4Material* copperMaterial,
    G4double meshScale,
    const std::vector<bool>& copperMeshes) {
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

    G4Material* material = defaultMaterial;
    G4Colour color(0.6, 0.6, 0.7);
    if (m < copperMeshes.size() && copperMeshes[m] && copperMaterial) {
      material = copperMaterial;
      color = G4Colour(0.8, 0.4, 0.2);
    }

    auto* lv = new G4LogicalVolume(tess, material, meshName + "_lv");
    lv->SetVisAttributes(MakeVisAttributes(color));
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
  auto thickness_from = [](const G4LogicalVolume* lv) -> G4double {
    if (!lv) {
      return 0.0;
    }
    if (const auto* box = dynamic_cast<const G4Box*>(lv->GetSolid())) {
      return 2.0 * box->GetZHalfLength();
    }
    return 0.0;
  };
  if (fCCDOverlayLogical) {
    return thickness_from(fCCDOverlayLogical);
  }
  return thickness_from(fPrimaryScoring);
}

void B02DetectorConstruction::PrintCCDInfo() const {
  auto print_box = [](const G4LogicalVolume* lv, const char* label) {
    if (!lv) {
      return;
    }
    const auto* solid = lv->GetSolid();
    if (!solid) {
      return;
    }
    const auto* box = dynamic_cast<const G4Box*>(solid);
    if (box) {
      const auto full_x = 2.0 * box->GetXHalfLength();
      const auto full_y = 2.0 * box->GetYHalfLength();
      const auto full_z = 2.0 * box->GetZHalfLength();
      G4cout << "[CCD] " << label << " solid=" << solid->GetName()
             << " size=(" << full_x / mm << ", " << full_y / mm << ", " << full_z / mm
             << ") mm thickness=" << full_z / mm << " mm" << G4endl;
    } else {
      G4cout << "[CCD] " << label << " solid=" << solid->GetName()
             << " (non-box, thickness unavailable)" << G4endl;
    }
  };

  G4cout << "[CCD] geometryMode=" << fOptions.geometryMode << G4endl;
  if (fPrimaryScoring) {
    print_box(fPrimaryScoring, "primary");
  }
  if (fCCDOverlayLogical) {
    print_box(fCCDOverlayLogical, "overlay");
  }
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
    candidates.emplace_back(fs::path(assetsDir) / "cad" / "new_assembly.dae");
    candidates.emplace_back(fs::path(assetsDir) / simccd::kCadRelativePath);
  }
  if (simccd::kBinaryAssetsDir && std::strlen(simccd::kBinaryAssetsDir) > 0) {
    candidates.emplace_back(fs::path(simccd::kBinaryAssetsDir) / "cad" /
                            "new_assembly.dae");
    candidates.emplace_back(
        fs::path(simccd::kBinaryAssetsDir) / simccd::kCadRelativePath);
  }
  if (simccd::kSourceAssetsDir && std::strlen(simccd::kSourceAssetsDir) > 0) {
    candidates.emplace_back(fs::path(simccd::kSourceAssetsDir) / "cad" /
                            "new_assembly.dae");
    candidates.emplace_back(
        fs::path(simccd::kSourceAssetsDir) / simccd::kCadRelativePath);
  }
  candidates.emplace_back(fs::current_path() / "build" / "assets" / "cad" /
                          "new_assembly.dae");
  candidates.emplace_back(fs::current_path() / "assets" / "cad" / "new_assembly.dae");
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
