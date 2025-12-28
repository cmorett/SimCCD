// 
// --------------------------------------------------------------
//      GEANT4 - B02 with CAD/primitive geometry switches
// --------------------------------------------------------------

#include "B02DetectorConstruction.hh"
#include "B02PhysicsList.hh"
#include "B02PrimaryGeneratorAction.hh"
#include "B02RunAction.hh"
#include "B02EventAction.hh"
#include "B02SteppingAction.hh"
#include "B02SteppingVerbose.hh"
#include "B02ActionInitialization.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4SteppingVerbose.hh"
#include "G4Version.hh"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

std::string ToLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

struct AppOptions {
  GeometryOptions geom;
  std::string macroFile;
  bool noVis = false;
  bool showHelp = false;
};

void PrintUsage() {
  std::cout << "Usage: b02_executable [options] [macro.mac]\n"
            << "Options:\n"
            << "  --geometry <cad|primitive>   Select geometry mode (default: cad)\n"
            << "  --cad-mode <merged|parts|tessellated|none>    CAD import mode (default: merged)\n"
            << "  --cad-file <path>            Path to assembly.dae\n"
            << "  --assets-dir <dir>           Path to assets directory containing cad/\n"
            << "  --no-vis                     Disable visualization manager\n"
            << "  --help                       Show this message\n";
}

void ApplyEnvOverrides(GeometryOptions& geom) {
  if (const char* envGeom = std::getenv("SIMCCD_GEOMETRY_MODE")) {
    geom.geometryMode = envGeom;
  }
  if (const char* envCadMode = std::getenv("SIMCCD_CAD_MODE")) {
    geom.cadMode = envCadMode;
  }
  if (const char* envCadFile = std::getenv("SIMCCD_CAD_FILE")) {
    if (geom.cadFile.empty()) {
      geom.cadFile = envCadFile;
    }
  }
  if (const char* envAssets = std::getenv("SIMCCD_ASSETS_DIR")) {
    if (geom.assetsDir.empty()) {
      geom.assetsDir = envAssets;
    }
  }
}

AppOptions ParseOptions(int argc, char** argv) {
  AppOptions opts;

  auto nextArg = [&](int& i) -> std::string {
    if (i + 1 < argc) {
      return std::string(argv[++i]);
    }
    return {};
  };

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--help" || arg == "-h") {
      opts.showHelp = true;
      break;
    } else if (arg == "--geometry") {
      opts.geom.geometryMode = nextArg(i);
    } else if (arg.rfind("--geometry=", 0) == 0) {
      opts.geom.geometryMode = arg.substr(std::string("--geometry=").size());
    } else if (arg == "--cad-mode") {
      opts.geom.cadMode = nextArg(i);
    } else if (arg.rfind("--cad-mode=", 0) == 0) {
      opts.geom.cadMode = arg.substr(std::string("--cad-mode=").size());
    } else if (arg == "--cad-file") {
      opts.geom.cadFile = nextArg(i);
    } else if (arg.rfind("--cad-file=", 0) == 0) {
      opts.geom.cadFile = arg.substr(std::string("--cad-file=").size());
    } else if (arg == "--assets-dir") {
      opts.geom.assetsDir = nextArg(i);
    } else if (arg.rfind("--assets-dir=", 0) == 0) {
      opts.geom.assetsDir = arg.substr(std::string("--assets-dir=").size());
    } else if (arg == "--no-vis") {
      opts.noVis = true;
    } else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << "\n";
    } else if (opts.macroFile.empty()) {
      opts.macroFile = arg;
    }
  }

  ApplyEnvOverrides(opts.geom);
  opts.geom.geometryMode = ToLower(opts.geom.geometryMode);
  opts.geom.cadMode = ToLower(opts.geom.cadMode);

  if (opts.geom.geometryMode != "cad" && opts.geom.geometryMode != "primitive") {
    std::cerr << "Unrecognized geometry mode '" << opts.geom.geometryMode
              << "', falling back to 'cad'.\n";
    opts.geom.geometryMode = "cad";
  }
  if (opts.geom.cadMode != "parts" && opts.geom.cadMode != "merged" &&
      opts.geom.cadMode != "tessellated" && opts.geom.cadMode != "none") {
    std::cerr << "Unrecognized CAD mode '" << opts.geom.cadMode
              << "', falling back to 'merged'.\n";
    opts.geom.cadMode = "merged";
  }

  return opts;
}

}  // namespace

int main(int argc, char** argv)
{
  auto options = ParseOptions(argc, argv);
  if (options.showHelp) {
    PrintUsage();
    return 0;
  }

  G4UIExecutive* ui = nullptr;
  if (options.macroFile.empty()) {
    ui = new G4UIExecutive(static_cast<G4int>(argc), argv);
  }
  
  // Choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // User Verbose output class
  G4VSteppingVerbose* verbosity = new B02SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // The Run manager
  G4RunManager* runManager = new G4RunManager;

  auto* detector = new B02DetectorConstruction(options.geom);
  runManager->SetUserInitialization(detector); 

  G4VUserPhysicsList* physics = new B02PhysicsList;
  runManager->SetUserInitialization(physics);

  runManager->SetUserInitialization(new B02ActionInitialization(detector));
  
  runManager->Initialize();
    
  // Visualization
  G4VisManager* visManager = nullptr;
  if (!options.noVis) {
#if defined(G4VERSION_NUMBER) && (G4VERSION_NUMBER >= 1110)
    visManager = new G4VisExecutive();
#else
    visManager = new G4VisExecutive(argc, argv);
#endif
    visManager->Initialize();
  }

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if (!ui) {
    if (!options.macroFile.empty()) {
      G4String command = "/control/execute ";
      G4String fileName = options.macroFile.c_str();
      UImanager->ApplyCommand(command + fileName);
    }
  }
  else {
    if (!options.noVis) {
      UImanager->ApplyCommand("/control/execute init_vis.mac");
    }
    ui->SessionStart();
    delete ui;
  }
  
  delete visManager;
  delete runManager;
 
}
