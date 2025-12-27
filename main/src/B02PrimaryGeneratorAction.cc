//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "B02PrimaryGeneratorAction.hh"

#include "B02DetectorConstruction.hh"
#include "B02PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4GenericMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <algorithm>
#include <cmath>
#include <utility>

namespace {

enum class MuonMode {
  ForcedFootprint,
  TierBPlaneFlux,
};

enum class MuonImpactMode {
  Unbiased,
  Targeted,
};

enum class MuonChargeMode {
  Equal,
  FixedRatio,
  MinusOnly,
  PlusOnly,
};

MuonMode ParseMuonMode(const G4String& value) {
  G4String lower = value;
  lower.toLower();
  if (lower == "tierb_plane_flux" || lower == "tierb" || lower == "plane_flux") {
    return MuonMode::TierBPlaneFlux;
  }
  return MuonMode::ForcedFootprint;
}

MuonImpactMode ParseMuonImpactMode(const G4String& value) {
  G4String lower = value;
  lower.toLower();
  if (lower == "targeted" || lower == "target") {
    return MuonImpactMode::Targeted;
  }
  return MuonImpactMode::Unbiased;
}

MuonChargeMode ParseChargeMode(const G4String& value) {
  G4String lower = value;
  lower.toLower();
  if (lower == "fixedratio" || lower == "fixed_ratio" || lower == "ratio") {
    return MuonChargeMode::FixedRatio;
  }
  if (lower == "minusonly" || lower == "mu-" || lower == "mu- only") {
    return MuonChargeMode::MinusOnly;
  }
  if (lower == "plusonly" || lower == "mu+" || lower == "mu+ only" || lower == "positive") {
    return MuonChargeMode::PlusOnly;
  }
  return MuonChargeMode::Equal;
}

G4String ToString(MuonChargeMode mode) {
  switch (mode) {
    case MuonChargeMode::FixedRatio:
      return "fixedRatio";
    case MuonChargeMode::MinusOnly:
      return "minusOnly";
    case MuonChargeMode::PlusOnly:
      return "plusOnly";
    case MuonChargeMode::Equal:
    default:
      return "equal";
  }
}

G4String ToString(MuonImpactMode mode) {
  switch (mode) {
    case MuonImpactMode::Targeted:
      return "targeted";
    case MuonImpactMode::Unbiased:
    default:
      return "unbiased";
  }
}

G4String ToString(MuonMode mode) {
  switch (mode) {
    case MuonMode::TierBPlaneFlux:
      return "tierB_plane_flux";
    case MuonMode::ForcedFootprint:
    default:
      return "forced_footprint";
  }
}

double Clamp(double x, double lo, double hi) {
  return std::max(lo, std::min(hi, x));
}

}  // namespace

G4int B02PrimaryGeneratorAction::GetMuonModeCode() const {
  return (ParseMuonMode(fMuonModeString) == MuonMode::TierBPlaneFlux) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02PrimaryGeneratorAction::B02PrimaryGeneratorAction(
    B02DetectorConstruction* myDC)
    : myDetector(myDC), rndmFlag("off"), fUseCosmicMuons(true) {
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  // create a messenger for this class
  gunMessenger = new B02PrimaryGeneratorMessenger(this);

  // default particle
  G4ParticleDefinition* particle =
      G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  particleGun->SetParticleDefinition(particle);

  DefineCommands();
  fMuonModeString = ToString(MuonMode::ForcedFootprint);
  fMuonImpactModeString = ToString(MuonImpactMode::Unbiased);
  fRequestedEminGeV = fEminGeV;
  fRequestedEmaxGeV = fEmaxGeV;
  fEffectiveEminGeV = fEminGeV;
  fEffectiveEmaxGeV = fEmaxGeV;
  MarkSamplerDirty();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02PrimaryGeneratorAction::~B02PrimaryGeneratorAction() {
  delete particleGun;
  delete gunMessenger;
  delete fMessenger;
  delete fMuonMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  SelectMuonDefinition();

  fEventWeight = 0.0;
  fEventLivetime = 0.0;
  fSampledEnergyGeV = fMuonEnergyGeV;
  fMuonCosTheta = 1.0;
  fGeomIntersectsCCD = false;
  fMuonIsTargeted = false;

  if (fUseCosmicMuons) {
    const MuonMode mode = ParseMuonMode(fMuonModeString);
    if (mode == MuonMode::TierBPlaneFlux) {
      GenerateTierBFlux(anEvent);
    } else {
      GenerateForcedFootprint(anEvent);
    }
    return;
  }

  // Non-cosmic mode: respect fixed energy if requested and let the configured
  // particleGun settings drive direction/position.
  if (fUseFixedEnergy) {
    particleGun->SetParticleEnergy(fMuonEnergyGeV * GeV);
  }
  particleGun->GeneratePrimaryVertex(anEvent);
}

void B02PrimaryGeneratorAction::GenerateForcedFootprint(G4Event* anEvent) {
  // Legacy biased mode: back-project through the CCD footprint.
  const G4double thetaUpper = std::min(fThetaMaxRad, 0.5 * pi);
  G4double theta = 0.0;
  while (true) {
    const G4double p = G4RandFlat::shoot(0.0, 1.0);
    theta = G4RandFlat::shoot(0.0, thetaUpper);
    const G4double q = std::pow(std::cos(theta), 2) * std::sin(theta);
    if (p <= q) {
      break;
    }
  }
  const G4double phi = G4RandFlat::shoot(0.0, twopi);

  const G4double cosTheta = std::cos(theta);
  const G4double sinTheta = std::sin(theta);
  const G4double X = R * sinTheta * std::cos(phi);
  const G4double Y = R * sinTheta * std::sin(phi);
  const G4double Z = R * cosTheta;

  const G4double u = G4RandFlat::shoot(-px / 2.0, px / 2.0);
  const G4double v = G4RandFlat::shoot(-py / 2.0, py / 2.0);
  const G4double x0_cm = X + u * cosTheta * std::cos(phi) - v * std::sin(phi);
  const G4double y0_cm = Y + u * cosTheta * std::sin(phi) + v * std::cos(phi);
  const G4double z0_cm = Z - u * std::sin(theta);

  const G4double targetZ_cm = fZImpactPlane / cm;
  const G4double deltaZ_cm = z0_cm - targetZ_cm;
  const G4double tanTheta = (cosTheta != 0.0) ? sinTheta / cosTheta : 0.0;
  const G4double xImp_cm = x0_cm - deltaZ_cm * tanTheta * std::cos(phi);
  const G4double yImp_cm = y0_cm - deltaZ_cm * tanTheta * std::sin(phi);

  fMuonX0 = x0_cm * cm;
  fMuonY0 = y0_cm * cm;
  fMuonZ0 = z0_cm * cm;
  fMuonXImp = xImp_cm * cm;
  fMuonYImp = yImp_cm * cm;
  fMuonZImp = fZImpactPlane;
  fMuonTheta = theta;
  fMuonPhi = phi;
  fMuonCosTheta = cosTheta;

  particleGun->SetParticlePosition(
      G4ThreeVector(fMuonX0, fMuonY0, fMuonZ0));
  particleGun->SetParticleMomentumDirection(
      G4ThreeVector(-sinTheta * std::cos(phi), -sinTheta * std::sin(phi),
                    -cosTheta));

  G4double kineticEnergy = fMuonEnergyGeV * GeV;
  if (!fUseFixedEnergy) {
    // Smith & Duller energy spectrum
    double Emin = -1;
    double Emax = 5;

    const int ee = 10000;
    double ES[ee];
    double dE_log = (Emax - Emin) / ee;

    double Eu;            // Variable de energia cinetica
    double Au = 2e9;                      // Parametros de la funcion de Smith
    double gu = 2.645;                    // ...
    double ru = 0.76;                     // ...
    double au = 2.5;
    double y0u = 1000.0;
    double bmu = 0.80;
    double cu = 299792458.0e2;
    double mmu = 105.7 / pow(cu, 2);
    double t0mu = 2.2e-6;
    double r0u = 0.00129;
    double Epu;
    double Bmu = bmu * mmu * y0u * cu / (t0mu * r0u);
    double Pmu;
    double lpu = 120.0;
    double bu = 0.771;
    double mpu = 139.6 / pow(cu, 2);
    double t0pu = 2.6e-8;
    double jpu = mpu * y0u * cu / (t0pu * r0u);

    for (int j = 0; j < ee; j++) {    // Construye la funcion de Smith en un arreglo
      Eu = pow(10, Emin + j * dE_log);
      Epu = (Eu + au * y0u * (1.0 / cosTheta - 0.100)) / ru;
      Pmu = pow(0.100 * cosTheta *
                    (1 - (au * (y0u / cosTheta - 100) / (ru * Epu))),
                (Bmu / ((ru * Epu + 100 * au) * cosTheta)));
      ES[j] = Au * (pow(Epu, -gu)) * Pmu * lpu * bu * jpu /
              (Epu * cosTheta + bu * jpu);
    }

    int nbins = ee;
    G4RandGeneral GenDist(ES, nbins);          // Distribucion de energias
    double E = pow(10, Emin + (GenDist.shoot()) * (Emax - Emin));   // Sampleo de la energia
    kineticEnergy = E * MeV;
  }

  particleGun->SetParticleEnergy(kineticEnergy);
  particleGun->GeneratePrimaryVertex(anEvent);
  fSampledEnergyGeV = kineticEnergy / GeV;
}

void B02PrimaryGeneratorAction::GenerateTierBFlux(G4Event* anEvent) {
  EnsureSamplerReady();
  const auto plane = ComputeSourcePlaneSize();
  const G4double planeLx = plane.first;
  const G4double planeLy = plane.second;
  if (!fPrintedFluxInfo) {
    const double planeArea_cm2 = (planeLx / cm) * (planeLy / cm);
    const double flux_quad = fFluxSampler.FluxIntegralQuadrature();  // cm^-2 s^-1 over phi
    const double flux_mc = fFluxSampler.FluxIntegralMonteCarlo(200000);
    const double relDiff = (flux_quad > 0.0) ? std::abs(flux_quad - flux_mc) / flux_quad : 0.0;
    const double rate_hz = flux_quad * planeArea_cm2;
    G4cout << "[TierB] Flux integral (phi integrated, quadrature) = " << flux_quad
           << " cm^-2 s^-1, plane area = " << planeArea_cm2 << " cm^2, expected rate = "
           << rate_hz << " Hz" << G4endl;
    G4cout << "[TierB] Flux integral (phi integrated, MonteCarlo) = " << flux_mc
           << " cm^-2 s^-1, rel diff vs quadrature = " << relDiff * 100.0 << " %" << G4endl;
    if (relDiff > 0.02) {
      G4cout << "[TierB][WARN] Flux integral MC vs quadrature differ by more than 2%. Check 2pi factors."
             << G4endl;
    }
    fPrintedFluxInfo = true;
  }

  auto sampled = fFluxSampler.SampleCosThetaAndEnergy();
  const double cosTheta = Clamp(sampled.second, 0.0, 1.0);
  const double theta = std::acos(cosTheta);
  const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
  const double phi = G4RandFlat::shoot(0.0, twopi);

  const G4double z0 = fSourcePlaneZ;

  const G4ThreeVector dir(sinTheta * std::cos(phi),
                          sinTheta * std::sin(phi),
                          -cosTheta);
  const G4double zImp = fZImpactPlane;
  const G4double tToImpact =
      (std::abs(dir.z()) > 0.0) ? (zImp - z0) / dir.z() : 0.0;

  const MuonImpactMode impactMode = ParseMuonImpactMode(fMuonImpactModeString);
  fMuonIsTargeted = (impactMode == MuonImpactMode::Targeted);
  G4double x0 = 0.0;
  G4double y0 = 0.0;
  G4double xImp = 0.0;
  G4double yImp = 0.0;
  if (fMuonIsTargeted) {
    const G4double halfTargetX = 0.5 * px * cm + fTargetMargin;
    const G4double halfTargetY = 0.5 * py * cm + fTargetMargin;
    xImp = G4RandFlat::shoot(-halfTargetX, halfTargetX);
    yImp = G4RandFlat::shoot(-halfTargetY, halfTargetY);
    x0 = xImp - tToImpact * dir.x();
    y0 = yImp - tToImpact * dir.y();
  } else {
    x0 = G4RandFlat::shoot(-0.5 * planeLx, 0.5 * planeLx);
    y0 = G4RandFlat::shoot(-0.5 * planeLy, 0.5 * planeLy);
    xImp = x0 + tToImpact * dir.x();
    yImp = y0 + tToImpact * dir.y();
  }

  fMuonX0 = x0;
  fMuonY0 = y0;
  fMuonZ0 = z0;
  fMuonXImp = xImp;
  fMuonYImp = yImp;
  fMuonZImp = zImp;
  fMuonTheta = theta;
  fMuonPhi = phi;
  fMuonCosTheta = cosTheta;
  fSampledEnergyGeV = fUseFixedEnergy ? fMuonEnergyGeV : sampled.first;
  const G4double kineticEnergy = fSampledEnergyGeV * GeV;
  const G4double halfPx = 0.5 * px * cm;
  const G4double halfPy = 0.5 * py * cm;
  fGeomIntersectsCCD = (std::abs(fMuonXImp) <= halfPx) && (std::abs(fMuonYImp) <= halfPy);

  particleGun->SetParticlePosition(G4ThreeVector(fMuonX0, fMuonY0, fMuonZ0));
  particleGun->SetParticleMomentumDirection(dir);
  particleGun->SetParticleEnergy(kineticEnergy);
  particleGun->GeneratePrimaryVertex(anEvent);

  if (fUseFixedEnergy || fMuonIsTargeted) {
    fEventWeight = 0.0;
    fEventLivetime = 0.0;
  } else {
    const double planeArea_cm2 = (planeLx / cm) * (planeLy / cm);
    const double flux = fFluxSampler.FluxIntegral();  // cm^-2 s^-1 over phase space
    if (flux > 0.0 && planeArea_cm2 > 0.0) {
      fEventWeight = 1.0 / (flux * planeArea_cm2);
      fEventLivetime = fEventWeight;
    } else {
      fEventWeight = 0.0;
      fEventLivetime = 0.0;
    }
  }
}

void B02PrimaryGeneratorAction::SelectMuonDefinition() {
  const MuonChargeMode mode = ParseChargeMode(fMuonChargeMode);
  double probPlus = 0.5;
  switch (mode) {
    case MuonChargeMode::FixedRatio:
      if (fMuonPlusToMinusRatio <= 0.0) {
        probPlus = 0.0;
      } else {
        probPlus = fMuonPlusToMinusRatio / (1.0 + fMuonPlusToMinusRatio);
      }
      break;
    case MuonChargeMode::MinusOnly:
      probPlus = 0.0;
      break;
    case MuonChargeMode::PlusOnly:
      probPlus = 1.0;
      break;
    case MuonChargeMode::Equal:
    default:
      probPlus = 0.5;
      break;
  }

  const bool choosePlus = (G4UniformRand() < probPlus);
  const char* particleName = choosePlus ? "mu+" : "mu-";
  G4ParticleDefinition* particle =
      G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if (!particle) {
    particle = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  }
  particleGun->SetParticleDefinition(particle);
  fMuonPDGCode = choosePlus ? -13 : 13;
  fMuonChargeSign = choosePlus ? +1 : -1;
}

std::pair<G4double, G4double> B02PrimaryGeneratorAction::ComputeSourcePlaneSize() const {
  if (!fSourcePlaneAutoSize) {
    return {std::max(1e-6, fSourcePlaneLx), std::max(1e-6, fSourcePlaneLy)};
  }

  const double theta = std::min(fThetaMaxRad, 89.9 * deg);
  const double tanMax = std::tan(theta);
  const G4double margin = std::max(0.0, fSourcePlaneMargin);
  const G4double targetX = px * cm;
  const G4double targetY = py * cm;
  const G4double lx = targetX + 2.0 * fSourcePlaneZ * tanMax + margin;
  const G4double ly = targetY + 2.0 * fSourcePlaneZ * tanMax + margin;
  return {std::max(lx, targetX), std::max(ly, targetY)};
}

void B02PrimaryGeneratorAction::EnsureSamplerReady() {
  const auto plane = ComputeSourcePlaneSize();
  if (std::abs(fThetaMaxRad - fCachedThetaMax) > 1e-9 ||
      std::abs(fEminGeV - fCachedEminGeV) > 1e-12 ||
      std::abs(fEmaxGeV - fCachedEmaxGeV) > 1e-12 ||
      fFluxModel != fCachedFluxModel ||
      std::abs(plane.first - fCachedPlaneLx) > 1e-6 ||
      std::abs(plane.second - fCachedPlaneLy) > 1e-6 ||
      std::abs(fSourcePlaneZ - fCachedPlaneZ) > 1e-6) {
    fSamplerDirty = true;
  }

  if (!fSamplerDirty) {
    return;
  }

  fFluxSampler.SetFluxModel(fFluxModel);
  fFluxSampler.SetEnergyRange(fEminGeV, fEmaxGeV);
  fFluxSampler.SetThetaMax(std::min(fThetaMaxRad, 89.9 * deg));
  fFluxSampler.RebuildIfNeeded();

  fCachedThetaMax = fThetaMaxRad;
  fCachedEminGeV = fEminGeV;
  fCachedEmaxGeV = fEmaxGeV;
  fCachedFluxModel = fFluxModel;
  fCachedPlaneZ = fSourcePlaneZ;
  fCachedPlaneLx = plane.first;
  fCachedPlaneLy = plane.second;
  fSamplerDirty = false;
}

void B02PrimaryGeneratorAction::MarkSamplerDirty() {
  fSamplerDirty = true;
  fPrintedFluxInfo = false;
}

void B02PrimaryGeneratorAction::SetMuonMode(const G4String& mode) {
  fMuonModeString = mode;
}

void B02PrimaryGeneratorAction::SetMuonImpactMode(const G4String& mode) {
  fMuonImpactModeString = ToString(ParseMuonImpactMode(mode));
}

void B02PrimaryGeneratorAction::SetFluxModel(const G4String& model) {
  fFluxModel = FluxModelFromString(model);
  MarkSamplerDirty();
}

void B02PrimaryGeneratorAction::SetSourcePlaneZ(G4double value) {
  fSourcePlaneZ = std::max(0.0, value);
  MarkSamplerDirty();
}

void B02PrimaryGeneratorAction::SetSourcePlaneLx(G4double value) {
  fSourcePlaneLx = std::max(0.0, value);
  MarkSamplerDirty();
}

void B02PrimaryGeneratorAction::SetSourcePlaneLy(G4double value) {
  fSourcePlaneLy = std::max(0.0, value);
  MarkSamplerDirty();
}

void B02PrimaryGeneratorAction::SetSourcePlaneAutoSize(G4bool value) {
  fSourcePlaneAutoSize = value;
  MarkSamplerDirty();
}

void B02PrimaryGeneratorAction::SetSourcePlaneMargin(G4double value) {
  fSourcePlaneMargin = std::max(0.0, value);
  MarkSamplerDirty();
}

void B02PrimaryGeneratorAction::SetTargetMargin(G4double value) {
  fTargetMargin = std::max(0.0, value);
}

void B02PrimaryGeneratorAction::SetEminGeV(G4double value) {
  fRequestedEminGeV = value;
  fEminGeV = std::max(1e-4, value);
  if (fEminGeV > fEmaxGeV) {
    fEmaxGeV = fEminGeV * 10.0;
  }
  fEffectiveEminGeV = fEminGeV;
  if (fEffectiveEminGeV != fRequestedEminGeV) {
    G4cout << "[TierB] Requested Emin=" << fRequestedEminGeV << " GeV but using "
           << fEffectiveEminGeV << " GeV (clamped)" << G4endl;
  }
  MarkSamplerDirty();
}

void B02PrimaryGeneratorAction::SetEmaxGeV(G4double value) {
  fRequestedEmaxGeV = value;
  fEmaxGeV = std::max(fEminGeV * 1.0001, value);
  fEffectiveEmaxGeV = fEmaxGeV;
  if (fEffectiveEmaxGeV != fRequestedEmaxGeV) {
    G4cout << "[TierB] Requested Emax=" << fRequestedEmaxGeV << " GeV but using "
           << fEffectiveEmaxGeV << " GeV (adjusted to exceed Emin)" << G4endl;
  }
  MarkSamplerDirty();
}

void B02PrimaryGeneratorAction::SetMuonChargeMode(const G4String& mode) {
  fMuonChargeMode = mode;
}

void B02PrimaryGeneratorAction::SetMuonChargeRatio(G4double value) {
  const G4double clamped = (value < 0.0) ? 0.0 : value;
  fMuonPlusToMinusRatio = clamped;
}

void B02PrimaryGeneratorAction::DefineCommands() {
  // Define /generator command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/generator/",
                                      "Primary generator control");

  fMuonMessenger = new G4GenericMessenger(this, "/sim/muon/",
                                          "Muon sampling controls");
  fMuonMessenger->DeclareMethod("mode", &B02PrimaryGeneratorAction::SetMuonImpactMode,
                                "Impact sampling mode: unbiased or targeted");
  fMuonMessenger->DeclareMethodWithUnit(
      "targetMargin", "cm", &B02PrimaryGeneratorAction::SetTargetMargin,
      "Additional margin around the CCD footprint when sampling targeted impact points.");

  fMessenger->DeclareProperty("radius", R, "Radius of the hemisphere (cm).");
  fMessenger->DeclareProperty("px", px,
                              "x-direction tangent plane size (cm, legacy).");
  fMessenger->DeclareProperty("py", py,
                              "y-direction tangent plane size (cm, legacy).");
  fMessenger->DeclarePropertyWithUnit("impactPlaneZ", "cm", fZImpactPlane,
                                      "Z position of the impact plane");
  fMessenger->DeclareProperty("SmithActivation", fUseCosmicMuons,
                              "Enable or disable cosmic muon primaries");
  fMessenger->DeclarePropertyWithUnit(
      "thetaMax", "deg", fThetaMaxRad,
      "Maximum zenith angle for cosmic muons (used by Tier B)");
  fMessenger->DeclareProperty("muonEnergyGeV", fMuonEnergyGeV,
                              "Fixed muon kinetic energy in GeV");
  fMessenger->DeclareProperty("useFixedEnergy", fUseFixedEnergy,
                              "Use fixed muon kinetic energy instead of sampling");

  fMessenger->DeclareMethod("muonMode", &B02PrimaryGeneratorAction::SetMuonMode,
                            "Muon source mode: forced_footprint or tierB_plane_flux");
  fMessenger->DeclareMethod("fluxModel",
                            &B02PrimaryGeneratorAction::SetFluxModel,
                            "Flux model: guan2015 (default) or pdg_gaisser");
  fMessenger->DeclareMethodWithUnit(
      "sourcePlaneZ", "cm", &B02PrimaryGeneratorAction::SetSourcePlaneZ,
      "Z position of the source plane (above the impact plane).");
  fMessenger->DeclareMethodWithUnit(
      "sourcePlaneLx", "cm", &B02PrimaryGeneratorAction::SetSourcePlaneLx,
      "Manual x-extent of the source plane (used when autoSize is false).");
  fMessenger->DeclareMethodWithUnit(
      "sourcePlaneLy", "cm", &B02PrimaryGeneratorAction::SetSourcePlaneLy,
      "Manual y-extent of the source plane (used when autoSize is false).");
  fMessenger->DeclareMethod(
      "sourcePlaneAutoSize",
      &B02PrimaryGeneratorAction::SetSourcePlaneAutoSize,
      "Enable auto-sizing of the source plane using thetaMax and margins.");
  fMessenger->DeclareMethodWithUnit(
      "sourcePlaneMargin", "cm",
      &B02PrimaryGeneratorAction::SetSourcePlaneMargin,
      "Additional margin applied when auto sizing the source plane.");
  fMessenger->DeclareMethod("EminGeV", &B02PrimaryGeneratorAction::SetEminGeV,
                            "Minimum energy (GeV) for flux sampling");
  fMessenger->DeclareMethod("EmaxGeV", &B02PrimaryGeneratorAction::SetEmaxGeV,
                            "Maximum energy (GeV) for flux sampling");
  fMessenger->DeclareMethod("muonChargeMode",
                            &B02PrimaryGeneratorAction::SetMuonChargeMode,
                            "Muon charge selection: equal (default), fixedRatio, plusOnly, minusOnly");
  fMessenger->DeclareMethod("muonPlusToMinusRatio",
                            &B02PrimaryGeneratorAction::SetMuonChargeRatio,
                            "Muon charge ratio N(mu+)/N(mu-) used when muonChargeMode=fixedRatio (default 1.25)");
}
