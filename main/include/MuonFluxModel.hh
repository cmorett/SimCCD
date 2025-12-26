// Utility for cosmic muon flux evaluation and sampling.
// Implements the Modified Gaisser (Guan et al. 2015) parameterisation and
// provides a lightweight CDF-based sampler over (cos(theta), energy).
#ifndef MUONFLUXMODEL_HH
#define MUONFLUXMODEL_HH

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

enum class MuonFluxModelType {
  Guan2015,
  PDGGaisser,
};

MuonFluxModelType FluxModelFromString(const std::string& name);
std::string ToString(MuonFluxModelType model);

// Effective zenith angle that accounts for Earth curvature (Guan et al. Eq. 2).
double cosThetaStar(double cosTheta);

// Differential intensity (cm^-2 s^-1 sr^-1 GeV^-1) using Modified Gaisser.
double dIdEdOmega_Guan2015(double E_GeV, double cosTheta);

// Standard PDG Gaisser form (no low-energy or curvature corrections).
double dIdEdOmega_PDG_Gaisser(double E_GeV, double cosTheta);

// Tabulated sampler for (cos(theta), E) following I(E,theta) * cos(theta).
class MuonFluxSampler {
 public:
  MuonFluxSampler();

  void SetFluxModel(MuonFluxModelType model);
  void SetEnergyRange(double eminGeV, double emaxGeV);
  void SetThetaMax(double thetaMaxRad);
  void SetGridSize(std::size_t nCos, std::size_t nEnergy);

  // Build tables if dirty; idempotent otherwise.
  void RebuildIfNeeded();

  bool IsReady() const { return fReady; }

  // Returns {energy_GeV, cosTheta}. Uses G4UniformRand for reproducibility.
  std::pair<double, double> SampleCosThetaAndEnergy();

  // Integrated flux over phi in the configured phase space:
  // 2*pi * ∫_{cosMin}^{1} cosθ dcosθ ∫_{Emin}^{Emax} I(E,θ) dE
  // Units: cm^-2 s^-1.
  double FluxIntegral() const { return fFluxIntegral; }
  double FluxIntegralQuadrature() const;
  double FluxIntegralMonteCarlo(std::size_t nSamples = 200000) const;
  double FluxIntegralQuadratureNoPhi() const;

  double MinCosTheta() const { return fCosMin; }
  double MaxCosTheta() const { return 1.0; }
  double Emin() const { return fEminGeV; }
  double Emax() const { return fEmaxGeV; }

 private:
  void BuildTables();
  double EvaluateIntensity(double E_GeV, double cosTheta) const;
  double SampleEnergyForBin(std::size_t cosIdx) const;

  MuonFluxModelType fModel;
  double fEminGeV;
  double fEmaxGeV;
  double fThetaMaxRad;
  std::size_t fNCos;
  std::size_t fNEnergy;
  bool fDirty;
  bool fReady;

  double fCosMin;
  double fFluxIntegral;  // Includes 2*pi over phi

  std::vector<double> fCosGrid;
  std::vector<double> fEnergyGrid;
  std::vector<double> fCosCDF;  // Normalized 0-1
  std::vector<double> fEnergyIntegral;
  std::vector<std::vector<double>> fEnergyCDF;  // [cosIdx][eIdx]
};

#endif  // MUONFLUXMODEL_HH
