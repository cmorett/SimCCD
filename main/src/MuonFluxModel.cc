// Cosmic muon flux utilities (Modified Gaisser) and fast sampling tables.
#include "MuonFluxModel.hh"

#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

double Clamp01(double v) {
  return std::max(0.0, std::min(1.0, v));
}

}  // namespace

MuonFluxModelType FluxModelFromString(const std::string& name) {
  std::string lower(name);
  std::transform(lower.begin(), lower.end(), lower.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lower == "pdg_gaisser" || lower == "pdg" || lower == "gaisser") {
    return MuonFluxModelType::PDGGaisser;
  }
  // Default
  return MuonFluxModelType::Guan2015;
}

std::string ToString(MuonFluxModelType model) {
  switch (model) {
    case MuonFluxModelType::PDGGaisser:
      return "pdg_gaisser";
    case MuonFluxModelType::Guan2015:
    default:
      return "guan2015";
  }
}

double cosThetaStar(double cosTheta) {
  // Guan et al. 2015, Eq. (2). Parameters P1..P5.
  constexpr double P1 = 0.102573;
  constexpr double P2 = -0.068287;
  constexpr double P3 = 0.958633;
  constexpr double P4 = 0.0407253;
  constexpr double P5 = 0.817285;

  const double c = Clamp01(cosTheta);
  const double term = (c * c) + (P1 * P1) + P2 * std::pow(c, P3) + P4 * std::pow(c, P5);
  return std::sqrt(std::max(0.0, term));
}

double dIdEdOmega_Guan2015(double E_GeV, double cosTheta) {
  if (E_GeV <= 0.0) {
    return 0.0;
  }
  const double cStar = cosThetaStar(cosTheta);
  if (cStar <= 0.0) {
    return 0.0;
  }
  const double energyTerm = std::pow(E_GeV, -2.7);
  const double corr = 1.0 / (1.0 + 3.64 / (E_GeV * std::pow(cStar, 1.29)));
  const double denom1 = 1.0 + 1.1 * E_GeV * cStar / 115.0;
  const double denom2 = 1.0 + 1.1 * E_GeV * cStar / 850.0;
  return 0.14 * energyTerm * (1.0 / denom1 + 0.054 / denom2) * corr;
}

double dIdEdOmega_PDG_Gaisser(double E_GeV, double cosTheta) {
  if (E_GeV <= 0.0) {
    return 0.0;
  }
  const double c = Clamp01(cosTheta);
  const double energyTerm = std::pow(E_GeV, -2.7);
  const double denom1 = 1.0 + 1.1 * E_GeV * c / 115.0;
  const double denom2 = 1.0 + 1.1 * E_GeV * c / 850.0;
  return 0.14 * energyTerm * (1.0 / denom1 + 0.054 / denom2);
}

MuonFluxSampler::MuonFluxSampler()
    : fModel(MuonFluxModelType::Guan2015),
      fEminGeV(0.1),
      fEmaxGeV(1.0e4),
      fThetaMaxRad(CLHEP::halfpi),
      fNCos(120),
      fNEnergy(200),
      fDirty(true),
      fReady(false),
      fCosMin(std::cos(CLHEP::halfpi)),
      fFluxIntegral(0.0) {}

void MuonFluxSampler::SetFluxModel(MuonFluxModelType model) {
  if (fModel != model) {
    fModel = model;
    fDirty = true;
  }
}

void MuonFluxSampler::SetEnergyRange(double eminGeV, double emaxGeV) {
  const double emin = std::max(1e-4, eminGeV);
  const double emax = std::max(emin * 1.0001, emaxGeV);
  if (emin != fEminGeV || emax != fEmaxGeV) {
    fEminGeV = emin;
    fEmaxGeV = emax;
    fDirty = true;
  }
}

void MuonFluxSampler::SetThetaMax(double thetaMaxRad) {
  const double clamped = std::max(0.0, std::min(thetaMaxRad, static_cast<double>(CLHEP::halfpi)));
  if (std::abs(clamped - fThetaMaxRad) > 1e-9) {
    fThetaMaxRad = clamped;
    fCosMin = std::cos(fThetaMaxRad);
    fDirty = true;
  }
}

void MuonFluxSampler::SetGridSize(std::size_t nCos, std::size_t nEnergy) {
  const std::size_t cosBins = std::max<std::size_t>(2, nCos);
  const std::size_t eBins = std::max<std::size_t>(2, nEnergy);
  if (cosBins != fNCos || eBins != fNEnergy) {
    fNCos = cosBins;
    fNEnergy = eBins;
    fDirty = true;
  }
}

void MuonFluxSampler::RebuildIfNeeded() {
  if (!fDirty) {
    return;
  }
  BuildTables();
}

double MuonFluxSampler::EvaluateIntensity(double E_GeV, double cosTheta) const {
  switch (fModel) {
    case MuonFluxModelType::PDGGaisser:
      return dIdEdOmega_PDG_Gaisser(E_GeV, cosTheta);
    case MuonFluxModelType::Guan2015:
    default:
      return dIdEdOmega_Guan2015(E_GeV, cosTheta);
  }
}

void MuonFluxSampler::BuildTables() {
  fReady = false;
  fFluxIntegral = 0.0;

  const std::size_t nCos = std::max<std::size_t>(2, fNCos);
  const std::size_t nE = std::max<std::size_t>(2, fNEnergy);

  fCosGrid.assign(nCos, 0.0);
  fEnergyGrid.assign(nE, 0.0);
  fCosCDF.assign(nCos, 0.0);
  fEnergyIntegral.assign(nCos, 0.0);
  fEnergyCDF.assign(nCos, std::vector<double>(nE, 0.0));

  const double cosMin = Clamp01(fCosMin);
  const double dcos = (1.0 - cosMin) / static_cast<double>(nCos - 1);
  for (std::size_t i = 0; i < nCos; ++i) {
    fCosGrid[i] = cosMin + dcos * static_cast<double>(i);
  }

  const double logEmin = std::log10(fEminGeV);
  const double logEmax = std::log10(fEmaxGeV);
  const double dlogE = (logEmax - logEmin) / static_cast<double>(nE - 1);
  for (std::size_t j = 0; j < nE; ++j) {
    fEnergyGrid[j] = std::pow(10.0, logEmin + dlogE * static_cast<double>(j));
  }

  // Build conditional CDFs over energy for each cos bin.
  for (std::size_t ic = 0; ic < nCos; ++ic) {
    double cumulative = 0.0;
    fEnergyCDF[ic][0] = 0.0;
    for (std::size_t j = 1; j < nE; ++j) {
      const double e0 = fEnergyGrid[j - 1];
      const double e1 = fEnergyGrid[j];
      const double i0 = EvaluateIntensity(e0, fCosGrid[ic]);
      const double i1 = EvaluateIntensity(e1, fCosGrid[ic]);
      cumulative += 0.5 * (i0 + i1) * (e1 - e0);
      fEnergyCDF[ic][j] = cumulative;
    }
    fEnergyIntegral[ic] = cumulative;
    const double norm = cumulative;
    if (norm > 0.0) {
      for (double& v : fEnergyCDF[ic]) {
        v /= norm;
      }
      fEnergyCDF[ic].back() = 1.0;
    } else {
      // Fallback to log-uniform if the intensity is zero in this bin.
      for (std::size_t j = 0; j < nE; ++j) {
        fEnergyCDF[ic][j] = static_cast<double>(j) / static_cast<double>(nE - 1);
      }
      fEnergyIntegral[ic] = 0.0;
    }
  }

  // Marginal over cos(theta) using the integral over energy times cos(theta).
  double cosCumulative = 0.0;
  fCosCDF[0] = 0.0;
  for (std::size_t ic = 1; ic < nCos; ++ic) {
    const double u0 = fCosGrid[ic - 1];
    const double u1 = fCosGrid[ic];
    const double f0 = fEnergyIntegral[ic - 1] * u0;
    const double f1 = fEnergyIntegral[ic] * u1;
    cosCumulative += 0.5 * (f0 + f1) * (u1 - u0);
    fCosCDF[ic] = cosCumulative;
  }

  const double total = fCosCDF.back();
  if (total > 0.0) {
    for (double& v : fCosCDF) {
      v /= total;
    }
    fCosCDF.back() = 1.0;
    // Prefer the quadrature result for reporting; total is used for CDF normalization.
    fFluxIntegral = FluxIntegralQuadrature();
    fReady = true;
  } else {
    fFluxIntegral = 0.0;
    fReady = false;
  }

  fDirty = false;
}

double MuonFluxSampler::SampleEnergyForBin(std::size_t cosIdx) const {
  const std::size_t idx = std::min(cosIdx, fEnergyCDF.size() - 1);
  const std::vector<double>& cdf = fEnergyCDF[idx];
  if (cdf.empty()) {
    return fEminGeV;
  }

  const double r = G4UniformRand();
  auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
  std::size_t j = static_cast<std::size_t>(std::distance(cdf.begin(), it));
  if (j == 0) {
    j = 1;
  }
  if (j >= cdf.size()) {
    j = cdf.size() - 1;
  }

  const double c0 = cdf[j - 1];
  const double c1 = cdf[j];
  const double e0 = fEnergyGrid[j - 1];
  const double e1 = fEnergyGrid[j];
  const double denom = (c1 > c0) ? (c1 - c0) : 1.0;
  const double alpha = (denom > 0.0) ? (r - c0) / denom : 0.0;
  return e0 + alpha * (e1 - e0);
}

std::pair<double, double> MuonFluxSampler::SampleCosThetaAndEnergy() {
  RebuildIfNeeded();
  if (!fReady || fCosCDF.size() < 2) {
    return {fEminGeV, 1.0};
  }

  const double rCos = G4UniformRand();
  auto it = std::lower_bound(fCosCDF.begin(), fCosCDF.end(), rCos);
  std::size_t idx = static_cast<std::size_t>(std::distance(fCosCDF.begin(), it));
  if (idx == 0) {
    idx = 1;
  }
  if (idx >= fCosCDF.size()) {
    idx = fCosCDF.size() - 1;
  }

  const double cdf0 = fCosCDF[idx - 1];
  const double cdf1 = fCosCDF[idx];
  const double cos0 = fCosGrid[idx - 1];
  const double cos1 = fCosGrid[idx];
  const double denom = (cdf1 > cdf0) ? (cdf1 - cdf0) : 1.0;
  const double alpha = (denom > 0.0) ? (rCos - cdf0) / denom : 0.0;
  const double cosTheta = cos0 + alpha * (cos1 - cos0);

  const double energyGeV = SampleEnergyForBin(idx);
  return {energyGeV, cosTheta};
}

double MuonFluxSampler::FluxIntegralQuadratureNoPhi() const {
  // ∫_{cosMin}^{1} cosθ dcosθ ∫_{Emin}^{Emax} I(E,θ) dE (no phi factor).
  const std::size_t nCos = std::max<std::size_t>(120, fNCos);
  const std::size_t nE = std::max<std::size_t>(200, fNEnergy);
  const double cosMin = Clamp01(fCosMin);
  const double dcos = (1.0 - cosMin) / static_cast<double>(nCos - 1);
  const double logEmin = std::log10(fEminGeV);
  const double logEmax = std::log10(fEmaxGeV);
  const double dlogE = (logEmax - logEmin) / static_cast<double>(nE - 1);

  double integral = 0.0;
  for (std::size_t ic = 0; ic < nCos; ++ic) {
    const double c = cosMin + dcos * static_cast<double>(ic);
    const double weightCos = (ic == 0 || ic == nCos - 1) ? 0.5 : 1.0;
    for (std::size_t ie = 0; ie < nE; ++ie) {
      const double logE = logEmin + dlogE * static_cast<double>(ie);
      const double E = std::pow(10.0, logE);
      const double weightE = (ie == 0 || ie == nE - 1) ? 0.5 : 1.0;
      const double jacE = E * std::log(10.0);  // dE when integrating in log10E
      const double integrand = EvaluateIntensity(E, c) * c * jacE;
      integral += integrand * weightCos * weightE;
    }
  }

  integral *= dcos * dlogE;
  return integral;
}

double MuonFluxSampler::FluxIntegralQuadrature() const {
  return 2.0 * CLHEP::pi * FluxIntegralQuadratureNoPhi();
}

double MuonFluxSampler::FluxIntegralMonteCarlo(std::size_t nSamples) const {
  const double cosMin = Clamp01(fCosMin);
  const double logEmin = std::log10(fEminGeV);
  const double logEmax = std::log10(fEmaxGeV);
  const double dcos = (1.0 - cosMin);
  const double dlogE = (logEmax - logEmin);
  if (nSamples == 0) {
    return 0.0;
  }

  double sum = 0.0;
  for (std::size_t i = 0; i < nSamples; ++i) {
    const double ucos = G4UniformRand();
    const double cosTheta = cosMin + ucos * dcos;
    const double logE = logEmin + G4UniformRand() * dlogE;
    const double E = std::pow(10.0, logE);
    const double integrand = EvaluateIntensity(E, cosTheta) * cosTheta * (E * std::log(10.0));
    sum += integrand;
  }

  const double domainVolume = dcos * dlogE * (2.0 * CLHEP::pi);
  return domainVolume * (sum / static_cast<double>(nSamples));
}
