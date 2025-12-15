//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "B02PrimaryGeneratorAction.hh"

#include "B02DetectorConstruction.hh"
#include "B02PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"
#include "G4GeneralParticleSource.hh"
#include "G4PhysicalConstants.hh"

#include <cmath>
#include <algorithm>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



B02PrimaryGeneratorAction::B02PrimaryGeneratorAction( 
							   B02DetectorConstruction* myDC)
  :myDetector(myDC), rndmFlag("off"), fUseCosmicMuons(true)
 
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new B02PrimaryGeneratorMessenger(this);

// default particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("mu-");    
  particleGun->SetParticleDefinition(particle);

// define commands for this class
  DefineCommands();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02PrimaryGeneratorAction::~B02PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("mu-"));

     if (fUseCosmicMuons) {
	// Sample cos(theta) with pdf ~ mu^2 on [muMin, 1] via inverse CDF.
	const G4double muMin = std::cos(fThetaMaxRad);
	const G4double muMin3 = muMin * muMin * muMin;
	const G4double mu = std::cbrt(G4UniformRand() * (1.0 - muMin3) + muMin3);
	const G4double cosTheta = mu;
	const G4double phi = G4RandFlat::shoot(0.0, twopi);
	const G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
	const G4double tanTheta = (mu > 0.0) ? sinTheta / mu : 0.0;

	const G4double xImp = G4RandFlat::shoot(-px / 2.0, px / 2.0);
	const G4double yImp = G4RandFlat::shoot(-py / 2.0, py / 2.0);
	const G4double zImp = fZImpactPlane;
	const G4double z0 = zImp + R;
	const G4double deltaZ = z0 - zImp;
	const G4double cosPhi = std::cos(phi);
	const G4double sinPhi = std::sin(phi);
	const G4double x0 = xImp - deltaZ * tanTheta * cosPhi;
	const G4double y0 = yImp - deltaZ * tanTheta * sinPhi;

	fMuonX0 = x0 * cm;
	fMuonY0 = y0 * cm;
	fMuonZ0 = z0 * cm;
	fMuonXImp = xImp * cm;
	fMuonYImp = yImp * cm;
	fMuonZImp = zImp * cm;
	fMuonTheta = std::acos(mu);
	fMuonPhi = phi;

	particleGun->SetParticlePosition(G4ThreeVector(fMuonX0, fMuonY0, fMuonZ0));
	particleGun->SetParticleMomentumDirection(G4ThreeVector(sinTheta * cosPhi,
	                                                       sinTheta * sinPhi,
	                                                       -mu));


// *****   Muons   *****

	G4double kineticEnergy = fMuonEnergyGeV * GeV;
	if (!fUseFixedEnergy) {
	  // Smith & Duller energy spectrum
	  double Emin = -1;
	  double Emax = 5;

	  const int ee = 10000;
	  double ES[ee];
	  double dE_log = (Emax-Emin)/ee;

	  double Eu;            // Variable de energia cinetica
	  double Au = 2e9;                      // Parametros de la funcion de Smith
	  double gu = 2.645;                    // ...
	  double ru = 0.76;                     // ...
	  double au = 2.5;
	  double y0u = 1000.0;
	  double bmu = 0.80;
	  double cu = 299792458.0e2;
	  double mmu = 105.7/pow(cu,2);
	  double t0mu = 2.2e-6;
	  double r0u = 0.00129;
	  double Epu;
	  double Bmu = bmu*mmu*y0u*cu/(t0mu*r0u);
	  double Pmu;
	  double lpu = 120.0;
	  double bu = 0.771;
	  double mpu = 139.6/pow(cu,2);
	  double t0pu = 2.6e-8;
	  double jpu = mpu*y0u*cu/(t0pu*r0u);
	
	  for (int j=0; j<ee; j++) {    // Construye la funcion de Smith en un arreglo
		Eu = pow(10, Emin+j*dE_log);
		Epu = (Eu+au*y0u*(1.0/cosTheta-0.100))/ru;
		Pmu = pow(0.100*cosTheta*(1-(au*(y0u/cosTheta-100)/(ru*Epu))),(Bmu/((ru*Epu+100*au)*cosTheta)));
		ES[j] = Au*(pow(Epu,-gu))*Pmu*lpu*bu*jpu/(Epu*cosTheta+bu*jpu);
	  }

	  int nbins = ee;
	  G4RandGeneral GenDist(ES,nbins);          // Distribucion de energias
	  double E = pow(10, Emin + (GenDist.shoot())*(Emax-Emin));   // Sampleo de la energia
	  kineticEnergy = E*MeV;
	}

       particleGun->SetParticleEnergy(kineticEnergy);
       particleGun->GeneratePrimaryVertex(anEvent);

      } // if fUseCosmicMuons


    else {
    
       if (fUseFixedEnergy) {
         particleGun->SetParticleEnergy(fMuonEnergyGeV * GeV);
       }
       particleGun->GeneratePrimaryVertex(anEvent);   
       //particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.*cm));
       //particleGun->SetParticleEnergy(10.0*MeV);
       //G4double zposition = -0.5*(myDetector->GetWorldFullLength());
       //particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm, 120*cm));
 
    }


} 

void B02PrimaryGeneratorAction::DefineCommands()
{
  // Define /generator command directory using generic messenger class
  fMessenger
    = new G4GenericMessenger(this,
                             "/generator/",
                             "Primary generator control");
  fMessenger->DeclareProperty("radius", R, "Radius of the hemisphere");
  fMessenger->DeclareProperty("px", px, "x-direction tangent plane");
  fMessenger->DeclareProperty("py", py, "y-direction tangent plane");
  fMessenger->DeclareProperty("impactPlaneZ", fZImpactPlane, "Z position of the impact plane");
  fMessenger->DeclareProperty("SmithActivation", fUseCosmicMuons, "Enable or disable cosmic muon primaries");
  fMessenger->DeclarePropertyWithUnit("thetaMax", "deg", fThetaMaxRad, "Maximum zenith angle for cosmic muons");
  fMessenger->DeclareProperty("muonEnergyGeV", fMuonEnergyGeV, "Fixed muon kinetic energy in GeV");
  fMessenger->DeclareProperty("useFixedEnergy", fUseFixedEnergy, "Use fixed muon kinetic energy instead of sampling");
  
}
