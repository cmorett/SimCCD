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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



B02PrimaryGeneratorAction::B02PrimaryGeneratorAction( 
							   B02DetectorConstruction* myDC)
  :myDetector(myDC), rndmFlag("off")
 
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

     if (doSmith) {
	
	//double R = 550.0;

	//double px = 480.0;
	//double py = 480.0;

	double theta;
	
	double h = 61.5;

	bool a = true;
	while (a) {
		double p = G4RandFlat::shoot(0.0, 1.0);
		theta = G4RandFlat::shoot(0.0, 3.1416/2);
		double q = pow(cos(theta), 2)*sin(theta);   // Muons
//		double q = sin(theta)*(0.0545+pow(cos(theta), 2.5))/1.0545;   // Electrons
//		double q = sin(theta)*(0.0266+pow(cos(theta), 2.6))/1.0266;   // Positrons
//		double q = sin(theta)*(0.0402+pow(cos(theta), 3.0))/1.0402;   // Photons
//		double q = sin(theta)*(0.276+pow(cos(theta), 2.4))/1.276;     // Neutrons
		if (p<=q) a = false;
		}

	double phi = G4RandFlat::shoot(0.0, 2*3.1416);
	// new R
	
	double X = R*sin(theta)*cos(phi);
	double Y = R*sin(theta)*sin(phi);
	double Z = R*cos(theta);
	// down the hemisphere
	double Z1 = R*cos(theta)-h;
	
	double u = G4RandFlat::shoot(-px/2, px/2);
	double v = G4RandFlat::shoot(-py/2, py/2);
	double x0 = X+u*cos(theta)*cos(phi)-v*sin(phi);
	double y0 = Y+u*cos(theta)*sin(phi)+v*cos(phi);
	double z0 = Z-u*sin(theta);
	//double z0 = Z;
	
        particleGun->SetParticlePosition(G4ThreeVector(x0*cm,y0*cm,z0*cm));
        particleGun->SetParticleMomentumDirection(G4ThreeVector(-sin(theta)*cos(phi),-sin(theta)*sin(phi),-cos(theta)));


// *****   Muons   *****

	// Smith & Duller

	// double Emin = 10;
	// double Emax = pow(10,5);

	double Emin = -1;
	double Emax = 5;

	const int ee = 10000;
	double ES[ee];
	// double dE_log = (log10(Emax)-log10(Emin))/ee;
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
		// Eu = pow(10, log10(Emin)+j*dE_log);
		Eu = pow(10, Emin+j*dE_log);
		Epu = (Eu+au*y0u*(1.0/cos(theta)-0.100))/ru;
		Pmu = pow(0.100*cos(theta)*(1-(au*(y0u/cos(theta)-100)/(ru*Epu))),(Bmu/((ru*Epu+100*au)*cos(theta))));
		ES[j] = Au*(pow(Epu,-gu))*Pmu*lpu*bu*jpu/(Epu*cos(theta)+bu*jpu);
		}

	int nbins = ee;
	G4RandGeneral GenDist(ES,nbins);          // Distribucion de energias
	// double E = pow(10, log10(Emin) + (GenDist.shoot())*(log10(Emax)-log10(Emin)));   // Sampleo de la energia
	double E = pow(10, Emin + (GenDist.shoot())*(Emax-Emin));   // Sampleo de la energia

       particleGun->SetParticleEnergy(E*MeV);
       particleGun->GeneratePrimaryVertex(anEvent);

      } // if doSmith


    else {
    
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
  fMessenger->DeclareProperty("SmithActivation", doSmith, "on or off smith");
  
}










