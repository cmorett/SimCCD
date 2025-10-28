#ifdef G4ANALYSIS_USE

#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4HCofThisEvent.hh"
#include "G4PrimaryVertex.hh"

 #include <AIDA/IAnalysisFactory.h>
 #include <AIDA/ITreeFactory.h>
 #include <AIDA/ITupleFactory.h>
 #include <AIDA/IHistogramFactory.h>
 #include <AIDA/ITree.h>
 #include <AIDA/IHistogram1D.h>
 #include <AIDA/ITuple.h>

#include "B02BarHit.hh"
#include "B02AnalysisManager.hh"

 B02AnalysisManager::B02AnalysisManager(AIDA::IAnalysisFactory* aAIDA)
:fBarCollID(-1)
,fAIDA(aAIDA)
,fTree(0)
,fEhitBar(0)
,fXhitBar(0)
,fHitTuple(0)
,fEvtTuple(0)
{

  // Could fail if no AIDA implementation found :
  if(!fAIDA) {
    G4cout << "AIDA analysis factory not found." << G4endl;
    return;
  }

  AIDA::ITreeFactory* treeFactory = fAIDA->createTreeFactory();
  if(!treeFactory) return;

  // Create a tree-like container to handle histograms.
  // This tree is associated to a B02.root file.  
  //std::string opts = "compress=yes";
  //std::string opts = "compress=no";
  //fTree = treeFactory->create("B02.aida","xml",false,true,opts);
 
 
  std::string opts = "export=root";
  char const* coutputfilename = getenv( "OUTROOTFILE" );
  std::string soutfilename ( coutputfilename );
  fTree = treeFactory->create(soutfilename,"ROOT",false,true,opts);
//  fTree = treeFactory->create("output/muonsconcrete_340.root","ROOT",false,true,opts);
//  fTree = treeFactory->create("output/muonsconcreteEXPACS_180.root","ROOT",false,true,opts);
//  fTree = treeFactory->create("output/muonsconcreterig_003.root","ROOT",false,true,opts);
//  fTree = treeFactory->create("output/neutronsconcreteEXPACS_031.root","ROOT",false,true,opts);
//  fTree = treeFactory->create("output/test.root","ROOT",false,true,opts);

  // Factories are not "managed" by an AIDA analysis system.
  // They must be deleted by the AIDA user code.
  delete treeFactory; 

  if(!fTree) return;

  //fTree->mkdir("histograms");
  //fTree->cd("histograms");

  // Create an histo factory that will create histo in the tree :
  AIDA::IHistogramFactory* histoFactory = 
    fAIDA->createHistogramFactory(*fTree);

  if(histoFactory) {
    fEpriHis = histoFactory->createHistogram1D("Epri",100,0,100);
    fThHis = histoFactory->createHistogram1D("ThAng",100,0,2);
    fEhitBar = histoFactory->createHistogram1D("EBar",100,0,100);
    fXhitBar = histoFactory->createHistogram1D("XBar",100,-1500,1500);
    delete histoFactory;
  }

  //fTree->cd("..");
  //fTree->mkdir("tuples");
  //fTree->cd("tuples");

  G4cout <<"dentro ..."<<G4endl; 

  // Get a tuple factory :
 
  AIDA::ITupleFactory* tupleFactory = fAIDA->createTupleFactory(*fTree);
  if(tupleFactory) {

    // Create a Hit tuple :
    fHitTuple = tupleFactory->create("B02Hit","B02Hit",
     "int evtID,double XhitBar,double YhitBar,double ZhitBar,double RhitBar,double EhitBar,double WhitBar");
    // Create a Evt tuple :
    fEvtTuple = tupleFactory->create("B02Evt","B02Evt",
    "int evtID,double EevtBar,double WevtBar, double GevtBar, double ADCevtBar, double EevtPri, double thetaPri, double phiPri, int nHitBar");
		//fThTuple = tupleFactory->create("B02Th", "B02Th", "double ThAng");

    delete tupleFactory;
  }

}

B02AnalysisManager::~B02AnalysisManager() {

}

void B02AnalysisManager::BeginOfRun(const G4Run* aRun){
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

}

void B02AnalysisManager::EndOfRun(const G4Run*){
  if(fTree) fTree->commit();
  if(fEhitBar) {
    G4cout<<"Histo: EpriHis (MeV): mean: "<<fEpriHis->mean()/MeV << " rms: "<<fEpriHis->rms()/MeV <<G4endl;
		G4cout<<"Histo: ThAng: mean: "<<fThHis->mean() << " rms: "<<fThHis->rms() <<G4endl;
    G4cout<<"Histo: EhitBar (MeV): mean: "<<fEhitBar->mean()/MeV << " rms: "<<fEhitBar->rms()/MeV <<G4endl;
    G4cout<<"Histo: XhitBar (cm) : mean: "<<fXhitBar->mean()/cm  << " rms: "<<fXhitBar->rms()/cm  <<G4endl;
  }
}

void B02AnalysisManager::BeginOfEvent(const G4Event* aEvent){
  if(fBarCollID==-1) {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    fBarCollID = SDman->GetCollectionID("barCollection");
  } 

  if(!fEpriHis)   return; // No histo booked !

  G4PrimaryVertex* primaryVertex = aEvent->GetPrimaryVertex();
  G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
  fEpri = primaryParticle->GetKineticEnergy();
	fpx = primaryParticle->GetPx();
	fpy = primaryParticle->GetPy();
	fpz = primaryParticle->GetPz();
	fp = sqrt(fpx*fpx+fpy*fpy+fpz*fpz);
	fth = acos(-fpz/fp);
	fph = atan2(fpy,fpx);
	fThHis->fill(cos(fth));
  fEpriHis->fill(fEpri/GeV);

}

// ***** Gaussian function *****
/*
double gaussfunction(double xg) {
	double sig = 0.058*sqrt(20*WevtBar);
	double gauss = exp(-pow(xg,2)/pow(sig,2));
	return xg+10;
}
*/
// ***** ***** ***** ***** *****

void B02AnalysisManager::EndOfEvent(const G4Event* aEvent){

  if(!fEhitBar) return; // No histo booked !
  if(!fHitTuple) return; // No tuple booked !
  if(!fEvtTuple) return; // No tuple booked !

  G4int event_id = aEvent->GetEventID();  
  //
  // get number of stored trajectories
  //
  G4TrajectoryContainer* trajectoryContainer = aEvent->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  //
  // periodic printing
  //
  if (event_id < 100 || event_id%100 == 0) {
    G4cout << ">>> Event " << aEvent->GetEventID() << G4endl;
    G4cout << "    " << n_trajectories 
	   << " trajectories stored in this event." << G4endl;
  }

  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  B02BarHitsCollection* BHC = 
    HCE ? (B02BarHitsCollection*)(HCE->GetHC(fBarCollID)) : 0;

  if(BHC) {
    G4int n_hit = BHC->entries();
    G4double EevtBar = 0.;
    G4double WevtBar = 0.;
    G4double GevtBar = 0.;
    G4double ADCevtBar = 0.;
    for (G4int i=0;i<n_hit;i++) {
      G4double EhitBar = (*BHC)[i]->GetEdep();
      G4double WhitBar = (*BHC)[i]->GetEvis();
      G4double XhitBar = (*BHC)[i]->GetPos().getX();
      G4double YhitBar = (*BHC)[i]->GetPos().getY();
      G4double ZhitBar = (*BHC)[i]->GetPos().getZ();
      G4double RhitBar = sqrt(XhitBar*XhitBar+YhitBar*YhitBar+ZhitBar*ZhitBar);
      fEhitBar->fill(EhitBar/MeV);
      fXhitBar->fill(XhitBar/cm);

      fHitTuple->fill(0,event_id);
      fHitTuple->fill(1,XhitBar/cm);
      fHitTuple->fill(2,YhitBar/cm);
      fHitTuple->fill(3,ZhitBar/cm);
      fHitTuple->fill(4,RhitBar/cm);
      fHitTuple->fill(5,EhitBar/MeV);
      fHitTuple->fill(6,WhitBar/MeV);
      fHitTuple->addRow();

      EevtBar += EhitBar;
      WevtBar += WhitBar;
    }
    if (n_hit==0) EevtBar = -5.;

		//double res = 0.05;                      //
		//double sig = res*sqrt(550.0*WevtBar);   //
		//double xg = G4RandGauss::shoot(0, sig); //
		//GevtBar = WevtBar+xg;

		//double p1 =  0.0;      // No offset
		//double p2 = 7.3;      // 15 PE/MEV
		//double p3 = 1.0/2000; // NL turns on at ~550

		//ADCevtBar = p1+p2*GevtBar/(1+p3*GevtBar);

          fEvtTuple->fill(0,event_id);
	  fEvtTuple->fill(1,EevtBar/MeV);
	  fEvtTuple->fill(2,WevtBar/MeV);
	  //fEvtTuple->fill(3,GevtBar/MeV);
	  //fEvtTuple->fill(4,ADCevtBar);
          fEvtTuple->fill(5,fEpri/GeV);
	  fEvtTuple->fill(6,fth);
	  fEvtTuple->fill(7,fph);
	  fEvtTuple->fill(8, n_hit);
    fEvtTuple->addRow();

  }
  
}

void B02AnalysisManager::Step(const G4Step*){}

#endif
