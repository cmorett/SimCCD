#ifndef B02AnalysisManager_h
#define B02AnalysisManager_h 1

//Required explicitly since G4 version 10.0
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#ifdef G4ANALYSIS_USE

class G4Run;
class G4Event;
class G4Step;

namespace AIDA {
  class IAnalysisFactory;
  class ITree;
  class IHistogram1D;
  class ITuple;
}

class B02AnalysisManager {
public:
  B02AnalysisManager(AIDA::IAnalysisFactory*);
  virtual ~B02AnalysisManager();
public:
  virtual void BeginOfRun(const G4Run*); 
  virtual void EndOfRun(const G4Run*); 
  virtual void BeginOfEvent(const G4Event*); 
  virtual void EndOfEvent(const G4Event*); 
  virtual void Step(const G4Step*);
private:
  int fBarCollID;                
  double fEpri;
	double fpx, fpy, fpz, fp, fth, fph;
  AIDA::IAnalysisFactory* fAIDA;
  AIDA::ITree* fTree;
	AIDA::IHistogram1D* fThHis;
  AIDA::IHistogram1D* fEpriHis;
  AIDA::IHistogram1D* fEhitBar;
  AIDA::IHistogram1D* fXhitBar;
  AIDA::ITuple* fHitTuple;
  AIDA::ITuple* fThTuple;	
  AIDA::ITuple* fEvtTuple;
};

#else

class B02AnalysisManager;

#endif

#endif
