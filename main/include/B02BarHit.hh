//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef B02BarHit_h
#define B02BarHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class B02BarHit : public G4VHit
{
  public:

      B02BarHit();
     ~B02BarHit();
      B02BarHit(const B02BarHit&);
      const B02BarHit& operator=(const B02BarHit&);
      G4int operator==(const B02BarHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  
      void SetTrackID  (G4int track)      { trackID = track; };
      void SetEdep     (G4double de)      { edep = de; };
      void SetEvis     (G4double ve)      { evis = ve; };
      void SetPos      (G4ThreeVector xyz){ pos = xyz; };
      
      G4int GetTrackID()    { return trackID; };
      G4double GetEdep()    { return edep; };      
      G4double GetEvis()    { return evis; };      
      G4ThreeVector GetPos(){ return pos; };
      
  private:
  
      G4int         trackID;
      G4double      edep;
      G4double      evis;
      G4ThreeVector pos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<B02BarHit> B02BarHitsCollection;

extern G4Allocator<B02BarHit> B02BarHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* B02BarHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) B02BarHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void B02BarHit::operator delete(void *aHit)
{
  B02BarHitAllocator.FreeSingle((B02BarHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
