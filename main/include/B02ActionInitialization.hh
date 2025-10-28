//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B02/include/ActionInitialization.hh
/// \brief Definition of the B02::ActionInitialization class

#ifndef B02ActionInitialization_h
#define B02ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "B02DetectorConstruction.hh"
#include "B02PrimaryGeneratorAction.hh"
#include "B02RunAction.hh"
#include "B02EventAction.hh"
#include "B02SteppingAction.hh"
#include "B02ActionInitialization.hh"

/// Action initialization class.


class B02ActionInitialization : public G4VUserActionInitialization
{
  public:
      B02ActionInitialization();
     ~B02ActionInitialization();
    
  virtual void BuildForMaster() const ;
  virtual void Build() const;
  // void Build() const override;   
};






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
