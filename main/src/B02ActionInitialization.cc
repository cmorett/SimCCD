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
/// \file B5/src/ActionInitialization.cc
/// \brief Implementation of the B02::ActionInitialization class

#include "B02ActionInitialization.hh"
#include "B02PrimaryGeneratorAction.hh"
#include "B02RunAction.hh"
#include "B02EventAction.hh"
#include "B02DetectorConstruction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B02ActionInitialization::B02ActionInitialization()
{}


B02ActionInitialization::~B02ActionInitialization()
{}

void B02ActionInitialization::BuildForMaster() const
{
   B02EventAction *eventAction = new B02EventAction();
   SetUserAction(new B02RunAction(eventAction));
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B02ActionInitialization::Build() const
{

B02DetectorConstruction *detConstruction = new B02DetectorConstruction();

  B02PrimaryGeneratorAction *genAction = new B02PrimaryGeneratorAction(detConstruction);
  SetUserAction(genAction); 
    
  B02EventAction *eventAction = new B02EventAction();
  SetUserAction(eventAction);
  
  B02RunAction *runAction = new B02RunAction(eventAction);
  SetUserAction(runAction);
 
  B02SteppingAction *steppingAction = new B02SteppingAction(eventAction);
  SetUserAction(steppingAction);     
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


