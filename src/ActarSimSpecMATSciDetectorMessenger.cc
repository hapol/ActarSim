// - AUTHOR: Hector Alvarez-Pol 06/2017
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/
//////////////////////////////////////////////////////////////////
/// \class ActarSimSciDetectorMessenger
/// Messenger of the Scintillator detector
/////////////////////////////////////////////////////////////////

#include "ActarSimSpecMATSciDetectorMessenger.hh"

#include "G4RunManager.hh"

#include "ActarSimSpecMATSciDetectorConstruction.hh"
#include "ActarSimPrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//////////////////////////////////////////////////////////////////
/// Constructor with complete functionality
ActarSimSpecMATSciDetectorMessenger::
ActarSimSpecMATSciDetectorMessenger(ActarSimSpecMATSciDetectorConstruction* ActarSimSpecMATSciDet)
  :ActarSimSciDetector(ActarSimSpecMATSciDet) {

  detDir = new G4UIdirectory("/ActarSim/det/SpecMAT/");
  detDir->SetGuidance("SpecMAT scintillator detector control");

  printCmd = new G4UIcmdWithoutParameter("/ActarSim/det/SpecMAT/print",this);
  printCmd->SetGuidance("Prints geometry.");
  printCmd->AvailableForStates(G4State_Idle);

}

//////////////////////////////////////////////////////////////////
/// Destructor
ActarSimSpecMATSciDetectorMessenger::~ActarSimSpecMATSciDetectorMessenger() {
  delete detDir;
  delete printCmd;
}

//////////////////////////////////////////////////////////////////
/// Setting the new values and connecting to actions
void ActarSimSpecMATSciDetectorMessenger::SetNewValue(G4UIcommand* command,
						      G4String newValue) {
  if( command == printCmd )
    ActarSimSciDetector->PrintDetectorParameters();

}
