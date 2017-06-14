// - AUTHOR: Hector Alvarez-Pol 06/2017
/******************************************************************
 * Copyright (C) 2005-2017, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/

#ifndef ActarSimSpecMATSciDetectorMessenger_h
#define ActarSimSpecMATSciDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ActarSimSpecMATSciDetectorConstruction;
class ActarSimPrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;

class ActarSimSpecMATSciDetectorMessenger: public G4UImessenger {
private:
  ActarSimSpecMATSciDetectorConstruction* ActarSimSciDetector; ///< Pointer to main sci detector class

  G4UIdirectory*             detDir;                    ///< Directory in messenger structure
  G4UIcmdWithoutParameter*   printCmd;                  ///< Prints geometry

public:
  ActarSimSpecMATSciDetectorMessenger(ActarSimSpecMATSciDetectorConstruction* );
  ~ActarSimSpecMATSciDetectorMessenger();

  void SetNewValue(G4UIcommand*, G4String);
  //G4String GetCurrentValue(G4UIcommand*);
};
#endif
