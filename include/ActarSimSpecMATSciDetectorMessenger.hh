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

class ActarSimSpecMATSciDetectorMessenger: public G4UImessenger {
private:
  ActarSimSpecMATSciDetectorConstruction* ActarSimSpecMATSciDetector; ///< Pointer to main sci detector class
  
  G4UIdirectory*           detDir;                      ///< Directory in messenger structure
  G4UIcmdWithoutParameter* printCmd;                    ///< Prints geometry
  
  G4UIcmdWithAString* sciSpecMATMaterialCmd;            ///< Select Material the scintillators are made of
  G4UIcmdWithAString* reflectorSpecMATMaterialCmd;      ///< Select Material the reflectors are made of
  G4UIcmdWithAString* windowSpecMATMaterialCmd;         ///< Select Material the windows are made of
  G4UIcmdWithAString* housingSpecMATMaterialCmd;        ///< Select Material the housing is made of
  G4UIcmdWithAString* vacuumFlangeSpecMATMaterialCmd;   ///< Select Material the vacuum flange is made of
  G4UIcmdWithAString* insulationTubeSpecMATMaterialCmd; ///< Select Material the insulation tube is made of

  G4UIcmdWithAString* vacuumChamberIncludedFlagCmd;     ///< Includes the vacuum chamber and associated elements (default on).
  G4UIcmdWithAString* insulationTubeIncludedFlagCmd;    ///< Includes the insulation tube and associated elements (default on).

  G4UIcmdWithAnInteger* nbSegmentsCmd;           ///< Sets the number of scintillator segments
  G4UIcmdWithAnInteger* nbCrystInSegmentRowCmd;  ///< Sets the number of scintillator crystals in a segment row
  G4UIcmdWithAnInteger* nbCrystInSegmentColumnCmd;  ///< Sets the number of scintillator crystals in a segment column

  //G4UIcmdWithADoubleAndUnit* chamberSizeXCmd;     ///< Sets the half-length X dimension of the Gas Chamber.
  //G4UIcmdWithADoubleAndUnit* yBoxHalfLengthCmd;   ///< Sets the y half length of the sci detectors box

public:
  ActarSimSpecMATSciDetectorMessenger(ActarSimSpecMATSciDetectorConstruction* );
  ~ActarSimSpecMATSciDetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  //G4String GetCurrentValue(G4UIcommand*);
};
#endif
