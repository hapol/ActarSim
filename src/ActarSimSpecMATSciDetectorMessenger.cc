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
/// Constructor 
/// command included in this SpecMATSciDetectorMessenger:
/// - /ActarSim/det/SpecMAT/print
/// - /ActarSim/det/SpecMAT/setSciMat
/// - /ActarSim/det/SpecMAT/setReflectorMat
/// - /ActarSim/det/SpecMAT/setWindowMat
/// - /ActarSim/det/SpecMAT/setHousingMat
/// - /ActarSim/det/SpecMAT/setVacuumFlangeMat
/// - /ActarSim/det/SpecMAT/setInsulationTubeMat
/// - /ActarSim/det/SpecMAT/vacuumChamberIncludedFlag
/// - /ActarSim/det/SpecMAT/insulationTubeIncludedFlagCmd
/// - /ActarSim/det/SpecMAT/nbSegmentsCmd
/// - /ActarSim/det/SpecMAT/nbCrystInSegmentRowCmd
/// - /ActarSim/det/SpecMAT/nbCrystInSegmentColumnCmd
ActarSimSpecMATSciDetectorMessenger::
ActarSimSpecMATSciDetectorMessenger(ActarSimSpecMATSciDetectorConstruction* ActarSimSpecMATSciDet)
  :ActarSimSpecMATSciDetector(ActarSimSpecMATSciDet) {

  detDir = new G4UIdirectory("/ActarSim/det/SpecMAT/");
  detDir->SetGuidance("SpecMAT detector control");

  printCmd = new G4UIcmdWithoutParameter("/ActarSim/det/SpecMAT/print",this);
  printCmd->SetGuidance("Prints geometry.");
  printCmd->AvailableForStates(G4State_Idle);

  sciSpecMATMaterialCmd = new G4UIcmdWithAString("/ActarSim/det/SpecMAT/setSciMat",this);
  sciSpecMATMaterialCmd->SetGuidance("Select Material the scintillators are made of.");
  sciSpecMATMaterialCmd->SetParameterName("sciMat",false);
  sciSpecMATMaterialCmd->SetDefaultValue("CeBr3");
  sciSpecMATMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  reflectorSpecMATMaterialCmd = new G4UIcmdWithAString("/ActarSim/det/SpecMAT/setReflectorMat",this);
  reflectorSpecMATMaterialCmd->SetGuidance("Select Material the reflectors are made of.");
  reflectorSpecMATMaterialCmd->SetParameterName("reflectorMat",false);
  reflectorSpecMATMaterialCmd->SetDefaultValue("TiO2");
  reflectorSpecMATMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  windowSpecMATMaterialCmd = new G4UIcmdWithAString("/ActarSim/det/SpecMAT/setWindowMat",this);
  windowSpecMATMaterialCmd->SetGuidance("Select Material the windows are made of.");
  windowSpecMATMaterialCmd->SetParameterName("windowMat",false);
  windowSpecMATMaterialCmd->SetDefaultValue("Quartz");
  windowSpecMATMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  housingSpecMATMaterialCmd = new G4UIcmdWithAString("/ActarSim/det/SpecMAT/setHousingMat",this);
  housingSpecMATMaterialCmd->SetGuidance("Select Material the housing is made of.");
  housingSpecMATMaterialCmd->SetParameterName("housingMat",false);
  housingSpecMATMaterialCmd->SetDefaultValue("AluminiumMat");
  housingSpecMATMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  vacuumFlangeSpecMATMaterialCmd = new G4UIcmdWithAString("/ActarSim/det/SpecMAT/setVacuumFlangeMat",this);
  vacuumFlangeSpecMATMaterialCmd->SetGuidance("Select Material the vacuum flange is made of.");
  vacuumFlangeSpecMATMaterialCmd->SetParameterName("vacuumFlangeMat",false);
  vacuumFlangeSpecMATMaterialCmd->SetDefaultValue("AluminiumMat");
  vacuumFlangeSpecMATMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  insulationTubeSpecMATMaterialCmd = new G4UIcmdWithAString("/ActarSim/det/SpecMAT/setInsulationTubeMat",this);
  insulationTubeSpecMATMaterialCmd->SetGuidance("Select Material the insulation tube is made of.");
  insulationTubeSpecMATMaterialCmd->SetParameterName("insulationTubeMat",false);
  insulationTubeSpecMATMaterialCmd->SetDefaultValue("Ceramic_Al2O3");
  insulationTubeSpecMATMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  vacuumChamberIncludedFlagCmd = new G4UIcmdWithAString("/ActarSim/det/SpecMAT/vacuumChamberIncludedFlag",this);
  vacuumChamberIncludedFlagCmd->SetGuidance("Includes the vacuum chamber and associated elements (default on).");
  vacuumChamberIncludedFlagCmd->SetGuidance("  Choice : on(default), off");
  vacuumChamberIncludedFlagCmd->SetParameterName("choice",true);
  vacuumChamberIncludedFlagCmd->SetDefaultValue("on");
  vacuumChamberIncludedFlagCmd->SetCandidates("on off");
  vacuumChamberIncludedFlagCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  insulationTubeIncludedFlagCmd = new G4UIcmdWithAString("/ActarSim/det/SpecMAT/insulationTubeIncludedFlag",this);
  insulationTubeIncludedFlagCmd->SetGuidance("Includes the insulation tube and associated elements (default on).");
  insulationTubeIncludedFlagCmd->SetGuidance("  Choice : on(default), off");
  insulationTubeIncludedFlagCmd->SetParameterName("choice",true);
  insulationTubeIncludedFlagCmd->SetDefaultValue("on");
  insulationTubeIncludedFlagCmd->SetCandidates("on off");
  insulationTubeIncludedFlagCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nbSegmentsCmd = new G4UIcmdWithAnInteger("/ActarSim/det/SpecMAT/nbSegmentsCmd",this);
  nbSegmentsCmd->SetGuidance("Sets the number of scintillator segments");
  nbSegmentsCmd->SetParameterName("nbSegments",false);
  nbSegmentsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nbCrystInSegmentRowCmd = new G4UIcmdWithAnInteger("/ActarSim/det/SpecMAT/nbCrystInSegmentRowCmd",this);
  nbCrystInSegmentRowCmd->SetGuidance("Sets the number of scintillator crystals in a segment row");
  nbCrystInSegmentRowCmd->SetParameterName("nbCrystInSegmentRow",false);
  nbCrystInSegmentRowCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nbCrystInSegmentColumnCmd = new G4UIcmdWithAnInteger("/ActarSim/det/SpecMAT/nbCrystInSegmentColumnCmd",this);
  nbCrystInSegmentColumnCmd->SetGuidance("Sets the number of scintillator crystals in a segment column");
  nbCrystInSegmentColumnCmd->SetParameterName("nbCrystInSegmentColumn",false);
  nbCrystInSegmentColumnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//////////////////////////////////////////////////////////////////
/// Destructor
ActarSimSpecMATSciDetectorMessenger::~ActarSimSpecMATSciDetectorMessenger() {
  delete detDir;
  delete printCmd;
  delete sciSpecMATMaterialCmd;
  delete reflectorSpecMATMaterialCmd;
  delete windowSpecMATMaterialCmd;
  delete housingSpecMATMaterialCmd;
  delete vacuumFlangeSpecMATMaterialCmd;
  delete insulationTubeSpecMATMaterialCmd;
  delete vacuumChamberIncludedFlagCmd;
  delete insulationTubeIncludedFlagCmd;
  delete nbSegmentsCmd;
  delete nbCrystInSegmentRowCmd;
  delete nbCrystInSegmentColumnCmd;
}

//////////////////////////////////////////////////////////////////
/// Setting the new values and connecting to actions
void ActarSimSpecMATSciDetectorMessenger::SetNewValue(G4UIcommand* command,
						      G4String newValue) {
  if( command == printCmd )
    ActarSimSpecMATSciDetector->PrintDetectorParameters();

  if(command == sciSpecMATMaterialCmd)
    ActarSimSpecMATSciDetector->SetSciSpecMATMaterial(newValue);

  if(command == reflectorSpecMATMaterialCmd)
    ActarSimSpecMATSciDetector->SetReflectorSpecMATMaterial(newValue);

  if(command == windowSpecMATMaterialCmd)
    ActarSimSpecMATSciDetector->SetWindowSpecMATMaterial(newValue);

  if(command == housingSpecMATMaterialCmd)
    ActarSimSpecMATSciDetector->SetHousingSpecMATMaterial(newValue);

  if(command == vacuumFlangeSpecMATMaterialCmd)
    ActarSimSpecMATSciDetector->SetVacuumFlangeSpecMATMaterial(newValue);

  if(command == insulationTubeSpecMATMaterialCmd)
    ActarSimSpecMATSciDetector->SetInsulationTubeSpecMATMaterial(newValue);

  if(command == vacuumChamberIncludedFlagCmd)
    ActarSimSpecMATSciDetector->SetVacuumChamberIncludedFlag(newValue);

  if(command == insulationTubeIncludedFlagCmd)
    ActarSimSpecMATSciDetector->SetInsulationTubeIncludedFlag(newValue);

  if( command == nbSegmentsCmd)
    ActarSimSpecMATSciDetector->SetNbSegments(nbSegmentsCmd->GetNewIntValue(newValue));

  if( command == nbCrystInSegmentRowCmd)
    ActarSimSpecMATSciDetector->SetNbCrystInSegmentRow(nbCrystInSegmentRowCmd->GetNewIntValue(newValue));

  if( command == nbCrystInSegmentColumnCmd)
    ActarSimSpecMATSciDetector->SetNbCrystInSegmentColumn(nbCrystInSegmentColumnCmd->GetNewIntValue(newValue));

}
