// - AUTHOR: Hector Alvarez-Pol 03/2017
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/
//////////////////////////////////////////////////////////////////
/// \class ActarSimSpecMATSciDetectorConstruction
/// Scintillator detector description for SpecMAT
/////////////////////////////////////////////////////////////////

#include "ActarSimSpecMATSciDetectorConstruction.hh"
#include "ActarSimDetectorConstruction.hh"
#include "ActarSimGasDetectorConstruction.hh"
#include "ActarSimSpecMATSciDetectorMessenger.hh"
#include "ActarSimROOTAnalysis.hh"
#include "ActarSimSpecMATSciSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4Transform3D.hh"
#include "G4SubtractionSolid.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "globals.hh"

//////////////////////////////////////////////////////////////////
/// Constructor. Sets the material and the pointer to the Messenger
ActarSimSpecMATSciDetectorConstruction::
ActarSimSpecMATSciDetectorConstruction(ActarSimDetectorConstruction* det)
  :	sciBulkMaterial(0),detConstruction(det) {
    SetSpecMATSciMaterial("CeBr3");
    SetReflectorSpecMATMaterial("TiO2");
    SetWindowSpecMATMaterial("Quartz");
    SetHousingSpecMATMaterial("Al_Alloy");

//// FROM THE ORIGINAL SPECMAT CODE
  nbSegments = 15;
  nbCrystInSegmentRow = 3;        //# of rings
  nbCrystInSegmentColumn = 1;     //# of crystals

  vacuumChamber = "yes"; //"yes"/"no"
  vacuumFlangeSizeX = 150*mm;
  vacuumFlangeSizeY = 19*mm;
  vacuumFlangeSizeZ = 10*mm;
  vacuumFlangeThickFrontOfScint = 2*mm;

  insulationTube = "yes"; //"yes"/"no"
  insulationTubeThickness = 10*mm;

  dPhi = twopi/nbSegments;
  half_dPhi = 0.5*dPhi;
  tandPhi = std::tan(half_dPhi);

  sciCrystSizeX = 24.*mm;								//Size and position of all components depends on Crystal size and position.
  sciCrystSizeY = 24.*mm;
  sciCrystSizeZ = 24.*mm;

  // Position of the crystal
  sciCrystPosX = 0;									//Position of the Crystal along the X axis
  sciCrystPosY = 0;									//Position of the Crystal along the Y axis
  sciCrystPosZ = 0; 			 						//Position of the Crystal along the Z axis

  // create commands for interactive definition of the calorimeter
  sciMessenger = new ActarSimSciDetectorMessenger(this);
}

//////////////////////////////////////////////////////////////////
/// Destructor.
ActarSimSciDetectorConstruction::~ActarSimSciDetectorConstruction(){
  delete sciMessenger;
}

//////////////////////////////////////////////////////////////////
///  Wrap for the construction function
G4VPhysicalVolume* ActarSimSciDetectorConstruction::Construct(G4LogicalVolume* worldLog) {
  return ConstructSci(worldLog);
}

//////////////////////////////////////////////////////////////////
/// Constructs the scintillator detector elements
///
/// -  Details should follow!!!!
G4VPhysicalVolume* ActarSimSciDetectorConstruction::ConstructSci(G4LogicalVolume* worldLog) {
  // Position of the elements
    G4ThreeVector sciCrystPos = G4ThreeVector(sciCrystPosX,sciCrystPosY,sciCrystPosZ);
    G4ThreeVector sciWindPos =

  // Define box for Crystal
  G4VSolid* sciCrystSolid = new G4Box("sciCrystSolid",sciCrystSizeX,sciCrystSizeY,sciCrystSizeZ);
  G4LogicalVolume* sciCrystLog = new G4LogicalVolume(sciCrystSolid,sciCrystMat,"crystal");

  // Thickness of reflector walls
  sciReflWallThickX = 0.5*mm;
  sciReflWallThickY = 0.5*mm;
  sciReflWindThick = 0.5*mm;

  // Outer dimensions of the reflector relative to the crystal size
  sciReflSizeX = sciCrystSizeX + sciReflWallThickX;
  sciReflSizeY = sciCrystSizeY + sciReflWallThickY;
  sciReflSizeZ = sciCrystSizeZ + sciReflWindThick/2;

  // Position of the reflector relative to the crystal position
  sciReflPosX = sciCrystPosX;
  sciReflPosY = sciCrystPosY;
  sciReflPosZ = sciCrystPosZ - sciReflWindThick/2;					//Position of the Reflector relative to the Al Housing along the Z axis

  G4ThreeVector sciReflPos = G4ThreeVector(sciReflPosX,sciReflPosY,sciReflPosZ);

  // Define box for Reflector
  G4VSolid* reflBoxSolid = new G4Box("reflBoxSolid",sciReflSizeX,sciReflSizeY,sciReflSizeZ);

  // Subtracts Crystal box from Reflector box
  G4VSolid* sciReflSolid = new G4SubtractionSolid("sciReflSolid",reflBoxSolid,sciCrystSolid,0,
                                                  G4ThreeVector(sciCrystPosX, sciCrystPosY, sciReflWindThick/2));

  // Define Logical Volume for Reflector//
  G4LogicalVolume* sciReflLog = new G4LogicalVolume(sciReflSolid, sciReflMat, "sciReflLog");





  //------------------------------------------------
  // Sensitive detectors
  //------------------------------------------------
  sciLog->SetSensitiveDetector( detConstruction->GetSciSD() );

  //------------------------------------------------------------------
  // Visualization attributes
  //------------------------------------------------------------------
  G4VisAttributes* sciCrystVisAtt =
	  new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));	//Instantiation of visualization attributes with blue colour
  sciCrystVisAtt->SetVisibility(true);						//Pass this object to Visualization Manager for visualization
  sciCrystVisAtt->SetForceWireframe(true);			  //I still believe that it might make Crystal transparent
  sciCrystLog->SetVisAttributes(sciCrystVisAtt); //Assignment of visualization attributes to the logical volume of the Crystal

  G4VisAttributes* sciReflVisAtt =
    new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));	//Instantiation of visualization attributes with yellow colour
  sciReflVisAtt->SetVisibility(true);							//Pass this object to Visualization Manager for visualization
  sciReflLog->SetVisAttributes(sciReflVisAtt);		//Assignment of visualization attributes to the logical volume of the Reflector

  // Printing the final settings...
  G4cout << "##################################################################"
   << G4endl
   << "#########  ActarSimSpecMATSciDetectorConstruction::ConstructSci()  #######"
   << G4endl
  G4cout << "##################################################################"
   << G4endl;

  return sciPhys;
}

//////////////////////////////////////////////////////////////////
/// Set the material the scintillator bulk is made of
void ActarSimSciDetectorConstruction::SetSciBulkMaterial (G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) sciBulkMaterial = pttoMaterial;
}

//////////////////////////////////////////////////////////////////
/// Updates Scintillator detector
void ActarSimSciDetectorConstruction::UpdateGeometry() {
  Construct(detConstruction->GetWorldLogicalVolume());
  G4RunManager::GetRunManager()->
    DefineWorldVolume(detConstruction->GetWorldPhysicalVolume());
}

//////////////////////////////////////////////////////////////////
/// Prints Scintillator detector parameters. TODO: To be filled
void ActarSimSciDetectorConstruction::PrintDetectorParameters() {
  G4cout << "##################################################################"
	 << G4endl
	 << "####  ActarSimSciDetectorConstruction::PrintDetectorParameters() ####"
	 << G4endl;
  G4cout << "##################################################################"
	 << G4endl;
}
