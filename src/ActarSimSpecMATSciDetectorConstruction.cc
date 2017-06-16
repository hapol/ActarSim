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
#include "ActarSimSciSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polyhedra.hh"
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
  : detConstruction(det) {
  SetSciSpecMATMaterial("CeBr3");
  SetReflectorSpecMATMaterial("TiO2");
  SetWindowSpecMATMaterial("Quartz");
  SetHousingSpecMATMaterial("AluminiumMat");
  SetVacuumFlangeSpecMATMaterial("AluminiumMat");
  SetInsulationTubeSpecMATMaterial("Ceramic_Al2O3");
  
  // FROM THE ORIGINAL SPECMAT CODE
  nbSegments = 15;            //the number of scintillator segments
  nbCrystInSegmentRow = 3;    //the number of scintillator crystals in a segment row
  nbCrystInSegmentColumn = 1; //the number of scintillator crystals in a segment column
  
  SetVacuumChamberIncludedFlag("on"); //default value "on" in the SpecMAT constructor
  vacuumFlangeSizeX = 150*mm;
  vacuumFlangeSizeY = 19*mm;
  vacuumFlangeSizeZ = 10*mm;
  vacuumFlangeThickFrontOfScint = 2*mm;
  
  SetInsulationTubeIncludedFlag("on"); //default value "on" in the SpecMAT constructor
  insulationTubeThickness = 10*mm;
  
  //****************************************************************************//
  //**************** CeBr3 cubic scintillator 1.5"x1.5"x1.5" *******************//
  //****************************************************************************//
  
  //--------------------------------------------------------//
  //***************** Scintillation crystal ****************//
  //--------------------------------------------------------//
  // Dimensions of the crystal
  sciCrystSizeX = 24.*mm;  //Size and position of all components depends on Crystal size and position.
  sciCrystSizeY = 24.*mm;
  sciCrystSizeZ = 24.*mm;
  
  // Position of the crystal
  sciCrystPosX = 0;       //Position of the Crystal along the X axis
  sciCrystPosY = 0;       //Position of the Crystal along the Y axis
  sciCrystPosZ = 0;       //Position of the Crystal along the Z axis
  
  // Thickness of reflector walls
  sciReflWallThickX = 0.5*mm;
  sciReflWallThickY = 0.5*mm;
  sciReflWindThick = 0.5*mm;
  
  //--------------------------------------------------------//
  //******************** Aluminum Housing ******************//
  //--------------------------------------------------------//
  // Dimensions of Housing (half-side)
  sciHousWallThickX = 3.0*mm;
  sciHousWallThickY = 3.0*mm;
  sciHousWindThick = 1.0*mm;

  // Outer dimensions of the housing relative to the crystal size and to the thickness of the reflector
  sciHousSizeX = sciCrystSizeX + sciReflWallThickX + sciHousWallThickX;
  sciHousSizeY = sciCrystSizeY + sciReflWallThickY + sciHousWallThickY;
  sciHousSizeZ = sciCrystSizeZ + sciReflWindThick/2 + sciHousWindThick/2;
 
  // create commands for interactive definition of the calorimeter
  sciMessenger = new ActarSimSpecMATSciDetectorMessenger(this);
}

//////////////////////////////////////////////////////////////////
/// Destructor.
ActarSimSpecMATSciDetectorConstruction::~ActarSimSpecMATSciDetectorConstruction(){
  delete sciMessenger;
}

//////////////////////////////////////////////////////////////////
///  Wrap for the construction function
G4VPhysicalVolume* ActarSimSpecMATSciDetectorConstruction::Construct(G4LogicalVolume* worldLog) {
  return ConstructSci(worldLog);
}

//////////////////////////////////////////////////////////////////
/// Constructs the scintillator detector elements
///
/// -  Details should follow!!!!
G4VPhysicalVolume* ActarSimSpecMATSciDetectorConstruction::ConstructSci(G4LogicalVolume* worldLog) {
  // Position of the elements
  G4ThreeVector sciCrystPos = G4ThreeVector(sciCrystPosX,sciCrystPosY,sciCrystPosZ);
  
  // Define box for Crystals
  G4VSolid* sciCrystSolid = new G4Box("sciCrystSolid",sciCrystSizeX,sciCrystSizeY,sciCrystSizeZ);
  G4LogicalVolume* sciCrystLog = new G4LogicalVolume(sciCrystSolid,sciSpecMATMaterial,"crystal");
  G4VPhysicalVolume* sciCrysPhys;

  // Outer dimensions of the reflector relative to the crystal size
  G4double sciReflSizeX = sciCrystSizeX + sciReflWallThickX;
  G4double sciReflSizeY = sciCrystSizeY + sciReflWallThickY;
  G4double sciReflSizeZ = sciCrystSizeZ + sciReflWindThick/2;
  
  // Position of the reflector relative to the crystal position
  G4double sciReflPosX = sciCrystPosX;
  G4double sciReflPosY = sciCrystPosY;
  G4double sciReflPosZ = sciCrystPosZ - sciReflWindThick/2; //Position of the Reflector relative to the Al Housing along the Z axis
  
  G4ThreeVector sciReflPos = G4ThreeVector(sciReflPosX,sciReflPosY,sciReflPosZ);
  
  // Define box for Reflector
  G4VSolid* reflBoxSolid = new G4Box("reflBoxSolid",sciReflSizeX,sciReflSizeY,sciReflSizeZ);
  
  // Subtracts Crystal box from Reflector box 
  G4VSolid* sciReflSolid = new G4SubtractionSolid("sciReflSolid",reflBoxSolid,sciCrystSolid,0,
                                                  G4ThreeVector(sciCrystPosX, sciCrystPosY, sciReflWindThick/2));
  
  // Define Logical Volume for Reflector
  G4LogicalVolume* sciReflLog = new G4LogicalVolume(sciReflSolid, reflectorSpecMATMaterial, "sciReflLog");
  
  //--------------------------------------------------------//
  //******************** Aluminum Housing ******************//
  //--------------------------------------------------------//
  // Outer dimensions of the housing relative to the crystal size and to the thickness of the reflector
  sciHousSizeX = sciCrystSizeX + sciReflWallThickX + sciHousWallThickX;
  sciHousSizeY = sciCrystSizeY + sciReflWallThickY + sciHousWallThickY;
  sciHousSizeZ = sciCrystSizeZ + sciReflWindThick/2 + sciHousWindThick/2;
  
  // Position of the housing relative to the crystal position
  G4double sciHousPosX = sciCrystPosX;
  G4double sciHousPosY = sciCrystPosY;
  G4double sciHousPosZ = sciCrystPosZ - (sciReflWindThick/2 + sciHousWindThick/2);
  
  G4ThreeVector sciHousPos = G4ThreeVector(sciHousPosX, sciHousPosY, sciHousPosZ);
  
  // Define box for Housing
  G4VSolid* housBoxASolid = new G4Box("housBoxASolid", sciHousSizeX, sciHousSizeY, sciHousSizeZ);
  
  // Subtracts Reflector box from Housing box
  G4VSolid* sciHousSolid = new G4SubtractionSolid("housBoxBSolid", housBoxASolid, reflBoxSolid,
						  0, G4ThreeVector(sciReflPosX, sciReflPosY, sciHousWindThick/2));
  
  // Define Logical Volume for Housing
  G4LogicalVolume* sciHousLog =
    new G4LogicalVolume(sciHousSolid, 	    	//Housing solid shape
			housingSpecMATMaterial,           	//Housing material
			"sciCaseLog");        	//Housing logic volume name
  
  //--------------------------------------------------------//
  //******************** Quartz window *********************//
  //--------------------------------------------------------//
  // Dimensions of the Window (half-side)
  G4double sciWindSizeX = sciCrystSizeX + sciReflWallThickX + sciHousWallThickX;	//X half-size of the Window
  G4double sciWindSizeY = sciCrystSizeY + sciReflWallThickY + sciHousWallThickY; //Y half-size of the Window
  G4double sciWindSizeZ = 1.*mm;						        //Z half-size of the Window
  
  // Position of the window relative to the crystal
  G4double sciWindPosX = sciCrystPosX ;				      //Position of the Window along the X axis
  G4double sciWindPosY = sciCrystPosY ;				      //Position of the Window along the Y axis
  G4double sciWindPosZ = sciCrystPosZ + sciCrystSizeZ + sciWindSizeZ;  //Position of the Window relative to the Al Housing along the Z axis

  //Position of the Window in space relative to the Al Housing
  G4ThreeVector sciWindPos = G4ThreeVector(sciWindPosX, 
					   sciWindPosY,
					   sciWindPosZ);
  
  // Define solid for the Window
  G4VSolid* sciWindSolid = 	        	//Define object for the Window's box
    new G4Box("sciWindSolid",		        //Name of the Window's box
	      sciWindSizeX, 	        	//X half_size of the box
	      sciWindSizeY, 	        	//Y half_size of the box
	      sciWindSizeZ);	        	//Z half_size of the box
  
  // Define Logical Volume for Window
  G4LogicalVolume* sciWindLog = new G4LogicalVolume(sciWindSolid, windowSpecMATMaterial, "sciWindLog");
  
  G4double circleRadius = ComputeCircleR1();
  
  // Define segment which will contain crystals
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* segment_mat = nist->FindOrBuildMaterial("G4_Galactic", false);
  G4VSolid* segmentBox = new G4Box("segmentBox",
				   sciHousSizeX*nbCrystInSegmentRow,
				   sciHousSizeY*nbCrystInSegmentColumn,
				   sciHousSizeZ+sciWindSizeZ);
  
  // Checking if the flange dimensions are not smaller than the segment dimensions
  if (vacuumFlangeSizeY<sciHousSizeY*nbCrystInSegmentColumn) vacuumFlangeSizeY=sciHousSizeY*nbCrystInSegmentColumn;
  
  //Define the vacuum chamber flange
  G4LogicalVolume* vacuumFlangeBoxLog;
  if (vacuumChamberIncludedFlag == "on") {
    G4VSolid* vacuumFlangeBox = 
      new G4Box("vacuumFlangeBox",
		vacuumFlangeSizeX,
		vacuumFlangeSizeY,
		vacuumFlangeSizeZ);
    // Subtracts Reflector box from Housing box
    G4VSolid* vacuumFlangeSolid =
      new G4SubtractionSolid("vacuumFlangeSolid",
			     vacuumFlangeBox,
			     segmentBox,
			     0,
			     G4ThreeVector(0, 0, (sciHousSizeZ+sciWindSizeZ)+vacuumFlangeSizeZ-(2*vacuumFlangeSizeZ-vacuumFlangeThickFrontOfScint)));
    
    vacuumFlangeBoxLog = 
      new G4LogicalVolume(vacuumFlangeSolid,
			  vacuumFlangeSpecMATMaterial,
			  "vacuumFlangeBoxLog");
    
    //modify if needed... for the moment same material
    G4Material* vacuumSideFlangeMat = vacuumFlangeSpecMATMaterial;
 
    G4RotationMatrix rotSideFlnge  = G4RotationMatrix();
    rotSideFlnge.rotateZ(twopi/nbSegments/2);
    G4ThreeVector positionSideFlange1 = G4ThreeVector(0, 0, vacuumFlangeSizeX);
    G4Transform3D transformSideFlange1 = G4Transform3D(rotSideFlnge, positionSideFlange1);
    G4ThreeVector positionSideFlange2 = G4ThreeVector(0, 0, -vacuumFlangeSizeX-2*vacuumFlangeSizeZ);
    G4Transform3D transformSideFlange2 = G4Transform3D(rotSideFlnge, positionSideFlange2);
    G4double vacuumChamberSideFlangeThickness[] = {0, 2*vacuumFlangeSizeZ};
    G4double vacuumChamberSideFlangeInnerR[] = {0, 0};
    G4double vacuumChamberSideFlangeOuterR[] = {circleRadius+2*vacuumFlangeSizeZ, circleRadius+2*vacuumFlangeSizeZ};
    
    G4VSolid* vacuumChamberSideFlange = 
      new G4Polyhedra("vacuumChamberSideFlange",
		      0,
		      2*3.1415926535897932384626433,
		      nbSegments,
		      2,
		      vacuumChamberSideFlangeThickness,
		      vacuumChamberSideFlangeInnerR,
		      vacuumChamberSideFlangeOuterR);
    
    G4LogicalVolume* vacuumChamberSideFlangeLog = 
      new G4LogicalVolume(vacuumChamberSideFlange,
			  vacuumSideFlangeMat,
			  "vacuumChamberSideFlangeLog");
    
    G4VPhysicalVolume* vacuumChamberSideFlange1Phys = 
      new G4PVPlacement(transformSideFlange1,
			vacuumChamberSideFlangeLog,                //its logical volume
			"VacuumChamberSideFlangeLog",              //its name
			worldLog,                                //its mother  volume
			false,                                     //no boolean operation
			1);                                         //copy number
    //fCheckOverlaps);                           // checking overlaps
    
    G4VPhysicalVolume* vacuumChamberSideFlange2Phys = 
      new G4PVPlacement(transformSideFlange2,
			vacuumChamberSideFlangeLog,                //its logical volume
			"VacuumChamberSideFlangeLog",              //its name
			worldLog,                                //its mother  volume
			false,                                     //no boolean operation
			2);                                         //copy number
    //fCheckOverlaps);                           // checking overlaps
    if(vacuumChamberSideFlange1Phys)  
      ;
    if(vacuumChamberSideFlange2Phys)  
      ;
  }
  
  
  // Defines insulation tube between the field cage and the vacuum chamber which 
  // might be used for preventing sparks in the real setup
  // And its stopping power should be simulated
  if (vacuumChamberIncludedFlag == "on" && insulationTubeIncludedFlag == "on") {
    //Geometry of the insulation Tube
    G4double insulationTubeInnerRadius = circleRadius-insulationTubeThickness;
    G4double insulationTubeOuterRadius = circleRadius;

  G4cout << "insulationTubeInnerRadius: " <<  circleRadius 
	 << " - " << insulationTubeThickness<< G4endl;

    G4VSolid* insulationTubeSolid = 
      new G4Tubs("insulationTubeSolid",
		 insulationTubeInnerRadius,
		 insulationTubeOuterRadius,
		 vacuumFlangeSizeX,
		 0*deg,
		 360*deg);
    
    G4LogicalVolume* insulationTubeLog = 
      new G4LogicalVolume(insulationTubeSolid,
			  insulationTubeSpecMATMaterial,
			  "insulationTubeLog");
    
    G4VPhysicalVolume* insulationTubePhys = 
      new G4PVPlacement(0,
			G4ThreeVector(0,0,0),
			insulationTubeLog,              //its logical volume
			"insulationTubePhys",           //its name
			worldLog,                     //its mother  volume
			false,                          //no boolean operation
			1);                              //copy number
      //fCheckOverlaps);                // checking overlaps
    
    if(insulationTubePhys)
      ;

    // Visualization attributes for the insulation tube
    G4VisAttributes* insulationTubeVisAtt =
      new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));	//Instantiation of visualization attributes with cyan colour
    insulationTubeVisAtt->SetVisibility(true);	      	//Pass this object to Visualization Manager for visualization
    insulationTubeVisAtt->SetForceWireframe(true);    	//I believe that it might make Window transparent
    insulationTubeLog->SetVisAttributes(insulationTubeVisAtt);//Assignment of visualization attributes to the logical volume of the Window
  }
  
  //Positioning of segments and crystals in the segment
  
  // In TotalCrystNb array will be stored coordinates of the all crystals, 
  // which could be used for further Doppler correction
  G4int TotalCrystNb = nbCrystInSegmentRow*nbCrystInSegmentColumn*nbSegments;  //Dimension of the dynamic the array
  G4ThreeVector *crystalPositionsArray = new G4ThreeVector[TotalCrystNb];      //Dinamic mamory allocation for the array
  for (int i=0; i<TotalCrystNb; i++) {
    crystalPositionsArray[i] = G4ThreeVector(0.,0.,0.); // Initialize all elements of the array to zero.
  }
  
  G4int i = 0;          //counter for reconstruction of crystal positions
  G4int crysNb = 1;     //crystal counter
  for (G4int iseg = 0; iseg < nbSegments ; iseg++) {
    G4double phi = iseg*twopi/nbSegments;
    G4RotationMatrix rotm  = G4RotationMatrix();     //** rotation matrix for positioning segments
    rotm.rotateY(90*deg);                            //** rotation matrix for positioning segments
    rotm.rotateZ(phi);                               //** rotation matrix for positioning segments
    
    G4RotationMatrix rotm2  = G4RotationMatrix();    //### rotation matrix for reconstruction of crystal positions
    rotm2.rotateX(360*deg - phi);                    //### rotation matrix for reconstruction of crystal positions
    G4RotationMatrix rotm3  = G4RotationMatrix();    //### rotation matrix for reconstruction of crystal positions
    rotm3.rotateY(90*deg);                           //### rotation matrix for reconstruction of crystal positions
    
    G4ThreeVector uz = G4ThreeVector(std::cos(phi), std::sin(phi), 0.); //coeficient which will be used for preliminary rotation of the segments and crystals
    
    G4LogicalVolume* segmentBoxLog = 
      new G4LogicalVolume(segmentBox,
			  segment_mat,
			  "segmentBoxLog");

    G4ThreeVector positionInSegment = G4ThreeVector(-(nbCrystInSegmentRow*sciHousSizeX-sciHousSizeX), 
						    -(nbCrystInSegmentColumn*sciHousSizeY-sciHousSizeY), 
						    (sciHousSizeZ-sciCrystSizeZ-sciWindSizeZ));
    
    for (G4int icrystRow = 0; icrystRow < nbCrystInSegmentColumn; icrystRow++) {
      for (G4int icrystCol = 0; icrystCol < nbCrystInSegmentRow; icrystCol++) {
	G4RotationMatrix rotm1  = G4RotationMatrix();
	
	G4ThreeVector positionCryst = (G4ThreeVector(0., 0., sciCrystPosZ) + positionInSegment);
	G4ThreeVector positionWind = (G4ThreeVector(0., 0., sciWindPosZ) + positionInSegment);
	G4ThreeVector positionRefl = (G4ThreeVector(0., 0., sciReflPosZ) + positionInSegment);
	G4ThreeVector positionHous = (G4ThreeVector(0., 0., sciHousPosZ) + positionInSegment);
	
	crystalPositionsArray[crysNb - 1] = positionCryst; //assigning initial crystal positions in a segment into array
	
	G4Transform3D transformCryst = G4Transform3D(rotm1,positionCryst);
	G4Transform3D transformWind = G4Transform3D(rotm1,positionWind);
	G4Transform3D transformRefl = G4Transform3D(rotm1,positionRefl);
	G4Transform3D transformHous = G4Transform3D(rotm1,positionHous);
	
	// Crystal position
	sciCrysPhys=
	  new G4PVPlacement(transformCryst,			  //no rotation here rotm1 is empty, position
			    sciCrystLog,                //its logical volume
			    "sciCrystPl",               //its name
			    segmentBoxLog,              //its mother  volume
			    false,                      //no boolean operation
			    crysNb);                     //crystal unique number will
	//fCheckOverlaps);            // checking overlaps
	
	new G4PVPlacement(transformWind,
			  sciWindLog,
			  "sciWindPl",
			  segmentBoxLog,
			  false,
			  crysNb);
	//fCheckOverlaps);
	
	new G4PVPlacement(transformRefl,
			  sciReflLog,
			  "sciReflPl",
			  segmentBoxLog,
			  false,
			  crysNb);
	//fCheckOverlaps);
	
	new G4PVPlacement(transformHous,
			  sciHousLog,
			  "sciHousPl",
			  segmentBoxLog,
			  false,
			  crysNb);
	//fCheckOverlaps);
	crysNb += 1;
	positionInSegment += G4ThreeVector(sciHousSizeX*2, 0., 0.);
      }
      positionInSegment -= G4ThreeVector(nbCrystInSegmentRow*sciHousSizeX*2, 0., 0.);
      positionInSegment += G4ThreeVector(0., sciHousSizeY*2, 0.);
    }

    //segment and flange positioning
    
    if (vacuumChamberIncludedFlag == "on") {
      //Flange positioning
      G4ThreeVector positionVacuumFlange = (circleRadius+vacuumFlangeSizeZ)*uz;
      G4Transform3D transformVacuumFlange = G4Transform3D(rotm, positionVacuumFlange);
      new G4PVPlacement(transformVacuumFlange, //position
			vacuumFlangeBoxLog,                //its logical volume
			"VacuumFlange",                    //its name
			worldLog,                         //its mother  volume
			false,                              //no boolean operation
			iseg);                               //copy number
      //fCheckOverlaps);                    // checking overlaps
      //Segment positioning
      G4ThreeVector positionSegment = (circleRadius+2*vacuumFlangeSizeZ+(sciHousSizeZ+sciWindSizeZ)-(2*vacuumFlangeSizeZ-vacuumFlangeThickFrontOfScint))*uz;
      G4Transform3D transformSegment = G4Transform3D(rotm, positionSegment);
      new G4PVPlacement(transformSegment, //position
			segmentBoxLog,                //its logical volume
			"Segment",                    //its name
			worldLog,                   //its mother  volume
			false,                        //no boolean operation
			iseg);                         //copy number
			//fCheckOverlaps);              // checking overlaps
      //Saving crystal positions in the crystalPositionsArray array
      for (; i < crysNb-1; i++) {
	G4AffineTransform TransformCrystPos1;
	TransformCrystPos1.SetNetRotation(rotm2); //rotates the crystal centers (in one segment) by angle phi around X
	crystalPositionsArray[i] = TransformCrystPos1.TransformPoint(crystalPositionsArray[i]);
	
	G4AffineTransform TransformCrystPos;
	TransformCrystPos.SetNetRotation(rotm3); //rotates the crystal centers (in one segment) by 90deg around Y
	TransformCrystPos.SetNetTranslation(positionSegment);
	crystalPositionsArray[i] = TransformCrystPos.TransformPoint(crystalPositionsArray[i]);
      }
    }
    //segment position in case vacuumChamberIncludedFlag is "off"
    else {
      //Segment positioning
      G4ThreeVector positionSegment = (circleRadius+(sciHousSizeZ+sciWindSizeZ))*uz;
      G4Transform3D transformSegment = G4Transform3D(rotm, positionSegment);
      new G4PVPlacement(transformSegment, //position
			segmentBoxLog,                //its logical volume
			"Segment",                    //its name
			worldLog,                   //its mother  volume
			false,                        //no boolean operation
			iseg);                         //copy number
			//fCheckOverlaps);              // checking overlaps
      //Saving crystal positions in the crystalPositionsArray array
      for (; i < crysNb-1; i++) {
	G4AffineTransform TransformCrystPos1;
	TransformCrystPos1.SetNetRotation(rotm2); //rotates the crystal centers (in one segment) by angle phi around X
	crystalPositionsArray[i] = TransformCrystPos1.TransformPoint(crystalPositionsArray[i]);
	
	G4AffineTransform TransformCrystPos;
	TransformCrystPos.SetNetRotation(rotm3); //rotates the crystal centers (in one segment) by 90deg around Y
	TransformCrystPos.SetNetTranslation(positionSegment);
	crystalPositionsArray[i] = TransformCrystPos.TransformPoint(crystalPositionsArray[i]);
      }
    }
  }
  
  //------------------------------------------------
  // Sensitive detectors
  //------------------------------------------------
  sciCrystLog->SetSensitiveDetector( detConstruction->GetSciSD() );
  
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
  
  G4VisAttributes* sciHousVisAtt =
    new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); //Instantiation of visualization attributes with grey colour
  sciHousVisAtt->SetVisibility(true);						  //Pass this object to Visualization Manager for visualization
  sciHousLog->SetVisAttributes(sciHousVisAtt);	  //Assignment of visualization attributes to the logical volume of the Housing
  
  G4VisAttributes* sciWindVisAtt =
    new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));					//Instantiation of visualization attributes with cyan colour
  sciWindVisAtt->SetVisibility(true);							//Pass this object to Visualization Manager for visualization
  sciWindVisAtt->SetForceWireframe(true);						//I believe that it might make Window transparent
  sciWindLog->SetVisAttributes(sciWindVisAtt);						//Assignment of visualization attributes to the logical volume of the Window
  
  // Printing the final settings...
  G4cout << "##################################################################"
	 << G4endl
	 << "#########  ActarSimSpecMATSciDetectorConstruction::ConstructSci()  #######"
	 << G4endl;
    G4cout << "##################################################################"
	 << G4endl;
  
  return sciCrysPhys;
}

//////////////////////////////////////////////////////////////////
/// Set the material the scintillator bulk is made of
void ActarSimSpecMATSciDetectorConstruction::SetSciSpecMATMaterial (G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) sciSpecMATMaterial = pttoMaterial;
}

//////////////////////////////////////////////////////////////////
/// Set the material the reflector is made of
void ActarSimSpecMATSciDetectorConstruction::SetReflectorSpecMATMaterial (G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) reflectorSpecMATMaterial = pttoMaterial;
}

//////////////////////////////////////////////////////////////////
/// Set the material the window is made of
void ActarSimSpecMATSciDetectorConstruction::SetWindowSpecMATMaterial (G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) windowSpecMATMaterial = pttoMaterial;
}

//////////////////////////////////////////////////////////////////
/// Set the material the housing is made of
void ActarSimSpecMATSciDetectorConstruction::SetHousingSpecMATMaterial (G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) housingSpecMATMaterial = pttoMaterial;
}

//////////////////////////////////////////////////////////////////
/// Set the material the vacuum flange is made of
void ActarSimSpecMATSciDetectorConstruction::SetVacuumFlangeSpecMATMaterial (G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) vacuumFlangeSpecMATMaterial = pttoMaterial;
}

//////////////////////////////////////////////////////////////////
/// Set the material the insulationTube is made of
void ActarSimSpecMATSciDetectorConstruction::SetInsulationTubeSpecMATMaterial (G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) insulationTubeSpecMATMaterial = pttoMaterial;
}



//////////////////////////////////////////////////////////////////
/// Updates Scintillator detector
void ActarSimSpecMATSciDetectorConstruction::UpdateGeometry() {
  Construct(detConstruction->GetWorldLogicalVolume());
  G4RunManager::GetRunManager()->
    DefineWorldVolume(detConstruction->GetWorldPhysicalVolume());
}

//////////////////////////////////////////////////////////////////
/// Prints Scintillator detector parameters. TODO: To be filled
void ActarSimSpecMATSciDetectorConstruction::PrintDetectorParameters() {
  G4cout << "##################################################################"
	 << G4endl
	 << "####  ActarSimSpecMATSciDetectorConstruction::PrintDetectorParameters() ####"
	 << G4endl;
  G4cout << "##################################################################"
	 << G4endl;
}

//////////////////////////////////////////////////////////////////
/// Calculates Auxiliar Radius for the location of the elements
G4double ActarSimSpecMATSciDetectorConstruction::ComputeCircleR1() {
  G4double circleR1;
  if (nbSegments == 1) {
    circleR1 = 0;
  }
  else if (nbSegments == 2) {
    circleR1 = 100;
  }
  else {
    if (vacuumChamberIncludedFlag == "on") {
      if (vacuumFlangeSizeY>sciHousSizeY*nbCrystInSegmentColumn) {
	circleR1 = vacuumFlangeSizeY/(std::tan(0.5*twopi/nbSegments));
      }
      else {
	circleR1 = sciHousSizeY*nbCrystInSegmentColumn/(std::tan(0.5*twopi/nbSegments));
      }
    }
    else {
      circleR1 = sciHousSizeY*nbCrystInSegmentColumn/(std::tan(0.5*twopi/nbSegments));
    }
  }
  return circleR1;
}
