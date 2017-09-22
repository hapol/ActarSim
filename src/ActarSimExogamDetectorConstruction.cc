// - AUTHOR: Pablo Cabanelas 12/2016
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/
//////////////////////////////////////////////////////////////////
/// \class ActarSimExogamDetectorConstruction
/// Exogam (array of HPGe detectors) detector description
/////////////////////////////////////////////////////////////////

#include "ActarSimExogamDetectorConstruction.hh"
#include "ActarSimDetectorConstruction.hh"
#include "ActarSimGasDetectorConstruction.hh"
#include "ActarSimROOTAnalysis.hh"
#include "ActarSimExogamSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4Transform3D.hh"

#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Para.hh"
#include "G4ThreeVector.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SDManager.hh"
#include <math.h>
#include "G4UserLimits.hh"
#include "G4ios.hh"

#include "globals.hh"

#include "HodoParametrisation.hh"
#include "G4PVParameterised.hh"

// The exogam array parts definitions:
//    #define PLAQALU
#define REPLICAS_CLOVER
#define REPLICAS

//////////////////////////////////////////////////////////////////
/// Constructor. Sets the material and the pointer to the Messenger
ActarSimExogamDetectorConstruction::
ActarSimExogamDetectorConstruction(ActarSimDetectorConstruction* det)
  :	exogamBulkMaterial(0),detConstruction(det) {

  SetExogamBulkMaterial("Germanium");

  // create commands for interactive definition of the detector
  //exogamMessenger = new ActarSimPlaDetectorMessenger(this);
}

//////////////////////////////////////////////////////////////////
/// Destructor
ActarSimExogamDetectorConstruction::~ActarSimExogamDetectorConstruction(){
  // delete exogamMessenger;
}

//////////////////////////////////////////////////////////////////
/// Wrap for the construction functions within the TOF
G4VPhysicalVolume* ActarSimExogamDetectorConstruction::Construct(G4LogicalVolume* worldLog) {
  //Introduce here other constructors for materials around the TOF (windows, frames...)
  //which can be controlled by the calMessenger
  return ConstructExogam(worldLog);
}

//////////////////////////////////////////////////////////////////
/// Construct Exogam Detector
G4VPhysicalVolume* ActarSimExogamDetectorConstruction::ConstructExogam(G4LogicalVolume* worldLog) {

  //Gas chamber position inside the chamber
  //ActarSimGasDetectorConstruction* gasDet = detConstruction->GetGasDetector();
  //G4double zGasBoxPosition=gasDet->GetGasBoxCenterZ();

  // Printing the final settings...
  G4cout << "##################################################################"
    << G4endl
    << "###########  ActarSimExogamDetectorConstruction::ConstructExogam() ####"
    << G4endl
    << " Detector material: "
    << exogamBulkMaterial
    << G4endl;
  G4cout << "##################################################################"
    << G4endl;


// =====================================================================

  //////////////////////////////////////////////////////////////////////////////////////////
  //
  //  EXOGAM - 4 clovers
  //
  //

  G4double xHalfLengthAlu1 = 6.2*CLHEP::cm;
  G4double yHalfLengthAlu1 = 6.2*CLHEP::cm;
  G4double zHalfLengthAlu1 = 5.50*CLHEP::cm;  // 5.35*cm;

  //G4double xHalfLengthPlAlu = 4.15*CLHEP::cm;
  //G4double yHalfLengthPlAlu = 4.15*CLHEP::cm;
  G4double zHalfLengthPlAlu = (0.15/2)*CLHEP::cm;     // 1.5 mm thick

  G4double xHalfLengthVac1 = (2.725+0.3/5)*CLHEP::cm;
  G4double yHalfLengthVac1 = (2.725+0.3/5)*CLHEP::cm;
  G4double zHalfLengthVac1 = 5.20*CLHEP::cm ; //4.55*cm;

  G4double pDzGer1 = 4.5*CLHEP::cm;
  G4double pThetaGer1 = 0.*CLHEP::deg;
  G4double pPhiGer1 = 0.*CLHEP::deg;
  G4double pDy1Ger1 = 1.104898186*CLHEP::cm;
  G4double pDx1Ger1 = 1.1385*CLHEP::cm;
  G4double pDx2Ger1 = 0.000000001*CLHEP::cm;
  G4double pAlp1Ger1 = 27.25778821*CLHEP::deg;
  G4double pDy2Ger1 = 1.104898186*CLHEP::cm;
  G4double pDx3Ger1 = 1.1385*CLHEP::cm;
  G4double pDx4Ger1 = 0.000000001*CLHEP::cm;
  G4double pAlp2Ger1 = 27.25778821*CLHEP::deg;

  G4double minRadiusGer2 = 0.*CLHEP::cm;
  G4double maxRadiusGer2 = 3.173*CLHEP::cm;
  G4double HalfLengthGer2 = 4.5*CLHEP::cm;
  G4double startPhiGer2 = -45.85811705*CLHEP::deg;
  G4double deltaPhiGer2 = 45.858117053*CLHEP::deg;

  G4double minRadiusGer3 = 0.*CLHEP::cm;
  G4double maxRadiusGer3 = 3.173*CLHEP::cm;
  G4double HalfLengthGer3 = 4.5*CLHEP::cm;
  G4double startPhiGer3 = 0*CLHEP::deg;
  G4double deltaPhiGer3 = 90*CLHEP::deg;

  G4double minRadiusGer4 = 0.*CLHEP::cm;
  G4double maxRadiusGer4 = 3.173*CLHEP::cm;
  G4double HalfLengthGer4 = 4.5*CLHEP::cm;
  G4double startPhiGer4 = 90*CLHEP::deg;
  G4double deltaPhiGer4 = 45.85811705*CLHEP::deg;

  G4double pDzGer5 = 4.5*CLHEP::cm;
  G4double pThetaGer5 = 0.*CLHEP::deg;
  G4double pPhiGer5 = 0.*CLHEP::deg;
  G4double pDy1Ger5 = 1.104898186*CLHEP::cm;
  G4double pDx1Ger5 = 1.1385*CLHEP::cm;
  G4double pDx2Ger5 = 0.000000001*CLHEP::cm;
  G4double pAlp1Ger5 = 27.25778821*CLHEP::deg;
  G4double pDy2Ger5 = 1.104898186*CLHEP::cm;
  G4double pDx3Ger5 = 1.1385*CLHEP::cm;
  G4double pDx4Ger5 = 0.000000001*CLHEP::cm;
  G4double pAlp2Ger5 = 27.25778821*CLHEP::deg;

  G4double pDzGer6 = 4.5*CLHEP::cm;
  G4double pThetaGer6 = 0.*CLHEP::deg;
  G4double pPhiGer6 = 0.*CLHEP::deg;
  G4double pDy1Ger6 = 0.837735041*CLHEP::cm;
  G4double pDx1Ger6 = 1.1385*CLHEP::cm;
  G4double pDx2Ger6 = 0.000000001*CLHEP::cm;
  G4double pAlp1Ger6 = 34.19653088*CLHEP::deg;
  G4double pDy2Ger6 = 0.837735041*CLHEP::cm;
  G4double pDx3Ger6 = 1.1385*CLHEP::cm;
  G4double pDx4Ger6 = 0.000000001*CLHEP::cm;
  G4double pAlp2Ger6 = 34.19653088*CLHEP::deg;

  G4double minRadiusGer7 = 0.*CLHEP::cm;
  G4double maxRadiusGer7 = 2.827*CLHEP::cm;
  G4double HalfLengthGer7 = 4.5*CLHEP::cm;
  G4double startPhiGer7 = 216.3465078*CLHEP::deg;
  G4double deltaPhiGer7 = 17.30698432*CLHEP::deg;

  G4double pDzGer8 = 4.5*CLHEP::cm;
  G4double pThetaGer8 = 0.*CLHEP::deg;
  G4double pPhiGer8 = 0.*CLHEP::deg;
  G4double pDy1Ger8 = 0.837735041*CLHEP::cm;
  G4double pDx1Ger8 = 1.1385*CLHEP::cm;
  G4double pDx2Ger8 = 0.000000001*CLHEP::cm;
  G4double pAlp1Ger8 = 34.19653088*CLHEP::deg;
  G4double pDy2Ger8 = 0.837735041*CLHEP::cm;
  G4double pDx3Ger8 = 1.1385*CLHEP::cm;
  G4double pDx4Ger8 = 0.000000001*CLHEP::cm;
  G4double pAlp2Ger8 = 34.19653088*CLHEP::deg;

//  G4double w1dx = (0.2)*CLHEP::cm;
//  G4double w1dy = (0.14)*CLHEP::cm;
//  G4double w1dz = (0.4)*CLHEP::cm;
  G4double w2dx1 = 1.243469*CLHEP::cm;
  G4double w2dx2 = (0.+0.001)*CLHEP::cm;
  G4double w2dy1 = 1.139*CLHEP::cm;
  G4double w2dy2 = 1.139*CLHEP::cm;
//  G4double w2dy1 = 1.1385*CLHEP::cm;
//  G4double w2dy2 = 1.1385*CLHEP::cm;
  G4double w2dz = (1.5+0.001)*CLHEP::cm;
  G4double w3dx1 = 1.243469*CLHEP::cm;
  G4double w3dx2 = (0.+0.001)*CLHEP::cm;
  G4double w3dy1 = (0.9652+1.)*CLHEP::cm;
  G4double w3dy2 = (0.9652+1.)*CLHEP::cm;
//  G4double w3dy1 = 0.965179656*CLHEP::cm;
//  G4double w3dy2 = 0.965179656*CLHEP::cm;
  G4double w3dz = (1.5+0.001)*CLHEP::cm;
  G4double w3bdx1 = 1.243469*CLHEP::cm;
  G4double w3bdx2 = (0.+0.001)*CLHEP::cm;
  G4double w3bdy1 = (0.9652+1.)*CLHEP::cm;
  G4double w3bdy2 = (0.9652+1.)*CLHEP::cm;
//  G4double w3bdy1 = 0.965179656*CLHEP::cm;
//  G4double w3bdy2 = 0.965179656*CLHEP::cm;
  G4double w3bdz = (1.5+0.001)*CLHEP::cm;
  G4double w4dx1 = 1.243469*CLHEP::cm;
  G4double w4dx2 = (0.+0.001)*CLHEP::cm;
  G4double w4dy1 = 1.139*CLHEP::cm;
  G4double w4dy2 = 1.139*CLHEP::cm;
//  G4double w4dy1 = 1.1385*CLHEP::cm;
//  G4double w4dy2 = 1.1385*CLHEP::cm;
  G4double w4dz = (1.5+0.001)*CLHEP::cm;


  ////////////////////////////////////////////////////////
  //
  //    Position Vectors

//Alu1 = Clover #0
  G4ThreeVector positionAlu1 = G4ThreeVector(-(zHalfLengthAlu1+13.25*CLHEP::cm),0.*CLHEP::cm,-5.5*CLHEP::cm);

#ifdef REPLICAS_CLOVER
  G4ThreeVector positionAlu2 = G4ThreeVector(-(zHalfLengthAlu1+13.25*CLHEP::cm),0.*CLHEP::cm, 5.5*CLHEP::cm);

  G4ThreeVector positionAlu3 = G4ThreeVector( (zHalfLengthAlu1+13.25*CLHEP::cm),0.*CLHEP::cm,-5.5*CLHEP::cm);
  G4ThreeVector positionAlu4 = G4ThreeVector( (zHalfLengthAlu1+13.25*CLHEP::cm),0.*CLHEP::cm, 5.5*CLHEP::cm);

  G4ThreeVector positionAlu5 = G4ThreeVector( -5.5*CLHEP::cm,0.*CLHEP::cm,(zHalfLengthAlu1+13.25*CLHEP::cm));
  G4ThreeVector positionAlu6 = G4ThreeVector(  5.5*CLHEP::cm,0.*CLHEP::cm,(zHalfLengthAlu1+13.25*CLHEP::cm));
#endif

//PlAlu
  G4ThreeVector positionPlAlu = G4ThreeVector(0.*CLHEP::cm,0.*CLHEP::cm,zHalfLengthPlAlu-zHalfLengthVac1);
// Vac1
  G4ThreeVector positionVac1 = G4ThreeVector(xHalfLengthVac1,yHalfLengthVac1,
                                                                       0.*CLHEP::cm);
// vac2
  G4ThreeVector positionVac2 = G4ThreeVector(-xHalfLengthVac1,yHalfLengthVac1,
                                                                     -0.*CLHEP::cm);
// vac3
  G4ThreeVector positionVac3 = G4ThreeVector(-xHalfLengthVac1,-yHalfLengthVac1,
                                                                     -0.*CLHEP::cm);
// vac4
  G4ThreeVector positionVac4 = G4ThreeVector(xHalfLengthVac1,-yHalfLengthVac1,
                                                                     -0.*CLHEP::cm);


//Ger1 & Wed1
  G4ThreeVector positionGer1 = G4ThreeVector((0.656898185-0.3/5)*CLHEP::cm,(-2.15575-0.3/5)*CLHEP::cm,0.*CLHEP::cm);
//  G4ThreeVector positionGer1 = G4ThreeVector(pDy1Ger1,-pDx1Ger1*1.5,0.*CLHEP::cm);
  G4ThreeVector positionG1 = G4ThreeVector(0.*CLHEP::cm,0.*CLHEP::cm,0.*CLHEP::cm);
  G4ThreeVector positionWed1 = G4ThreeVector((pDx1Ger1/2-0.15*CLHEP::cm),(pDy1Ger1),(-pDzGer1+0.25*CLHEP::cm));
//  G4RotationMatrix rm12;  // Wed1
//  rm12.rotateX(M_PI*(-22.5/180.));

// Ger2 & Wed2
  G4ThreeVector positionGer2 = G4ThreeVector((-0.448-0.3/5)*CLHEP::cm,(-0.448-0.3/5)*CLHEP::cm,0.*CLHEP::cm);

//  G4ThreeVector positionGer2 = G4ThreeVector(0.*CLHEP::cm,0.*CLHEP::cm,0.*CLHEP::cm);
  G4ThreeVector positionWed2 = G4ThreeVector(3.173*CLHEP::cm,-1.1385*CLHEP::cm,(-3.0-0.001)*CLHEP::cm);

// Ger3
  G4ThreeVector positionGer3 = G4ThreeVector((-0.448-0.3/5)*CLHEP::cm,(-0.448-0.3/5)*CLHEP::cm,0.*CLHEP::cm);
//   G4ThreeVector positionGer3 = G4ThreeVector(0.*CLHEP::cm,0.*CLHEP::cm,0.*CLHEP::cm);
   G4ThreeVector positionWed3 = G4ThreeVector(3.173*CLHEP::cm,(w3dy1-0.001*CLHEP::cm),(-3.0-0.001)*CLHEP::cm);
   G4ThreeVector positionWed3b = G4ThreeVector((w3bdy1-0.001*CLHEP::cm),3.173*CLHEP::cm,(-3.0-0.001)*CLHEP::cm);

// Ger4 & Wed4
  G4ThreeVector positionGer4 = G4ThreeVector((-0.448-0.3/5)*CLHEP::cm,(-0.448-0.3/5)*CLHEP::cm,0.*CLHEP::cm);
//  G4ThreeVector positionGer4 = G4ThreeVector(0.*CLHEP::cm,0.*CLHEP::cm,0.*CLHEP::cm);
  G4ThreeVector positionWed4 = G4ThreeVector(-1.1385*CLHEP::cm,3.173*CLHEP::cm,(-3.0-0.001)*CLHEP::cm);

// Ger5
  G4ThreeVector positionGer5 = G4ThreeVector((-2.15575-0.3/5)*CLHEP::cm,(0.656898185-0.3/5)*CLHEP::cm,0.*CLHEP::cm);
//  G4ThreeVector positionGer5 = G4ThreeVector(-pDx1Ger5*1.5,pDy1Ger5,0.*CLHEP::cm);
  G4RotationMatrix rm51;  // Ger5
  rm51.rotateY(M_PI*(1.));

// Ger6
  G4ThreeVector positionGer6 = G4ThreeVector((-2.15575-0.3/5)*CLHEP::cm,(-1.285735041-0.3/5)*CLHEP::cm,0.*CLHEP::cm);
//  G4ThreeVector positionGer6 = G4ThreeVector(-pDx1Ger6*1.5,-pDy1Ger6,0.*CLHEP::cm);
  G4RotationMatrix rm61;  // Ger6
  rm61.rotateZ(M_PI*(1.));

// Ger7
  G4ThreeVector positionGer7 = G4ThreeVector((-0.448-0.3/5)*CLHEP::cm,(-0.448-0.3/5)*CLHEP::cm,0.*CLHEP::cm);
//  G4ThreeVector positionGer7 = G4ThreeVector(0.*CLHEP::cm,0.*CLHEP::cm,0.*CLHEP::cm);

// Ger8
  G4ThreeVector positionGer8 = G4ThreeVector((-1.285735041-0.3/5)*CLHEP::cm,(-2.15575-0.3/5)*CLHEP::cm,0.*CLHEP::cm);
//  G4ThreeVector positionGer8 = G4ThreeVector(-pDy1Ger8,-pDx1Ger8*1.5,0.*CLHEP::cm);

/*
  G4RotationMatrix rmAlu1y; // Alu1
  rmAlu1y.rotateY(M_PI*(1./2.));
  G4RotationMatrix rmAlu1x; // Alu1
  rmAlu1x.rotateX(M_PI*(-1./4.));
  G4RotationMatrix rmAlu1; // Alu1
  rmAlu1 = rmAlu1y*rmAlu1x;
*/


  G4RotationMatrix rmAlu1; // Alu1
  rmAlu1.rotateY(M_PI*(-1./2.));


#ifdef REPLICAS_CLOVER
  G4RotationMatrix rmAlu3; // Alu3
  rmAlu3.rotateY(M_PI*(1./2.));

  G4RotationMatrix rmAluZ2;
  rmAluZ2.rotateZ(M_PI*(1./2.));

  G4RotationMatrix rmAluZ3;
  rmAluZ3.rotateZ(M_PI*(1./1.));

  G4RotationMatrix rmAluZ4;
  rmAluZ4.rotateZ(M_PI*(3./2.));

  G4RotationMatrix rmAlu2; // Alu2
  rmAlu2 = rmAluZ2*rmAlu3 ;

  G4RotationMatrix rmAlu4; // Alu4
  rmAlu4.rotateZ(M_PI*(-1./2.));

  G4RotationMatrix rmAlu5; // Alu5
  rmAlu5.rotateZ(M_PI*(-1./2.));

#endif

  G4RotationMatrix rm;
  rm.rotateZ(M_PI*(1./2.));

  G4RotationMatrix rm1;
  rm1.rotateZ(M_PI*(1./1.));

  G4RotationMatrix rm2;
  rm2.rotateZ(M_PI*(3./2.));

  G4RotationMatrix rm11;  // Ger1
  rm11.rotateZ(M_PI*(-1./2.));

  G4RotationMatrix rm81;
  rm81.rotateZ(M_PI*(3./2.));

  G4RotationMatrix rm82;
  rm82.rotateX(M_PI*(1./1.));

  G4RotationMatrix rm8; // Ger8
  rm8 = rm81*rm82 ;


  ////////////////////////////////////////////////////////
  //
  //   Al capsule ( "Alu1" )

  G4double density, a, z;
  G4Material* Vacuum =
    new G4Material("Galactic", z=1., a=1.01*CLHEP::g/CLHEP::mole,density= CLHEP::universe_mean_density,
                   kStateGas, 3.e-18*CLHEP::pascal, 2.73*CLHEP::kelvin);
  Vacuum->SetChemicalFormula("NOMATTER");

  G4Box* solidAlu1 = new G4Box("Alu1",xHalfLengthAlu1,yHalfLengthAlu1,
                                                     zHalfLengthAlu1);
  G4LogicalVolume* logicAlu1 = new G4LogicalVolume(solidAlu1,Vacuum,"Alu1");
  G4VPhysicalVolume* physiAlu1 = new G4PVPlacement(G4Transform3D(rmAlu1,positionAlu1),logicAlu1,"Alu1",worldLog,false,0); //coupled to Alu2

#ifdef REPLICAS_CLOVER
   new G4PVPlacement(G4Transform3D(rmAlu1,positionAlu2),logicAlu1,"Alu2",worldLog,false,1); //coupled to Alu1

   new G4PVPlacement(G4Transform3D(rmAlu3,positionAlu3),logicAlu1,"Alu3",worldLog,false,2); //coupled to Alu4
   new G4PVPlacement(G4Transform3D(rmAlu3,positionAlu4),logicAlu1,"Alu4",worldLog,false,3); //coupled to Alu3

   new G4PVPlacement(G4Transform3D(rmAlu5,positionAlu5),logicAlu1,"Alu5",worldLog,false,4); //coupled to Alu6
   new G4PVPlacement(G4Transform3D(rmAlu5,positionAlu6),logicAlu1,"Alu6",worldLog,false,5); //coupled to Alu5
#endif



  // Vacuum Inside The Al Capsule ( "Vac1" )
  G4Box* solidVac1 = new G4Box("Vac1",xHalfLengthVac1,yHalfLengthVac1,
                                                     zHalfLengthVac1);

  G4LogicalVolume* logicVac1 = new G4LogicalVolume(solidVac1,Vacuum,"Vac1");
  G4VPhysicalVolume* physiVac1 = new G4PVPlacement(0,positionVac1,"Vac1",
                                            logicVac1,physiAlu1,false,0);


  // Vacuum Inside The Al Capsule ( "Vac2" )
  G4Box* solidVac2 = new G4Box("Vac2",xHalfLengthVac1,yHalfLengthVac1,
                                                     zHalfLengthVac1);
  G4LogicalVolume* logicVac2 = new G4LogicalVolume(solidVac2,Vacuum,"Vac2");
  G4VPhysicalVolume* physiVac2 = new G4PVPlacement(G4Transform3D(rm,
                  positionVac2),"Vac2",logicVac2,physiAlu1,false,0);

  // Vacuum Inside The Al Capsule ( "Vac3" )
  G4Box* solidVac3 = new G4Box("Vac3",xHalfLengthVac1,yHalfLengthVac1,
                                                     zHalfLengthVac1);
  G4LogicalVolume* logicVac3 = new G4LogicalVolume(solidVac3,Vacuum,"Vac3");
  G4VPhysicalVolume* physiVac3 = new G4PVPlacement(G4Transform3D(rm1,
                   positionVac3),"Vac3",logicVac3,physiAlu1,false,0);

  // Vacuum Inside The Al Capsule ( "Vac4" )
  G4Box* solidVac4 = new G4Box("Vac4",xHalfLengthVac1,yHalfLengthVac1,
                                                     zHalfLengthVac1);
  G4LogicalVolume* logicVac4 = new G4LogicalVolume(solidVac4,Vacuum,"Vac4");
  G4VPhysicalVolume* physiVac4 = new G4PVPlacement(G4Transform3D(rm2,
                   positionVac4),"Vac4",logicVac4,physiAlu1,false,0);

#ifdef PLAQALU
  // Plaque Aluminium in front Germanium detector
  G4Box* solidPlAlu = new G4Box("PlAlu",xHalfLengthPlAlu,yHalfLengthPlAlu,
                                                     zHalfLengthPlAlu);

  G4LogicalVolume* logicPlAlu = new G4LogicalVolume(solidPlAlu,Aluminium,"PlAlu");

  G4VPhysicalVolume* physiPlAlu = new  G4PVPlacement(0,positionPlAlu,"PlAlu",logicPlAlu,physiVac1,false,0);

#ifdef REPLICAS
  physiPlAlu = new  G4PVPlacement(0,positionPlAlu,"PlAlu",logicPlAlu,physiVac2,false,1);
  physiPlAlu = new  G4PVPlacement(0,positionPlAlu,"PlAlu",logicPlAlu,physiVac3,false,2);
  physiPlAlu = new  G4PVPlacement(0,positionPlAlu,"PlAlu",logicPlAlu,physiVac4,false,3);
#endif

#endif

// Germanium Detector Component Number One ( "Ger1" )
   G4Trap* solidGer1 = new G4Trap("Ger1",pDzGer1,pThetaGer1,pPhiGer1,
                                         pDy1Ger1,pDx1Ger1,pDx2Ger1,pAlp1Ger1,
                                         pDy2Ger1,pDx3Ger1,pDx4Ger1,pAlp2Ger1);
   G4LogicalVolume* logicGer1 = new G4LogicalVolume(solidGer1,exogamBulkMaterial,"Ger1");
   G4VPhysicalVolume* physiGer1 = new G4PVPlacement(G4Transform3D(rm11,
                     positionGer1),"Ger1",logicGer1,physiVac1,false,0);

#ifdef REPLICAS

   G4LogicalVolume* logicGer12 = new G4LogicalVolume(solidGer1,exogamBulkMaterial,"Ger1");
   G4VPhysicalVolume* physiGer12 = new G4PVPlacement(G4Transform3D(rm11,
                     positionGer1),"Ger1",logicGer12,physiVac2,false,1);

   G4LogicalVolume* logicGer22 = new G4LogicalVolume(solidGer1,exogamBulkMaterial,"Ger1");
   G4VPhysicalVolume* physiGer22 = new G4PVPlacement(G4Transform3D(rm11,
                     positionGer1),"Ger1",logicGer22,physiVac3,false,2);

   G4LogicalVolume* logicGer31 = new G4LogicalVolume(solidGer1,exogamBulkMaterial,"Ger1");
   G4VPhysicalVolume* physiGer31 = new G4PVPlacement(G4Transform3D(rm11,
                     positionGer1),"Ger1",logicGer31,physiVac4,false,3);

#endif


// exogamBulkMaterial Detector Component Number Two ( "Ger2" )

   G4Tubs* solidG2 = new G4Tubs("G2",minRadiusGer2,maxRadiusGer2,
                          HalfLengthGer2,startPhiGer2,deltaPhiGer2);

// We now substract the wedge (Wed2)

   G4Trd* solidWed2 = new G4Trd("Wed2",w2dx1,w2dx2,w2dy1,w2dy2,w2dz);


   G4SubtractionSolid* G2minusW2=new G4SubtractionSolid("Ger2",solidG2,solidWed2,0,positionWed2);
   G4LogicalVolume* logicGer2 = new G4LogicalVolume(G2minusW2,exogamBulkMaterial,"Ger2");
   G4VPhysicalVolume* physiGer2 = new G4PVPlacement(0,
                                   positionGer2,"Ger2",
                                   logicGer2,physiVac1,false,0);

#ifdef REPLICAS
//
// REPLICAS
//
   G4LogicalVolume* logicGer13 = new G4LogicalVolume(G2minusW2,exogamBulkMaterial,"Ger2");
   G4VPhysicalVolume* physiGer13 = new G4PVPlacement(0,positionGer2,"Ger2",
                                               logicGer13,physiVac2,false,1);

   G4LogicalVolume* logicGer23 = new G4LogicalVolume(G2minusW2,exogamBulkMaterial,"Ger2");
   G4VPhysicalVolume* physiGer23 = new G4PVPlacement(0,positionGer2,"Ger2",
                                               logicGer23,physiVac3,false,2);

   G4LogicalVolume* logicGer32 = new G4LogicalVolume(G2minusW2,exogamBulkMaterial,"Ger2");
   G4VPhysicalVolume* physiGer32 = new G4PVPlacement(0,positionGer2,"Ger2",
                                               logicGer32,physiVac4,false,3);

#endif


// exogamBulkMaterial Detector Component Number Three ( "Ger3" )

   G4Tubs* solidG3 = new G4Tubs("G3",minRadiusGer3,maxRadiusGer3,
                          HalfLengthGer3,startPhiGer3,deltaPhiGer3);

// We now substract Wed3 + Wed3b

   G4Trd* solidWed3 = new G4Trd("Wed3",w3dx1,w3dx2,w3dy1,w3dy2,w3dz);
   G4Trd* solidWed3b = new G4Trd("Wed3b",w3bdy1,w3bdy2,w3bdx1,w3bdx2,w3bdz);

   G4SubtractionSolid* G3minusW3=new G4SubtractionSolid("Ger3-Wed3",solidG3,solidWed3,0,positionWed3);
   G4SubtractionSolid* G3minusW3b=new G4SubtractionSolid("Ger3",G3minusW3,solidWed3b,0,positionWed3b);
   G4LogicalVolume* logicGer3 = new G4LogicalVolume(G3minusW3b,exogamBulkMaterial,"Ger3");
   G4VPhysicalVolume* physiGer3 = new G4PVPlacement(0,
                                   positionGer3,"Ger3",
                                   logicGer3,physiVac1,false,0);

#ifdef REPLICAS

//   // REPLICAS

   G4LogicalVolume* logicGer14 = new G4LogicalVolume(G3minusW3b,exogamBulkMaterial,"Ger3");
   G4VPhysicalVolume* physiGer14 = new G4PVPlacement(0,positionGer3,"Ger3",
                                               logicGer14,physiVac2,false,1);

   G4LogicalVolume* logicGer24 = new G4LogicalVolume(G3minusW3b,exogamBulkMaterial,"Ger3");
   G4VPhysicalVolume* physiGer24 = new G4PVPlacement(0,positionGer3,"Ger3",
                                               logicGer24,physiVac3,false,2);

   G4LogicalVolume* logicGer25 = new G4LogicalVolume(G3minusW3b,exogamBulkMaterial,"Ger3");
   G4VPhysicalVolume* physiGer25 = new G4PVPlacement(0,positionGer3,"Ger3",
                                               logicGer25,physiVac4,false,3);

#endif


// exogamBulkMaterial Detector Component Number Two ( "Ger4" )

   G4Tubs* solidG4 = new G4Tubs("G4",minRadiusGer4,maxRadiusGer4,
                          HalfLengthGer4,startPhiGer4,deltaPhiGer4);

// wed4

   G4Trd* solidWed4 = new G4Trd("Wed4",w4dy1,w4dy2,w4dx1,w4dx2,w4dz);

   G4SubtractionSolid* G4minusW4=new G4SubtractionSolid("Ger4",solidG4,solidWed4,0,positionWed4);
   G4LogicalVolume* logicGer4 = new G4LogicalVolume(G4minusW4,exogamBulkMaterial,
                                 "Ger4");
   G4VPhysicalVolume* physiGer4 = new G4PVPlacement(0,
                                   positionGer4,"Ger4",
                                   logicGer4,physiVac1,false,0);

#ifdef REPLICAS
//
// REPLICAS

   G4LogicalVolume* logicGer15 = new G4LogicalVolume(G4minusW4,exogamBulkMaterial,"Ger4");
   G4VPhysicalVolume* physiGer15 = new G4PVPlacement(0,positionGer4,"Ger4",
                                               logicGer15,physiVac2,false,1);

   G4LogicalVolume* logicGer17 = new G4LogicalVolume(G4minusW4,exogamBulkMaterial,"Ger4");
   G4VPhysicalVolume* physiGer17 = new G4PVPlacement(0,positionGer4,"Ger4",
                                               logicGer17,physiVac3,false,2);

   G4LogicalVolume* logicGer26 = new G4LogicalVolume(G4minusW4,exogamBulkMaterial,"Ger4");
   G4VPhysicalVolume* physiGer26 = new G4PVPlacement(0,positionGer4,"Ger4",
                                               logicGer26,physiVac4,false,3);

#endif


// exogamBulkMaterial Detector Component Number Five ( "Ger5" )

   G4Trap* solidGer5 = new G4Trap("Ger5",pDzGer5,pThetaGer5,pPhiGer5,
   pDy1Ger5,pDx1Ger5,pDx2Ger5,pAlp1Ger5,
                              pDy2Ger5,pDx3Ger5,pDx4Ger5,pAlp2Ger5);
   G4LogicalVolume* logicGer5 = new G4LogicalVolume(solidGer5,exogamBulkMaterial,"Ger5");
   G4VPhysicalVolume* physiGer5 = new G4PVPlacement(G4Transform3D(rm51,
                                                  positionGer5),"Ger5",
                                                  logicGer5,physiVac1,false,0);

#ifdef REPLICAS

//
//     REPLICAS
//
   G4LogicalVolume* logicGer16 = new G4LogicalVolume(solidGer5,exogamBulkMaterial,"Ger5");
   G4VPhysicalVolume* physiGer16 = new G4PVPlacement(G4Transform3D(rm51,
                                                  positionGer5),"Ger5",
                                                  logicGer16,physiVac2,false,1);

   G4LogicalVolume* logicGer18 = new G4LogicalVolume(solidGer5,exogamBulkMaterial,"Ger5");
   G4VPhysicalVolume* physiGer18 = new G4PVPlacement(G4Transform3D(rm51,
                                                  positionGer5),"Ger5",
                                                  logicGer18,physiVac3,false,2);

   G4LogicalVolume* logicGer27 = new G4LogicalVolume(solidGer5,exogamBulkMaterial,"Ger5");
   G4VPhysicalVolume* physiGer27 = new G4PVPlacement(G4Transform3D(rm51,
                                                  positionGer5),"Ger5",
                                                  logicGer27,physiVac4,false,3);

#endif

// exogamBulkMaterial Detector Component Number Six ( "Ger6" )

   G4Trap* solidGer6 = new G4Trap("Ger6",pDzGer6,pThetaGer6,pPhiGer6,
   pDy1Ger6,pDx1Ger6,pDx2Ger6,pAlp1Ger6,
                              pDy2Ger6,pDx3Ger6,pDx4Ger6,pAlp2Ger6);
   G4LogicalVolume* logicGer6 = new G4LogicalVolume(solidGer6,exogamBulkMaterial,"Ger6");
   G4VPhysicalVolume* physiGer6 = new G4PVPlacement(G4Transform3D(rm61,
                             positionGer6),"Ger6",logicGer6,physiVac1,false,0);

#ifdef REPLICAS

//
// REPLICAS
//
   G4LogicalVolume* logicGer9 = new G4LogicalVolume(solidGer6,exogamBulkMaterial,"Ger6");
   G4VPhysicalVolume* physiGer9 = new G4PVPlacement(G4Transform3D(rm61,
                             positionGer6),"Ger6",logicGer9,physiVac2,false,1);

   G4LogicalVolume* logicGer19 = new G4LogicalVolume(solidGer6,exogamBulkMaterial,"Ger6");
   G4VPhysicalVolume* physiGer19 = new G4PVPlacement(G4Transform3D(rm61,
                             positionGer6),"Ger6",logicGer19,physiVac3,false,2);

   G4LogicalVolume* logicGer28 = new G4LogicalVolume(solidGer6,exogamBulkMaterial,"Ger6");
   G4VPhysicalVolume* physiGer28 = new G4PVPlacement(G4Transform3D(rm61,
                             positionGer6),"Ger6",logicGer28,physiVac4,false,3);

#endif


// exogamBulkMaterial Detector Component Number Seven ( "Ger7" )

   G4Tubs* solidGer7 = new G4Tubs("Ger7",minRadiusGer7,maxRadiusGer7,
                          HalfLengthGer7,startPhiGer7,deltaPhiGer7);
   G4LogicalVolume* logicGer7 = new G4LogicalVolume(solidGer7,exogamBulkMaterial,"Ger7");
   G4VPhysicalVolume* physiGer7 = new G4PVPlacement(0,positionGer7,"Ger7",
                                           logicGer7, physiVac1,false,0);

#ifdef REPLICAS

//
// REPLICAS
//
   G4LogicalVolume* logicGer10 = new G4LogicalVolume(solidGer7,exogamBulkMaterial,"Ger7");
   G4VPhysicalVolume* physiGer10 = new G4PVPlacement(0,positionGer7,"Ger7",
                                               logicGer10,physiVac2,false,1);

   G4LogicalVolume* logicGer20 = new G4LogicalVolume(solidGer7,exogamBulkMaterial,"Ger7");
   G4VPhysicalVolume* physiGer20 = new G4PVPlacement(0,positionGer7,"Ger7",
                                               logicGer20,physiVac3,false,2);

   G4LogicalVolume* logicGer29 = new G4LogicalVolume(solidGer7,exogamBulkMaterial,"Ger7");
   G4VPhysicalVolume* physiGer29 = new G4PVPlacement(0,positionGer7,"Ger7",
                                               logicGer29,physiVac4,false,3);

#endif

// exogamBulkMaterial Detector Component Number Eight ( "Ger8" )

   G4Trap* solidGer8 = new G4Trap("Ger8",pDzGer8,pThetaGer8,pPhiGer8,
   pDy1Ger8,pDx1Ger8,pDx2Ger8,pAlp1Ger8,
                              pDy2Ger8,pDx3Ger8,pDx4Ger8,pAlp2Ger8);

   G4LogicalVolume* logicGer8 = new G4LogicalVolume(solidGer8,exogamBulkMaterial,"Ger8");
   G4VPhysicalVolume* physiGer8 = new G4PVPlacement(G4Transform3D(rm8,
                    positionGer8),"Ger8",logicGer8,physiVac1,false,0);

#ifdef REPLICAS

//
// REPLICAS
//
   G4LogicalVolume* logicGer11 = new G4LogicalVolume(solidGer8,exogamBulkMaterial,"Ger8");
   G4VPhysicalVolume* physiGer11 = new G4PVPlacement(G4Transform3D(rm8,
                                                 positionGer8),"Ger8",
                                                  logicGer11,physiVac2,false,1);

   G4LogicalVolume* logicGer21 = new G4LogicalVolume(solidGer8,exogamBulkMaterial,"Ger8");
   G4VPhysicalVolume* physiGer21 = new G4PVPlacement(G4Transform3D(rm8,
                                                 positionGer8),"Ger8",
                                                  logicGer21,physiVac3,false,2);

   G4LogicalVolume* logicGer30 = new G4LogicalVolume(solidGer8,exogamBulkMaterial,"Ger8");
   G4VPhysicalVolume* physiGer30 = new G4PVPlacement(G4Transform3D(rm8,
                                                 positionGer8),"Ger8",
                                                  logicGer30,physiVac4,false,3);

#endif

// Visualization attributes

  G4VisAttributes* AluVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0));  // Yellow
  AluVisAtt ->SetForceWireframe(true);
//logicAlu1 ->SetVisAttributes(AluVisAtt);
  logicAlu1 ->SetVisAttributes(G4VisAttributes::Invisible);


#ifdef PLAQALU
  G4VisAttributes* PlAluVisAtt= new G4VisAttributes(G4Colour(0.0,0.8,0.8));  // Turquoise
  PlAluVisAtt ->SetForceSolid(true);
  logicPlAlu ->SetVisAttributes(PlAluVisAtt);
#endif

  G4VisAttributes* VacVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0));  // Magenta
  VacVisAtt ->SetForceWireframe(true);
  logicVac1 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicVac2 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicVac3 ->SetVisAttributes(G4VisAttributes::Invisible);
  logicVac4 ->SetVisAttributes(G4VisAttributes::Invisible);


  G4VisAttributes* Crystal1VisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));  // Red (Crystal A)
  Crystal1VisAtt ->SetForceSolid(true);
  logicGer1 ->SetVisAttributes(Crystal1VisAtt);
  logicGer2 ->SetVisAttributes(Crystal1VisAtt);
  logicGer3 ->SetVisAttributes(Crystal1VisAtt);
  logicGer4 ->SetVisAttributes(Crystal1VisAtt);
  logicGer5 ->SetVisAttributes(Crystal1VisAtt);
  logicGer6 ->SetVisAttributes(Crystal1VisAtt);
  logicGer7 ->SetVisAttributes(Crystal1VisAtt);
  logicGer8 ->SetVisAttributes(Crystal1VisAtt);


#ifdef REPLICAS
  G4VisAttributes* Crystal2VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0));  // Green (Crystal B)

//  Crystal2VisAtt ->SetForceWireframe(true);
  Crystal2VisAtt ->SetForceSolid(true);
  logicGer9 ->SetVisAttributes(Crystal2VisAtt);
  logicGer10 ->SetVisAttributes(Crystal2VisAtt);
  logicGer11 ->SetVisAttributes(Crystal2VisAtt);
  logicGer12 ->SetVisAttributes(Crystal2VisAtt);
  logicGer13 ->SetVisAttributes(Crystal2VisAtt);
  logicGer14 ->SetVisAttributes(Crystal2VisAtt);
  logicGer15 ->SetVisAttributes(Crystal2VisAtt);
  logicGer16 ->SetVisAttributes(Crystal2VisAtt);

  G4VisAttributes* Crystal3VisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5));  // Gray (Crystal C)

//  Crystal3VisAtt ->SetForceWireframe(true);
  Crystal3VisAtt ->SetForceSolid(true);
  logicGer17 ->SetVisAttributes(Crystal3VisAtt);
  logicGer18 ->SetVisAttributes(Crystal3VisAtt);
  logicGer19 ->SetVisAttributes(Crystal3VisAtt);
  logicGer20 ->SetVisAttributes(Crystal3VisAtt);
  logicGer21 ->SetVisAttributes(Crystal3VisAtt);
  logicGer22 ->SetVisAttributes(Crystal3VisAtt);
  logicGer23 ->SetVisAttributes(Crystal3VisAtt);
  logicGer24 ->SetVisAttributes(Crystal3VisAtt);

  G4VisAttributes* Crystal4VisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));  // Blue (Crystal D)

//  Crystal4VisAtt ->SetForceWireframe(true);
  Crystal4VisAtt ->SetForceSolid(true);
  logicGer25 ->SetVisAttributes(Crystal4VisAtt);
  logicGer26 ->SetVisAttributes(Crystal4VisAtt);
  logicGer27 ->SetVisAttributes(Crystal4VisAtt);
  logicGer28 ->SetVisAttributes(Crystal4VisAtt);
  logicGer29 ->SetVisAttributes(Crystal4VisAtt);
  logicGer30 ->SetVisAttributes(Crystal4VisAtt);
  logicGer31 ->SetVisAttributes(Crystal4VisAtt);
  logicGer32 ->SetVisAttributes(Crystal4VisAtt);
#endif  // REPLICAS


// =====================================================================
  //------------------------------------------------
  // Sensitive detectors
  //------------------------------------------------
  logicGer1 ->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer2 ->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer3 ->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer4 ->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer5 ->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer6 ->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer7 ->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer8 ->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer9 ->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer10->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer11->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer12->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer13->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer14->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer15->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer16->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer17->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer18->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer19->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer20->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer21->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer22->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer23->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer24->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer25->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer26->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer27->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer28->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer29->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer20->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer31->SetSensitiveDetector( detConstruction->GetExogamSD() );
  logicGer32->SetSensitiveDetector( detConstruction->GetExogamSD() );


  return physiGer1;
  return physiGer2;
  return physiGer3;
  return physiGer4;
  return physiGer5;
  return physiGer6;
  return physiGer7;
  return physiGer8;
  return physiGer9;
  return physiGer10;
  return physiGer11;
  return physiGer12;
  return physiGer13;
  return physiGer14;
  return physiGer15;
  return physiGer16;
  return physiGer17;
  return physiGer18;
  return physiGer19;
  return physiGer20;
  return physiGer21;
  return physiGer22;
  return physiGer23;
  return physiGer24;
  return physiGer25;
  return physiGer26;
  return physiGer27;
  return physiGer28;
  return physiGer29;
  return physiGer30;
  return physiGer31;
  return physiGer32;

}

//////////////////////////////////////////////////////////////////
/// Set the material the scintillator bulk is made of
void ActarSimExogamDetectorConstruction::SetExogamBulkMaterial (G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) exogamBulkMaterial = pttoMaterial;
}

//////////////////////////////////////////////////////////////////
/// Updates Exogam detector
void ActarSimExogamDetectorConstruction::UpdateGeometry() {
  Construct(detConstruction->GetWorldLogicalVolume());
  G4RunManager::GetRunManager()->
    DefineWorldVolume(detConstruction->GetWorldPhysicalVolume());
}

//////////////////////////////////////////////////////////////////
/// Prints Exogam detector parameters. To be filled
void ActarSimExogamDetectorConstruction::PrintDetectorParameters() {
  G4cout << "########################################################################"
	 << G4endl
	 << "####  ActarSimExogamDetectorConstruction::PrintDetectorParameters() ####"
	 << G4endl;
  G4cout << "########################################################################"
	 << G4endl;
}
