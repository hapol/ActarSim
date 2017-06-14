// - AUTHOR: Hector Alvarez-Pol 03/20017
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/

#ifndef ActarSimSpecMATSciDetectorConstruction_h
#define ActarSimSpecMATSciDetectorConstruction_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class ActarSimSpecMATSciDetectorMessenger;
class ActarSimDetectorConstruction;

class ActarSimSpecMATSciDetectorConstruction {
private:
  G4Material* sciSpecMATMaterial;      ///< Pointer to the material the SpecMAT scintillator is made of
  G4Material* windowSpecMATMaterial;   ///< Pointer to the material the SpecMAT scintillator window is made of
  G4Material* reflectorSpecMATMaterial;///< Pointer to the material the SpecMAT scintillator reflector is made of
  G4Material* housingSpecMATMaterial;  ///< Pointer to the material the SpecMAT scintillator housing is made of
  G4Material* vacuumFlangeSpecMATMaterial;  ///< Pointer to the material the SpecMAT scintillator vacuum flange is made of
  G4Material* insulationTubeSpecMATMaterial; ///< Pointer to the material the SpecMAT scintillator insulator is made of

  ActarSimSpecMATSciDetectorMessenger* sciMessenger;    ///< Pointer to the Messenger
  ActarSimDetectorConstruction* detConstruction; ///< Pointer to the global detector

  /////////////COPIED - CHECK THEM
  G4String vacuumChamber;
  G4String insulationTube;

  G4int nbSegments;
  G4int nbCrystInSegmentRow;
  G4int nbCrystInSegmentColumn;

  G4double dPhi;
  G4double half_dPhi;
  G4double tandPhi;

  G4double circleR1;

  G4double sciCrystSizeX;
  G4double sciCrystSizeY;
  G4double sciCrystSizeZ;

  G4double sciCrystPosX;
  G4double sciCrystPosY;
  G4double sciCrystPosZ;

  G4double sciWindSizeX;
  G4double sciWindSizeY;
  G4double sciWindSizeZ;

  G4double sciWindPosX;
  G4double sciWindPosY;
  G4double sciWindPosZ;

  G4double sciReflWallThickX;
  G4double sciReflWallThickY;
  G4double sciReflWindThick;

  G4double sciReflSizeX;
  G4double sciReflSizeY;
  G4double sciReflSizeZ;

  G4double sciReflPosX;
  G4double sciReflPosY;
  G4double sciReflPosZ;

  G4double sciHousWallThickX;
  G4double sciHousWallThickY;
  G4double sciHousWindThick;

  G4double sciHousSizeX;
  G4double sciHousSizeY;
  G4double sciHousSizeZ;

  G4double sciHousPosX;
  G4double sciHousPosY;
  G4double sciHousPosZ;

  G4double insulationTubeThickness;

  G4double vacuumFlangeSizeX;
  G4double vacuumFlangeSizeY;
  G4double vacuumFlangeSizeZ;

  G4double vacuumFlangeThickFrontOfScint;

  //////////////////////////

  G4VPhysicalVolume* ConstructSci(G4LogicalVolume*);

public:
  ActarSimSpecMATSciDetectorConstruction(ActarSimDetectorConstruction*);
  ~ActarSimSpecMATSciDetectorConstruction();

  G4VPhysicalVolume* Construct(G4LogicalVolume*);

  void SetSciSpecMATMaterial (G4String);
  void SetReflectorSpecMATMaterial(G4String);
  void SetWindowSpecMATMaterial(G4String);
  void SetHousingSpecMATMaterial(G4String);
  void SetVacuumFlangeSpecMATMaterial(G4String);
  void SetInsulationTubeSpecMATMaterial(G4String);

  void SetSciCrystSizeX(G4double val){sciCrystSizeX = val;}
  void SetSciCrystSizeY(G4double val){sciCrystSizeY = val;}
  void SetSciCrystSizeZ(G4double val){sciCrystSizeZ = val;}

  void SetNbSegments(G4int val){nbSegments = val;}
  void SetNbCrystInSegmentRow(G4int val){nbCrystInSegmentRow = val;}
  void SetNbCrystInSegmentColumn(G4int val){nbCrystInSegmentColumn = val;}

  G4Material* GetSpecMATSciMaterial() {return sciSpecMATMaterial;}
  G4Material* GetReflectorSpecMATMaterial() {return reflectorSpecMATMaterial;}
  G4Material* GetWindowSpecMATMaterial() {return windowSpecMATMaterial;}
  G4Material* GetHousingSpecMATMaterial() {return housingSpecMATMaterial;}
  G4Material* GetVacuumFlangeSpecMATMaterial() {return vacuumFlangeSpecMATMaterial;}
  G4Material* GetInsulationTubeSpecMATMaterial() {return insulationTubeSpecMATMaterial;}

  G4double GetSciCrystSizeX(void){return sciCrystSizeX;}
  G4double GetSciCrystSizeY(void){return sciCrystSizeY;}
  G4double GetSciCrystSizeZ(void){return sciCrystSizeZ;}

  G4double GetNbSegments(void){return nbSegments;}
  G4double GetNbCrystInSegmentRow(void){return nbCrystInSegmentRow;}
  G4double GetNbCrystInSegmentColumn(void){return nbCrystInSegmentColumn;}

  G4double GetInsulationTubeThickness(void){return insulationTubeThickness;}
  G4double GetVacuumFlangeSizeX(void){return vacuumFlangeSizeX;}

  G4double ComputeCircleR1();

  void UpdateGeometry();
  void PrintDetectorParameters();
};
#endif
