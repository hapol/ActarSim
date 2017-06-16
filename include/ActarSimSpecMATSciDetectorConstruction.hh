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
  G4Material* sciSpecMATMaterial;            ///< Pointer to the material the SpecMAT scintillator is made of
  G4Material* windowSpecMATMaterial;         ///< Pointer to the material the SpecMAT scintillator window is made of
  G4Material* reflectorSpecMATMaterial;      ///< Pointer to the material the SpecMAT scintillator reflector is made of
  G4Material* housingSpecMATMaterial;        ///< Pointer to the material the SpecMAT scintillator housing is made of
  G4Material* vacuumFlangeSpecMATMaterial;   ///< Pointer to the material the SpecMAT scintillator vacuum flange is made of
  G4Material* insulationTubeSpecMATMaterial; ///< Pointer to the material the SpecMAT scintillator insulator is made of
  
  ActarSimSpecMATSciDetectorMessenger* sciMessenger; ///< Pointer to the Messenger
  ActarSimDetectorConstruction* detConstruction;     ///< Pointer to the global detector
  
  G4String vacuumChamberIncludedFlag;   ///< Includes the vacuum chamber and associated elements (default on).
  G4String insulationTubeIncludedFlag;  ///< Includes the insulation tube and associated elements (default on).
  
  G4int nbSegments;                 ///< Number of scintillator segments
  G4int nbCrystInSegmentRow;        ///< Number of scintillator crystals in a segment row
  G4int nbCrystInSegmentColumn;     ///< Number of scintillator crystals in a segment column
  
  G4double vacuumFlangeSizeX;              ///< Vacuum flange Half Length in X direction
  G4double vacuumFlangeSizeY;              ///< Vacuum flange Half Length in Y direction
  G4double vacuumFlangeSizeZ;              ///< Vacuum flange Half Length in Z direction
  G4double vacuumFlangeThickFrontOfScint;  ///< Vacuum flange thickness in front of the scintillator

  G4double insulationTubeThickness; ///< Insulation tube thickness
    
  G4double sciCrystSizeX;   ///< Half-Size of sci crystal in X
  G4double sciCrystSizeY;   ///< Half-Size of sci crystal in Y
  G4double sciCrystSizeZ;   ///< Half-Size of sci crystal in Z
  
  G4double sciCrystPosX;   ///< Position of sci crystal in X
  G4double sciCrystPosY;   ///< Half-Size of sci crystal in X
  G4double sciCrystPosZ;   ///< Half-Size of sci crystal in X
 
  G4double sciHousSizeX;   ///< Outer dimensions of the housing X 
  G4double sciHousSizeY;   ///< Outer dimensions of the housing Y 
  G4double sciHousSizeZ;   ///< Outer dimensions of the housing Z
  
  G4double sciReflWallThickX;   ///< Thickness of reflector walls X
  G4double sciReflWallThickY;   ///< Thickness of reflector walls Y
  G4double sciReflWindThick;    ///< Thickness of reflector walls

  G4double sciHousWallThickX;  ///< Dimensions of Housing (half-side X)
  G4double sciHousWallThickY;  ///< Dimensions of Housing (half-side Y)
  G4double sciHousWindThick;   ///< Dimensions of Housing 

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
  
  void SetVacuumChamberIncludedFlag(G4String val){vacuumChamberIncludedFlag=val;}
  void SetInsulationTubeIncludedFlag(G4String val){insulationTubeIncludedFlag=val;}

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
  
  G4String GetVacuumChamberIncludedFlag(void){return vacuumChamberIncludedFlag;}
  G4String GetInsulationTubeIncludedFlag(void){return insulationTubeIncludedFlag;}

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
