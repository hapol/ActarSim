// - AUTHOR: Pablo Cabanelas 12/2016
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/

#ifndef ActarSimExogamDetectorConstruction_h
#define ActarSimExogamDetectorConstruction_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
//class ActarSimExogamDetectorMessenger;
class ActarSimDetectorConstruction;

class ActarSimExogamDetectorConstruction {
private:
  G4Material* exogamBulkMaterial;                   ///< Materials
  //ActarSimExogamDetectorMessenger* exogamMessenger;    //pointer to the Messenger
  ActarSimDetectorConstruction* detConstruction; ///< pointer to the global detector

  G4VPhysicalVolume* ConstructExogam(G4LogicalVolume*);

public:
  ActarSimExogamDetectorConstruction(ActarSimDetectorConstruction*);
  ~ActarSimExogamDetectorConstruction();

  G4VPhysicalVolume* Construct(G4LogicalVolume*);

  void SetExogamBulkMaterial (G4String);

  G4Material* GetExogamBulkMaterial() {return exogamBulkMaterial;}

  void UpdateGeometry();
  void PrintDetectorParameters();
};
#endif
