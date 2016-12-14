// - AUTHOR: Pablo Cabanelas 12/2016
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/

#ifndef ActarSimExogamSD_h
#define ActarSimExogamSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ActarSimExogamGeantHit.hh"

class G4Step;
class G4HCofThisEvent;

class ActarSimExogamSD : public G4VSensitiveDetector {
private:
  ActarSimExogamGeantHitsCollection* hitsCollection; ///< Geant step-like hits collect.

public:
  ActarSimExogamSD(G4String);
  ~ActarSimExogamSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
};
#endif
