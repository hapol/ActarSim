// - AUTHOR: Pablo Cabanelas 12/2016
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/
//////////////////////////////////////////////////////////////////
/// \class ActarSimExogamHit
/// An Exogam (germanium) hit.
/////////////////////////////////////////////////////////////////

#include "ActarSimExogamHit.hh"
#include "G4ios.hh"
#include "globals.hh"

ClassImp(ActarSimExogamHit)

//////////////////////////////////////////////////////////////////
/// Constructor with initialization to zero
ActarSimExogamHit::ActarSimExogamHit(){
  detectorID = 0;
  xpos=0;
  ypos=0;
  zpos=0;
  time = 0.;
  energy = 0.;
  eBeforeExogam = 0.;
  eAfterExogam = 0.;
  eventID = 0;
  runID = 0;
  trackID=0;
  particleID=0;
  particleCharge=0.;
  particleMass=0.;
  stepsContributing = 0;
}

//////////////////////////////////////////////////////////////////
/// Destructor, nothing to do
ActarSimExogamHit::~ActarSimExogamHit(){
}

//////////////////////////////////////////////////////////////////
/// Printing data information (NOT IMPLEMENTED!)
void ActarSimExogamHit::print(void){
}
