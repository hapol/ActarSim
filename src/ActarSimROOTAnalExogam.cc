// - AUTHOR: Hector Alvarez-Pol 05/2008
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/
//////////////////////////////////////////////////////////////////
/// \class ActarSimROOTAnalExogam
/// The ACTAR Scintillator detectorpart of the ROOT Analysis
/////////////////////////////////////////////////////////////////

#include "ActarSimROOTAnalExogam.hh"
#include "ActarSimExogamHit.hh"
#include "ActarSimExogamGeantHit.hh"

#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Track.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4TrackStatus.hh"
#include "G4Step.hh"
#include "G4Types.hh"

//#include "G4PhysicalConstants.hh"
//#include "G4SystemOfUnits.hh"

#include "TROOT.h"
//#include "TApplication.h"
//#include "TSystem.h"
//#include "TH1.h"
//#include "TH2.h"
//#include "TPad.h"
//#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"

//////////////////////////////////////////////////////////////////
/// Constructor
ActarSimROOTAnalExogam::ActarSimROOTAnalExogam() {

  //The simulation file
  simFile = ((ActarSimROOTAnalysis*)gActarSimROOTAnalysis)->GetSimFile();
  simFile->cd();

  //The tree
  eventTree = ((ActarSimROOTAnalysis*)gActarSimROOTAnalysis)->GetEventTree();

  //TODO->Implement SimHits/Hits duality as R3B...

  //The clones array of Sci Hits
  exogamHitCA = new TClonesArray("ActarSimExogamHit",2);

  exogamHitsBranch = eventTree->Branch("exogamHits",&exogamHitCA);
}

//////////////////////////////////////////////////////////////////
/// Destructor. Makes nothing
ActarSimROOTAnalExogam::~ActarSimROOTAnalExogam() {
}

//////////////////////////////////////////////////////////////////
/// Actions to perform in the scintillator anal when generating the primaries
void ActarSimROOTAnalExogam::GeneratePrimaries(const G4Event *anEvent){
  if(anEvent){;} // to quiet the compiler
}

//////////////////////////////////////////////////////////////////
/// Actions to perform in the scintillator anal at the begining of the run
void ActarSimROOTAnalExogam::BeginOfRunAction(const G4Run *aRun) {
  if (aRun){;} /* keep the compiler "quiet" */

  //Storing the runID
  SetTheRunID(aRun->GetRunID());

  char newDirName[255];
  sprintf(newDirName,"%s%i","Histos",aRun->GetRunID());
  simFile->cd(newDirName);

  dirName = new char[255];

  sprintf(dirName,"%s","exogam");
  gDirectory->mkdir(dirName,dirName);
  gDirectory->cd(dirName);

  simFile->cd();
}

//////////////////////////////////////////////////////////////////
/// Actions to perform in the scintillator anal at the begining of the event
void ActarSimROOTAnalExogam::BeginOfEventAction(const G4Event *anEvent) {
  SetTheEventID(anEvent->GetEventID());
}

//////////////////////////////////////////////////////////////////
/// Actions to perform in the scintillator anal at the beginning of the run
void ActarSimROOTAnalExogam::EndOfEventAction(const G4Event *anEvent) {
  FillingHits(anEvent);
}

//////////////////////////////////////////////////////////////////
/// Actions to perform in the scintillator detector analysis after each step
void ActarSimROOTAnalExogam::UserSteppingAction(const G4Step *aStep){
  if(aStep){;} // to quiet the compiler
}

//////////////////////////////////////////////////////////////////
/// Defining the ActarSimExogamHits from the ActarSimExogamGeantHits
void ActarSimROOTAnalExogam::FillingHits(const G4Event *anEvent) {

  //Hit Container ID for ActarSimExogamGeantHit
  G4int hitsCollectionID =
    G4SDManager::GetSDMpointer()->GetCollectionID("ExogamCollection");

  G4HCofThisEvent* HCofEvent = anEvent->GetHCofThisEvent();

  ActarSimExogamGeantHitsCollection* hitsCollection =
    (ActarSimExogamGeantHitsCollection*) HCofEvent->GetHC(hitsCollectionID);

  //Number of R3BCalGeantHit (or steps) in the hitsCollection
  G4int NbHits = hitsCollection->entries();
  G4int NbHitsWithSomeEnergy = NbHits;
  //G4cout << " NbHits: " << NbHits << G4endl;

  //We accept edep=0 GeantHits, we have to remove them
  //from the total number of GeantHits accepted for creating a real CrystalHit
  for (G4int i=0;i<NbHits;i++)
    if((*hitsCollection)[i]->GetEdep()==0.)
      NbHitsWithSomeEnergy--;
  //G4cout << " NbHitsWithSomeEnergy: " << NbHitsWithSomeEnergy << G4endl;

  G4int NbCrystalsWithHit=NbHitsWithSomeEnergy;
  G4String* nameWithSomeEnergy;
  G4int* detIDWithSomeEnergy;
  detIDWithSomeEnergy = new G4int[NbHitsWithSomeEnergy];
  nameWithSomeEnergy = new G4String[NbHitsWithSomeEnergy];

  //keep the crystal identifiers of the GeantHits with some energy
  G4int hitsWithEnergyCounter=0;
  for (G4int i=0;i<NbHits;i++) {
    if((*hitsCollection)[i]->GetEdep()>0.){
      nameWithSomeEnergy[hitsWithEnergyCounter] = (*hitsCollection)[i]->GetDetName();
      detIDWithSomeEnergy[hitsWithEnergyCounter] = (*hitsCollection)[i]->GetDetID();
      //G4cout << "with some energy:   name:" << nameWithSomeEnergy[hitsWithEnergyCounter]
      //     << " detID:"<<detIDWithSomeEnergy[hitsWithEnergyCounter] << G4endl;
      hitsWithEnergyCounter++;
    }
  }

  //if(hitsWithEnergyCounter != NbHitsWithSomeEnergy) {
  //  G4cout << "ERROR in ActarSimROOTAnalExogam::FillingHits(): hitsWithEnergyCounter!=NbHitsWithSomeEnergy" << G4endl;
  //  G4cout << " hitsWithEnergyCounter = " << hitsWithEnergyCounter << " while NbHitsWithSomeEnergy = " << NbHitsWithSomeEnergy << G4endl;
  //}

  //We ask for the number of crystals with signal and
  //this loop calculates them from the number of ActarSimExogamGeantHits with energy
  for (G4int i=0;i<NbHitsWithSomeEnergy;i++) {
    for(G4int j=0;j<i;j++){
      if(nameWithSomeEnergy[i] == nameWithSomeEnergy[j] && detIDWithSomeEnergy[i] == detIDWithSomeEnergy[j]){
	NbCrystalsWithHit--;
	break;  //break stops the second for() sentence as soon as one repetition is found
      }
    }
  }
  //G4cout << " NbCrystalsWithHit: " << NbCrystalsWithHit << G4endl;

  G4String* nameWithHit;
  G4int* detIDWithHit;
  detIDWithHit = new G4int[NbCrystalsWithHit];
  nameWithHit = new G4String[NbCrystalsWithHit];
  hitsWithEnergyCounter=0;

  //keep the crystal identifiers of the final crystalHits
  G4int hitsCounter=0;
  for (G4int i=0;i<NbHits;i++) {
    if ((*hitsCollection)[i]->GetEdep()>0.) {
      if(hitsCounter==0) { //the first geantHit with energy always creates a newHit
	nameWithHit[hitsCounter] = (*hitsCollection)[i]->GetDetName();
	detIDWithHit[hitsCounter] = (*hitsCollection)[i]->GetDetID();
	//G4cout << "with hit:   name:" << nameWithHit[hitsCounter]
	//       << " detID:"<<detIDWithHit[hitsCounter] << G4endl;
	hitsCounter++;
      }
      else {
	for (G4int j=0;j<hitsCounter;j++) {
	  if ( (*hitsCollection)[i]->GetDetName() == nameWithHit[j] &&
	       (*hitsCollection)[i]->GetDetID() == detIDWithHit[j] ) {
	    break;  //break stops the second for() sentence as soon as one repetition is found
	  }
	  if(j==hitsCounter-1){ //when the previous if was never true, there is no repetition
	    nameWithHit[hitsCounter] = (*hitsCollection)[i]->GetDetName();
	    detIDWithHit[hitsCounter] = (*hitsCollection)[i]->GetDetID();
	    //G4cout << "with hit:   name:" << nameWithHit[hitsCounter]
	    //   << " detID:"<<detIDWithHit[hitsCounter] << G4endl;
	    hitsCounter++;
	  }
	}
      }
    }
  }
  //DELETE ME AFTER THE PRELIMINARY TESTS!
  //if(hitsCounter != NbCrystalsWithHit) {
  //  G4cout << "ERROR in R3BROOTAnalCal::FillingHits(): hitsCounter!=NbCrystalsWithHit" << G4endl;
  //  G4cout << " hitsCounter = " << hitsCounter << " while NbCrystalsWithHit = " << NbCrystalsWithHit << G4endl;
  //}

  //Let us create as many R3BCalCrystalHit as NbCrystalsWithHit
  //TODO->Recover here the simhit/hit duality if needed!!
  //if(((ActarSimROOTAnalysis*)gActarSimROOTAnalysis)->GetUseCrystalHitSim()==0){
  theExogamHit = new ActarSimExogamHit*[NbCrystalsWithHit];
  for (G4int i=0;i<NbCrystalsWithHit;i++)
    theExogamHit[i] = new ActarSimExogamHit();
  //}
  //else{  //In case that the simulated hits were used!
  //  theExogamHit = new ActarSimExogamHit*[NbCrystalsWithHit];
  //  for (G4int i=0;i<NbCrystalsWithHit;i++)
  //    theExogamHit[i] = new ActarSimExogamHitSim();
  //}

  //Clear the ClonesArray before filling it
  exogamHitCA->Clear();

  G4bool counted = 0;
  hitsCounter = 0; //now this variable is going to count the added GeantHits
  G4String* name;
  G4int* detID;
  detID = new G4int[NbCrystalsWithHit];
  name = new G4String[NbCrystalsWithHit];

  G4bool shouldThisBeStored = 0;
  //Filling the ActarSimExogamHit from the R3BCalGeantHit (or step)
  for (G4int i=0;i<NbHits;i++) {
    //do not accept GeantHits with edep=0
    //if there is no other GeantHit with edep>0 in the same crystal!
    shouldThisBeStored=0;
    counted =0;
    //G4cout << "ADDING THE HITS. GeantHit with name:" << (*hitsCollection)[i]->GetDetName()
    //   << " detID:"<< (*hitsCollection)[i]->GetDetID() << " edep:"<< (*hitsCollection)[i]->GetEdep() <<" under consideration"<< G4endl;
    for (G4int j=0;j<NbCrystalsWithHit;j++) {
      if( (*hitsCollection)[i]->GetDetName() == nameWithHit[j] &&
	  (*hitsCollection)[i]->GetDetID() == detIDWithHit[j] ) {
	shouldThisBeStored=1;
	break;  //break stops the for() sentence as soon as one "energetic" partner is found
      }
    }
    if(!shouldThisBeStored) continue; //so, continue with next, if there is no "energetic partner"
    if(hitsCounter==0){ //only for the first geantHit accepted for storage
      name[hitsCounter] = (*hitsCollection)[i]->GetDetName();
      detID[hitsCounter] = (*hitsCollection)[i]->GetDetID();
      AddCalCrystalHit(theExogamHit[hitsCounter],(*hitsCollection)[i],0);
      //G4cout << "ADD hit:   name:" << name[hitsCounter]
      //     << " detID:"<<detID[hitsCounter]  << " with code 0 (initial)" << G4endl;
      hitsCounter++;
    }
    else {   // from the second in advance, compare if a R3BCalCrystalHit
             // in the same volume already exists
      for (G4int j=0;j<hitsCounter;j++) {
	if( (*hitsCollection)[i]->GetDetName() == name[j] &&
	    (*hitsCollection)[i]->GetDetID() == detID[j]     ){ //identical Hit already present
	  AddCalCrystalHit(theExogamHit[j],(*hitsCollection)[i],1);
	  //G4cout << "ADD hit:   name:" << name[j]
	  // << " detID:"<< detID[j]  << " with code 1" << G4endl;
	  counted = 1;
	  break;   // stops the for() loop; the info is written in the first identical Hit
	}
      }
      if(counted==0) {	//No identical Hit present.
	name[hitsCounter] = (*hitsCollection)[i]->GetDetName();
	detID[hitsCounter] = (*hitsCollection)[i]->GetDetID();
	AddCalCrystalHit(theExogamHit[hitsCounter],(*hitsCollection)[i],0);
	//G4cout << "ADD hit:   name:" << name[hitsCounter]
	//     << " detID:"<<detID[hitsCounter]  << " with code 0" << G4endl;
	hitsCounter++;
      }
    }
  }

  //BORRAME
  //G4cout << G4endl <<"************************  END OF DEBUGGING *****************************************" << G4endl;
  //G4cout <<"************************************************************************************" << G4endl;

  //at the end, fill the ClonesArray
  //TODO-> Recover here the simhit/hit duality if needed!!
  //if(((ActarSimROOTAnalysis*)gActarSimROOTAnalysis)->GetUseCrystalHitSim()==0){
  for (G4int i=0;i<NbCrystalsWithHit;i++)
    new((*exogamHitCA)[i])ActarSimExogamHit(*theExogamHit[i]);
  //}
  //else{  //In case that the simulated hits were used!
  //  ActarSimExogamHitSim* testSim;
  //  for (G4int i=0;i<NbCrystalsWithHit;i++) {
  //    testSim =  (ActarSimExogamHitSim*) (theExogamHit[i]);
  //    new((*exogamHitCA)[i])ActarSimExogamHitSim(*testSim);
  //  }
  //}

  // all branches at the same time!
  //G4cout << " #@BITCOUNT "<< crystalHitsBranch->Fill() << G4endl;
  //G4cout << " #@BITCOUNT "<< theR3BTree->Fill() << G4endl;

  delete [] detIDWithSomeEnergy;
  delete [] nameWithSomeEnergy;
  delete [] detIDWithHit;
  delete [] nameWithHit;
  delete [] detID;
  delete [] name;
  for (G4int i=0;i<NbCrystalsWithHit;i++) delete theExogamHit[i];
  delete [] theExogamHit;
}

//////////////////////////////////////////////////////////////////
/// Function to move the information from the ActarSimExogamGeantHit (a step hit)
/// to ActarSimExogamHit (an event hit) for the Darmstadt-Heidelberg Crystall Ball.
/// Two modes are possible:
/// - mode == 0 : creation; the ActarSimExogamHit is void and is
///             filled by the data from the ActarSimExogamGeantHit
/// - mode == 1 : addition; the ActarSimExogamHit was already created
///             by other ActarSimExogamGeantHit and some data members are updated
void ActarSimROOTAnalExogam::AddCalCrystalHit(ActarSimExogamHit* cHit,
					   ActarSimExogamGeantHit* gHit,
					   G4int mode) {

  if(mode == 0) { //creation

    G4int copy = 0;

    cHit->SetEnergy(gHit->GetEdep()/ CLHEP::MeV);
    cHit->SetTime(gHit->GetToF() / CLHEP::ns);

    cHit->SetXPos(gHit->GetPos().x()/CLHEP::mm);
    cHit->SetYPos(gHit->GetPos().y()/CLHEP::mm);
    cHit->SetZPos(gHit->GetPos().z()/CLHEP::mm);

    cHit->SetEventID(GetTheEventID());
    cHit->SetRunID(GetTheRunID());

    cHit->SetParticleID(gHit->GetParticleID());
    cHit->SetParticleCharge(gHit->GetParticleCharge());
    cHit->SetParticleMass(gHit->GetParticleMass());

    cHit->SetStepsContributing(1);

    if(cHit->GetXPos()<0 && cHit->GetZPos()<0)                                  copy=0;
    if(cHit->GetXPos()<0 && cHit->GetZPos()>0 && cHit->GetZPos()<20.*CLHEP::cm) copy=1;
    if(cHit->GetXPos()<0 && cHit->GetZPos()>20.*CLHEP::cm)                      copy=2;
    if(cHit->GetXPos()>0 && cHit->GetZPos()>20.*CLHEP::cm)                      copy=3;
    if(cHit->GetXPos()>0 && cHit->GetZPos()>0 && cHit->GetZPos()<20.*CLHEP::cm) copy=4;
    if(cHit->GetXPos()>0 && cHit->GetZPos()<0)                                  copy=5;
    cHit->SetDetectorID(gHit->GetDetID() + 4*copy);

  } else if(mode==1){ //addition

    cHit->SetEnergy(cHit->GetEnergy() + gHit->GetEdep()/ CLHEP::MeV);
    if(gHit->GetToF()<cHit->GetTime()) cHit->SetTime(gHit->GetToF()/ CLHEP::ns);

    cHit->SetStepsContributing(cHit->GetStepsContributing()+1);

    cHit->SetXPos(cHit->GetXPos() +
        ((gHit->GetPos().x()/CLHEP::mm)-cHit->GetXPos())/((G4double)cHit->GetStepsContributing()));
    cHit->SetYPos(cHit->GetYPos() +
        ((gHit->GetPos().y()/CLHEP::mm)-cHit->GetYPos())/((G4double)cHit->GetStepsContributing()));
    cHit->SetZPos(cHit->GetZPos() +
        ((gHit->GetPos().z()/CLHEP::mm)-cHit->GetZPos())/((G4double)cHit->GetStepsContributing()));

  }
}
