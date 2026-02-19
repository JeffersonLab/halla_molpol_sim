
#include "MolPolEventAction.hh"
#include "MolPolDetectorHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#include "MolPolIO.hh"
//#include <iostream>
//using namespace std;

MolPolEventAction::MolPolEventAction() {
  fRecordHitsOnly = false;
}

MolPolEventAction::~MolPolEventAction(){
}


void MolPolEventAction::BeginOfEventAction(const G4Event*ev) {
    // Pretty ongoing status with flush
    if( (ev->GetEventID() % 100) == 0 ){
	printf("Event %8d\r", ev->GetEventID() );
	fflush(stdout);
    }

    return;
}

void MolPolEventAction::EndOfEventAction(const G4Event* evt ) {
  //cout << "End of event action" << endl;
  //G4SDManager   *SDman = G4SDManager::GetSDMpointer();
  G4HCofThisEvent *HCE = evt->GetHCofThisEvent();

  G4VHitsCollection *thiscol;

  G4bool hasHits = false;
  G4int totalHits = 0;
  // Count total number o hits in the HC
  for (int hcidx = 0; hcidx < HCE->GetCapacity(); hcidx++){
    thiscol = HCE->GetHC(hcidx);
    if(thiscol){
      totalHits += thiscol->GetSize();
    }
  }
  if(totalHits > 0){
    hasHits = true;
  }

  // Traverse all hit collections, sort by output type
  for( int hcidx = 0; hcidx < HCE->GetCapacity(); hcidx++ ){
      thiscol = HCE->GetHC(hcidx);
      if(thiscol){ // This is NULL if nothing is stored

        //// Detector Hits ///////////////////////////////////
        if( MolPolDetectorHitsCollection *thiscast = 
            dynamic_cast<MolPolDetectorHitsCollection *>(thiscol)){
            for( unsigned int hidx = 0; hidx < thiscast->GetSize(); hidx++ ){
            fIO->AddDetectorHit(
                (MolPolDetectorHit *) thiscast->GetHit(hidx) );
            }
        }
        
      }
  }
  
  //cout << "Flushing" << endl;
  // Fill tree and reset buffers
  if(fRecordHitsOnly){
    // Only record events with hits
    if(hasHits){
      fIO->FillTree();
    }
  } else {
    // Record all events
    fIO->FillTree();
  }

  fIO->Flush();

  return;
}



