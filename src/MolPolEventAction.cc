
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

  
  // Traverse all hit collections, sort by output type
  for( int hcidx = 0; hcidx < HCE->GetCapacity(); hcidx++ ){
      thiscol = HCE->GetHC(hcidx);
      if(thiscol){ // This is NULL if nothing is stored
	  // Dyanmic cast to test types, process however see fit and feed to IO
	  
	  ////  Detector Hits ///////////////////////////////////
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
  fIO->FillTree();
  fIO->Flush();

  return;
}



