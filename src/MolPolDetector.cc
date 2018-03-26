#include "MolPolDetector.hh"
#include "G4SDManager.hh"

MolPolDetector::MolPolDetector( G4String name, G4int detnum ) : G4VSensitiveDetector(name){
    char colname[255];

    fDetNo = detnum;
    assert( fDetNo > 0 );

    fHCID = -1;

    fTrackSecondaries = false;
    //    fTrackSecondaries = true;

    sprintf(colname, "genhit_%s_%d", name.data(), detnum);
    collectionName.insert(G4String(colname));

}

MolPolDetector::~MolPolDetector(){
}

void MolPolDetector::Initialize(G4HCofThisEvent *){

    fHitColl = new MolPolDetectorHitsCollection( SensitiveDetectorName, collectionName[0] );
}

///////////////////////////////////////////////////////////////////////

G4bool MolPolDetector::ProcessHits( G4Step *step, G4TouchableHistory *){
    G4bool badedep = false;
    G4bool badhit  = false;

    // Get touchable volume info
    G4TouchableHistory *hist = 
	(G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
    G4int  copyID = hist->GetReplicaNumber();

    G4StepPoint *prestep = step->GetPreStepPoint();
    G4Track     *track   = step->GetTrack();

    G4TouchableHandle theTouchable = prestep->GetTouchableHandle();

//    G4Material* material = track->GetMaterial();

//    printf("Standard detector %d hit by %s!\n", fDetNo, track->GetParticleDefinition()->GetParticleName().data());

//    G4double edep = step->GetTotalEnergyDeposit();

    //  Make pointer to new hit if it's a valid track
    MolPolDetectorHit *thishit;
    if( !badhit ){
	thishit = new MolPolDetectorHit(fDetNo, copyID);
	fHitColl->insert( thishit );
    }

    if( !badhit ){
	// Hit
	thishit->f3X = prestep->GetPosition();
	thishit->f3lX = theTouchable->GetHistory()->GetTopTransform().TransformPoint(prestep->GetPosition());
	thishit->f3V = track->GetVertexPosition();
	thishit->f3D = track->GetVertexMomentumDirection();
	thishit->f3P = track->GetMomentum();

	thishit->fP = track->GetMomentum().mag();
	thishit->fE = track->GetTotalEnergy();
	thishit->fM = track->GetDefinition()->GetPDGMass();

	thishit->fTrID  = track->GetTrackID();
	thishit->fmTrID = track->GetParentID();
	thishit->fPID   = track->GetDefinition()->GetPDGEncoding();

	// FIXME - Enumerate encodings
	thishit->fGen   = (long int) track->GetCreatorProcess();
    }

    if( fDetNo == 9 )
      track->SetTrackStatus(fStopAndKill);

    return !badedep && !badhit;
}

///////////////////////////////////////////////////////////////////////

void MolPolDetector::EndOfEvent(G4HCofThisEvent*HCE) {
    G4SDManager *sdman = G4SDManager::GetSDMpointer();

    if(fHCID<0){ fHCID = sdman->GetCollectionID(collectionName[0]); }

    HCE->AddHitsCollection( fHCID, fHitColl );
    return;
}


