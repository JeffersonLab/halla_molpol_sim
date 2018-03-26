#ifndef __MOLPOLDETECTOR_HH
#define __MOLPOLDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "MolPolDetectorHit.hh"

#include <map>

/*! 
      Default detector class.  This will record information on:

      - Primary generated hit information

*/

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class MolPolDetector : public G4VSensitiveDetector {
    public:
	MolPolDetector( G4String name, G4int detnum );
	virtual ~MolPolDetector();

	virtual void Initialize(G4HCofThisEvent*);
	virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	virtual void EndOfEvent(G4HCofThisEvent*);

    private:
	MolPolDetectorHitsCollection *fHitColl;
	G4int fHCID;

	G4bool fTrackSecondaries;
	G4int fDetNo;

	G4double fDontRecordThresh;
};

#endif//__MOLPOLDETECTOR_HH
