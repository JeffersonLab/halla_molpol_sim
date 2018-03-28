#ifndef __MOLPOLDETECTORHIT_HH
#define __MOLPOLDETECTORHIT_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class MolPolDetectorHit : public G4VHit {
    public:
	MolPolDetectorHit(G4int, G4int);
	~MolPolDetectorHit();

	MolPolDetectorHit(const MolPolDetectorHit &right);
	const MolPolDetectorHit& operator=(const MolPolDetectorHit &right);
	G4int operator==(const MolPolDetectorHit &right) const;

	inline void *operator new(size_t);
	inline void operator delete(void *aHit);
	void *operator new(size_t,void*p){return p;}

    private:

    public:
	G4int fDetID;
	G4int fCopyID;

	// Position and momentum in lab coordinates
	G4ThreeVector f3X;
	G4ThreeVector f3lX;
	G4ThreeVector f3P;
	// Total momentum, energy, mass
	G4double fP, fE, fM;
	// Origin
	G4ThreeVector f3V;
	G4ThreeVector f3D;
	// Geant4 track ID, particle type, and mother ID
	G4int    fTrID, fPID, fmTrID;
	// Process generator type
	G4int    fGen;
};


typedef G4THitsCollection<MolPolDetectorHit> MolPolDetectorHitsCollection;

extern G4Allocator<MolPolDetectorHit> MolPolDetectorHitAllocator;

inline void* MolPolDetectorHit::operator new(size_t){
    void *aHit;
    aHit = (void *) MolPolDetectorHitAllocator.MallocSingle();
    return aHit;
}

inline void MolPolDetectorHit::operator delete(void *aHit){
    MolPolDetectorHitAllocator.FreeSingle( (MolPolDetectorHit*) aHit);
}

#endif//__MOLPOLDETECTORHIT_HH
