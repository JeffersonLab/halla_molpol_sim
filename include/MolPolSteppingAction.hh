#ifndef __MOLPOLSTEPPINGACTION_HH
#define __MOLPOLSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class MolPolSteppingAction : public G4UserSteppingAction
{
  public:
    MolPolSteppingAction();

    virtual ~MolPolSteppingAction(){};
    virtual void UserSteppingAction(const G4Step*);

    void    SetStepActKryptEdge(G4bool val){fStepActKryptEdge = val;G4cout << "fStepActKryptEdge set to: " << fStepActKryptEdge << G4endl;};
    void    SetMollerTracksOnly(G4bool val){fTrackMollersOnly = val;G4cout << "fTrackMollersOnly set to: " << fTrackMollersOnly << G4endl;};

  private:
    G4bool  fTrackMollersOnly;
    G4bool  fStepActKryptEdge;

};

#endif//__MOLPOLSTEPPINGACTION_HH
