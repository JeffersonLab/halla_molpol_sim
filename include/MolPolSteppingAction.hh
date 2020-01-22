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

    void    SetStepActKryptEdge(G4bool val);
    void    SetMollerTracksOnly(G4bool val);

  private:
    G4bool  fTrackMollersOnly;
    G4bool  fStepActKryptEdge;

};

#endif//__MOLPOLSTEPPINGACTION_HH
