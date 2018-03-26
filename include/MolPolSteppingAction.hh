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

  private:
    G4bool drawFlag;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };
};

#endif//__MOLPOLSTEPPINGACTION_HH
