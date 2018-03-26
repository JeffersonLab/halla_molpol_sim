
#ifndef MolPolEventAction_h
#define MolPolEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class MolPolIO;

class MolPolEventAction : public G4UserEventAction
{
  public:
    MolPolEventAction();
    virtual ~MolPolEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void SetIO( MolPolIO *io ){ fIO = io; }

  private:

    MolPolIO *fIO;

  public:
};

#endif

    
