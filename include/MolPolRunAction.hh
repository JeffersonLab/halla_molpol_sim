
#ifndef MolPolRunAction_h
#define MolPolRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;
class MolPolIO;

class MolPolRunAction : public G4UserRunAction
{
  public:
    MolPolRunAction();
    ~MolPolRunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

    void SetIO( MolPolIO *io ){ fIO = io; }

  private:
    G4Timer* timer;

    MolPolIO *fIO;
};

#endif

