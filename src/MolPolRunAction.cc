
// Make this appear first!
#include "G4Timer.hh"

#include "MolPolRunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "MolPolIO.hh"

MolPolRunAction::MolPolRunAction()
{
  timer = new G4Timer;
}

MolPolRunAction::~MolPolRunAction()
{
  delete timer;
}

void MolPolRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  //  timer->Start();
  fIO->InitializeTree();

}

void MolPolRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  //       << " " << *timer << G4endl;

  fIO->WriteTree();
}

