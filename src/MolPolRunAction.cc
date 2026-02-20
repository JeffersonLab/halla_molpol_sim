
// Make this appear first!
#include "G4Timer.hh"

#include "MolPolRunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "MolPolIO.hh"
#include "MolPolPrimaryGeneratorAction.hh"

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
  timer->Start();
  fIO->InitializeTree();

}

void MolPolRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "Total Number of Events = " << aRun->GetNumberOfEvent()  << G4endl;
  G4cout << "     User Elapsed Time = " << std::fixed << std::setprecision(2) << timer->GetUserElapsed()   << " seconds" << G4endl;
  G4cout << "   + Syst Elapsed Time = " << std::fixed << std::setprecision(2) << timer->GetSystemElapsed() << " seconds" << G4endl;
  G4cout << "   = Real Elapsed Time = " << std::fixed << std::setprecision(2) << timer->GetRealElapsed()   << " seconds" << G4endl;

  fIO->WriteTree();
}

