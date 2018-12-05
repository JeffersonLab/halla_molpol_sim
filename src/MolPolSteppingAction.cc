#include "MolPolSteppingAction.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"
#include "G4String.hh"
#include <string>

MolPolSteppingAction::MolPolSteppingAction()
:drawFlag(false)
{

}

void MolPolSteppingAction::UserSteppingAction(const G4Step *aStep) {
  G4Track * aTrack = aStep->GetTrack();

  //max step num 1024, if more than this number, do not pass to root tree
  //G4cout << "test" << G4endl;
  /*G4cout<<"Stepping Action:\n\tPDGid: "<<aTrack->GetDefinition()->GetPDGEncoding()
	<<"\tTid: "<<aTrack->GetTrackID()
	<<"\tPos: "<<aStep->GetPostStepPoint()->GetPosition()
	<<"\tMom "<<aStep->GetPostStepPoint()->GetMomentum()<<G4endl;
	std::cin.ignore();*/

  G4String strPhysVolName = aTrack->GetVolume()->GetName();

  std::size_t foundVB = strPhysVolName.find("virtualBoundaryPhys_det");
  if( foundVB != G4String::npos ){
    //G4cout << strPhysVolName << G4endl;
    //G4cout << "det!" << G4endl;
    aTrack->SetTrackStatus(fStopAndKill);
  }

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // These can be removed at any time or perhaps have macro made. Kills all tracks
  // which are greater than 2 and kills in any material other than Vacuum and Target.
  if( aTrack->GetParentID() > 0 ){
      aTrack->SetTrackStatus(fStopAndKill);
      return;
  }
  if( aTrack->GetMaterial()->GetName() != "Vacuum" &&
      aTrack->GetVolume()->GetName() != "Target"  &&
      aTrack->GetVolume()->GetName() != "DipoleExitWindowR" &&
      aTrack->GetVolume()->GetName() != "DipoleExitWindowL")
  {
    aTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////

}
