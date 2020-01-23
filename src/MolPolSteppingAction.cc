#include "MolPolSteppingAction.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"
#include "G4String.hh"
#include <string>

MolPolSteppingAction::MolPolSteppingAction()
:fTrackMollersOnly(true),  //Default option, tracking all will substantially increase file sizes
 fStepActKryptEdge(false)
{

}

void MolPolSteppingAction::UserSteppingAction(const G4Step *aStep) {
  G4Track * aTrack = aStep->GetTrack();

  G4String strPhysVolName = aTrack->GetVolume()->GetName();

  std::size_t foundVB = strPhysVolName.find("virtualBoundaryPhys_det");
  if( foundVB != G4String::npos ){aTrack->SetTrackStatus(fStopAndKill);}

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Condition to track mollers only.
  if(fTrackMollersOnly){
    if( aTrack->GetParentID() > 0 ){ aTrack->SetTrackStatus(fStopAndKill); return; }
  }

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Condition to make kryptonite-like edges without kryptonite.
  if(fStepActKryptEdge){
    if(aTrack->GetMaterial()->GetName() != "MP_Vacuum" &&
       aTrack->GetMaterial()->GetName() != "MP_Air" &&
       aTrack->GetMaterial()->GetName() != "MP_Scint" &&
       aTrack->GetVolume()->GetName() != "Target"  &&
       aTrack->GetVolume()->GetName() != "DipoleExitWindowR" &&
       aTrack->GetVolume()->GetName() != "DipoleExitWindowL")
    {
      aTrack->SetTrackStatus(fStopAndKill);
      return;
    }
  }

}
