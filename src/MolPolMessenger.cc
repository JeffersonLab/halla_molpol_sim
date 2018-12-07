#include "MolPolMessenger.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"

#include "MolPolDetectorConstruction.hh"
#include "MolPolIO.hh"
#include "MolPolEventAction.hh"
#include "MolPolPrimaryGeneratorAction.hh"
#include "MolPolSteppingAction.hh"
//#include "MolPolEMFieldSetup.hh"

#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4RunManager.hh"

#include "G4VPhysicalVolume.hh"

MolPolMessenger::MolPolMessenger(){
    /*  Initialize all the things it talks to to NULL */

    fIO           = NULL;
    fdetcon       = NULL;
    fevact        = NULL;
    fprigen       = NULL;
    fStepAct      = NULL;
    //    fFieldSet     = NULL;

    fMolPolDir = new G4UIdirectory("/MolPol/");
    fMolPolDir->SetGuidance("MolPol UI commands for anything other than field controls.");

    fileCmd = new G4UIcmdWithAString("/MolPol/filename",this);
    fileCmd->SetGuidance("Sets the output ROOT file name.");
    fileCmd->SetParameterName("filename", false);

    seedCmd = new G4UIcmdWithAnInteger("/MolPol/seed",this);
    seedCmd->SetGuidance("Sets a specific seed for the engine.");
    seedCmd->SetParameterName("seed", false);

    genSelectCmd = new G4UIcmdWithAString("/MolPol/gen",this);
    genSelectCmd->SetGuidance("Sets the event generator for the simulation.");
    genSelectCmd->SetGuidance("Setting this to 'moller' will set the Moller Generator (of course).");
    genSelectCmd->SetParameterName("generator", false);

    fXminCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/xmin", this);
    fXminCmd->SetGuidance("Set x range minimum");
    fXminCmd->SetParameterName("xmin", false);

    fXmaxCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/xmax", this);
    fXmaxCmd->SetGuidance("Set x range maximum");
    fXmaxCmd->SetParameterName("xmax", false);

    fYminCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/ymin", this);
    fYminCmd->SetGuidance("Set y range minimum");
    fYminCmd->SetParameterName("ymin", false);

    fYmaxCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/ymax", this);
    fYmaxCmd->SetGuidance("Set y range maximum");
    fYmaxCmd->SetParameterName("ymax", false);

    fBeamECmd = new G4UIcmdWithADoubleAndUnit("/MolPol/beamE", this);
    fBeamECmd->SetGuidance("Sets the beam energy in GeV.");
    fBeamECmd->SetParameterName("beamE", false);

    fEminCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/emin", this);
    fEminCmd->SetGuidance("Set energy range minimum");
    fEminCmd->SetParameterName("emin", false);

    fEmaxCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/emax", this);
    fEmaxCmd->SetGuidance("Set Energy range maximum");
    fEmaxCmd->SetParameterName("emax", false);

    fthetaComMinCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/thcommin", this);
    fthetaComMinCmd->SetGuidance("Sets the theta_com range minimum");
    fthetaComMinCmd->SetParameterName("thcommin", false);

    fthetaComMaxCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/thcommax", this);
    fthetaComMaxCmd->SetGuidance("Sets the theta_com range maximum");
    fthetaComMaxCmd->SetParameterName("thcommax", false);

    fthetaMinCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/thetamin", this);
    fthetaMinCmd->SetGuidance("Set theta range minimum");
    fthetaMinCmd->SetParameterName("thetamin", false);

    fthetaMaxCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/thetamax", this);
    fthetaMaxCmd->SetGuidance("Set theta range maximum");
    fthetaMaxCmd->SetParameterName("thetamax", false);

    fphiMinCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/phimin", this);
    fphiMinCmd->SetGuidance("Set phi range minimum");
    fphiMinCmd->SetParameterName("phimin", false);

    fphiMaxCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/phimax", this);
    fphiMaxCmd->SetGuidance("Set phi range maximum");
    fphiMaxCmd->SetParameterName("phimax", false);

    fLevchukEffectCmd = new G4UIcmdWithABool("/MolPol/calculateLevchuk", this);
    fLevchukEffectCmd->SetGuidance("Set Levchuck Effect On:True Off:False");
    fLevchukEffectCmd->SetParameterName("calculateLevchuk",false);

    fZCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/fz", this);
    fZCmd->SetGuidance("Set particle z");
    fZCmd->SetParameterName("fz", false);

    fRadCorrCmd = new G4UIcmdWithABool("/MolPol/radCorrections",this);
    fRadCorrCmd->SetGuidance("Radiative corrections? True:On False:Off");
    fRadCorrCmd->SetParameterName("radCorrections",false);

    fRemollMSFlagCmd = new G4UIcmdWithABool("/MolPol/remollMS",this);
    fRemollMSFlagCmd->SetGuidance("Remoll Multiple Scattering? True:On False:Off");
    fRemollMSFlagCmd->SetParameterName("remollMS",false);

    fTargPolCmd = new G4UIcmdWithADouble("/MolPol/targetPolPct",this);
    fTargPolCmd->SetGuidance("Target polarization percentage? (Between 0 and 1)");
    fTargPolCmd->SetParameterName("targetPolPct",false);
}

MolPolMessenger::~MolPolMessenger(){
}


void MolPolMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue){
  if( cmd == fTargPolCmd ){
    G4double x = fTargPolCmd->GetNewDoubleValue(newValue);
    if( x >= 0.0 && x <= 1.0) fprigen->fTargPol = x;
    else  G4Exception("MolPolMessenger.cc","",RunMustBeAborted,"targetPolPct set outside of allowable bounds; default value being used.");
  }
  if( cmd == fRadCorrCmd ){
    G4bool flag = fRadCorrCmd->GetNewBoolValue(newValue);
    fprigen->fRadCorrFlag = flag;
  }
  if( cmd == fRemollMSFlagCmd ){
    G4bool flag = fRemollMSFlagCmd->GetNewBoolValue(newValue);
    fprigen->fRemollMSFlag = flag;
  }
  if( cmd == fLevchukEffectCmd ){
    G4bool flag = fLevchukEffectCmd->GetNewBoolValue(newValue);
    fprigen->fLevchukFlag = flag;
  }
  if( cmd == fileCmd ){
    fIO->SetFilename(newValue);
  }
  if( cmd == seedCmd ){
    G4int seed = seedCmd->GetNewIntValue(newValue);
    CLHEP::HepRandom::setTheSeed(seed);
  }

  // POSSCAN
  /*
  if (cmd == fDetPosXCmd ) {
    G4double x = fDetPosXCmd->GetNewDoubleValue(newValue);
    fdetcon->fDetPosX = x;
  }

  if (cmd == fDetPosYCmd ) {
    G4double x = fDetPosYCmd->GetNewDoubleValue(newValue);
    fdetcon->fDetPosY = x;
  }
  */

  if( cmd == genSelectCmd ){
    fprigen->SetGenerator(newValue);
  }
  if( cmd == fXminCmd ){
    G4double x = fXminCmd->GetNewDoubleValue(newValue);
    fprigen->fXmin = x;
  }
  if( cmd == fXmaxCmd ){
    G4double x = fXmaxCmd->GetNewDoubleValue(newValue);
    fprigen->fXmax = x;
  }
  if( cmd == fYminCmd ){
    G4double x = fYminCmd->GetNewDoubleValue(newValue);
    fprigen->fYmin = x;
  }
  if( cmd == fYmaxCmd ){
    G4double x = fYmaxCmd->GetNewDoubleValue(newValue);
    fprigen->fYmax = x;
  }
  if( cmd == fBeamECmd ){
    G4double x = fBeamECmd->GetNewDoubleValue(newValue);
    fprigen->fBeamE = x;
  }
  if( cmd == fEminCmd ){
    G4double x = fEminCmd->GetNewDoubleValue(newValue);
    fprigen->fEmin = x;
  }
  if( cmd == fEmaxCmd ){
    G4double x = fEmaxCmd->GetNewDoubleValue(newValue);
    fprigen->fEmax = x;
  }
  if( cmd == fthetaComMinCmd ){
    G4double x = fthetaComMinCmd->GetNewDoubleValue(newValue);
    fprigen->fthetaComMin = x;
  }
  if( cmd == fthetaComMaxCmd ){
    G4double x = fthetaComMaxCmd->GetNewDoubleValue(newValue);
    fprigen->fthetaComMax = x;
  }
  if( cmd == fthetaMinCmd ){
    G4double x = fthetaMinCmd->GetNewDoubleValue(newValue);
    fprigen->fthetaMin = x;
  }
  if( cmd == fthetaMaxCmd ){
    G4double x = fthetaMaxCmd->GetNewDoubleValue(newValue);
    fprigen->fthetaMax = x;
  }
  if( cmd == fphiMinCmd ){
    G4double x = fphiMinCmd->GetNewDoubleValue(newValue);
    fprigen->fphiMin = x;
  }
  if( cmd == fphiMaxCmd ){
    G4double x = fphiMaxCmd->GetNewDoubleValue(newValue);
    fprigen->fphiMax = x;
  }
  if( cmd == fZCmd ){
    G4double x = fZCmd->GetNewDoubleValue(newValue);
    fprigen->fZ = x;
  }

}
