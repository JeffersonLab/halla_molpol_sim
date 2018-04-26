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
    fMolPolDir->SetGuidance("UI commands of this code");

    fileCmd = new G4UIcmdWithAString("/MolPol/filename",this);
    fileCmd->SetGuidance("Output filename");
    fileCmd->SetParameterName("filename", false);

    seedCmd = new G4UIcmdWithAnInteger("/MolPol/seed",this);
    seedCmd->SetGuidance("Set random engine seed");
    seedCmd->SetParameterName("seed", false);

    genSelectCmd = new G4UIcmdWithAString("/MolPol/gen",this);
    genSelectCmd->SetGuidance("Select generator");
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
    fBeamECmd->SetGuidance("Set beam energy");
    fBeamECmd->SetParameterName("beamE", false);

    fEminCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/emin", this);
    fEminCmd->SetGuidance("Set energy range minimum");
    fEminCmd->SetParameterName("emin", false);

    fEmaxCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/emax", this);
    fEmaxCmd->SetGuidance("Set Energy range maximum");
    fEmaxCmd->SetParameterName("emax", false);

    fthetaComMinCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/thcommin", this);
    fthetaComMinCmd->SetGuidance("Set thcom range minimum");
    fthetaComMinCmd->SetParameterName("thcommin", false);

    fthetaComMaxCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/thcommax", this);
    fthetaComMaxCmd->SetGuidance("Set thcom range maximum");
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

    fBeamRadCorrCmd = new G4UIcmdWithABool("/MolPol/beamRadCorrections",this);
    fBeamRadCorrCmd->SetGuidance("Radiative corrections for beam? True:On False:Off");
    fBeamRadCorrCmd->SetParameterName("beamRadCorrections",false);

    fElectronsRadCorrCmd = new G4UIcmdWithABool("/MolPol/electronsRadCorrections",this);
    fElectronsRadCorrCmd->SetGuidance("Radiative corrections for electrons? True:On False:Off");
    fElectronsRadCorrCmd->SetParameterName("electronsRadCorrections",false);

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
  if( cmd == fBeamRadCorrCmd ){
    G4bool flag = fBeamRadCorrCmd->GetNewBoolValue(newValue);
    fprigen->fBeamRadCorrFlag = flag;
  }
  if( cmd == fElectronsRadCorrCmd ){
    G4bool flag = fElectronsRadCorrCmd->GetNewBoolValue(newValue);
    fprigen->fElectronsRadCorrFlag = flag;
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
  /*
  if( cmd == fMagSourceCmd ){
    G4double x = fMagSourceCmd->GetNewIntValue(newValue);
    fFieldSet->fMagSourceMode = x;
  }
  if( cmd == fQ1ACmd ){
    G4double x = fQ1ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ1A = x;
  }
  if( cmd == fQ2ACmd ){
    G4double x = fQ2ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ2A = x;
  }
  if( cmd == fQ3ACmd ){
    G4double x = fQ3ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ3A = x;
  }
  if( cmd == fQ4ACmd ){
    G4double x = fQ4ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ4A = x;
  }
  if( cmd == fQ5ACmd ){
    G4double x = fQ5ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ5A = x;
  }
  if( cmd == fQ1ACmd ){
    G4double x = fQ1ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ1T = x;
  }
  if( cmd == fQ2ACmd ){
    G4double x = fQ2ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ2T = x;
  }
  if( cmd == fQ3ACmd ){
    G4double x = fQ3ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ3T = x;
  }
  if( cmd == fQ4ACmd ){
    G4double x = fQ4ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ4T = x;
  }
  if( cmd == fQ5ACmd ){
    G4double x = fQ5ACmd->GetNewDoubleValue(newValue);
    fFieldSet->fQ5T = x;
  }
  */
  if( cmd == fZCmd ){
    G4double x = fZCmd->GetNewDoubleValue(newValue);
    fprigen->fZ = x;
  }
    
}
