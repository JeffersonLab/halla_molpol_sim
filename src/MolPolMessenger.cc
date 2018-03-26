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

    /*
    fMagSourceCmd = new G4UIcmdWithAnInteger("/MolPol/MagSourceMode", this);
    fMagSourceCmd->SetGuidance("Set source mode for mag field setting");
    fMagSourceCmd->SetParameterName("magsource", false);

    fQ1ACmd = new G4UIcmdWithADouble("/MolPol/Q1A", this);
    fQ1ACmd->SetGuidance("Set Q1 current");
    fQ1ACmd->SetParameterName("Q1A", false);

    fQ2ACmd = new G4UIcmdWithADouble("/MolPol/Q2A", this);
    fQ2ACmd->SetGuidance("Set Q2 current");
    fQ2ACmd->SetParameterName("Q2A", false);

    fQ3ACmd = new G4UIcmdWithADouble("/MolPol/Q3A", this);
    fQ3ACmd->SetGuidance("Set Q3 current");
    fQ3ACmd->SetParameterName("Q3A", false);

    fQ4ACmd = new G4UIcmdWithADouble("/MolPol/Q4A", this);
    fQ4ACmd->SetGuidance("Set Q4 current");
    fQ4ACmd->SetParameterName("Q4A", false);

    fQ5ACmd = new G4UIcmdWithADouble("/MolPol/Q5A", this);
    fQ5ACmd->SetGuidance("Set Dipole current");
    fQ5ACmd->SetParameterName("Q5A", false);

    fQ1TCmd = new G4UIcmdWithADouble("/MolPol/Q1T", this);
    fQ1TCmd->SetGuidance("Set Q1 field in Tesla");
    fQ1TCmd->SetParameterName("Q1T", false);

    fQ2TCmd = new G4UIcmdWithADouble("/MolPol/Q2T", this);
    fQ2TCmd->SetGuidance("Set Q2 field in Tesla" );
    fQ2TCmd->SetParameterName("Q2T", false);

    fQ3TCmd = new G4UIcmdWithADouble("/MolPol/Q3T", this);
    fQ3TCmd->SetGuidance("Set Q3 field in Tesla");
    fQ3TCmd->SetParameterName("Q3T", false);

    fQ4TCmd = new G4UIcmdWithADouble("/MolPol/Q4T", this);
    fQ4TCmd->SetGuidance("Set Q4 field in Tesla");
    fQ4TCmd->SetParameterName("Q4T", false);

    fQ5TCmd = new G4UIcmdWithADouble("/MolPol/Q5T", this);
    fQ5TCmd->SetGuidance("Set Dipole field in Tesla");
    fQ5TCmd->SetParameterName("Q5T", false);
    */

    fZCmd = new G4UIcmdWithADoubleAndUnit("/MolPol/fz", this);
    fZCmd->SetGuidance("Set particle z");
    fZCmd->SetParameterName("fz", false);

}

MolPolMessenger::~MolPolMessenger(){
}


void MolPolMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue){
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
