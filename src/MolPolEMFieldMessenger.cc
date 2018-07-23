#include "MolPolEMFieldMessenger.hh"

#include "MolPolEMFieldSetup.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
MolPolEMFieldMessenger::MolPolEMFieldMessenger(MolPolEMFieldSetup* fieldSetup)
  : G4UImessenger(),
    fEMfieldSetup(fieldSetup),
    fFieldDir(0)
{
  fFieldDir = new G4UIdirectory("/field/");
  fFieldDir->SetGuidance("MolPolEM field tracking control.");

  fMagSourceCmd = new G4UIcmdWithAnInteger("/field/MagSourceMode", this);
  fMagSourceCmd->SetGuidance("Set source mode for mag field setting");
  fMagSourceCmd->SetParameterName("magsource", false);
  fMagSourceCmd->SetDefaultValue(1);

  fQ1ACmd = new G4UIcmdWithADouble("/field/setQ1A", this);
  fQ1ACmd->SetGuidance("Set Q1 current");
  fQ1ACmd->SetParameterName("Q1A", false);

  fQ2ACmd = new G4UIcmdWithADouble("/field/setQ2A", this);
  fQ2ACmd->SetGuidance("Set Q2 current");
  fQ2ACmd->SetParameterName("Q2A", false);

  fQ3ACmd = new G4UIcmdWithADouble("/field/setQ3A", this);
  fQ3ACmd->SetGuidance("Set Q3 current");
  fQ3ACmd->SetParameterName("Q3A", false);

  fQ4ACmd = new G4UIcmdWithADouble("/field/setQ4A", this);
  fQ4ACmd->SetGuidance("Set Q4 current");
  fQ4ACmd->SetParameterName("Q4A", false);

  fQ5ACmd = new G4UIcmdWithADouble("/field/setQ5A", this);
  fQ5ACmd->SetGuidance("Set Dipole current");
  fQ5ACmd->SetParameterName("Q5A", false);

  fQ6ACmd = new G4UIcmdWithADouble("/field/setQ6A", this);
  fQ6ACmd->SetGuidance("Set Dipole current");
  fQ6ACmd->SetParameterName("Q6A", false);

  fQ1TCmd = new G4UIcmdWithADouble("/field/setQ1T", this);
  fQ1TCmd->SetGuidance("Set Q1 field in Tesla");
  fQ1TCmd->SetParameterName("Q1T", false);

  fQ2TCmd = new G4UIcmdWithADouble("/field/setQ2T", this);
  fQ2TCmd->SetGuidance("Set Q2 field in Tesla" );
  fQ2TCmd->SetParameterName("Q2T", false);

  fQ3TCmd = new G4UIcmdWithADouble("/field/setQ3T", this);
  fQ3TCmd->SetGuidance("Set Q3 field in Tesla");
  fQ3TCmd->SetParameterName("Q3T", false);

  fQ4TCmd = new G4UIcmdWithADouble("/field/setQ4T", this);
  fQ4TCmd->SetGuidance("Set Q4 field in Tesla");
  fQ4TCmd->SetParameterName("Q4T", false);

  fQ5TCmd = new G4UIcmdWithADouble("/field/setQ5T", this);
  fQ5TCmd->SetGuidance("Set Dipole field in Tesla");
  fQ5TCmd->SetParameterName("Q5T", false);

  fQ6TCmd = new G4UIcmdWithADouble("/field/setQ6T", this);
  fQ6TCmd->SetGuidance("Set Holding field in Tesla");
  fQ6TCmd->SetParameterName("Q6T", false);

  // This can probably be done with just a single command since the global field
  // doesn't actually care about this...
  fToscaQ1Cmd = new G4UIcmdWithAString("/field/setToscaQ1",this);
  fToscaQ1Cmd->SetGuidance("Q1 TOSCA Map [relative file location] [field scale] [beamline offset of field in centimeters]");
  fToscaQ1Cmd->SetParameterName("ToscaQ1",false);

  fToscaQ2Cmd = new G4UIcmdWithAString("/field/setToscaQ2",this);
  fToscaQ2Cmd->SetGuidance("Q2 TOSCA Map [relative file location] [field scale] [beamline offset of field in centimeters]");
  fToscaQ2Cmd->SetParameterName("ToscaQ2",false);

  fToscaQ3Cmd = new G4UIcmdWithAString("/field/setToscaQ3",this);
  fToscaQ3Cmd->SetGuidance("Q3 TOSCA Map [relative file location] [field scale] [beamline offset of field in centimeters]");
  fToscaQ3Cmd->SetParameterName("ToscaQ3",false);

  fToscaQ4Cmd = new G4UIcmdWithAString("/field/setToscaQ4",this);
  fToscaQ4Cmd->SetGuidance("Q4 TOSCA Map [relative file location] [field scale] [beamline offset of field in centimeters]");
  fToscaQ4Cmd->SetParameterName("ToscaQ4",false);

  fToscaQ5Cmd = new G4UIcmdWithAString("/field/setToscaQ5",this);
  fToscaQ5Cmd->SetGuidance("Q5 TOSCA Map [relative file location] [field scale] [beamline offset of field in centimeters]");
  fToscaQ5Cmd->SetParameterName("ToscaQ5",false);

  fToscaQ6Cmd = new G4UIcmdWithAString("/field/setToscaQ6",this);
  fToscaQ6Cmd->SetGuidance("Q6 TOSCA Map [relative file location] [field scale] [beamline offset of field in centimeters]");
  fToscaQ6Cmd->SetParameterName("ToscaQ6",false);


  fUpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  fUpdateCmd->SetGuidance("This command MUST be applied after setting field values ");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MolPolEMFieldMessenger::~MolPolEMFieldMessenger()
{
  delete fFieldDir;
  delete fUpdateCmd;
  delete fMagSourceCmd;
  delete fQ1ACmd;
  delete fQ2ACmd;
  delete fQ3ACmd;
  delete fQ4ACmd;
  delete fQ5ACmd;
  delete fQ1TCmd;
  delete fQ2TCmd;
  delete fQ3TCmd;
  delete fQ4TCmd;
  delete fQ5TCmd;
  delete fToscaQ1Cmd;
  delete fToscaQ2Cmd;
  delete fToscaQ3Cmd;
  delete fToscaQ4Cmd;
  delete fToscaQ5Cmd;
  delete fToscaQ6Cmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MolPolEMFieldMessenger::SetNewValue( G4UIcommand* cmd, G4String newValue)
{

  if( cmd == fMagSourceCmd ){
    G4double x = fMagSourceCmd->GetNewIntValue(newValue);
    fEMfieldSetup->fMagSourceMode = x;
  }else if( cmd == fQ1ACmd ){
    G4double x = fQ1ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ1A = x;
  }else if( cmd == fQ2ACmd ){
    G4double x = fQ2ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ2A = x;
  }else if( cmd == fQ3ACmd ){
    G4double x = fQ3ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ3A = x;
  }else if( cmd == fQ4ACmd ){
    G4double x = fQ4ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ4A = x;
  }else if( cmd == fQ5ACmd ){
    G4double x = fQ5ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ5A = x;
  }else if( cmd == fQ6ACmd ){
    G4double x = fQ6ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ6A = x;
  }else if( cmd == fQ1TCmd ){
    G4double x = fQ1TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ1T = x;
  }else if( cmd == fQ2TCmd ){
    G4double x = fQ2TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ2T = x;
  }else if( cmd == fQ3TCmd ){
    G4double x = fQ3TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ3T = x;
  }else if( cmd == fQ4TCmd ){
    G4double x = fQ4TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ4T = x;
  }else if( cmd == fQ5TCmd ){
    G4double x = fQ5TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ5T = x;
  }else if( cmd == fQ6TCmd ){
    G4double x = fQ6TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fQ6T = x;
  }else if( cmd == fToscaQ1Cmd ){
    G4String fname;
    G4double scale;
    G4double offset;
    if(std::istringstream(newValue) >> fname >> scale >> offset || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q1." << G4endl;
      } else {
        G4cout << "TOSCA File Q1 name: " << fname << G4endl
               << "Scale: " << scale << G4endl
               << "Offset: " << offset << G4endl;
      }
      fEMfieldSetup->fToscaFields[0] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ1." << G4endl;
    }
  }else if( cmd == fToscaQ2Cmd ){
    G4String fname;
    G4double scale;
    G4double offset;
    if(std::istringstream(newValue) >> fname >> scale >> offset || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q2." << G4endl;
      } else {
        G4cout << "TOSCA File Q2 name: " << fname << G4endl
               << "Scale: " << scale << G4endl
               << "Offset: " << offset << G4endl;
      }
      fEMfieldSetup->fToscaFields[1] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ2y." << G4endl;
    }
  }else if( cmd == fToscaQ3Cmd ){
    G4String fname;
    G4double scale;
    G4double offset;
    if(std::istringstream(newValue) >> fname >> scale >> offset || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q3." << G4endl;
      } else {
        G4cout << "TOSCA File Q3 name: " << fname << G4endl
               << "Scale: " << scale << G4endl
               << "Offset: " << offset << G4endl;
      }
      fEMfieldSetup->fToscaFields[2] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ3." << G4endl;
    }
  }else if( cmd == fToscaQ4Cmd ){
    G4String fname;
    G4double scale;
    G4double offset;
    if(std::istringstream(newValue) >> fname >> scale >> offset || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q4." << G4endl;
      } else {
        G4cout << "TOSCA File Q4 name: " << fname << G4endl
               << "Scale: " << scale << G4endl
               << "Offset: " << offset << G4endl;
      }
      fEMfieldSetup->fToscaFields[3] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ4." << G4endl;
    }
  }else if( cmd == fToscaQ5Cmd ){
    G4String fname;
    G4double scale;
    G4double offset;
    if(std::istringstream(newValue) >> fname >> scale >> offset || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q5." << G4endl;
      } else {
        G4cout << "TOSCA File Q5 name: " << fname << G4endl
               << "Scale: " << scale << G4endl
               << "Offset: " << offset << G4endl;
      }
      fEMfieldSetup->fToscaFields[4] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ5." << G4endl;
    }
  }else if( cmd == fToscaQ6Cmd ){
    G4String fname;
    G4double scale;
    G4double offset;
    if(std::istringstream(newValue) >> fname >> scale >> offset || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q6." << G4endl;
      } else {
        G4cout << "TOSCA File Q6 name: " << fname << G4endl
               << "Scale: " << scale << G4endl
               << "Offset: " << offset << G4endl;
      }
      fEMfieldSetup->fToscaFields[5] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ6." << G4endl;
    }
  }else if( cmd == fUpdateCmd ){
    G4cout << "Updating magnetic field configuration... " << G4endl;
    fEMfieldSetup->UpdateConfiguration();
  }else{
    G4cout<<__PRETTY_FUNCTION__<<" at line "<<__LINE__<<G4endl;
    G4cerr <<"Don't know this command :"<<cmd<<G4endl;
  }

  /*
    if( command == fStepperCmd )
    fEMfieldSetup->SetStepperType(fStepperCmd->GetNewIntValue(newValue));
    if( command == fUpdateCmd )
    fEMfieldSetup->CreateStepperAndChordFinder();
    if( command == fMagFieldCmd )
    fEMfieldSetup->SetFieldValue(fMagFieldCmd->GetNewDoubleValue(newValue));
    if( command == fMinStepCmd )
    fEMfieldSetup->SetMinStep(fMinStepCmd->GetNewDoubleValue(newValue));
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
