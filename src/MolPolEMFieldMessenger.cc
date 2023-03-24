#include "MolPolEMFieldMessenger.hh"

#include "MolPolEMFieldSetup.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

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

  fQ1XposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ1XOffset",this);
  fQ1XposCmd->SetGuidance("Set Q1 X-Position Offset with unit");
  fQ1XposCmd->SetParameterName("Q1XOffset", false);

  fQ2XposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ2XOffset",this);
  fQ2XposCmd->SetGuidance("Set Q2 X-Position Offset with unit");
  fQ2XposCmd->SetParameterName("Q2XOffset", false);

  fQ3XposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ3XOffset",this);
  fQ3XposCmd->SetGuidance("Set Q3 X-Position Offset with unit");
  fQ3XposCmd->SetParameterName("Q3XOffset", false);

  fQ4XposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ4XOffset",this);
  fQ4XposCmd->SetGuidance("Set Q4 X-Position Offset with unit");
  fQ4XposCmd->SetParameterName("Q4XOffset", false);

  fQ5XposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ45Offset",this);
  fQ5XposCmd->SetGuidance("Set Q5 X-Position Offset with unit");
  fQ5XposCmd->SetParameterName("Q5XOffset", false);

  fQ6XposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ6XOffset",this);
  fQ6XposCmd->SetGuidance("Set Q6 X-Position Offset with unit");
  fQ6XposCmd->SetParameterName("Q6XOffset", false);

  fQ1YposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ1YOffset",this);
  fQ1YposCmd->SetGuidance("Set Q1 Y-Position Offset with unit");
  fQ1YposCmd->SetParameterName("Q1YOffset", false);

  fQ2YposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ2YOffset",this);
  fQ2YposCmd->SetGuidance("Set Q2 Y-Position Offset with unit");
  fQ2YposCmd->SetParameterName("Q2YOffset", false);

  fQ3YposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ3YOffset",this);
  fQ3YposCmd->SetGuidance("Set Q3 Y-Position Offset with unit");
  fQ3YposCmd->SetParameterName("Q3YOffset", false);

  fQ4YposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ4YOffset",this);
  fQ4YposCmd->SetGuidance("Set Q4 Y-Position Offset with unit");
  fQ4YposCmd->SetParameterName("Q4YOffset", false);

  fQ5YposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ5YOffset",this);
  fQ5YposCmd->SetGuidance("Set Q5 Y-Position Offset with unit");
  fQ5YposCmd->SetParameterName("Q5YOffset", false);

  fQ6YposCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ6YOffset",this);
  fQ6YposCmd->SetGuidance("Set Q6 Y-Position Offset with unit");
  fQ6YposCmd->SetParameterName("Q6YOffset", false);

  fQ6XrotCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ6XRot",this);
  fQ6XrotCmd->SetGuidance("Set Solenoid X Rotation with unit");
  fQ6XrotCmd->SetParameterName("setQ6XRot", false);

  fQ6YrotCmd = new G4UIcmdWithADoubleAndUnit("/field/setQ6YRot",this);
  fQ6YrotCmd->SetGuidance("Set Solenoid Y Rotation with unit");
  fQ6YrotCmd->SetParameterName("setQ6YRot", false);
  
  // This can probably be done with just a single command since the global field
  // doesn't actually care about this...
  fToscaQ1Cmd = new G4UIcmdWithAString("/field/setToscaQ1",this);
  fToscaQ1Cmd->SetGuidance("Q1 TOSCA Map [relative file location] [TOSCA map pole tip strength] [desired pole tip strength]");
  fToscaQ1Cmd->SetParameterName("ToscaQ1",false);

  fToscaQ2Cmd = new G4UIcmdWithAString("/field/setToscaQ2",this);
  fToscaQ2Cmd->SetGuidance("Q2 TOSCA Map [relative file location] [TOSCA map pole tip strength] [desired pole tip strength]");
  fToscaQ2Cmd->SetParameterName("ToscaQ2",false);

  fToscaQ3Cmd = new G4UIcmdWithAString("/field/setToscaQ3",this);
  fToscaQ3Cmd->SetGuidance("Q3 TOSCA Map [relative file location] [TOSCA map pole tip strength] [desired pole tip strength]");
  fToscaQ3Cmd->SetParameterName("ToscaQ3",false);

  fToscaQ4Cmd = new G4UIcmdWithAString("/field/setToscaQ4",this);
  fToscaQ4Cmd->SetGuidance("Q4 TOSCA Map [relative file location] [TOSCA map pole tip strength] [desired pole tip strength]");
  fToscaQ4Cmd->SetParameterName("ToscaQ4",false);

  fToscaQ5Cmd = new G4UIcmdWithAString("/field/setToscaQ5",this);
  fToscaQ5Cmd->SetGuidance("Q5 TOSCA Map [relative file location] [TOSCA map pole tip strength] [desired pole tip strength]");
  fToscaQ5Cmd->SetParameterName("ToscaQ5",false);

  fToscaQ6Cmd = new G4UIcmdWithAString("/field/setToscaQ6",this);
  fToscaQ6Cmd->SetGuidance("Q6 TOSCA Map [relative file location] [TOSCA map pole tip strength] [desired pole tip strength]");
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
  }else if( cmd == fQ1XposCmd ){
    G4double x = fQ1XposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fXoffsetQ1 = x;
  }else if( cmd == fQ2XposCmd ){
    G4double x = fQ2XposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fXoffsetQ2 = x;
  }else if( cmd == fQ3XposCmd ){
    G4double x = fQ3XposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fXoffsetQ3 = x;
  }else if( cmd == fQ4XposCmd ){
    G4double x = fQ4XposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fXoffsetQ4 = x;
  }else if( cmd == fQ5XposCmd ){
    G4double x = fQ5XposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fXoffsetQ5 = x;
  }else if( cmd == fQ6XposCmd ){
    G4double x = fQ6XposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fXoffsetQ6 = x;
  }else if( cmd == fQ1YposCmd ){
    G4double x = fQ1YposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fYoffsetQ1 = x;
  }else if( cmd == fQ2YposCmd ){
    G4double x = fQ2YposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fYoffsetQ2 = x;
  }else if( cmd == fQ3YposCmd ){
    G4double x = fQ3YposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fYoffsetQ3 = x;
  }else if( cmd == fQ4YposCmd ){
    G4double x = fQ4YposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fYoffsetQ4 = x;
  }else if( cmd == fQ5YposCmd ){
    G4double x = fQ5YposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fYoffsetQ5 = x;
  }else if( cmd == fQ6YposCmd ){
    G4double x = fQ6YposCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fYoffsetQ6 = x;
  }else if( cmd == fQ6XrotCmd ){
    G4double x = fQ6XrotCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fSolRotX = x;
  }else if( cmd == fQ6YrotCmd ){
    G4double x = fQ6YrotCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->fSolRotY = x;
    
  }else if( cmd == fToscaQ1Cmd ){
    G4String fname;
    G4double mapPoleTipStr;
    G4double desPolTipStr;
    if(std::istringstream(newValue) >> fname >> mapPoleTipStr >> desPolTipStr || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q1." << G4endl;
      } else {
        G4cout << "  TOSCA File Q1 name: " << fname << G4endl
               << "    Map Pole Tip (G): " <<  mapPoleTipStr << G4endl
               << "Desired Pole Tip (G): " << desPolTipStr << G4endl;
      }
      fEMfieldSetup->fToscaFields[0] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ1." << G4endl;
    }
  }else if( cmd == fToscaQ2Cmd ){
    G4String fname;
    G4double mapPoleTipStr;
    G4double desPolTipStr;
    if(std::istringstream(newValue) >> fname >> mapPoleTipStr >> desPolTipStr || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q2." << G4endl;
      } else {
        G4cout << "  TOSCA File Q2 name: " << fname << G4endl
               << "    Map Pole Tip (G): " <<  mapPoleTipStr << G4endl
               << "Desired Pole Tip (G): " << desPolTipStr << G4endl;
      }
      fEMfieldSetup->fToscaFields[1] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ2." << G4endl;
    }
  }else if( cmd == fToscaQ3Cmd ){
    G4String fname;
    G4double mapPoleTipStr;
    G4double desPolTipStr;
    if(std::istringstream(newValue) >> fname >> mapPoleTipStr >> desPolTipStr || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q3." << G4endl;
      } else {
        G4cout << "TOSCA File Q3 name: " << fname << G4endl
               << "    Map Pole Tip (G): " <<  mapPoleTipStr << G4endl
               << "Desired Pole Tip (G): " << desPolTipStr << G4endl;
      }
      fEMfieldSetup->fToscaFields[2] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ3." << G4endl;
    }
  }else if( cmd == fToscaQ4Cmd ){
    G4String fname;
    G4double mapPoleTipStr;
    G4double desPolTipStr;
    if(std::istringstream(newValue) >> fname >> mapPoleTipStr >> desPolTipStr || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q4." << G4endl;
      } else {
        G4cout << "TOSCA File Q4 name: " << fname << G4endl
               << "    Map Pole Tip (G): " <<  mapPoleTipStr << G4endl
               << "Desired Pole Tip (G): " << desPolTipStr << G4endl;
      }
      fEMfieldSetup->fToscaFields[3] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ4." << G4endl;
    }
  }else if( cmd == fToscaQ5Cmd ){
    G4String fname;
    G4double mapPoleTipStr;
    G4double desPolTipStr;
    if(std::istringstream(newValue) >> fname >> mapPoleTipStr >> desPolTipStr || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q5." << G4endl;
      } else {
        G4cout << "TOSCA File Q5 name: " << fname << G4endl
               << "    Map Pole Tip (G): " <<  mapPoleTipStr << G4endl
               << "Desired Pole Tip (G): " << desPolTipStr << G4endl;
      }
      fEMfieldSetup->fToscaFields[4] = newValue;
    } else {
      G4cout << "Incorrect macro entry for ToscaQ5." << G4endl;
    }
  }else if( cmd == fToscaQ6Cmd ){
    G4String fname;
    G4double mapPoleTipStr;
    G4double desPolTipStr;
    if(std::istringstream(newValue) >> fname >> mapPoleTipStr >> desPolTipStr || newValue == "none" ){
      if(newValue == "none"){
        G4cout << "No TOSCA File for Q6." << G4endl;
      } else {
        G4cout << "TOSCA File Q6 name: " << fname << G4endl
               << "    Map Pole Tip (G): " <<  mapPoleTipStr << G4endl
               << "Desired Pole Tip (G): " << desPolTipStr << G4endl;
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
