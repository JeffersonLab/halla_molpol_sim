// *************************************************************** (╯°□°）╯︵ ┻━┻
//
//	MolPolEMFieldMessenger.cc
//
//  Updated field messenger for new integrated field types.
//
//
//
//	Eric King - 2018-11-19
//
// *****************************************************************************

#include "MolPolEMFieldMessenger.hh"

#include "MolPolEMFieldSetup.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
MolPolEMFieldMessenger::MolPolEMFieldMessenger(MolPolEMFieldSetup* fieldSetup)
  : G4UImessenger(),
    fEMfieldSetup(fieldSetup),
    fFieldDir(0)
{
  fFieldDir = new G4UIdirectory("/field/");
  fFieldDir->SetGuidance("MolPolEM field tracking control.");

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // These are the commands to set with with field maps
  fQ1MCmd = new G4UIcmdWithAString("/field/setQ1M", this);
  fQ1MCmd->SetGuidance("Set Q1 TOSCA map usage. String [fileLoc] [mapStrength] [desiredStrength]");
  fQ1MCmd->SetParameterName("Q1M", false);

  fQ2MCmd = new G4UIcmdWithAString("/field/setQ2M", this);
  fQ2MCmd->SetGuidance("Set Q2 TOSCA map usage. String [fileLoc] [mapStrength] [desiredStrength]");
  fQ2MCmd->SetParameterName("Q2M", false);

  fQ3MCmd = new G4UIcmdWithAString("/field/setQ3M", this);
  fQ3MCmd->SetGuidance("Set Q3 TOSCA map usage. String [fileLoc] [mapStrength] [desiredStrength]");
  fQ3MCmd->SetParameterName("Q3M", false);

  fQ4MCmd = new G4UIcmdWithAString("/field/setQ4M", this);
  fQ4MCmd->SetGuidance("Set Q4 TOSCA map usage. String [fileLoc] [mapStrength] [desiredStrength]");
  fQ4MCmd->SetParameterName("Q4M", false);

  fDipMCmd = new G4UIcmdWithAString("/field/setDipM", this);
  fDipMCmd->SetGuidance("Set dipole TOSCA map usage. String [fileLoc] [mapStrength] [desiredStrength]");
  fDipMCmd->SetParameterName("DipM", false);

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // These are the commands to set with ideal field by current
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

  fDipACmd = new G4UIcmdWithADouble("/field/setDipA", this);
  fDipACmd->SetGuidance("Set Dipole current");
  fDipACmd->SetParameterName("DipA", false);

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // These are the commands to set with ideal field by pole value.
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

  fDipTCmd = new G4UIcmdWithADouble("/field/setDipT", this);
  fDipTCmd->SetGuidance("Set Dipole field in Tesla");
  fDipTCmd->SetParameterName("DipT", false);

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // This is the new Solenoid Command.  It's map only and you merely set the strength in Tesla.
  fSetSolenoidCmd = new G4UIcmdWithADouble("/field/setSolenoid", this);
  fSetSolenoidCmd->SetGuidance("Set the solenoid value in teslas. Map will scale appropriately.");
  fSetSolenoidCmd->SetParameterName("setSolenoid", false);

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // This is the update command.
  fUpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  fUpdateCmd->SetGuidance("This command MUST be applied after setting field values");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MolPolEMFieldMessenger::~MolPolEMFieldMessenger()
{
  delete fFieldDir;
  delete fUpdateCmd;
  delete fQ1ACmd;
  delete fQ2ACmd;
  delete fQ3ACmd;
  delete fQ4ACmd;
  delete fDipACmd;
  delete fQ1TCmd;
  delete fQ2TCmd;
  delete fQ3TCmd;
  delete fQ4TCmd;
  delete fDipTCmd;
  delete fQ1MCmd;
  delete fQ2MCmd;
  delete fQ3MCmd;
  delete fQ4MCmd;
  delete fDipMCmd;
  delete fSetSolenoidCmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MolPolEMFieldMessenger::SetNewValue( G4UIcommand* cmd, G4String newValue)
{
  if( cmd == fQ1TCmd ){
    G4double x = fQ1TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iQuad1Type = 1;
    fEMfieldSetup->dQuad1RelevantStr = x;
    G4cout << "Set Q1 Type '1' with field " << x << G4endl;
  }else if( cmd == fQ2TCmd ){
    G4double x = fQ2TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iQuad2Type = 1;
    fEMfieldSetup->dQuad2RelevantStr = x;
    G4cout << "Set Q2 Type '1' with field " << x << G4endl;
  }else if( cmd == fQ3TCmd ){
    G4double x = fQ3TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iQuad3Type = 1;
    fEMfieldSetup->dQuad3RelevantStr = x;
    G4cout << "Set Q3 Type '1' with field " << x << G4endl;
  }else if( cmd == fQ4TCmd ){
    G4double x = fQ4TCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iQuad4Type = 1;
    fEMfieldSetup->dQuad4RelevantStr = x;
    G4cout << "Set Q4 Type '1' with field " << x << G4endl;
  }else if( cmd == fDipTCmd ){
    G4double x = fDipTCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iDipoleType = 1;
    fEMfieldSetup->dDipRelevantStr = x;
    G4cout << "Set Dipole Type '1' with field " << x << G4endl;
  } else if( cmd == fQ1MCmd ){
    G4String fileLocation;
    G4double mapStrength;
    G4double desiredStrength;
    if(std::istringstream(newValue) >> fileLocation >> mapStrength >> desiredStrength){
      G4cout << "Set Q1 Type '2' with the following attributes," << G4endl
             << "TOSCA Q1 map location: " << fileLocation    << G4endl
             << "        Map Magnitude: " << mapStrength     << G4endl
             << "    Desired Magnitude: " << desiredStrength << G4endl;
      fEMfieldSetup->iQuad1Type = 2;
      fEMfieldSetup->strQuad1MapLoc = fileLocation;
      fEMfieldSetup->dQuad1ToscaMapStr = mapStrength;
      fEMfieldSetup->dQuad1RelevantStr = desiredStrength;
    } else {
      G4cout << "Incorrect macro entry for setQ1M. Must be string in format [fileLoc] [mapStrength] [desiredStrength]" << G4endl;
    }
  } else if( cmd == fQ2MCmd ){
    G4String fileLocation;
    G4double mapStrength;
    G4double desiredStrength;
    if(std::istringstream(newValue) >> fileLocation >> mapStrength >> desiredStrength){
      G4cout << "Set Q2 Type '2' with the following attributes," << G4endl
             << "TOSCA Q2 map location: " << fileLocation    << G4endl
             << "        Map Magnitude: " << mapStrength     << G4endl
             << "    Desired Magnitude: " << desiredStrength << G4endl;
      fEMfieldSetup->iQuad2Type = 2;
      fEMfieldSetup->strQuad2MapLoc = fileLocation;
      fEMfieldSetup->dQuad2ToscaMapStr = mapStrength;
      fEMfieldSetup->dQuad2RelevantStr = desiredStrength;
    } else {
      G4cout << "Incorrect macro entry for setQ2M. Must be string in format [fileLoc] [mapStrength] [desiredStrength]" << G4endl;
    }
  } else if( cmd == fQ3MCmd ){
    G4String fileLocation;
    G4double mapStrength;
    G4double desiredStrength;
    if(std::istringstream(newValue) >> fileLocation >> mapStrength >> desiredStrength){
      G4cout << "Set Q3 Type '2' with the following attributes," << G4endl
             << "TOSCA Q3 map location: " << fileLocation    << G4endl
             << "        Map Magnitude: " << mapStrength     << G4endl
             << "    Desired Magnitude: " << desiredStrength << G4endl;
      fEMfieldSetup->iQuad3Type = 2;
      fEMfieldSetup->strQuad3MapLoc = fileLocation;
      fEMfieldSetup->dQuad3ToscaMapStr = mapStrength;
      fEMfieldSetup->dQuad3RelevantStr = desiredStrength;
    } else {
      G4cout << "Incorrect macro entry for setQ3M. Must be string in format [fileLoc] [mapStrength] [desiredStrength]" << G4endl;
    }
  } else if( cmd == fQ4MCmd ){
    G4String fileLocation;
    G4double mapStrength;
    G4double desiredStrength;
    if(std::istringstream(newValue) >> fileLocation >> mapStrength >> desiredStrength){
      G4cout << "Set Q4 Type '2' with the following attributes," << G4endl
             << "TOSCA Q4 map location: " << fileLocation    << G4endl
             << "        Map Magnitude: " << mapStrength     << G4endl
             << "    Desired Magnitude: " << desiredStrength << G4endl;
      fEMfieldSetup->iQuad4Type = 2;
      fEMfieldSetup->strQuad4MapLoc = fileLocation;
      fEMfieldSetup->dQuad4ToscaMapStr = mapStrength;
      fEMfieldSetup->dQuad4RelevantStr = desiredStrength;
    } else {
      G4cout << "Incorrect macro entry for setQ4M. Must be string in format [fileLoc] [mapStrength] [desiredStrength]" << G4endl;
    }
  } else if( cmd == fDipMCmd ){
    G4String fileLocation;
    G4double mapStrength;
    G4double desiredStrength;
    if(std::istringstream(newValue) >> fileLocation >> mapStrength >> desiredStrength){
      G4cout << "Set Dipole to Type '2' with the following attributes," << G4endl
             << "TOSCA Dipole Map Location: " << fileLocation    << G4endl
             << "            Map Magnitude: " << mapStrength     << G4endl
             << "        Desired Magnitude: " << desiredStrength << G4endl;
      fEMfieldSetup->iDipoleType = 2;
      fEMfieldSetup->strDipoleMapLoc = fileLocation;
      fEMfieldSetup->dDipToscaMapStr = mapStrength;
      fEMfieldSetup->dDipRelevantStr = desiredStrength;
    } else {
      G4cout << "Incorrect macro entry for setDipM (Dipole). Must be string in format [fileLoc] [mapStrength] [desiredStrength]" << G4endl;
    }
  } else if( cmd == fSetSolenoidCmd){
    G4double x = fSetSolenoidCmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->dSolRelevantStr = x;
    G4cout << "Set dSolRelevantStr to: " << x << G4endl;
  } else if( cmd == fQ1ACmd ){
    G4double x = fQ1ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iQuad1Type = 0;
    fEMfieldSetup->dQuad1RelevantStr = x;
    G4cout << "Set Q1 Type '0' with current " << x << G4endl;
  } else if( cmd == fQ2ACmd ){
    G4double x = fQ2ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iQuad2Type = 0;
    fEMfieldSetup->dQuad2RelevantStr = x;
    G4cout << "Set Q2 Type '0' with current " << x << G4endl;
  } else if( cmd == fQ3ACmd ){
    G4double x = fQ3ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iQuad3Type = 0;
    fEMfieldSetup->dQuad3RelevantStr = x;
    G4cout << "Set Q3 Type '0' with current " << x << G4endl;
  } else if( cmd == fQ4ACmd ){
    G4double x = fQ4ACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iQuad4Type = 0;
    fEMfieldSetup->dQuad4RelevantStr = x;
    G4cout << "Set Q4 Type '0' with current " << x << G4endl;
  } else if( cmd == fDipACmd ){
    G4double x = fDipACmd->GetNewDoubleValue(newValue);
    fEMfieldSetup->iDipoleType = 0;
    fEMfieldSetup->dDipRelevantStr = x;
    G4cout << "Set Dipole Type '0' with current " << x << G4endl;
  } else if( cmd == fUpdateCmd ){
    G4cout << "Updating magnetic field configuration... " << G4endl;
    fEMfieldSetup->UpdateConfiguration();
  } else {
    G4cout << __PRETTY_FUNCTION__ << " at line " << __LINE__<<G4endl;
    G4cerr << "Don't know this command :" << cmd << G4endl;
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
