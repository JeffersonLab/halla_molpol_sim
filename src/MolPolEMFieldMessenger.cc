// *************************************************************** (╯°□°）╯︵ ┻━┻
//
//	MolPolEMFieldMessenger.cc
//
//  Updated field messenger for new integrated field types.
//
//  When changes are made:
//    (1) Please be as specific as necessary with SetGuidance()
//    (2) Please update fFieldDir->SetGuidance() line with new name and date.
//    (3) Generate new HTML via G4 command line with: /control/createHTML /field/
//    (4) Push the new field file to the repository.
//
//	Eric King - 2018-11-24
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
  fFieldDir->SetGuidance("MolPolEM field control. ");
  fFieldDir->SetGuidance("Interacts indirectly with MolPolEMField through MolPolEMFieldSetup class.");
  fFieldDir->SetGuidance("Last edited by Eric King on 11-24-2018");

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // These are the commands to set with with field maps
  fQ1MCmd = new G4UIcmdWithAString("/field/setQ1M", this);
  fQ1MCmd->SetGuidance("Sets Q1 TOSCA map usage.");
  fQ1MCmd->SetGuidance("Requires three arguments in string command.");
  fQ1MCmd->SetGuidance("[fileLoc] [mapStrength] [desiredStrength]");
  fQ1MCmd->SetGuidance("Map strength and desired strength can be arbitrarily chosen--center gradient, pole tip, etc. No units required.");
  fQ1MCmd->SetParameterName("Q1M", false);

  fQ2MCmd = new G4UIcmdWithAString("/field/setQ2M", this);
  fQ2MCmd->SetGuidance("Sets Q2 TOSCA map usage.");
  fQ2MCmd->SetGuidance("Requires three arguments in string command.");
  fQ2MCmd->SetGuidance("[fileLoc] [mapStrength] [desiredStrength]");
  fQ2MCmd->SetGuidance("Map strength and desired strength can be arbitrarily chosen--center gradient, pole tip, etc. No units required.");
  fQ2MCmd->SetParameterName("Q2M", false);

  fQ3MCmd = new G4UIcmdWithAString("/field/setQ3M", this);
  fQ3MCmd->SetGuidance("Sets Q3 TOSCA map usage.");
  fQ3MCmd->SetGuidance("Requires three arguments in string command.");
  fQ3MCmd->SetGuidance("[fileLoc] [mapStrength] [desiredStrength]");
  fQ3MCmd->SetGuidance("Map strength and desired strength can be arbitrarily chosen--center gradient, pole tip, etc. No units required.");
  fQ3MCmd->SetParameterName("Q3M", false);

  fQ4MCmd = new G4UIcmdWithAString("/field/setQ4M", this);
  fQ4MCmd->SetGuidance("Sets Q4 TOSCA map usage.");
  fQ4MCmd->SetGuidance("Requires three arguments in string command.");
  fQ4MCmd->SetGuidance("[fileLoc] [mapStrength] [desiredStrength]");
  fQ4MCmd->SetGuidance("Map strength and desired strength can be arbitrarily chosen--center gradient, pole tip, etc. No units required.");
  fQ4MCmd->SetParameterName("Q4M", false);

  fDipMCmd = new G4UIcmdWithAString("/field/setDipM", this);
  fDipMCmd->SetGuidance("Sets dipole TOSCA map usage.");
  fDipMCmd->SetGuidance("Requires three arguments in string command.");
  fDipMCmd->SetGuidance("[fileLoc] [mapStrength] [desiredStrength]");
  fDipMCmd->SetGuidance("Map strength and desired strength can be arbitrarily chosen--center gradient, pole tip, etc. No units required.");
  fDipMCmd->SetParameterName("DipM", false);

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // These are the commands to set with ideal field by current
  fQ1ACmd = new G4UIcmdWithADouble("/field/setQ1A", this);
  fQ1ACmd->SetGuidance("Sets an ideal Q1 field by current.");
  fQ1ACmd->SetGuidance("Current must be in Amps.");
  fQ1ACmd->SetGuidance("Calculates pole tip from Sasha G's GL data fit.");
  fQ1ACmd->SetParameterName("Q1A", false);

  fQ2ACmd = new G4UIcmdWithADouble("/field/setQ2A", this);
  fQ2ACmd->SetGuidance("Sets an ideal Q2 field by current.");
  fQ2ACmd->SetGuidance("Current must be in Amps.");
  fQ2ACmd->SetGuidance("Calculates pole tip from Sasha G's GL data fit.");
  fQ2ACmd->SetParameterName("Q2A", false);

  fQ3ACmd = new G4UIcmdWithADouble("/field/setQ3A", this);
  fQ3ACmd->SetGuidance("Sets an ideal Q3 field by current.");
  fQ3ACmd->SetGuidance("Current must be in Amps.");
  fQ3ACmd->SetGuidance("Calculates pole tip from Sasha G's GL data fit.");
  fQ3ACmd->SetParameterName("Q3A", false);

  fQ4ACmd = new G4UIcmdWithADouble("/field/setQ4A", this);
  fQ4ACmd->SetGuidance("Sets an ideal Q4 field by current.");
  fQ4ACmd->SetGuidance("Current must be in Amps.");
  fQ4ACmd->SetGuidance("Calculates pole tip from Sasha G's GL data fit.");
  fQ4ACmd->SetParameterName("Q4A", false);

  fDipACmd = new G4UIcmdWithADouble("/field/setDipA", this);
  fDipACmd->SetGuidance("Sets an ideal Dipole field by current.");
  fDipACmd->SetGuidance("Current must be in Amps.");
  fDipACmd->SetGuidance("Calculates pole tip from Sasha G's GL data fit.");
  fDipACmd->SetParameterName("DipA", false);

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // These are the commands to set with ideal field by pole value.
  fQ1TCmd = new G4UIcmdWithADouble("/field/setQ1T", this);
  fQ1TCmd->SetGuidance("Sets an ideal field for Q1");
  fQ1TCmd->SetGuidance("This value must be in teslas");
  fQ1TCmd->SetGuidance("Utilizing this command with a value of zero will turn off field.");
  fQ1TCmd->SetParameterName("Q1T", false);

  fQ2TCmd = new G4UIcmdWithADouble("/field/setQ2T", this);
  fQ2TCmd->SetGuidance("Sets an ideal field for Q2" );
  fQ2TCmd->SetGuidance("This value must be in teslas" );
  fQ2TCmd->SetGuidance("Utilizing this command with a value of zero will turn off field.");
  fQ2TCmd->SetParameterName("Q2T", false);

  fQ3TCmd = new G4UIcmdWithADouble("/field/setQ3T", this);
  fQ3TCmd->SetGuidance("Sets an ideal field for Q3");
  fQ3TCmd->SetGuidance("This value must be in teslas");
  fQ3TCmd->SetGuidance("Utilizing this command with a value of zero will turn off field.");
  fQ3TCmd->SetParameterName("Q3T", false);

  fQ4TCmd = new G4UIcmdWithADouble("/field/setQ4T", this);
  fQ4TCmd->SetGuidance("Sets an ideal field for Q4");
  fQ4TCmd->SetGuidance("This value must be in teslas");
  fQ4TCmd->SetGuidance("Utilizing this command with a value of zero will turn off field.");
  fQ4TCmd->SetParameterName("Q4T", false);

  fDipTCmd = new G4UIcmdWithADouble("/field/setDipT", this);
  fDipTCmd->SetGuidance("Sets an ideal Dipole field");
  fDipTCmd->SetGuidance("This value must be in teslas");
  fDipTCmd->SetGuidance("Utilizing this command with a value of zero will turn off field.");
  fDipTCmd->SetParameterName("DipT", false);

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // This is the new Solenoid Command.  It's map only and you merely Sets the strength in Tesla.
  fSetSolenoidCmd = new G4UIcmdWithADouble("/field/setSolenoid", this);
  fSetSolenoidCmd->SetGuidance("Set the solenoid value in teslas; this will allow MolPol to scale solenoid field appropriately.");
  fSetSolenoidCmd->SetGuidance("Setting the value to zero will turn off the field.");
  fSetSolenoidCmd->SetParameterName("setSolenoid", false);

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // This is the update command.
  fUpdateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  fUpdateCmd->SetGuidance("This command MUST be applied after setting field values. Executes UpdateConfiguration() in MolPolEMFieldSetup.");

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
