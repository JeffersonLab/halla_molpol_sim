#ifndef MolPolEMFieldMessenger_h
#define MolPolEMFieldMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"

class MolPolEMFieldSetup;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MolPolEMFieldMessenger: public G4UImessenger
{
  public:
    MolPolEMFieldMessenger(MolPolEMFieldSetup* );
    virtual ~MolPolEMFieldMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:

    MolPolEMFieldSetup*        fEMfieldSetup;

    G4UIdirectory*             fFieldDir;

    G4UIcmdWithAnInteger*      fMagSourceCmd;
    G4UIcmdWithADouble*        fQ1ACmd;
    G4UIcmdWithADouble*        fQ2ACmd;
    G4UIcmdWithADouble*        fQ3ACmd;
    G4UIcmdWithADouble*        fQ4ACmd;
    G4UIcmdWithADouble*        fQ5ACmd;
    G4UIcmdWithADouble*        fQ6ACmd;

    G4UIcmdWithADouble*        fQ1TCmd;
    G4UIcmdWithADouble*        fQ2TCmd;
    G4UIcmdWithADouble*        fQ3TCmd;
    G4UIcmdWithADouble*        fQ4TCmd;
    G4UIcmdWithADouble*        fQ5TCmd;
    G4UIcmdWithADouble*        fQ6TCmd;

    G4UIcmdWithADoubleAndUnit* fQ1XposCmd;
    G4UIcmdWithADoubleAndUnit* fQ2XposCmd;
    G4UIcmdWithADoubleAndUnit* fQ3XposCmd;
    G4UIcmdWithADoubleAndUnit* fQ4XposCmd;
    G4UIcmdWithADoubleAndUnit* fQ5XposCmd;
    G4UIcmdWithADoubleAndUnit* fQ6XposCmd;

    G4UIcmdWithADoubleAndUnit* fQ1YposCmd;
    G4UIcmdWithADoubleAndUnit* fQ2YposCmd;
    G4UIcmdWithADoubleAndUnit* fQ3YposCmd;
    G4UIcmdWithADoubleAndUnit* fQ4YposCmd;
    G4UIcmdWithADoubleAndUnit* fQ5YposCmd;
    G4UIcmdWithADoubleAndUnit* fQ6YposCmd;

    G4UIcmdWithADoubleAndUnit* fQ6XrotCmd;
    G4UIcmdWithADoubleAndUnit* fQ6YrotCmd;
    
    G4UIcmdWithAString*        fToscaQ1Cmd;
    G4UIcmdWithAString*        fToscaQ2Cmd;
    G4UIcmdWithAString*        fToscaQ3Cmd;
    G4UIcmdWithAString*        fToscaQ4Cmd;
    G4UIcmdWithAString*        fToscaQ5Cmd;
    G4UIcmdWithAString*        fToscaQ6Cmd;

    G4UIcmdWithoutParameter*   fUpdateCmd;

};

#endif
