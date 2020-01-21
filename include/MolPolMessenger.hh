#ifndef MolPolMessenger_HH
#define MolPolMessenger_HH

#include "globals.hh"
#include "MolPoltypes.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4VModularPhysicsList.hh"

/*!
 *   Global messenger class
 */

class MolPolIO;
class MolPolDetectorConstruction;
class MolPolEventAction;
class MolPolPrimaryGeneratorAction;
class MolPolSteppingAction;
//class MolPolEMFieldSetup;

class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIdirectory;

class MolPolMessenger : public G4UImessenger {
    public:
       	MolPolMessenger();
       	~MolPolMessenger();

	void SetIO( MolPolIO *io ){ fIO = io; }
	void SetPriGen( MolPolPrimaryGeneratorAction *pg ){ fprigen = pg; }
	void SetDetCon( MolPolDetectorConstruction *dc ){ fdetcon= dc; }
	void SetEvAct( MolPolEventAction *ev ){ fevact = ev; }
	void SetStepAct( MolPolSteppingAction *st ){ fStepAct = st; }
  //    void SetFieldSet( MolPolEMFieldSetup* fs ){ fFieldSet = fs; }

	void SetNewValue(G4UIcommand* cmd, G4String newValue);

    private:
	MolPolIO *fIO;
	MolPolDetectorConstruction *fdetcon;
	MolPolEventAction *fevact;
	MolPolPrimaryGeneratorAction *fprigen;
	MolPolSteppingAction *fStepAct;

        G4UIdirectory *fMolPolDir;

        //added in targpolcmd for testing/validation to avoid recompiling every time needed to change.
        G4UIcmdWithABool         *fLevchukEffectCmd;
        G4UIcmdWithADouble       *fTargPolCmd;
        G4UIcmdWithABool         *fRadCorrCmd;
        G4UIcmdWithABool         *fRemollMSFlagCmd;

	G4UIcmdWithAnInteger *seedCmd;
	G4UIcmdWithAString   *fileCmd;
        G4UIcmdWithAString   *genSelectCmd;

	G4UIcmdWithADoubleAndUnit *fXminCmd;
	G4UIcmdWithADoubleAndUnit *fXmaxCmd;
	G4UIcmdWithADoubleAndUnit *fYminCmd;
	G4UIcmdWithADoubleAndUnit *fYmaxCmd;

        // beam profile
	G4UIcmdWithADoubleAndUnit *fXsmearCmd;
	G4UIcmdWithADoubleAndUnit *fYsmearCmd;

	G4UIcmdWithADoubleAndUnit *fBeamECmd;
	G4UIcmdWithADoubleAndUnit *fEminCmd;
	G4UIcmdWithADoubleAndUnit *fEmaxCmd;

	G4UIcmdWithADoubleAndUnit *fthetaComMinCmd;
	G4UIcmdWithADoubleAndUnit *fthetaComMaxCmd;

	G4UIcmdWithADoubleAndUnit *fthetaMinCmd;
	G4UIcmdWithADoubleAndUnit *fthetaMaxCmd;
	G4UIcmdWithADoubleAndUnit *fphiMinCmd;
	G4UIcmdWithADoubleAndUnit *fphiMaxCmd;

	G4UIcmdWithADoubleAndUnit *fXCmd;
	G4UIcmdWithADoubleAndUnit *fYCmd;
	G4UIcmdWithADoubleAndUnit *fZCmd;

    G4UIcmdWithADoubleAndUnit *fBeamRotZXCmd;
    G4UIcmdWithADoubleAndUnit *fBeamRotZYCmd;

};

#endif//MolPolMessenger_HH
