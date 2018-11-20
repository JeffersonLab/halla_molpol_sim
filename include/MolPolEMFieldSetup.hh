#ifndef MolPolEMFieldSetup_H
#define MolPolEMFieldSetup_H 1

#include "MolPolEMField.hh"
#include "MolPolQuad.hh"
#include "MolPolDipole.hh"
#include "MolPolSolenoid.hh"

class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_UsualEqRhs;
class G4Mag_EqRhs;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
class G4MagInt_Driver;
class G4UniformMagField;
class MolPolEMFieldMessenger;

class MolPolEMFieldSetup
{
  //public:
  //Static method which returns the singleton pointer of this class.
  //    static MolPolEMFieldSetup* GetMolPolEMFieldSetup();
  /*
    private:
    static MolPolEMFieldSetup* fMolPolEMFieldSetup;
  */

public:
  MolPolEMFieldSetup() ;
  ~MolPolEMFieldSetup() ;

  void SetStepper();
  void UpdateField();

  void InitialseAll();

  inline  void SetStepperType( G4int val) { fStepperType = val ; }
  inline  G4int GetStepperType() {return fStepperType; }

  inline void SetMinStep(G4double val) { fMinStep = val ; }
  inline G4double GetMinStep() { return fMinStep ; }

  G4FieldManager* GetFieldManager(){return fFieldManager;}

  void UpdateConfiguration();

private:
  MolPolEMField*              fEMfield;
  G4FieldManager*             fFieldManager;
  G4ChordFinder*              fChordFinder ;
  G4EqMagElectricField*       fEquation ;
  G4MagIntegratorStepper*     fStepper ;
  G4MagInt_Driver*            fIntgrDriver;

  MolPolEMFieldMessenger*     fFieldMessenger;

  G4int                       fStepperType ;
  G4double                    fMinStep;

  // These don't change. I hate to add class variables. But I think they're better here
  G4double                    ORIGINQ1   =  75.19 * cm;
  G4double                    ORIGINQ2   = 140.46 * cm;
  G4double                    ORIGINQ3   = 209.08 * cm;
  G4double                    ORIGINQ4   = 274.59 * cm;
  G4double                    ORIGIND    = 423.4  * cm;
  G4double                    ORIGINQ6   = 6.9    * cm;
  G4double                    BORERADIUS = 5.08 * cm;

  G4RotationMatrix*           NOROT = new G4RotationMatrix;

public:
  // Types same as before 0: Ideal by amps, 1: Ideal by teslas, 2: TOSCA Map
  G4int                       iDipoleType = 1;
  G4int                       iQuad1Type  = 1;
  G4int                       iQuad2Type  = 1;
  G4int                       iQuad3Type  = 1;
  G4int                       iQuad4Type  = 1;

public:
  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Values which come from the macro. Must be public unless you want to write functions to fill.... :/
  //Dipole
  G4double                    dDipRelevantStr = 0.0;
  G4double                    dDipToscaMapStr = 1.0;
  G4String                    strDipoleMapLoc = " ";
  //Q4
  G4double                    dQuad4RelevantStr = 0.0;
  G4double                    dQuad4ToscaMapStr = 1.0;
  G4String                    strQuad4MapLoc = " ";
  //Q3
  G4double                    dQuad3RelevantStr = 0.0;
  G4double                    dQuad3ToscaMapStr = 1.0;
  G4String                    strQuad3MapLoc = " ";
  //Q2
  G4double                    dQuad2RelevantStr = 0.0;
  G4double                    dQuad2ToscaMapStr = 1.0;
  G4String                    strQuad2MapLoc = " ";
  //Q1
  G4double                    dQuad1RelevantStr = 0.0;
  G4double                    dQuad1ToscaMapStr = 1.0;
  G4String                    strQuad1MapLoc = " ";
  //SOLENOID
  G4double                    dSolRelevantStr = 0.0;

  G4double CalA2T(G4double current, G4int magnet);

};

#endif
