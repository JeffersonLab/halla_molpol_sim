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

  //Local field  FZB1
  void UpdateFieldFZB1();
  void SetBField3VFZB1(G4double fieldGradient);
  G4FieldManager* GetFieldManagerFZB1(){return fLocalFieldManagerFZB1;}
  //Local field  FZB2
  void UpdateFieldFZB2();
  void SetBField3VFZB2(G4double fieldGradient);
  G4FieldManager* GetFieldManagerFZB2(){return fLocalFieldManagerFZB2;}
  //Local field  FZB3
  void UpdateFieldFZB3();
  void SetBField3VFZB3(G4double fieldGradient);
  G4FieldManager* GetFieldManagerFZB3(){return fLocalFieldManagerFZB3;}
  //Local field  FZB4
  void UpdateFieldFZB4();
  void SetBField3VFZB4(G4double fieldGradient);
  G4FieldManager* GetFieldManagerFZB4(){return fLocalFieldManagerFZB4;}
  //Local field  FZB5
  void UpdateFieldFZB5();
  void SetBField3VFZB5(G4double fieldGradient);
  G4FieldManager* GetFieldManagerFZB5(){return fLocalFieldManagerFZB5;}

  //Local field  FZB6
  void UpdateFieldFZB6();
  void SetBField3VFZB6(G4double fieldGradient);
  G4FieldManager* GetFieldManagerFZB6(){return fLocalFieldManagerFZB6;}

  void UpdateConfiguration();

  G4int fMagSourceMode;

  //external input current values
  G4double                    fQ1A;
  G4double                    fQ2A;
  G4double                    fQ3A;
  G4double                    fQ4A;
  G4double                    fQ5A;
  G4double                    fQ6A;
  //external input field (pole tip) values
  G4double                    fQ1T;
  G4double                    fQ2T;
  G4double                    fQ3T;
  G4double                    fQ4T;
  G4double                    fQ5T;
  G4double                    fQ6T;

  G4double                    ORIGINQ1 =  75.19 * cm;;
  G4double                    ORIGINQ2 = 140.46 * cm;;
  G4double                    ORIGINQ3 = 209.08 * cm;;
  G4double                    ORIGINQ4 = 274.59 * cm;;
  G4double                    ORIGIND  = 423.4  * cm;;
  G4double                    ORIGINQ6 = 6.9    * cm;;
  G4double                    BORERADIUS = 5.08 * cm;

  G4RotationMatrix*           NOROT = new G4RotationMatrix;

  G4String                    fToscaFields[6];
  G4double                    fToscaOffset[6] = {75.19 * cm, 140.46 *cm, 209.08 * cm, 274.59 * cm, 423.4 * cm, 6.9 * cm}; 

private:
  MolPolEMField*              fEMfield;
  G4FieldManager*             fFieldManager;
  G4ChordFinder*              fChordFinder ;
  G4EqMagElectricField*       fEquation ;
  G4MagIntegratorStepper*     fStepper ;
  G4MagInt_Driver*            fIntgrDriver;
  G4MagInt_Driver*            fIntgrDriverFZB1;
  G4MagInt_Driver*            fIntgrDriverFZB2;
  G4MagInt_Driver*            fIntgrDriverFZB3;
  G4MagInt_Driver*            fIntgrDriverFZB4;
  G4MagInt_Driver*            fIntgrDriverFZB5;
  G4MagInt_Driver*            fIntgrDriverFZB6;

  MolPolEMFieldMessenger*     fFieldMessenger;

  G4int                       fStepperType ;
  G4double                    fMinStep ;

  std::vector<G4String>       fileNames;
  std::vector<G4double>       fileScales;
  std::vector<G4double>       fileOffsets;

  //for local field at FZB1 and FZB2
  MolPolQuad*                 fMagFieldFZB1 ;
  G4Mag_UsualEqRhs*           fEquationFZB1 ;
  G4ChordFinder*              fChordFinderFZB1 ;
  G4MagIntegratorStepper*     fStepperFZB1 ;
  G4FieldManager*             fLocalFieldManagerFZB1;
  MolPolQuad*                 fMagFieldFZB2 ;
  G4Mag_UsualEqRhs*           fEquationFZB2 ;
  G4ChordFinder*              fChordFinderFZB2 ;
  G4MagIntegratorStepper*     fStepperFZB2 ;
  G4FieldManager*             fLocalFieldManagerFZB2;
  MolPolQuad*                 fMagFieldFZB3 ;
  G4Mag_UsualEqRhs*           fEquationFZB3 ;
  G4ChordFinder*              fChordFinderFZB3 ;
  G4MagIntegratorStepper*     fStepperFZB3 ;
  G4FieldManager*             fLocalFieldManagerFZB3;
  MolPolQuad*                 fMagFieldFZB4 ;
  G4Mag_UsualEqRhs*           fEquationFZB4 ;
  G4ChordFinder*              fChordFinderFZB4 ;
  G4MagIntegratorStepper*     fStepperFZB4 ;
  G4FieldManager*             fLocalFieldManagerFZB4;
  MolPolDipole*               fMagFieldFZB5 ;
  G4Mag_UsualEqRhs*           fEquationFZB5 ;
  G4ChordFinder*              fChordFinderFZB5 ;
  G4MagIntegratorStepper*     fStepperFZB5 ;
  G4FieldManager*             fLocalFieldManagerFZB5;
  MolPolSolenoid*             fMagFieldFZB6 ;
  G4Mag_UsualEqRhs*           fEquationFZB6 ;
  G4ChordFinder*              fChordFinderFZB6 ;
  G4MagIntegratorStepper*     fStepperFZB6 ;
  G4FieldManager*             fLocalFieldManagerFZB6;

  G4double CalA2T(G4double current, G4int magnet);

};

#endif
