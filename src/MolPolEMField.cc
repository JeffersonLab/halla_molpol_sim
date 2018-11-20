// ********************************************************************
//
// $Id: MolPolEMField.hh,v 1.0, 2010/12/26   MolPol Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//   User Field class Setup implementation.
//
//
#include "MolPolEMField.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

MolPolEMField::MolPolEMField()
{
  EField3V.set(0,0,0);
  BField3V.set(0,0,0);
}

//////////////////////////////////////////////////////////////////////////
//
//  Deconstructors:
MolPolEMField::~MolPolEMField()
{

}


//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Get Field Value
inline void MolPolEMField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{
  // Member func is 'const' so everything used within must be declared within.
  G4double Bsum[3]  = {0.,0.,0.};
  G4double thisB[3] = {0.,0.,0.};
  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // B-Field
  Bfield[0]=BField3V.x();
  Bfield[1]=BField3V.y();
  Bfield[2]=BField3V.z();
  //G4cout << "CALCULATE THE FIELD VALUE......................................................................................." << G4endl;
  // Get solenoid field
  if(bBfieldInUse[0] == true){
    cSolenoidField->GetFieldValue(Point,thisB);
    //G4cout << "\nVolume Summary Solenoid @ P(" << Point[0]/cm << "," << Point[1]/cm << "," << Point[2]/cm << ") and B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[1] == true){
    cQuad1Field->GetFieldValue(Point,thisB);
    //G4cout << "Volume Summary Quad1 @ P(" << Point[0]/cm << "," << Point[1]/cm << "," << Point[2]/cm << ") and B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[2] == true){
    cQuad2Field->GetFieldValue(Point,thisB);
    //G4cout << "Volume Summary Quad2 @ P(" << Point[0]/cm << "," << Point[1]/cm << "," << Point[2]/cm << ") and B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[3] == true){
    cQuad3Field->GetFieldValue(Point,thisB);
    //G4cout << "Volume Summary Quad3 @ P(" << Point[0]/cm << "," << Point[1]/cm << "," << Point[2]/cm << ") and B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[4] == true){
    cQuad4Field->GetFieldValue(Point,thisB);
    //G4cout << "Volume Summary Quad4 @ P(" << Point[0]/cm << "," << Point[1]/cm << "," << Point[2]/cm << ") and B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[5] == true){
    cDipoleField->GetFieldValue(Point,thisB);
    //G4cout << "Volume Summary Dipole @ P(" << Point[0]/cm << "," << Point[1]/cm << "," << Point[2]/cm << ") and B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }

  // Put Bsum into BField array which will
  for (int i = 0; i < 3; i++) Bfield[i] = Bsum[i];

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // E-Field -- if ever so desired
  double  *Efield=&Bfield[3];
  Efield[0]=EField3V.x();
  Efield[1]=EField3V.y();
  Efield[2]=EField3V.z();
}
