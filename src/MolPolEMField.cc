// *************************************************************** (╯°□°）╯︵ ┻━┻
//
//	MolPolEMfield.cc
//
//  This class serves, in the integrated field version of G4MolPol as the global
//  field manager which manages all G4MagneticField derived class objects in the
//  simulation. These fields ARE DECOUPLED from the geometry.
//
//  Bfield and Efield initialize in constructors as (0,0,0) everywhere.
//
//  Objects in MolPolEMfield() are managed through MolPolEMfieldSetup(). this
//  class should STRICTLY serve as the field object.
//
//  thisB[] is reset after probing each field region.
//
//  ***There are probably ways to slightly speed up this portion of the simulation.
//  (1) Fields now have on-off boolean assigned.
//  (2) Perhaps limit boundaries for field checking by Zeff? Logical checks probably
//      longer than calculation of field.
//
//	Eric King - 2018-11-19
//
// *****************************************************************************
#include "MolPolEMField.hh"

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Constructors
MolPolEMField::MolPolEMField()
{
  EField3V.set(0,0,0);
  BField3V.set(0,0,0);
}

MolPolEMField::~MolPolEMField()
{

}

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Obligatory GetFieldValue() method for any G4MagneticField derived class.
inline void MolPolEMField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{
  // Member func is 'const' so everything used within must be declared within.
  G4double Bsum[3]  = {0.,0.,0.};
  G4double thisB[3] = {0.,0.,0.};
  ////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // B-Field
  Bfield[0]=BField3V.x();
  Bfield[1]=BField3V.y();
  Bfield[2]=BField3V.z();
  G4cout << "------------------------------------------------------------------------------------\nCALCULATE THE FIELD VALUES @ P(" << Point[0]/cm << "," << Point[1]/cm << "," << Point[2]/cm << ")" << G4endl;
  if(bBfieldInUse[0] == true){
    G4cout << "Solenoid" << G4endl;
    for(G4int i = 0; i < 3; i++) thisB[i] = 0.0;
    cSolenoidField->GetFieldValue(Point,thisB);
    G4cout << "  Solenoid reported B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[1] == true){
    G4cout << "Quad 1" << G4endl;
    for(G4int i = 0; i < 3; i++) thisB[i] = 0.0;
    cQuad1Field->GetFieldValue(Point,thisB);
    G4cout << "  Quad1 reported B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[2] == true){
    G4cout << "Quad 2" << G4endl;
    for(G4int i = 0; i < 3; i++) thisB[i] = 0.0;
    cQuad2Field->GetFieldValue(Point,thisB);
    G4cout << "  Quad2 reported B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[3] == true){
    G4cout << "Quad 3" << G4endl;
    for(G4int i = 0; i < 3; i++) thisB[i] = 0.0;
    cQuad3Field->GetFieldValue(Point,thisB);
    G4cout << "  Quad3 reported B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[4] == true){
    G4cout << "Quad 4" << G4endl;
    for(G4int i = 0; i < 3; i++) thisB[i] = 0.0;
    cQuad4Field->GetFieldValue(Point,thisB);
    G4cout << "  Quad4 reported B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }
  if(bBfieldInUse[5] == true){
    G4cout << "Dipole" << G4endl;
    for(G4int i = 0; i < 3; i++) thisB[i] = 0.0;
    cDipoleField->GetFieldValue(Point,thisB);
    G4cout << "  Dipole reported B(" << thisB[0]/tesla << "," << thisB[1]/tesla << "," << thisB[2]/tesla << ")" << G4endl;
    for (G4int i = 0; i < 3; i++) Bsum[i] += thisB[i];
  }

  // Put Bsum into BField array which will
  for (int i = 0; i < 3; i++) Bfield[i] = Bsum[i];

  G4cout << "TOTAL Bfield reported is (" << Bfield[0]/tesla << "," << Bfield[1]/tesla << "," << Bfield[2]/tesla << ")" << G4endl;

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // E-Field -- if ever so desired
  G4double *Efield = &Bfield[3];
  Efield[0] = EField3V.x();
  Efield[1] = EField3V.y();
  Efield[2] = EField3V.z();
}
