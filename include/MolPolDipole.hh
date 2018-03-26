#ifndef MOLPOLDIPOLE_HH
#define MOLPOLDIPOLE_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

class MolPolDipole : public G4MagneticField
{
public: // with description

  MolPolDipole(G4double          pBend);
  MolPolDipole(G4double          pBend,
               G4ThreeVector     pOrigin,
               G4RotationMatrix* pMatrix);
  ~MolPolDipole();

  void UpdateDipole(G4double, G4ThreeVector, G4RotationMatrix*);
  void GetFieldValue(const G4double yTrack[],
                     G4double B[]     ) const;
private:

  G4double          fBend;
  G4ThreeVector     fOrigin;
  G4RotationMatrix* fpMatrix;
};
#endif
