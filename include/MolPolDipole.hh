#ifndef MOLPOLDIPOLE_HH
#define MOLPOLDIPOLE_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

class MolPolDipole : public G4MagneticField
{
public: // with description

  MolPolDipole(G4double pBend, G4double pZeff);
  MolPolDipole(G4double pBend, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pZeff);
  ~MolPolDipole();

  void updateDipole(G4double, G4ThreeVector, G4RotationMatrix*);
  void GetFieldValue(const G4double yTrack[], G4double B[]) const;

private:
  G4double          fBend;
  G4ThreeVector     fOrigin = G4ThreeVector(0.,0.,423.4*cm);
  G4RotationMatrix* fMatrix;
  G4double          fZeff;

public:
  void setDipoleStrength(G4double pBend);
  void setDipoleZeff(G4double pZeff);
  void setRotationMatrix(G4RotationMatrix* pMatrix);
  void setDipoleOrigin(G4ThreeVector pOrigin);
  void updateDipole(G4double pBend, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pZeff);


};
#endif
