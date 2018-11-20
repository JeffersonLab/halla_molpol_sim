#ifndef MOLPOLQUAD_HH
#define MOLPOLQUAD_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
class MolPolQuad : public G4MagneticField
{
public: // with description
  MolPolQuad(G4double pGradient,G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius, G4double pZeff);
  ~MolPolQuad();
  void GetFieldValue(const G4double yTrack[], G4double B[] ) const;
  G4ThreeVector GetFringeField(G4ThreeVector) const;
  void updateQuad(G4double pGradient,G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius, G4double pZeff);
  void updateQuad(G4double pGradient);

private:
  G4double          fGradient;
  G4ThreeVector     fOrigin;
  G4RotationMatrix* fMatrix;
  G4double          fRadius;
  G4double          fZeff;

public:
  void setQuadGradient(G4double);
  void setQuadOrigin(G4ThreeVector);
  void setQuadRotation(G4RotationMatrix*);
  void setQuadBoreRadius(G4double);
  void setQuadZeff(G4double);

};
#endif
