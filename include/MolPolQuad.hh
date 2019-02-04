#ifndef MOLPOLQUAD_HH
#define MOLPOLQUAD_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
class MolPolQuad : public G4MagneticField
{
public: // with description

  MolPolQuad(G4double          pGradient,
             G4ThreeVector     pOrigin,
             G4RotationMatrix* pMatrix,
             G4double          pRadius,
	     G4double          pMultContrib
            );

  //Including second constructor with default 5% multipole contribution so we don't have
  //to change anything in field setup yet.
  MolPolQuad(G4double          pGradient,
             G4ThreeVector     pOrigin,
             G4RotationMatrix* pMatrix,
             G4double          pRadius
            );


  ~MolPolQuad();

  void GetFieldValue( const G4double yTrack[], G4double B[] ) const;
  G4ThreeVector GetFringeField(G4ThreeVector) const;
  void UpdateQuad(G4double, G4ThreeVector, G4RotationMatrix*, G4double);

private:
  G4double          fGradient;
  G4ThreeVector     fOrigin;
  G4RotationMatrix* fpMatrix;
  G4double          fRadius;
  G4double          fMultContrib;

};
#endif
