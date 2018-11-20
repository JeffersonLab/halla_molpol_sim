// *************************************************************** (╯°□°）╯︵ ┻━┻
//
//  MolPolQuad.cc
//
//  Updated MolPolQuad for use with integrated field types within a global
//  EMField manager. Relevenant changes are strictly limited to a change from
//  local coordinates to global coordinates in the calculation. Previously the
//  class was written for use in local field volumes.
//
//  There are some G4cout lines which were used while debugging. These can be
//  safely removed in the future.
//
//  Constructor now also takes the Zeff for any idealized quad field.
//
//  Radius could probably be removed from the constructors.
//
//
//  Eric King - $date
//
// *****************************************************************************

#include "MolPolQuad.hh"
#include "G4RotationMatrix.hh"

static G4RotationMatrix IdentityMatrix;

MolPolQuad::MolPolQuad(G4double pGradient, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius, G4double pZeff)
{
  fGradient    = pGradient ;
  fOrigin      = pOrigin ;
  fMatrix      = pMatrix ;
  fRadius      = pRadius;
  fZeff        = pZeff;
}


/////////////////////////////////////////////////////////////////////////

void MolPolQuad::updateQuad(G4double pGradient, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius, G4double pZeff){
  fGradient    = pGradient ;
  fOrigin      = pOrigin ;
  fMatrix      = pMatrix ;
  fRadius      = pRadius;
  fZeff        = pZeff;
}

void MolPolQuad::updateQuad(G4double pGradient){
  fGradient    = pGradient;
}

/////////////////////////////////////////////////////////////////////////
MolPolQuad::~MolPolQuad(){
}

////////////////////////////////////////////////////////////////////////
//  Allow displaced origin and rotation
//  Extensions by Björn Riese (GSI)

void MolPolQuad::GetFieldValue( const G4double y[4], G4double B[3]  ) const
{
  // note: now that the Quads are part of global field rather than being calculated
  // by local field manager in logical volume y[] are alread in global coordinates
  // so the calculation of local is unnecessary.

  //initialize B array to zero.
  B[0] = 0;
  B[1] = 0;
  B[2] = 0;

  G4ThreeVector position = G4ThreeVector(y[0],y[1],y[2]);
  //G4cout << "  Quad @ P(" << y[0]/10 << "," << y[1]/10 << "," << y[2]/10 << ") cm" << G4endl;
  //G4cout << "  Quad Boundary Z(" << (fOrigin.z() - 0.5 * fZeff)/10 << "," << (fOrigin.z() + 0.5 * fZeff)/10 << ") cm" << G4endl;

  G4double rsquared = (position.x() * position.x() + position.y() * position.y());

  // first, check the least likely thing that z is within allowed space for idealized quad field
  // then, check rsquared
  // then calculate and return field

  if( (position.z() > (fOrigin.z() - 0.5 * fZeff)) && (position.z() < (fOrigin.z() + 0.5 * fZeff)) ){
    if(rsquared < (fRadius * fRadius)){

      // idealized quad field Bx = kappa * y, By = kappa * x
      G4ThreeVector idealizedField = G4ThreeVector(fGradient * position.y(),
                                                   fGradient * position.x(),
                                                   0.0);
      // do we want quad field misaligned (leaving as option for Don)
      G4ThreeVector rotatedField = G4ThreeVector(fMatrix->rowX() * idealizedField,
                                                 fMatrix->rowY() * idealizedField,
                                                 fMatrix->rowZ() * idealizedField);
      B[0] = -1.0 * rotatedField.x() ;
      B[1] = -1.0 * rotatedField.y() ;
      B[2] = rotatedField.z() ;
    }
  } else {
    //G4cout << "  OUT OF QUAD BOUNDARY!" << G4endl;
  }

}
