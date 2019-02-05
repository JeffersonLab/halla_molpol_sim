#include "MolPolQuad.hh"
#include "G4RotationMatrix.hh"

static G4RotationMatrix IdentityMatrix;

MolPolQuad::MolPolQuad(G4double pGradient, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius, G4double pMultContrib)
{
  fGradient    = pGradient ;
  fOrigin      = pOrigin ;
  fpMatrix     = pMatrix ;
  fRadius      = pRadius;
  fMultContrib = pMultContrib;
}

//Including second constructor with default negative 5% multipole contribution so we don't have to
//currently change anything in the EMFieldSetup.
MolPolQuad::MolPolQuad(G4double pGradient, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius)
{
  fGradient    = pGradient ;
  fOrigin      = pOrigin ;
  fpMatrix     = pMatrix ;
  fRadius      = pRadius;
  fMultContrib = -0.05;
}


void MolPolQuad::UpdateQuad(G4double pGradient, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius)
{
  fGradient    = pGradient ;
  fOrigin      = pOrigin ;
  fpMatrix     = pMatrix ;
  fRadius      = pRadius;
}


MolPolQuad::~MolPolQuad()
{
}

////////////////////////////////////////////////////////////////////////
//  Allow displaced origin and rotation
//  Extensions by Björn Riese (GSI)
void MolPolQuad::GetFieldValue( const G4double y[4], G4double B[3]  ) const
{
  B[0] = 0;
  B[1] = 0;
  B[2] = 0;

  G4ThreeVector r_global= G4ThreeVector(
     y[0] - fOrigin.x(),
     y[1] - fOrigin.y(),
     y[2] - fOrigin.z()
  );

  G4ThreeVector r_local = G4ThreeVector(
     fpMatrix->colX() * r_global,
     fpMatrix->colY() * r_global,
     fpMatrix->colZ() * r_global
  );

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //CARTESIAN VERSION
  G4double phiquad_local = atan2( r_local.y() , r_local.x() );
  G4double rhoquad_local = sqrt( r_local.x()*r_local.x() + r_local.y()*r_local.y() );

  G4ThreeVector B_local = G4ThreeVector(
    fGradient * rhoquad_local * ( sin(phiquad_local) + fMultContrib * pow(rhoquad_local,4 ) / pow(fRadius,4 ) * sin(5*phiquad_local) ),
    fGradient * rhoquad_local * ( cos(phiquad_local) + fMultContrib * pow(rhoquad_local,4 ) / pow(fRadius,4 ) * cos(5*phiquad_local) ),
    0
  );
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  /*
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //CYLINDRICAL CONVERTED TO CARTESIAN
  //Get r and phi in local coordinates
  G4double phiquad_local = atan2( r_local.y() , r_local.x() );
  G4double rhoquad_local = sqrt(r_local.x() * r_local.x() + r_local.y() * r_local.y());

  //cylindrical field (rho, phi, z)
  G4ThreeVector BCyl2_local = G4ThreeVector(
    fGradient * rhoquad_local * sin( 2 * phiquad_local ),
    fGradient * rhoquad_local * cos( 2 * phiquad_local ),
    0
  );

  G4ThreeVector BCyl6_local = G4ThreeVector(
    fMultContrib * fGradient * rhoquad_local * sin( 6 * phiquad_local ) * pow( rhoquad_local , 4 ) / pow( fRadius , 4 ),
    fMultContrib * fGradient * rhoquad_local * cos( 6 * phiquad_local ) * pow( rhoquad_local , 4 ) / pow( fRadius , 4 ),
    0
  );
  
  G4ThreeVector BCylTotal_local = BCyl2_local + BCyl6_local;

  //Convert to cartesian field
  G4ThreeVector B_local = G4ThreeVector(
    BCylTotal_local.x() * cos( phiquad_local ) - BCylTotal_local.y() * sin( phiquad_local ),
    BCylTotal_local.x() * sin( phiquad_local ) + BCylTotal_local.y() * cos( phiquad_local ),
    0
  );
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  */


  G4ThreeVector B_global = G4ThreeVector(
     fpMatrix->rowX() * B_local,
     fpMatrix->rowY() * B_local,
     fpMatrix->rowZ() * B_local
  );

  G4double rquad = (r_global.x() * r_global.x() + r_global.y() * r_global.y());

  if(rquad < (fRadius * fRadius)){
    B[0] = -1.0 * B_global.x() ;
    B[1] = -1.0 * B_global.y() ;
    B[2] = B_global.z() ;
  }
}
