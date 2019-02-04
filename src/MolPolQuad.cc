
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: MolPolQuad.cc 69786 2013-05-15 09:38:51Z gcosmo $
//
// -------------------------------------------------------------------

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


/////////////////////////////////////////////////////////////////////////

void MolPolQuad::UpdateQuad(G4double pGradient, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius)
{
  fGradient    = pGradient ;
  fOrigin      = pOrigin ;
  fpMatrix     = pMatrix ;
  fRadius      = pRadius;
}

/////////////////////////////////////////////////////////////////////////
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

  G4ThreeVector r_global= G4ThreeVector
    (y[0] - fOrigin.x(),
     y[1] - fOrigin.y(),
     y[2] - fOrigin.z());

  G4ThreeVector r_local = G4ThreeVector
    (fpMatrix->colX() * r_global,
     fpMatrix->colY() * r_global,
     fpMatrix->colZ() * r_global);

  //G4ThreeVector B_local = G4ThreeVector
  //  (fGradient * r_local.y(),
  //   fGradient * r_local.x(),
  //   0);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //PLEASE SEE THIS FIRST
  //MY RECOLLECTION WAS THAT I COULD DO THIS AS FOLLOWS, BUT IT HAS BEEN A WHILE, AND IT SHOULD BE r^N-2 / R^N but one R get's absorbed to make the gradient...  
  //Get phi
  //G4double phiquad_local = atan2( r_local.y() , r_local.x() );
  //Bx and By
  //G4ThreeVector B_local = G4ThreeVector(
  //  fGradient * r_local.() * ( sin(phiquad_local) + fMultContrib * pow(rquad_local,4 ) / pow(fRadius,4 ) * sin(5*phiquad_local) ),
  //  fGradient * r_local.() * ( sin(phiquad_local) + fMultContrib * pow(rquad_local,4 ) / pow(fRadius,4 ) * cos(5*phiquad_local) ),
  //  0
  //);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //WILL TYPE UP CYLINDRICAL VERSION FOR WHICH I HAVE DOCUMENTATION ON HAND FOR
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Get r and phi in local coordinates
  G4double phiquad_local = atan2( r_local.y() , r_local.x() );
  G4double rquad_local = sqrt(r_local.x() * r_local.x() + r_local.y() * r_local.y());
  //Three vector for B field in R, Phi, Z for the quadrupole and 12-pole
  G4ThreeVector BCyl2_local = G4ThreeVector(
    fGradient * rquad_local * sin( 2 * phiquad_local ),
    fGradient * rquad_local * cos( 2 * phiquad_local ),
    0
  );
  G4ThreeVector BCyl6_local = G4ThreeVector(
    fMultContrib * fGradient * rquad_local * sin( 6 * phiquad_local ) * pow( rquad_local , 4 ) / pow( fRadius , 4 ),
    fMultContrib * fGradient * rquad_local * cos( 6 * phiquad_local ) * pow( rquad_local , 4 ) / pow( fRadius , 4 ),
    0
  );
  //Sum the vectors
  G4ThreeVector BCylTotal_local = BCyl2_local + BCyl6_local;
  //Convert to cartesian space - THIS REPLACES THE BLOCAL COMMENTED OUT ABOVE
  //AFTERWARDS IT PROCEEDS AS BEFORE
  G4ThreeVector B_local = G4ThreeVector(
    BCylTotal_local.x() * cos( phiquad_local ) - BCylTotal_local.y() * sin( phiquad_local ),
    BCylTotal_local.x() * sin( phiquad_local ) + BCylTotal_local.y() * cos( phiquad_local ),
    0
  );
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




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
