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
  //  G4cout << "  Quad @ P(" << y[0]/10 << "," << y[1]/10 << "," << y[2]/10 << ") cm" << G4endl;
  //  G4cout << "  Quad Boundary Z(" << (fOrigin.z() - 0.5 * fZeff)/10 << "," << (fOrigin.z() + 0.5 * fZeff)/10 << ") cm" << G4endl;


  // what is rsquared
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
