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

MolPolQuad::MolPolQuad(G4double pGradient, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius)
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

void MolPolQuad::GetFieldValue( const G4double y[4],
                                G4double B[3]  ) const
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

  G4ThreeVector B_local = G4ThreeVector
    (fGradient * r_local.y(),
     fGradient * r_local.x(),
     0);

  G4ThreeVector B_global = G4ThreeVector
    (fpMatrix->rowX() * B_local,
     fpMatrix->rowY() * B_local,
     fpMatrix->rowZ() * B_local);

  G4double rquad = (r_global.x() * r_global.x() + r_global.y() * r_global.y());
  if(rquad < (fRadius * fRadius))
    {
      B[0] = -1.0 * B_global.x() ;
      B[1] = -1.0 * B_global.y() ;
      B[2] = B_global.z() ;
    }
}
