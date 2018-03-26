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
// $Id: MolPolDipole.cc 69786 2013-05-15 09:38:51Z gcosmo $
//
// -------------------------------------------------------------------

#include "MolPolDipole.hh"
#include "G4RotationMatrix.hh"
#include <math.h>

//#define DEBUG_BFIELD_DIPOLE 0

extern"C" {
  void dfringe_(float* x,    float* y,  float* z,
                float* fact, float* a,
                float* fx,   float* fy, float* fz);
  void snakedipole_(float* x,  float* y,  float* z,
		    float* bx, float* by, float* bz,
		    float* b0);
}

static G4RotationMatrix IdentityMatrix; 

MolPolDipole::MolPolDipole(G4double pBend)
{
   fBend    = pBend ;
   fOrigin  = G4ThreeVector( 0.0, 0.0, 0.0) ;
   fpMatrix = &IdentityMatrix;
}

/////////////////////////////////////////////////////////////////////////

MolPolDipole::MolPolDipole(G4double pBend, G4ThreeVector
pOrigin, G4RotationMatrix* pMatrix)
{
   fBend    = pBend   ;
   fOrigin  = pOrigin ;
   fpMatrix = pMatrix ;
}

/////////////////////////////////////////////////////////////////////////

MolPolDipole::~MolPolDipole()
{
}

////////////////////////////////////////////////////////////////////////
//  Allow displaced origin and rotation 
//  Extensions by Björn Riese (GSI)

void MolPolDipole::GetFieldValue( const G4double y[7],
				 G4double B[3]  ) const  
{

  G4ThreeVector r_global= G4ThreeVector
    (y[0] - fOrigin.x(),
     y[1] - fOrigin.y(),
     y[2] - fOrigin.z());

  G4ThreeVector r_local = G4ThreeVector
    (fpMatrix->colX() * r_global,
     fpMatrix->colY() * r_global,
     fpMatrix->colZ() * r_global);

  //indexed dipole
  //B = ( r / r_0 ) ^ -n  * B_0
  //G4double B_0 = fBend;

  float B_0 = (float)fBend;

  //float snx = r_local.x();
  //float sny = r_local.y();
  //float snz = r_local.z();
  //float snbx,snby,snbz;

  //snakedipole_(&snx,  &sny,  &snz,
  //&snbx, &snby, &snbz,
  //       &B_0);

  //  snbx *= -tesla;
  //snby *= -tesla;
  //snbz *= -tesla;

  //G4ThreeVector B_local = G4ThreeVector( snbx, snby, snbz );
  
  //G4ThreeVector B_global = G4ThreeVector
  //(fpMatrix->rowX() * B_local,
  // fpMatrix->rowY() * B_local,
  // fpMatrix->rowZ() * B_local);


  //if( sqrt( snbx * snbx + snby * snby ) > 2. * cm ){
  //B[0] = B_global.x() ;
    //B[1] = B_global.y() ;
    //B[2] = B_global.z() ;
    //}else{
    //B[0] = 0.;

  B[0] = B_0;
  B[1] = 0.;
  B[2] = 0.;
    //}


}
