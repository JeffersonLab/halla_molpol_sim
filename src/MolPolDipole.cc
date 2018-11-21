// *************************************************************** (╯°□°）╯︵ ┻━┻
//
//	MolPolDipole.cc
//
//  Ideal dipole implementation for integrated field types. Object constructed
//  as part of global field managed from MolPolEMfield().
//
//  This now takes Zeff as an argument and should be updated to include bounds
//  for x and y as well.
//
//  Since this script was modified from one utilized in local field management
//  this was changed to utilize global field coordinates.
//
//  Will clean up file at later date. :/
//
//	Eric King - 2018-11-19
//
// *****************************************************************************

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

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Dipole Constructors -- First (standard), Second (specified position and rotation)
MolPolDipole::MolPolDipole(G4double pBend, G4double pZeff){
  fBend   = pBend;
  fZeff   = pZeff;
}

MolPolDipole::MolPolDipole(G4double pBend, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pZeff){
  fBend   = pBend;
  fOrigin = pOrigin;
  fMatrix = pMatrix;
  fZeff   = pZeff;
}

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Obligatory Deconstructor
MolPolDipole::~MolPolDipole()
{
}

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Member functions to update dipole -- maybe just move these into class library ... not sure if individuals are needed.  Probably not thinking about it.
void MolPolDipole::setDipoleStrength(G4double pBend){
  fBend    = pBend;
}
void MolPolDipole::setDipoleZeff(G4double pZeff){
  fZeff    = pZeff;
}
void MolPolDipole::setRotationMatrix(G4RotationMatrix* pMatrix){
  fMatrix = pMatrix;
}
void MolPolDipole::setDipoleOrigin(G4ThreeVector pOrigin){
  fOrigin = pOrigin;
}

void MolPolDipole::updateDipole(G4double pBend, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pZeff){
  fBend    = pBend;
  fOrigin  = pOrigin;
  fMatrix  = pMatrix;
  fZeff    = pZeff;
}

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Get Field Value ... Not sure who Bjorn Riese is... haven't touched any calculations below.
////////////////////////////////////////////////////////////////////////
//  Allow displaced origin and rotation
//  Extensions by Björn Riese (GSI)

void MolPolDipole::GetFieldValue( const G4double y[7], G4double B[3]  ) const
{

  G4ThreeVector position= G4ThreeVector(y[0],y[1],y[2]);

  //Don't need to turn this into local coordinates since EMfield is global.

  //G4ThreeVector r_local = G4ThreeVector
  //  (fMatrix->colX() * r_global,
  //   fMatrix->colY() * r_global,
  //   fMatrix->colZ() * r_global);

  //indexed dipole
  //B = ( r / r_0 ) ^ -n  * B_0
  //G4double B_0 = fBend;

  G4double B_0 = (G4double)fBend;

  //float snx = r_local.x();
  //float sny = r_local.y();
  //float snz = r_local.z();
  //float snbx,snby,snbz;

  //snakedipole_(&snx,  &sny,  &snz,
  //&snbx, &snby, &snbz,
  //       &B_0);

  //snbx *= -tesla;
  //snby *= -tesla;
  //snbz *= -tesla;

  //G4ThreeVector B_local = G4ThreeVector( snbx, snby, snbz );

  //G4ThreeVector B_global = G4ThreeVector
  //(fMatrix->rowX() * B_local,
  // fMatrix->rowY() * B_local,
  // fMatrix->rowZ() * B_local);


  //if( sqrt( snbx * snbx + snby * snby ) > 2. * cm ){
  //B[0] = B_global.x() ;
  //B[1] = B_global.y() ;
  //B[2] = B_global.z() ;
  //}else{
  //B[0] = 0.;

  G4cout << "Dipole Boundary Z(" << (fOrigin.z() - 0.5 * fZeff)/10 << "," << (fOrigin.z() + 0.5 * fZeff)/10 << ") cm" << G4endl;

  if( (position.z() > (fOrigin.z() - 0.5 * fZeff)) && (position.z() < (fOrigin.z() + 0.5 * fZeff)) ){
    //if( (position.y() > (fOrigin.y() - 300.)) && (position.y() < (fOrigin.y() + 300.)) ){
      //if( (position.x() > (fOrigin.x() - 300.)) && (position.x() < (fOrigin.x() + 300.)) ){
        B[0] = B_0;
        B[1] = 0.;
        B[2] = 0.;
      //}
    //}
  } else {
    G4cout << "OUT OF DIPOLE Z BOUNDARY!" << G4endl;
  }
  //}


}
