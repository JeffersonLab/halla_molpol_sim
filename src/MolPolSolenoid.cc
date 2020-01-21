// *************************************************************** (╯°□°）╯︵ ┻━┻
//
//	MolPolSolenoid.cc
//
//  Permenant implementation of solenoid map from Bill Henry. Essentially
//  copied MolPolTOSCAField to MolPol solenoid to act as dedicated class that
//  only reqires the desired field strength.
//
//  Constructor and UpdateSolenoid() both take the field strength in teslas
//  as their argument.
//
//
//
//
//
//	Eric King - 2018-11-19
//
// *****************************************************************************

#include "MolPolSolenoid.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>

MolPolSolenoid::MolPolSolenoid(G4double pFieldStrengthT, G4ThreeVector pCenter, G4double pXrot, G4double pYrot){
  fZoffset    = pCenter.z();
  fYoffset    = pCenter.y();
  fXoffset    = pCenter.x();
  fYrot       = pYrot;
  fXrot       = pXrot;
  fFilename   = "../solenoid.map";
  fInit       = false;
  fFieldScale = pFieldStrengthT / 5.0;
  ReadFieldMap();
}

void MolPolSolenoid::UpdateSolenoid(G4double pFieldStrengthT, G4ThreeVector pCenter, G4double pXrot, G4double pYrot){
  fFieldScale = pFieldStrengthT / 5.0;
  fZoffset    = pCenter.z();
  fYoffset    = pCenter.y();
  fXoffset    = pCenter.x();
  fYrot       = pYrot;
  fXrot       = pXrot;
  G4cout << "  New solenoid desired field strength: " << pFieldStrengthT << G4endl;
  G4cout << "  Scale factor required for 5.0 tesla map: " << fFieldScale << G4endl;
  G4cout << "  Solenoid Offset in X: " << fXoffset / mm << G4endl;
  G4cout << "  Solenoid Offset in Y: " << fYoffset / mm << G4endl;
  G4cout << "  Solenoid Offset in Z: " << fZoffset / mm << G4endl;
  G4cout << "  Solenoid Rotation in X: " << fXrot / deg << G4endl;
  G4cout << "  Solenoid Rotation in Y: " << fYrot / deg << G4endl;
}

MolPolSolenoid::~MolPolSolenoid(){
}

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Initialize Solenoid Grid -- Set up size of vector for field data
void MolPolSolenoid::InitializeGrid() {
  if( NX <= 0 || NY <= 0 || NZ <= 0 ){
    G4cerr << "Error " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ": grid size invalid.  Aborting" << G4endl;
    exit(1);
  }
  G4cout << "  Initializing field map grid for " << fFilename << G4endl;
  for( int i = 0; i < 3; i++ ){
    fBFieldData[i].clear();
    fBFieldData[i].resize(NX);
  	for( int x = 0; x < NX; x++) {
      fBFieldData[i][x].resize(NY);
      for( int y = 0; y < NY; y++) {
    		fBFieldData[i][x][y].resize(NZ);
    		for( int z = 0; z < NZ; z++) {
    	    fBFieldData[i][x][y][z] = 0.0;
        }
      }
    }
  }
  G4cout << "  Solenoid map " << fFilename     << " initialized" << G4endl;
  G4cout << "  With scale: "  << fFieldScale   << " set"         << G4endl;
  return;
}

void MolPolSolenoid::ReadFieldMap(){
  G4cout << "  " <<  __PRETTY_FUNCTION__ << " : " << fFilename << G4endl;
  G4int nX = 0, nY = 0, nZ = 0;
  G4int nlines = 0;
  G4double raw_xpos;
  G4double raw_ypos;
  G4double raw_zpos;
  G4double raw_B;
  G4double raw_Bx;
  G4double raw_By;
  G4double raw_Bz;

  std::ifstream inputfile;
  inputfile.open(fFilename.data());

  if (!inputfile.good() ){
    G4cerr << "Error " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ") File " << fFilename << " could not open.  Aborting" << G4endl;
    exit(1);
  }
  std::string inputline;
  // Get grid dimensions off of first line
  getline(inputfile,inputline);
  if (std::istringstream(inputline) >> NZ >> NY >> NX) {
    G4cout << "  GridSize[" << NX << "," << NY << "," << NZ << "]" << G4endl;
  } else {
    G4cerr << "(" << __PRETTY_FUNCTION__ << " line " << __LINE__ << ") Error! File "
           << fFilename << " contains unreadable header information. Aborting field mapping." << G4endl;
    exit(1);
  }

  InitializeGrid();

  std::vector<G4double> xpos;
  std::vector<G4double> ypos;
  std::vector<G4double> zpos;
  for(nX = 0; nX < NX; nX++){
    for(nY = 0; nY < NY; nY++){
      for(nZ = 0; nZ < NZ; nZ++){
        getline(inputfile,inputline);

        if (std::istringstream(inputline) >> raw_xpos >> raw_ypos >> raw_zpos >> raw_B >> raw_Bx >> raw_By >> raw_Bz) {
          nlines++;
        } else {
          G4cerr << "Error " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ") File " << fFilename << " contains invalid data on line " << nlines << " Aborting" << G4endl;
          exit(1);
        }
        xpos.push_back(raw_xpos*cm);
        ypos.push_back(raw_ypos*cm);
        zpos.push_back(raw_zpos*cm);

    		fBFieldData[0][nX][nY][nZ] = raw_Bx * gauss;
    		fBFieldData[1][nX][nY][nZ] = raw_By * gauss;
    		fBFieldData[2][nX][nY][nZ] = raw_Bz * gauss;
      }
    }
  }

  XMIN = *std::min_element(xpos.begin(),xpos.end());
  XMAX = *std::max_element(xpos.begin(),xpos.end());
  YMIN = *std::min_element(ypos.begin(),ypos.end());
  ZMIN = *std::min_element(zpos.begin(),zpos.end());
  YMAX = *std::max_element(ypos.begin(),ypos.end());
  ZMAX = *std::max_element(zpos.begin(),zpos.end());

  fInit = true;
  G4cout << "  ... read " << nlines << " lines." << G4endl;

  // The Ultimate Style Guide :P
  // https://twitter.com/ThePracticalDev/status/710156980535558144
}

void MolPolSolenoid::GetFieldValue(const G4double Point[4], G4double *Bfield ) const {

  G4double    fracVal[3], intVal[3];
  G4int       gridIndex[3];

  // ***** print checks - leaving in for the moment. will remove after independent validation from Don. ***** //
  G4bool   printCheck(false);  if(printCheck) G4cout << "===================================================================" << G4endl;
  if(printCheck) G4cout << "Global Point:" << G4endl;
  if(printCheck) G4cout << "    X: " << Point[0] / mm << "mm" << G4endl;
  if(printCheck) G4cout << "    Y: " << Point[1] / mm << "mm" << G4endl;
  if(printCheck) G4cout << "    Z: " << Point[2] / mm << "mm" << G4endl;

  G4ThreeVector rho_map_local = G4ThreeVector(
    Point[0] - fXoffset,
    Point[1] - fYoffset,
    Point[2] - fZoffset
  );

  //print checks - remove before pulling into anything else.
  if(printCheck) G4cout << "Offset Point:" << G4endl;
  if(printCheck) G4cout << "    x: " << rho_map_local.x() / mm << "mm" << G4endl;
  if(printCheck) G4cout << "    y: " << rho_map_local.y() / mm << "mm" << G4endl;
  if(printCheck) G4cout << "    z: " << rho_map_local.z() / mm << "mm" << G4endl;

  // Construct the G4 rotation matrix
  G4RotationMatrix solenoidRotation = G4RotationMatrix();
  solenoidRotation.rotateX( fXrot );
  solenoidRotation.rotateY( fYrot );
  solenoidRotation.rotateZ( 0. * deg);

  // Rotate rho (in inverse direction) and find poition on map.
  G4ThreeVector rho_map_rotated = G4ThreeVector(
    solenoidRotation.colX() * rho_map_local,
    solenoidRotation.colY() * rho_map_local,
    solenoidRotation.colZ() * rho_map_local
  );

  //print checks - remove before pulling into anything else.
  if(printCheck) G4cout << "Rotated Point:" << G4endl;
  if(printCheck) G4cout << "   x': " << rho_map_rotated.x() / mm << "mm" << G4endl;
  if(printCheck) G4cout << "   y': " << rho_map_rotated.y() / mm << "mm" << G4endl;
  if(printCheck) G4cout << "   z': " << rho_map_rotated.z() / mm << "mm" << G4endl;

  // Check that the rotated point is within the defined region before interpolation.  If it is outside then the field is zero
  if( rho_map_rotated.x()  < XMIN || rho_map_rotated.x() >= XMAX ||
      rho_map_rotated.y() >= YMAX || rho_map_rotated.y() <  YMIN ||
      rho_map_rotated.z() >= ZMAX || rho_map_rotated.z() <  ZMIN ){
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
    return;
  }

  // Ensures we're going to get our grid indices correct
  assert( XMIN <= rho_map_rotated.x() && rho_map_rotated.x() < XMAX );
  assert( YMIN <= rho_map_rotated.y() && rho_map_rotated.y() < YMAX );
  assert( ZMIN <= rho_map_rotated.z() && rho_map_rotated.z() < ZMAX );

  // Gets the fractional and integer values for interpolation
  fracVal[0] = modf( ( rho_map_rotated.x() - XMIN )*(NX-1)/( XMAX - XMIN ), &(intVal[0]) );
  fracVal[1] = modf( ( rho_map_rotated.y() - YMIN )*(NY-1)/( YMAX - YMIN ), &(intVal[1]) );
  fracVal[2] = modf( ( rho_map_rotated.z() - ZMIN )*(NZ-1)/( ZMAX - ZMIN ), &(intVal[2]) );

  // This was from the HallC TOSCA code and should be covered by the first check to see if the particle is in the map-defined region, but I'll leave it.
  for( G4int i = 0; i < 3; i++ ) gridIndex[i] = (G4int) intVal[i];
  assert( 0 <= gridIndex[0] && gridIndex[0] < NX );
  assert( 0 <= gridIndex[1] && gridIndex[1] < NY );
  assert( 0 <= gridIndex[2] && gridIndex[2] < NZ );
  for( G4int i = 0; i < 3; i++ ) gridIndex[i] = (G4int) intVal[i];

  // Interpolate the magnetic field from the map.
  G4ThreeVector Bint_local = G4ThreeVector(0.,0.,0.);
  G4double c00, c10, c01, c11, c0, c1;
  for( G4int i = 0; i < 3; i++ ){
    c00 = fBFieldData[i][ gridIndex[0] ][ gridIndex[1] ][ gridIndex[2] ] * (1.0-fracVal[0])
          + fBFieldData[i][ gridIndex[0]+1 ][ gridIndex[1] ][ gridIndex[2] ]*fracVal[0];
    c01 = fBFieldData[i][ gridIndex[0] ][ gridIndex[1] ][ gridIndex[2]+1 ] * (1.0-fracVal[0])
          + fBFieldData[i][ gridIndex[0]+1 ][ gridIndex[1] ][ gridIndex[2]+1 ]*fracVal[0];
    c10 = fBFieldData[i][ gridIndex[0] ][ gridIndex[1]+1 ][ gridIndex[2] ] * (1.0-fracVal[0])
          + fBFieldData[i][ gridIndex[0]+1 ][ gridIndex[1]+1 ][ gridIndex[2] ]*fracVal[0];
    c11 = fBFieldData[i][ gridIndex[0] ][ gridIndex[1]+1 ][ gridIndex[2]+1 ] * (1.0-fracVal[0])
          + fBFieldData[i][ gridIndex[0]+1 ][ gridIndex[1]+1 ][ gridIndex[2]+1 ]*fracVal[0];
    c0  = c00 * (1.0-fracVal[1]) + c10 * fracVal[1];
    c1  = c01 * (1.0-fracVal[1]) + c11 * fracVal[1];
    Bint_local[i] = c0 * (1.0-fracVal[2]) + c1 * fracVal[2];
  }

  //print checks - remove before pulling into anything else.
  if(printCheck) G4cout << "Interpolated BField at Rotated Point:" << G4endl;
  if(printCheck) G4cout << "   Bx: " << Bint_local[0] / tesla << "T" << G4endl;
  if(printCheck) G4cout << "   By: " << Bint_local[1] / tesla << "T" << G4endl;
  if(printCheck) G4cout << "   Bz: " << Bint_local[2] / tesla << "T" << G4endl;

  // Perform the proper rotation on the Bfield this gives us the global field from the local map position field.
  G4ThreeVector Bint_global = G4ThreeVector(
    solenoidRotation.rowX() * Bint_local,
    solenoidRotation.rowY() * Bint_local,
    solenoidRotation.rowZ() * Bint_local
  );

  //print checks - remove before pulling into anything else.
  if(printCheck) G4cout << "Rotated BField at Offset Point:" << G4endl;
  if(printCheck) G4cout << "  Bx': " << Bint_global.x() / tesla << "T" << G4endl;
  if(printCheck) G4cout << "  By': " << Bint_global.y() / tesla << "T" << G4endl;
  if(printCheck) G4cout << "  Bz': " << Bint_global.z() / tesla << "T" << G4endl;

  // Scale the field according to whatever fieldscale was set on constructing the solenoid class.
  Bfield[0] = Bint_global.x() * fFieldScale;
  Bfield[1] = Bint_global.y() * fFieldScale;
  Bfield[2] = Bint_global.z() * fFieldScale;

}
