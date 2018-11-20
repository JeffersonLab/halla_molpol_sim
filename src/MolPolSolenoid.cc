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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>

MolPolSolenoid::MolPolSolenoid(G4double desiredFieldStrengthT){
  fZoffset = 6.9*cm;
  fFilename = "../solenoid.map";
  fInit = false;
  fFieldScale = desiredFieldStrengthT / 5.0;
  ReadFieldMap();
}

void MolPolSolenoid::UpdateSolenoid(G4double desiredFieldStrengthT){
  fFieldScale = desiredFieldStrengthT / 5.0;
  G4cout << "  New solenoid desired field strength: " << desiredFieldStrengthT << G4endl;
  G4cout << "  Scale factor required for 5.0 tesla map: " << fFieldScale << G4endl;
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
  G4double    x, y, z, fracVal[3], intVal[3];
  G4int       gridIndex[3];

  x = Point[0];
  y = Point[1];
  z = Point[2] - fZoffset;

  // Check that the point is within the defined region before interpolation.  If it is outside then the field is zero
  if( x < XMIN || x >= XMAX || y >= YMAX || y < YMIN || z >= ZMAX || z < ZMIN ) {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
    return;
  }

  // Ensures we're going to get our grid indices correct
  assert( XMIN <= x && x < XMAX );
  assert( YMIN <= y && y < YMAX );
  assert( ZMIN <= z && z < ZMAX );

  fracVal[0] = modf( ( x - XMIN )*(NX-1)/( XMAX - XMIN ), &(intVal[0]) );
  fracVal[1] = modf( ( y - YMIN )*(NY-1)/( YMAX - YMIN ), &(intVal[1]) );
  fracVal[2] = modf( ( z - ZMIN )*(NZ-1)/( ZMAX - ZMIN ), &(intVal[2]) );

  for( G4int i = 0; i < 3; i++ ) gridIndex[i] = (G4int) intVal[i];
  assert( 0 <= gridIndex[0] && gridIndex[0] < NX );
  assert( 0 <= gridIndex[1] && gridIndex[1] < NY );
  assert( 0 <= gridIndex[2] && gridIndex[2] < NZ );
  for( G4int i = 0; i < 3; i++ ) gridIndex[i] = (G4int) intVal[i];

  G4double Bint[3] = {0,0,0};
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
    Bint[i] = c0 * (1.0-fracVal[2]) + c1 * fracVal[2];
  }

  Bfield[0] = Bint[0] * fFieldScale;
  Bfield[1] = Bint[1] * fFieldScale;
  Bfield[2] = Bint[2] * fFieldScale;
}
