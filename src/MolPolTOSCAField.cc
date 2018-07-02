#include "MolPolTOSCAField.hh"
//#include "G4UImanager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <iostream>
#include <iomanip>
#include <fstream>

#include <assert.h>
#include <math.h>

// Boost headers
#ifdef __USE_BOOST_IOSTREAMS
// This supports gzipped iostreams as magnetic field maps.
// Compile with -D __USE_BOOST_IOSTREAMS to use.
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#endif

MolPolTOSCAField::MolPolTOSCAField( G4String filename, G4double scale, G4double offset ){
  fZoffset = offset;
  fFilename = filename;
  fInit = false;
  fFieldScale = scale;
  ReadFieldMap();
}

MolPolTOSCAField::~MolPolTOSCAField(){
}

G4String MolPolTOSCAField::GetName(){
  if( !fInit ){
    G4cerr << "WARNING " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ": access uninitialized field." << G4endl;
    return G4String("");
  }
  return fFilename;
}

void MolPolTOSCAField::SetFieldScale(G4double s){
  fFieldScale = s;
  G4cout << fFilename << " scale set to " << s << G4endl;
  return;
}

void MolPolTOSCAField::SetMagnetCurrent(G4double s){
  if( fMagCurrent0 > 0.0 ) SetFieldScale(s/fMagCurrent0);
  else G4cerr << "Warning:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ": Field current not specified in map " << fFilename << " - Ignoring and proceeding " << G4endl;
  return;
}

void MolPolTOSCAField::InitializeGrid() {

  if( NX <= 0 || NY <= 0 || NZ <= 0 ){
    G4cerr << "Error " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ": grid size invalid.  Aborting" << G4endl;
    exit(1);
  }

  G4cout << "Initializing field map grid for " << fFilename << G4endl;

  for( int i = 0; i < 3; i++ ){
    fBFieldData[i].clear();
    fBFieldData[i].resize(NX);
  	for( int x = 0; x < NX; x++) {
      fBFieldData[i][x].resize(NY);
      for( int y = 0; y < NY; y++) {
    		fBFieldData[i][x][y].resize(NZ);
    		for( int z = 0; z < NZ; z++) {
    	    fBFieldData[i][x][y][z] = 0.0;
          //G4cout << "Setting fBFieldData[" << i << "][" << x << "][" << y << "][" << z << "] to 0.0" << G4endl;
        } // end of z
      } // end of y
    } // end of x
  } // end B-vector directions

  G4cout << "Map grid for " << fFilename << " initialized" << G4endl;
  G4cout << "With scale: " << fFieldScale << " set" << G4endl;
  G4cout << "With Offset: " << fZoffset / cm << " cm" << G4endl;

  return;

}

void MolPolTOSCAField::ReadFieldMap(){
  G4cout << __PRETTY_FUNCTION__ << " : " << fFilename << G4endl;

  G4int nX = 0, nY = 0, nZ = 0; //Grid point counters
  G4int nlines = 0; // # of lines

// Used in remoll for boost.
#ifdef __USE_BOOST_IOSTREAMS
  // Create Boost istream
  boost::iostreams::filtering_istream inputfile;
  // If the filename has .gz somewhere (hopefully the end)
  if (fFilename.find(".gz") != std::string::npos) {
    // Add gzip decompressor to stream
    inputfile.push(boost::iostreams::gzip_decompressor());
  }
  // Set file as source
  inputfile.push(boost::iostreams::file_source(fFilename));
#else
  // Create STL ifstream
  std::ifstream inputfile;
  // If the filename has .gz somewhere, fail ungracefully
  if (fFilename.find(".gz") != std::string::npos) {
    G4cerr << "Compressed input files not supported!" << G4endl;
    exit(1);
  }
  // Set file as source
  inputfile.open(fFilename.data());
#endif

  if (!inputfile.good() ){
    G4cerr << "Error " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ") File " << fFilename << " could not open.  Aborting" << G4endl;
    exit(1);
  }

  // Variable that will contain single lines
  std::string inputline;

  // Read in data about grid
  getline(inputfile,inputline);
  if (std::istringstream(inputline) >> NZ >> NY >> NX >> DUMMY) {
    G4cout << "GridSize[" << NX << "," << NY << "," << NZ << "]" << G4endl;
  } else {
    G4cerr << "(" << __PRETTY_FUNCTION__ << " line " << __LINE__ << ") Error! File "
           << fFilename << " contains unreadable header information. Aborting field mapping." << G4endl;
    exit(1);
  }

  // CURRENTLY SKIP THE REMAINING LINES OF THE HEADER
  for(int i = 0; i < 8; i++){
    getline(inputfile,inputline);
  }

  //////////////////////////////////////////////////////////////////////////////

  InitializeGrid();

  std::vector<G4double> xpos;
  std::vector<G4double> ypos;
  std::vector<G4double> zpos;
  for(nX = 0; nX < NX; nX++){
    for(nY = 0; nY < NY; nY++){
      for(nZ = 0; nZ < NZ; nZ++){
        getline(inputfile,inputline);

        // Read in field values and assign units
        if (std::istringstream(inputline) >> raw_xpos >> raw_ypos >> raw_zpos >> raw_B >> raw_Bx >> raw_By >> raw_Bz) {
          nlines++;
        } else {
          G4cerr << "Error " << __PRETTY_FUNCTION__ << " line " << __LINE__ << ") File " << fFilename << " contains invalid data.  Aborting" << G4endl;
          exit(1);
        }
        xpos.push_back(raw_xpos*cm);
        ypos.push_back(raw_ypos*cm);
        zpos.push_back(raw_zpos*cm);

        //G4cout << "Raw read in for field " << raw_Bz << G4endl;

        // FIXME: Double check this.
        // Convert Oersted to Teslas: Oe * 1e-4 ~ T
        raw_Bx *= 0.0001;
        raw_By *= 0.0001;
        raw_Bz *= 0.0001;

    		// Set the grid values to the values which have been read-in
    		fBFieldData[0][nX][nY][nZ] = raw_Bx * tesla;
    		fBFieldData[1][nX][nY][nZ] = raw_By * tesla;
    		fBFieldData[2][nX][nY][nZ] = raw_Bz * tesla;

        //FIXME: Delete me after bug checking... just sampling read in values.
        //if( nX%20==0 && nY%20==0 && nZ%20 == 0 ) G4cout << " Read-in B-Field(" << fBFieldData[0][nX][nY][nZ] / tesla << "," << fBFieldData[1][nX][nY][nZ] / tesla << "," << fBFieldData[2][nX][nY][nZ] / tesla << ")" << G4endl;

        //G4cout << "Read in for field " << raw_Bz / tesla << " T" << G4endl << G4endl;
      }
    }
  }

  XMIN = *std::min_element(xpos.begin(),xpos.end());
  XMAX = *std::max_element(xpos.begin(),xpos.end());
  YMIN = *std::min_element(ypos.begin(),ypos.end());
  ZMIN = *std::min_element(zpos.begin(),zpos.end());
  YMAX = *std::max_element(ypos.begin(),ypos.end());
  ZMAX = *std::max_element(zpos.begin(),zpos.end());

  //FIXME: CAN DELETE AFTER BUG CHECKING
  //G4cout << "ReadFieldMap() ... " << G4endl;
  //G4cout << "  NX: " << NX << G4endl;
  //G4cout << "  NY: " << NY << G4endl;
  //G4cout << "  NZ: " << NZ << G4endl;
  //G4cout << "FILE: " << fFilename << G4endl;
  //G4cout << "XMIN: " << XMIN / cm << G4endl;
  //G4cout << "XMAX: " << XMAX / cm << G4endl;
  //G4cout << "YMIN: " << YMIN / cm << G4endl;
  //G4cout << "YMAX: " << YMAX / cm << G4endl;
  //G4cout << "ZMIN: " << ZMIN / cm << G4endl;
  //G4cout << "ZMAX: " << ZMAX / cm << G4endl; // Looks like these produce just fine.


  fInit = true;
  G4cout << "... done reading " << nlines << " lines." << G4endl<< G4endl;

  // The Ultimate Style Guide
  // https://twitter.com/ThePracticalDev/status/710156980535558144

}

void MolPolTOSCAField::GetFieldValue(const G4double Point[4], G4double *Bfield ) const {

  //G4cout << G4endl << "Map received global point: (" << Point[0] / cm << "," << Point[1] / cm << "," << Point[2] / cm << "," << Point[3] / cm << ")." << G4endl;

  //FIXME: Detlet after debugging
  //G4cout << "TOSCA Received Point(" <<Point[0] / cm <<","<<Point[1] / cm <<","<<Point[2] / cm <<")"<<G4endl;

  //G4cout << "zOffset for this map is: " << fZoffset / cm << " cm " << G4endl;

  G4double    x, y, z;
  G4double    fracVal[3], intVal[3];
  G4int       gridIndex[3]; // THESE ARE THE GRID LOCATIONS

  x = Point[0];
  y = Point[1];
  z = Point[2] - fZoffset; // What point are we looking at in z (in terms of the TOSCA Map Coordinates)

  // Check that the point is within the defined region before interpolation.  If it is outside then the field is zero
  if( x < XMIN || x >= XMAX || y >= YMAX || y < YMIN || z >= ZMAX || z < ZMIN ) {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
    //FIXME: Delete after buh checking
    //G4cout << "Sending field: (" << Bfield[0] << "," << Bfield[1] << "," << Bfield[2] << ")" << G4endl;
    return;
  } else {
    //G4cout << ">>>>>>>>>>>>>>> point exists in TOSCA field map..." << G4endl;
  }

  //G4cout << "Which is map local point: (" << Point[0] / cm << "," << Point[1] / cm << "," << Point[2] / cm << "," << Point[3] / cm << ")." << G4endl;

  // Ensure we're going to get our grid indices correct
  assert( XMIN <= x && x < XMAX );
  //G4cout << "Checking that local point x = " << x / cm << " falls between " << XMIN / cm << " and " << XMAX / cm << G4endl;
  assert( YMIN <= y && y < YMAX );
  //G4cout << "Checking that local point y = " << y / cm << " falls between " << YMIN / cm << " and " << YMAX / cm << G4endl;;
  assert( ZMIN <= z && z < ZMAX );
  //G4cout << "Checking that local point z = " << z / cm << " falls between " << ZMIN / cm << " and " << ZMAX / cm << G4endl;;

  // Get interpolation variables; modf finds integer and fractional parts of input.
  // Also: the Nx-1 here is fencepost problem as noted in remoll code.
  fracVal[0] = modf( ( x - XMIN )*(NX-1)/( XMAX - XMIN ), &(intVal[0]) );
  //G4cout << "modf( ( " << x / cm << " - " << XMIN / cm << " )*(" << NX << "-1)/( " << XMAX / cm << " - " << XMIN / cm << " ), &(intval[0]) );" << G4endl;
  //G4cout << "Point x has fracVal: " << fracVal[0] << " , and intVal: " << intVal[0] << G4endl;
  fracVal[1] = modf( ( y - YMIN )*(NY-1)/( YMAX - YMIN ), &(intVal[1]) );
  //G4cout << "modf( ( " << y / cm << " - " << YMIN / cm << " )*(" << NY << "-1)/( " << YMAX / cm << " - " << YMIN / cm << " ), &(intval[1]) );" << G4endl;
  //G4cout << "Point y has fracVal: " << fracVal[1] << " , and intVal: " << intVal[1] << G4endl;
  fracVal[2] = modf( ( z - ZMIN )*(NZ-1)/( ZMAX - ZMIN ), &(intVal[2]) );
  //G4cout << "modf( ( " << z / cm << " - " << ZMIN nX/ cm << " )*(" << NZ << "-1)/( " << ZMAX / cm << " - " << ZMIN / cm << " ), &(intval[2]) );" << G4endl;
  //G4cout << "Point z has fracVal: " << fracVal[2] << " , and intVal: " << intVal[2] << G4endl;

  // Cast these to integers for indexing and check
  for( G4int i = 0; i < 3; i++ ) gridIndex[i] = (G4int) intVal[i];

  assert( 0 <= gridIndex[0] && gridIndex[0] < NX );
  assert( 0 <= gridIndex[1] && gridIndex[1] < NY );
  assert( 0 <= gridIndex[2] && gridIndex[2] < NZ );

  // NOTE: THIS IS IN THE REMOLL CODE BUT I AM UNSURE WHY IT IS REPEATED ...
  // Interpolate
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
    //FIXME: DELETE AFTER BUG CHECKING
    //G4cout << "TOSCA Point[" << i << "]:" << Point[i] << G4endl;
    //G4cout << " TOCSA Bint[" << i << "]:" << Bint[i] << G4endl;
  }

  Bfield[0] = Bint[0] * fFieldScale;
  Bfield[1] = Bint[1] * fFieldScale;
  Bfield[2] = Bint[2] * fFieldScale;
  //FIXME: Delete after bug checking
  //G4cout << std::setprecision(6) << std::scientific << "TOSCA Sending B-field: ( " << Bfield[0] /tesla << " , " << Bfield[1] / tesla << " , " << Bfield[2] / tesla << " ) Tesla"
  //       << std::setprecision(6) << std::fixed << " at Point (" << Point[0] / cm << "," << Point[1] /cm << "," << Point[2] / cm << ") cm" << G4endl;
}
