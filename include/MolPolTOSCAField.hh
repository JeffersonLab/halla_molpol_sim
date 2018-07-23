/* MolPolToscaField adopted from remoll */

#ifndef MolPolTOSCAField_H
#define MolPolTOSCAField_H

#include <vector>

#include "G4String.hh"
#include "G4MagneticField.hh"

class MolPolTOSCAField : public G4MagneticField {
  public:
    MolPolTOSCAField( G4String filename, G4double scale, G4double offset );
    virtual ~MolPolTOSCAField();
    void GetFieldValue( const   G4double Point[4], G4double *Bfield ) const;
    void InitializeGrid();
    void ReadFieldMap();
    void SetFieldScale(G4double s);
    void SetMagnetCurrent(G4double s);
    void SetZoffset(G4double z){ fZoffset = z; }
    G4bool IsInit(){ return fInit; }

  private:

    G4String fFilename;
    G4double XMIN;
    G4double XMAX;
    G4double YMIN;
    G4double YMAX;
    G4double ZMIN;
    G4double ZMAX;
    G4double fZoffset;
    G4int NX;
    G4int NY;
    G4int NZ;
    std::vector< std::vector< std::vector< G4double > > > fBFieldData[3];
    G4double fFieldScale; // Scale overall field by this amount
    G4bool fInit;
};

#endif
