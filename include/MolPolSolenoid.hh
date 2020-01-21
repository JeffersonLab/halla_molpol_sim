#ifndef MOLPOLSOLENOID_HH
#define MOLPOLSOLENOID_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

class MolPolSolenoid : public G4MagneticField
{
  public:
    MolPolSolenoid(G4double pFieldStrengthT, G4ThreeVector pCenter, G4double pXrot, G4double pYrot);
    virtual ~MolPolSolenoid();
    void UpdateSolenoid(G4double pFieldStrengthT, G4ThreeVector pCenter, G4double pXrot, G4double pYrot);
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
    G4double fYoffset;
    G4double fXoffset;
    G4double fYrot;
    G4double fXrot;
    G4int NX;
    G4int NY;
    G4int NZ;
    std::vector< std::vector< std::vector< G4double > > > fBFieldData[3];
    G4double fFieldScale;
    G4bool fInit;

};
#endif
