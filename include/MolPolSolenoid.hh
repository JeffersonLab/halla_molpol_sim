#ifndef MOLPOLSOLENOID_HH
#define MOLPOLSOLENOID_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

class MolPolSolenoid : public G4MagneticField
{
  public: // with description
    MolPolSolenoid(G4double Bz, G4double fz, G4ThreeVector pOrigin);
    ~MolPolSolenoid();
    
    void GetFieldValue(const G4double yTrack[], G4double B[] ) const;
    void UpdateSolenoid(G4double Bz, G4double fz, G4ThreeVector pOrigin);

  private:
    G4double          fBfield;
    G4double          fFringeZ;
    G4double          fFieldLength;
    G4double          fFieldRadius;
    G4ThreeVector     fOrigin;

    G4bool IsOutside(G4ThreeVector& local) const;
    G4bool IsWithin(G4ThreeVector& local) const;

};
#endif
