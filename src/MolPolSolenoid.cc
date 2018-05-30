#include "MolPolSolenoid.hh"
#include <math.h>
#include "G4ios.hh"

MolPolSolenoid::MolPolSolenoid(G4double Bz, G4double fz, G4ThreeVector pOrigin)
{
  fBfield = Bz;
  fFringeZ = fz;
  fOrigin = pOrigin;
  fFieldLength = 38.1 * cm;
  fFieldRadius = 25.4 * cm;
}

void MolPolSolenoid::UpdateSolenoid(G4double Bz, G4double fz, G4ThreeVector pOrigin)
{
  fBfield = Bz;
  fFringeZ = fz;
  fOrigin = pOrigin;
  fFieldLength = 38.1 * cm;
  fFieldRadius = 25.4 * cm;
}

MolPolSolenoid::~MolPolSolenoid()
{
}

void MolPolSolenoid::GetFieldValue(const G4double y[7], G4double B[3]) const
{
  G4ThreeVector global(y[0], y[1], y[2]);
  G4ThreeVector local
  (y[0] - fOrigin.x(),
   y[1] - fOrigin.y(),
   y[2] - fOrigin.z());

  if (IsOutside(local)) return;

  //  G4ThreeVector B(0.0,0.0,fBfield);
  //  B = fGlobal2local.Inverse().TransformAxis(B);

  B[0] = 0.;
  B[1] = 0.;
  B[2] = fBfield;

  //  field[0] += B[0];
  //  field[1] += B[1];
  //  field[2] += B[2];
}


G4bool MolPolSolenoid::IsOutside(G4ThreeVector& local) const
{
  //  EInside inside = tubs->Inside(local);
  //  return (inside == kOutside);
  G4double r = std::sqrt(local.x()*local.x()+local.y()*local.y());
  return (r > fFieldRadius || std::fabs(local.z()) > fFieldLength);
}


G4bool MolPolSolenoid::IsWithin(G4ThreeVector& local) const
{
  //  EInside inside = tubs->Inside(local);
  //  return (inside == kInside);
  G4double r = std::sqrt(local.x()*local.x()+local.y()*local.y());
  return (r < fFieldRadius && std::fabs(local.z()) < fFieldLength/2.0);
}
