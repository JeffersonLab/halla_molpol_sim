
#ifndef MolPolDetectorConstruction_h
#define MolPolDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MolPolEMFieldSetup;
class G4Material;

class MolPolDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MolPolDetectorConstruction();
   ~MolPolDetectorConstruction();

   public:
    G4VPhysicalVolume* Construct();

  G4Material* GetTarget() {return fTargetMaterial;}
  G4double GetTargetFullLength() {return fTargetFullLength;}

  private:

  MolPolEMFieldSetup* mEMFieldSetup;

  G4bool fCheckOverlaps;
  G4Material* fTargetMaterial;
  G4double fTargetFullLength;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*MolPolDetectorConstruction_h*/
