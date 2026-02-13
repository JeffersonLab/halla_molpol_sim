
#ifndef MolPolDetectorConstruction_h
#define MolPolDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class MolPolEMFieldSetup;
class G4Material;
class G4VPhysicalVolume;
class G4GenericMessenger;
class G4Tubs;

class MolPolDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MolPolDetectorConstruction();
   ~MolPolDetectorConstruction();

   public:
    G4VPhysicalVolume* Construct();
    void ConstructMaterials();

    void SetDipolePbJawsGap(G4double val);
    void SetTargetThickness(G4double val);
    void SetTargetZPosition(G4double val);
    void SetDPVPSensitiveDetectors();
    void SetDipoleInternalVPlanes(G4String val);

    void UpdateGeometry();

    G4Material * GetTarget()           {return fTargetMaterial;}
    G4double     GetTargetFullLength() {return fTargetFullLength;}

  private:
    void DefineGeometryCommands();

    MolPolEMFieldSetup*   mEMFieldSetup;

    G4bool                fCheckOverlaps;

    G4Material *          fTargetMaterial;     //Used for Sanghwa's update to PrimaryGenerator

    G4GenericMessenger *  fMessenger;

    G4bool                fEnableDipoleInternalVPlanes;

    G4double              fLeadJawGapWidth;
    G4double              fLeadJawsHLength;
    G4double              fLeadJawsZOrigin;
    G4double              fLeadJawsYOrigin;
    G4double              fLeadJawsXOrigin;
    G4VPhysicalVolume *   fLeadJawsPhysicalLT;
    G4VPhysicalVolume *   fLeadJawsPhysicalLB;
    G4VPhysicalVolume *   fLeadJawsPhysicalRT;
    G4VPhysicalVolume *   fLeadJawsPhysicalRB;

    G4LogicalVolume *DP0L = nullptr, *DP0R = nullptr,
                    *DP1L = nullptr, *DP1R = nullptr,
                    *DP2L = nullptr, *DP2R = nullptr,
                    *DP3L = nullptr, *DP3R = nullptr,
                    *DP4L = nullptr, *DP4R = nullptr,
                    *DP5L = nullptr, *DP5R = nullptr,
                    *DP6L = nullptr, *DP6R = nullptr,
                    *DP7L = nullptr, *DP7R = nullptr,
                    *DP8L = nullptr, *DP8R = nullptr,
                    *DP9L = nullptr, *DP9R = nullptr;

    G4VPhysicalVolume *   fTargetPhysVolume;
    G4double              fTargetFullLength;   //Used for Sanghwa's update to PrimaryGenerator and DetCon's Generic Messenger
    G4double              fTargetFullRadius;
    G4double              fTargetBeamlinePx;
    G4double              fTargetBeamlinePy;
    G4double              fTargetBeamlinePz;
    G4double              fHelmCoilBeamPosZ;   //Position for solenoid center, required for proper placing of the target.
    G4Tubs *              fTargetVSolidTubs;   //Solid volume for target
};

#endif
