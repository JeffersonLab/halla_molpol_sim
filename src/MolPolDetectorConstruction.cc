#include "MolPolDetectorConstruction.hh"
#include "MolPolEMFieldSetup.hh"
#include "MolPolDetector.hh"

#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

MolPolDetectorConstruction::MolPolDetectorConstruction():
  mEMFieldSetup(0),
  fCheckOverlaps(true),        //Check for overlaps, always.
  fTargetMaterial(0),
  fMessenger(nullptr),         //Prob not necessary since assigned in DefGeoCom() which is also run.
  fEnableDipoleInternalVPlanes(false),  //Default: disabled (planes created but not sensitive)
  fEnableFluxVPlanes(false),   //Default: disabled (planes created but not sensitive)
  fLeadJawGapWidth(3.6*cm),    //Default jaw gap fully open.
  fTargetFullLength(0.013*mm), //Default target width 13 microns.
  fTargetFullRadius(15.0*mm),
  fTargetBeamlinePx( 0.0*mm),
  fTargetBeamlinePy( 0.0*mm),
  fTargetBeamlinePz(67.4*mm)   //Default target-center position on beamline.
{
  DefineGeometryCommands();
}

MolPolDetectorConstruction::~MolPolDetectorConstruction(){
  if(mEMFieldSetup) delete mEMFieldSetup;
  if(fMessenger) delete fMessenger;
}

G4VPhysicalVolume* MolPolDetectorConstruction::Construct() {

  ConstructMaterials();

  G4cout << " fLeadJawGapWidth: " << fLeadJawGapWidth << G4endl;
  G4cout << "fTargetFullLength: " << fTargetFullLength << G4endl;
  G4cout << "fTargetBeamlinePz: " << fTargetBeamlinePz << G4endl;

  auto MolPol_Air       = G4Material::GetMaterial("MP_Air");
  auto MolPol_Aluminum  = G4Material::GetMaterial("MP_Aluminum");
  auto MolPol_Copper    = G4Material::GetMaterial("MP_Copper");
  auto MolPol_Iron      = G4Material::GetMaterial("MP_Iron");
  auto MolPol_Lead      = G4Material::GetMaterial("MP_Lead");
  auto MolPol_Scint     = G4Material::GetMaterial("MP_Scint");
  auto MolPol_SiliSteel = G4Material::GetMaterial("MP_SiliconSteel");
  auto MolPol_Stainless = G4Material::GetMaterial("MP_StainlessSteel304");
  auto MolPol_Titanium  = G4Material::GetMaterial("MP_Titanium");
  auto MolPol_Vacuum    = G4Material::GetMaterial("MP_Vacuum");


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  //// Define Visual Attributes
  G4double alphaVacuum = 0.15;
  G4double alphaMatStd = 0.50;
  G4double alphaTarget = 0.85;
  G4VisAttributes* IronVisAtt = new G4VisAttributes( G4Colour( 10./255., 10./255.,10./255.,alphaTarget) );
  G4VisAttributes* PbVisAtt   = new G4VisAttributes( G4Colour(149./255.,149./255.,100./255.,alphaMatStd) );//Use for wireframed lead
  G4VisAttributes* LeadVisAtt = new G4VisAttributes( G4Colour(149./255.,149./255.,100./255.,alphaMatStd) );
  G4VisAttributes* SteelVisAtt= new G4VisAttributes( G4Colour(  0./255., 80./255.,225./255.,alphaMatStd) );
  G4VisAttributes* AlumVisAtt = new G4VisAttributes( G4Colour(  0./255.,237./255.,  0./255.,alphaMatStd) );
  G4VisAttributes* BPVisAtt   = new G4VisAttributes( G4Colour(175./255.,175./255.,175./255.,alphaMatStd) );
  G4VisAttributes* VacVisAtt  = new G4VisAttributes( G4Colour(255./255.,255./255.,255./255.,alphaVacuum) );
  G4VisAttributes* CuVisAtt   = new G4VisAttributes( G4Colour(178./255.,102./255., 26./255.,alphaMatStd) );
  G4VisAttributes* ScintVisAtt= new G4VisAttributes( G4Colour(  0./255.,100./255.,100./255.,alphaMatStd) );
  G4VisAttributes* DipVisAtt  = new G4VisAttributes( G4Colour(  0./255., 80./255.,225./255.,alphaVacuum) );


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Build world
  G4double world_x = 10*m;  G4double world_y = 10*m;  G4double world_z = 10*m;
  G4Box* world_box = new G4Box("World",world_x,world_y,world_z);
  G4LogicalVolume* world_log = new G4LogicalVolume(world_box,MolPol_Vacuum,"World",0,0,0);
  world_log->SetVisAttributes(G4VisAttributes::GetInvisible() );
  G4VPhysicalVolume* world_phys = new G4PVPlacement(0,G4ThreeVector(),world_log,"World",0,false,fCheckOverlaps);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Field Setup
  mEMFieldSetup = new MolPolEMFieldSetup();  //setup the field,
  G4FieldManager* Q1FieldManager = mEMFieldSetup->GetFieldManagerFZB1();
  G4FieldManager* Q2FieldManager = mEMFieldSetup->GetFieldManagerFZB2();
  G4FieldManager* Q3FieldManager = mEMFieldSetup->GetFieldManagerFZB3();
  G4FieldManager* Q4FieldManager = mEMFieldSetup->GetFieldManagerFZB4();
  G4FieldManager* DFieldManager  = mEMFieldSetup->GetFieldManagerFZB5();
  G4FieldManager* Q6FieldManager = mEMFieldSetup->GetFieldManagerFZB6();
  G4bool allLocal = true;


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // HELMHOLTZ COIL Magnetic Volume
  G4double pQ6Rin  =  0    * cm;  G4double pQ6Rout =  25.4 * cm;  G4double pQ6HL   = 49.53 * cm;  fHelmCoilBeamPosZ = 6.9 * cm;//49.53*cm is half-length of solenoid map
  G4VSolid* Q6MagSolid = new G4Tubs( "Q6MagTubs", pQ6Rin, pQ6Rout, pQ6HL, 0.0, 360.0 * deg);
  G4LogicalVolume* Q6MagLogical = new G4LogicalVolume(Q6MagSolid, MolPol_Vacuum, "Q6Mag", 0,0,0);
  Q6MagLogical->SetFieldManager(Q6FieldManager, allLocal);
  Q6MagLogical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0, G4ThreeVector(0, 0, fHelmCoilBeamPosZ), Q6MagLogical, "SolenoidMag", world_log, 0, 0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Target Beam pipe - Alternative setup inside solenoid magnetic volume only.
  G4double pBPRin  = 4.78 * cm;
  G4double pBPRout = 5.08 * cm;
  G4VSolid * BPITAlum = new G4Tubs( "BPITAlum", 0., pBPRout, pQ6HL, 0.0 * deg, 360.0 * deg );
  G4VSolid * BPITVac  = new G4Tubs( "BPITVacm", 0., pBPRin,  pQ6HL, 0.0 * deg, 360.0 * deg );
  G4LogicalVolume * BPITLogical = new G4LogicalVolume(BPITAlum, MolPol_Aluminum, "BPalum_Targ", 0, 0, 0);
  BPITLogical->SetVisAttributes(BPVisAtt);
  G4LogicalVolume * BPVTLogical = new G4LogicalVolume(BPITVac,  MolPol_Vacuum,   "BPvacm_Targ",  0, 0, 0);
  BPVTLogical->SetVisAttributes(BPVisAtt);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),BPITLogical, "BeamPipeAlmnum_Targ", Q6MagLogical, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),BPVTLogical, "BeamPipeVacuum_Targ", BPITLogical, 0, 0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Target
  //G4double pMTATRin   = 0.0 * cm; G4double pMTATRout  = 1.5 * cm; G4double pMTATHLZ = 0.0062 * mm;/*Moved to method objects since modifiable*/
  //G4double pMTATPos_X = 0.0 * cm; G4double pMTATPos_Y = 0.0 * cm; G4double pMTATPos_Z = 6.9 * cm;/*Moved to method objects since modifiable*/
  fTargetVSolidTubs = new G4Tubs( "MTATTube", 0.0*cm, fTargetFullRadius, (fTargetFullLength / 2.0), 0.0*deg, 360.0 * deg );

  G4LogicalVolume* TargetLogical = new G4LogicalVolume(fTargetVSolidTubs, MolPol_Iron, "Target", 0, 0, 0);
  TargetLogical->SetVisAttributes(IronVisAtt);

  fTargetPhysVolume = new G4PVPlacement(0, G4ThreeVector(fTargetBeamlinePx,fTargetBeamlinePy,fTargetBeamlinePz - fHelmCoilBeamPosZ), TargetLogical, "Target", BPVTLogical, 0, 0, fCheckOverlaps);

  fTargetMaterial   = TargetLogical->GetMaterial();

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Helmholtz Coil 'Physical' Volume
  G4double pHLMZRin    = 5.10 *cm;   G4double pHLMZRout  = 15.0 * cm;   G4double pHLMZHLZ    = 5.0 * cm;
  G4double pHLMZ1Pos_X = 0.0 * cm;   G4double pHLMZ1Pos_Y = 0.0 * cm;   G4double pHLMZ1Pos_Z =-5.7 * cm;
  G4double pHLMZ2Pos_X = 0.0 * cm;   G4double pHLMZ2Pos_Y = 0.0 * cm;   G4double pHLMZ2Pos_Z =19.5 * cm;

  G4VSolid* HLMZ1Solid = new G4Tubs( "HLMZ1Tube", pHLMZRin, pHLMZRout, pHLMZHLZ, 0.0, 360.0 * deg );
  G4VSolid* HLMZ2Solid = new G4Tubs( "HLMZ2Tube", pHLMZRin, pHLMZRout, pHLMZHLZ, 0.0, 360.0 * deg );

  G4LogicalVolume* HLMZ1Logical = new G4LogicalVolume(HLMZ1Solid, MolPol_Copper, "Helmholtz1", 0, 0, 0);
  G4LogicalVolume* HLMZ2Logical = new G4LogicalVolume(HLMZ2Solid, MolPol_Copper, "Helmholtz2", 0, 0, 0);
  HLMZ1Logical->SetVisAttributes(CuVisAtt);
  HLMZ2Logical->SetVisAttributes(CuVisAtt);

  new G4PVPlacement(0, G4ThreeVector(pHLMZ1Pos_X, pHLMZ1Pos_Y, pHLMZ1Pos_Z - fHelmCoilBeamPosZ), HLMZ1Logical, "Helmholtz1", Q6MagLogical, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(pHLMZ2Pos_X, pHLMZ2Pos_Y, pHLMZ2Pos_Z - fHelmCoilBeamPosZ), HLMZ2Logical, "Helmholtz2", Q6MagLogical, 0, 0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // QUADRUPOLE MAGNETS
  G4double pQ1Rin  =  5.08 * cm;  G4double pQ1Rout = 20.00 * cm;  G4double pQ1HL   = 18.29 * cm;  G4double pQ1Pos_Z= 75.19 * cm;
  G4double pQ2Rin  =  5.08 * cm;  G4double pQ2Rout = 20.00 * cm;  G4double pQ2HL   = 22.30 * cm;  G4double pQ2Pos_Z=140.46 * cm;
  G4double pQ3Rin  =  5.08 * cm;  G4double pQ3Rout = 20.00 * cm;  G4double pQ3HL   = 18.37 * cm;  G4double pQ3Pos_Z=209.08 * cm;
  G4double pQ4Rin  =  5.08 * cm;  G4double pQ4Rout = 20.00 * cm;  G4double pQ4HL   = 18.37 * cm;  G4double pQ4Pos_Z=274.59 * cm;

  G4VSolid* Q1Solid = new G4Tubs( "Q1Tubs", pQ1Rin, pQ1Rout, pQ1HL, 0.0, 360.0 * deg );
  G4VSolid* Q2Solid = new G4Tubs( "Q2Tubs", pQ2Rin, pQ2Rout, pQ2HL, 0.0, 360.0 * deg );
  G4VSolid* Q3Solid = new G4Tubs( "Q3Tubs", pQ3Rin, pQ3Rout, pQ3HL, 0.0, 360.0 * deg );
  G4VSolid* Q4Solid = new G4Tubs( "Q4Tubs", pQ4Rin, pQ4Rout, pQ4HL, 0.0, 360.0 * deg );
  G4VSolid* Q1MagSolid = new G4Tubs( "Q1MagTubs", 0., pQ1Rin, pQ1HL, 0.0, 360.0 * deg );
  G4VSolid* Q2MagSolid = new G4Tubs( "Q2MagTubs", 0., pQ2Rin, pQ2HL, 0.0, 360.0 * deg );
  G4VSolid* Q3MagSolid = new G4Tubs( "Q3MagTubs", 0., pQ3Rin, pQ3HL, 0.0, 360.0 * deg );
  G4VSolid* Q4MagSolid = new G4Tubs( "Q4MagTubs", 0., pQ4Rin, pQ4HL, 0.0, 360.0 * deg );

  G4LogicalVolume* Q1Logical = new G4LogicalVolume(Q1Solid,MolPol_SiliSteel,"Q1Logical",0,0,0);
  G4LogicalVolume* Q2Logical = new G4LogicalVolume(Q2Solid,MolPol_SiliSteel,"Q2Logical",0,0,0);
  G4LogicalVolume* Q3Logical = new G4LogicalVolume(Q3Solid,MolPol_SiliSteel,"Q3Logical",0,0,0);
  G4LogicalVolume* Q4Logical = new G4LogicalVolume(Q4Solid,MolPol_SiliSteel,"Q4Logical",0,0,0);
  G4LogicalVolume* Q1MagLogical = new G4LogicalVolume(Q1MagSolid,MolPol_Vacuum,"Q1MagLogical",0,0,0);
  G4LogicalVolume* Q2MagLogical = new G4LogicalVolume(Q2MagSolid,MolPol_Vacuum,"Q2MagLogical",0,0,0);
  G4LogicalVolume* Q3MagLogical = new G4LogicalVolume(Q3MagSolid,MolPol_Vacuum,"Q3MagLogical",0,0,0);
  G4LogicalVolume* Q4MagLogical = new G4LogicalVolume(Q4MagSolid,MolPol_Vacuum,"Q4MagLogical",0,0,0);

  Q1Logical->SetVisAttributes(SteelVisAtt);
  Q2Logical->SetVisAttributes(SteelVisAtt);
  Q3Logical->SetVisAttributes(SteelVisAtt);
  Q4Logical->SetVisAttributes(SteelVisAtt);

  Q1MagLogical->SetVisAttributes(VacVisAtt);
  Q1MagLogical->SetFieldManager(Q1FieldManager,allLocal);

  Q2MagLogical->SetVisAttributes(VacVisAtt);
  Q2MagLogical->SetFieldManager(Q2FieldManager,allLocal);

  Q3MagLogical->SetVisAttributes(VacVisAtt);
  Q3MagLogical->SetFieldManager(Q3FieldManager,allLocal);

  Q4MagLogical->SetVisAttributes(VacVisAtt);
  Q4MagLogical->SetFieldManager(Q4FieldManager,allLocal);

  new G4PVPlacement(0,G4ThreeVector(0,0,pQ1Pos_Z),Q1Logical,"Q1Phys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ2Pos_Z),Q2Logical,"Q2Phys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ3Pos_Z),Q3Logical,"Q3Phys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ4Pos_Z),Q4Logical,"Q4Phys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ1Pos_Z),Q1MagLogical,"Q1MagPhys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ2Pos_Z),Q2MagLogical,"Q2MagPhys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ3Pos_Z),Q3MagLogical,"Q3MagPhys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ4Pos_Z),Q4MagLogical,"Q4MagPhys",world_log,0,0,fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Dipole magnetic field volume
  G4double pDMagHLX   =  8.0   * cm;  G4double pDMagHLY   = 30.0 * cm;  G4double pDMagHLZ   = 82.25 * cm;
  G4double pDMagPos_X =  0.0   * cm;  G4double pDMagPos_Y =  0.0 * cm;  G4double pDMagPos_Z =  423.4 * cm;

  G4VSolid* DMagSolid = new G4Box ( "DMagBox" , pDMagHLX , pDMagHLY , pDMagHLZ );

  /*NEEDS SUB4 VOLUME SUBTRACTED FROM THIS SPACE
   *G4LogicalVolume* DLogical = new G4LogicalVolume ( DMagSolid, MolPol_Vacuum, "DipoleMag", 0, 0, 0);
   *DLogical->SetFieldManager(DFieldManager,allLocal);
   *DLogical->SetVisAttributes(VacVisAtt);
   *new G4PVPlacement(0,G4ThreeVector(pDMagPos_X, pDMagPos_Y - 9*cm, pDMagPos_Z), DLogical,"DipoleMag",world_log,0,0,fCheckOverlaps);*/


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DIPOLE BOX
  G4double pDBI1HLX =  6.0   * cm;  G4double pDBI1HLY = 16.5 * cm;  G4double pDBI1HLZ = 98.5  * cm;
  G4double pDBV1HLX =  5.295 * cm;  G4double pDBV1HLY = 16.0 * cm;  G4double pDBV1HLZ = 98.5  * cm;
  G4double pDBW1HLX = 11.66  * cm;  G4double pDBW1HLY = 21.4 * cm;  G4double pDBW1HLZ =  0.64 * cm;
  G4double pDBW2Rin =  0.0   * cm;  G4double pDBW2Rout=  3.0 * cm;  G4double pDBW2HL  =  0.75 * cm;
  G4double pDBW3HLX =  1.18  * cm;  G4double pDBW3HLY =  8.0 * cm;  G4double pDBW3HLZ =  0.75 * cm;
  G4double pDBW4HLX =  1.18  * cm;  G4double pDBW4HLY =  8.0 * cm;  G4double pDBW4HLZ =  0.75 * cm;
  G4double pDBT3HLX =  1.18  * cm;  G4double pDBT3HLY =  8.0 * cm;  G4double pDBT3HLZ =  0.005* cm;
  G4double pDBT4HLX =  1.18  * cm;  G4double pDBT4HLY =  8.0 * cm;  G4double pDBT4HLZ =  0.005* cm;
  G4double pDMS1HLX =  3.0   * cm;  G4double pDMS1HLY = 15.0 * cm;  G4double pDMS1HLZ = 96.2  * cm;
  G4double pDMBHRin =  0.0   * cm;  G4double pDMBHRout=  2.0 * cm;  G4double pDMBHHL  = 96.2  * cm;
  G4double pDMFRRin =  1.27  * cm;  G4double pDMFRRout=  2.0 * cm;  G4double pDMFRHL  =  9.0  * cm;

  G4double pDBI1Pos_X =  0.0   * cm;  G4double pDBI1Pos_Y = -9.0 * cm;  G4double pDBI1Pos_Z =422.8  * cm;
  G4double pDBV1Pos_X =  0.0   * cm;  G4double pDBV1Pos_Y =  0.0 * cm;  G4double pDBV1Pos_Z =  0.0  * cm;
  G4double pDBW1Pos_X =  0.0   * cm;  G4double pDBW1Pos_Y = -9.0 * cm;  G4double pDBW1Pos_Z =521.94 * cm;
  G4double pDBW2Pos_X =  0.0   * cm;  G4double pDBW2Pos_Y =  9.0 * cm;  G4double pDBW2Pos_Z =  0.0  * cm;
  G4double pDBW3Pos_X = -4.13  * cm;  G4double pDBW3Pos_Y = -7.0 * cm;  G4double pDBW3Pos_Z =  0.0  * cm;
  G4double pDBW4Pos_X =  4.13  * cm;  G4double pDBW4Pos_Y = -7.0 * cm;  G4double pDBW4Pos_Z =  0.0  * cm;
  G4double pDMS1Pos_X =  0.0   * cm;  G4double pDMS1Pos_Y =  0.0 * cm;  G4double pDMS1Pos_Z =  0.0  * cm;
  G4double pDMBHPos_X =  0.0   * cm;  G4double pDMBHPos_Y =  9.0 * cm;  G4double pDMBHPos_Z =  0.0  * cm;
  G4double pDMFRPos_X =  0.0   * cm;  G4double pDMFRPos_Y =  0.0 * cm;  G4double pDMFRPos_Z =-87.2  * cm;

  G4VSolid* DBI1Solid = new G4Box ( "DBI1Box" , pDBI1HLX , pDBI1HLY  , pDBI1HLZ );
  G4VSolid* DBV1Solid = new G4Box ( "DBV1Box" , pDBV1HLX , pDBV1HLY  , pDBV1HLZ );
  G4VSolid* DBW1Solid = new G4Box ( "DBW1Box" , pDBW1HLX , pDBW1HLY  , pDBW1HLZ );
  G4VSolid* DBW2Solid = new G4Tubs( "DBW2Tubs", pDBW2Rin , pDBW2Rout , pDBW2HL  , 0.0, 360.0 * deg );
  G4VSolid* DBW3Solid = new G4Box ( "DBW3Box" , pDBW3HLX , pDBW3HLY  , pDBW3HLZ );
  G4VSolid* DBW4Solid = new G4Box ( "DBW4Box" , pDBW4HLX , pDBW4HLY  , pDBW4HLZ );
  G4VSolid* DBT3Solid = new G4Box ( "DBT3Box" , pDBT3HLX , pDBT3HLY  , pDBT3HLZ );
  G4VSolid* DBT4Solid = new G4Box ( "DBT4Box" , pDBT4HLX , pDBT4HLY  , pDBT4HLZ );
  G4VSolid* DMS1Solid = new G4Box ( "DMS1Box" , pDMS1HLX , pDMS1HLY  , pDMS1HLZ );
  G4VSolid* DMBHSolid = new G4Tubs( "DMBHTubs", pDMBHRin , pDMBHRout , pDMBHHL  , 0.0, 360.0 * deg );
  G4VSolid* DMFRSolid = new G4Tubs( "DMFRTubs", pDMFRRin , pDMFRRout , pDMFRHL  , 0.0, 360.0 * deg );

  G4SubtractionSolid* sub1 = new G4SubtractionSolid("sub1", DBI1Solid, DBV1Solid, 0,
                           G4ThreeVector(pDBV1Pos_X, pDBV1Pos_Y, pDBV1Pos_Z) );

  G4UnionSolid*       sub2 = new G4UnionSolid("sub2", sub1 , DMS1Solid, 0,
                           G4ThreeVector(pDMS1Pos_X + pDBV1Pos_X,
                           pDMS1Pos_Y + pDBV1Pos_Y, pDMS1Pos_Z + pDBV1Pos_Z) );

  G4SubtractionSolid* sub3 = new G4SubtractionSolid("sub3", sub2 , DMBHSolid, 0,
                           G4ThreeVector(pDMBHPos_X + pDMS1Pos_X + pDBV1Pos_X,
                           pDMBHPos_Y + pDMS1Pos_Y + pDBV1Pos_Y,
                           pDMBHPos_Z + pDMS1Pos_Z + pDBV1Pos_Z) );

  G4UnionSolid*       sub4 = new G4UnionSolid      ("sub4", sub3 , DMFRSolid, 0,
                           G4ThreeVector(pDMFRPos_X + pDMBHPos_X + pDMS1Pos_X + pDBV1Pos_X,
                           pDMFRPos_Y + pDMBHPos_Y + pDMS1Pos_Y + pDBV1Pos_Y,
                           pDMFRPos_Z + pDMBHPos_Z + pDMS1Pos_Z + pDBV1Pos_Z) );

  // Create Sub0 which subtracts the sub4 volume from the DMagSolid Volume, create DLogical from sub0, and place as before.
  G4SubtractionSolid* sub0 = new G4SubtractionSolid("sub0", DMagSolid , sub4 , 0 , G4ThreeVector (0,0,0));
  G4LogicalVolume* DLogical = new G4LogicalVolume ( sub0, MolPol_Vacuum, "DipoleMag", 0, 0, 0);
  DLogical->SetFieldManager(DFieldManager,allLocal);
  DLogical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0,G4ThreeVector(pDMagPos_X, pDMagPos_Y - 9*cm, pDMagPos_Z), DLogical,"DipoleMag",world_log,0,0,fCheckOverlaps);
  //Create the sub4 logical volume and place this in the WORLD_LOG not the DipoleMag volume.
  G4LogicalVolume* sub4Logical = new G4LogicalVolume ( sub4, MolPol_SiliSteel, "DipoleVacuumBox", 0, 0, 0);
  sub4Logical->SetVisAttributes(SteelVisAtt);
  new G4PVPlacement(0,G4ThreeVector(pDBI1Pos_X, pDBI1Pos_Y, pDBI1Pos_Z),sub4Logical,"DipoleVacuumBox",world_log,0,0,fCheckOverlaps);
  //This is the back dipole window. Start with the steel plate, take out the windows (DBW3/4), cut out the beampipe hole (DBW2)
  G4SubtractionSolid* sub9  = new G4SubtractionSolid("sub9" , DBW1Solid, DBW2Solid, 0, G4ThreeVector(pDBW2Pos_X, pDBW2Pos_Y, pDBW2Pos_Z) );
  G4SubtractionSolid* sub10 = new G4SubtractionSolid("sub10", sub9     , DBW3Solid, 0, G4ThreeVector(pDBW3Pos_X, pDBW3Pos_Y, pDBW3Pos_Z) );
  G4SubtractionSolid* sub11 = new G4SubtractionSolid("sub11", sub10    , DBW4Solid, 0, G4ThreeVector(pDBW4Pos_X, pDBW4Pos_Y, pDBW4Pos_Z) );
  G4LogicalVolume* sub11Logical = new G4LogicalVolume ( sub11, MolPol_SiliSteel, "sub13Logical", 0, 0, 0);
  sub11Logical->SetVisAttributes(SteelVisAtt);
  new G4PVPlacement(0,G4ThreeVector(pDBW1Pos_X, pDBW1Pos_Y, pDBW1Pos_Z),sub11Logical,"DipoleVacuumBox",world_log,0,0,fCheckOverlaps);
  //Thin titanium windows for the dipole exit windows, must be placed at z = pDBW1Pos_Z and place at the positions of the cutouts (pDBW3/4) plus the dipole vacuum box Y offset
  G4LogicalVolume* DipoleExitWindowLLogical = new G4LogicalVolume ( DBT3Solid, MolPol_Titanium, "DipoleExitWindowL", 0, 0, 0);
  G4LogicalVolume* DipoleExitWindowRLogical = new G4LogicalVolume ( DBT4Solid, MolPol_Titanium, "DipoleExitWindowR", 0, 0, 0);
  new G4PVPlacement(0,G4ThreeVector(pDBW3Pos_X, pDBW3Pos_Y + pDBI1Pos_Y, pDBW1Pos_Z),DipoleExitWindowLLogical,"DipoleExitWindowL",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(pDBW4Pos_X, pDBW4Pos_Y + pDBI1Pos_Y, pDBW1Pos_Z),DipoleExitWindowRLogical,"DipoleExitWindowR",world_log,0,0,fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DIPOLE LEAD COLLIMATOR JAWS
  G4double pDCOLHLX   =  1.00 * cm;  G4double pDCOLHLY   = 2.00 * cm;  G4double pDCOLHLZ   =  3.00  * cm;
  G4double pDCOLPos_X = -4.00 * cm;  G4double pDCOLPos_Y = 9.00 * cm;  G4double pDCOLPos_Z =-92.00  * cm;
  G4VSolid* DCOLSolid = new G4Box ( "DCOLBox" , pDCOLHLX , pDCOLHLY  , pDCOLHLZ );
  G4LogicalVolume * fLeadJawsLogicalLT = new G4LogicalVolume ( DCOLSolid, MolPol_Lead, "CollimatorLT", 0, 0, 0);
  G4LogicalVolume * fLeadJawsLogicalLB = new G4LogicalVolume ( DCOLSolid, MolPol_Lead, "CollimatorLB", 0, 0, 0);
  G4LogicalVolume * fLeadJawsLogicalRT = new G4LogicalVolume ( DCOLSolid, MolPol_Lead, "CollimatorRT", 0, 0, 0);
  G4LogicalVolume * fLeadJawsLogicalRB = new G4LogicalVolume ( DCOLSolid, MolPol_Lead, "CollimatorRB", 0, 0, 0);

  fLeadJawsLogicalLT ->SetVisAttributes(LeadVisAtt);
  fLeadJawsLogicalLB ->SetVisAttributes(LeadVisAtt);
  fLeadJawsLogicalRT ->SetVisAttributes(LeadVisAtt);
  fLeadJawsLogicalRB ->SetVisAttributes(LeadVisAtt);

  fLeadJawsHLength = pDCOLHLY;
  fLeadJawsZOrigin = pDCOLPos_Z + pDBV1Pos_Z + pDBI1Pos_Z;
  fLeadJawsYOrigin = pDCOLPos_Y + pDBV1Pos_Y + pDBI1Pos_Y;
  fLeadJawsXOrigin = pDCOLPos_X + pDBV1Pos_X + pDBI1Pos_X;

  fLeadJawsPhysicalLT = new G4PVPlacement(0,G4ThreeVector(-fLeadJawsXOrigin, (fLeadJawsYOrigin + (fLeadJawsHLength + (fLeadJawGapWidth/2.0))), fLeadJawsZOrigin),fLeadJawsLogicalLT,"CollimatorLT",world_log,0,0,fCheckOverlaps);
  fLeadJawsPhysicalLB = new G4PVPlacement(0,G4ThreeVector(-fLeadJawsXOrigin, (fLeadJawsYOrigin - (fLeadJawsHLength + (fLeadJawGapWidth/2.0))), fLeadJawsZOrigin),fLeadJawsLogicalLB,"CollimatorLB",world_log,0,0,fCheckOverlaps);
  fLeadJawsPhysicalRT = new G4PVPlacement(0,G4ThreeVector( fLeadJawsXOrigin, (fLeadJawsYOrigin + (fLeadJawsHLength + (fLeadJawGapWidth/2.0))), fLeadJawsZOrigin),fLeadJawsLogicalRT,"CollimatorRT",world_log,0,0,fCheckOverlaps);
  fLeadJawsPhysicalRB = new G4PVPlacement(0,G4ThreeVector( fLeadJawsXOrigin, (fLeadJawsYOrigin - (fLeadJawsHLength + (fLeadJawGapWidth/2.0))), fLeadJawsZOrigin),fLeadJawsLogicalRB,"CollimatorRB",world_log,0,0,fCheckOverlaps);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // UPSTREAM FLANGE ATTACHED TO DIPOLE BOX
  G4double pDUPFLRin = 0.5 * 10.16 * cm;  G4double pDUPFLRout = 0.5 * 17.145 * cm;  G4double pDUPFLHLZ = 0.5 * 2.1336 * cm;
  G4VSolid * DUpstreamFlange = new G4Tubs( "DUpstreamFlange", pDUPFLRin , pDUPFLRout , pDUPFLHLZ  , 0.0, 360.0 * deg );
  G4LogicalVolume * DUpstreamFlangeLogical = new G4LogicalVolume(DUpstreamFlange,MolPol_Stainless,"DUpstreamFlangeLogical",0,0,0);
  DUpstreamFlangeLogical->SetVisAttributes(LeadVisAtt);
  new G4PVPlacement(0,G4ThreeVector( 0 , 0 , pDBI1Pos_Z - pDBI1HLZ - pDUPFLHLZ ),DUpstreamFlangeLogical,"DipoleUpstreamFlange",world_log,0,0,fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Dipole Magnets Physical
  G4double pDHLX   = 10.00 * cm;  G4double pDHLY   = 40.00 * cm;  G4double pDHLZ   = 76.00 * cm;
  G4double pD1Pos_X=-16.00 * cm;  G4double pD2Pos_X= 16.00 * cm;  G4double pDPos_Y = -9.00 * cm;  G4double pDPos_Z =422.80 * cm;
  G4VSolid* DSolid  = new G4Box ( "DBox"  , pDHLX , pDHLY  , pDHLZ );
  G4LogicalVolume* DMagLogical  = new G4LogicalVolume(DSolid ,MolPol_SiliSteel, "DMagLogical" ,0,0,0);
  DMagLogical->SetVisAttributes(DipVisAtt);

  new G4PVPlacement(0,G4ThreeVector(pD1Pos_X,pDPos_Y,pDPos_Z),DMagLogical,"Dipole1",world_log,0,0,0);
  new G4PVPlacement(0,G4ThreeVector(pD2Pos_X,pDPos_Y,pDPos_Z),DMagLogical,"Dipole2",world_log,0,0,0);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // BEAM PIPE, BUILDS SEGMENTS FROM SOLENOID THROUGH QUADS TO G3 SPECIFIED END AT DIPOLE ENTRANCE
  G4double pBPpos[10] = {fHelmCoilBeamPosZ + pQ6HL,      //end of solenoid (world coord's)
                         pQ1Pos_Z - pQ1HL,      //begin Q1 (world coord's)
                         pQ1Pos_Z + pQ1HL,      //end Q1   (world coord's)
                         pQ2Pos_Z - pQ2HL,      //begin Q2 (world coord's)
                         pQ2Pos_Z + pQ2HL,      //end Q2   (world coord's)
                         pQ3Pos_Z - pQ3HL,      //begin Q3 (world coord's)
                         pQ3Pos_Z + pQ3HL,      //end Q3   (world coord's)
                         pQ4Pos_Z - pQ4HL,      //begin Q4 (world coord's)
                         pQ4Pos_Z + pQ4HL,      //end Q4   (world coord's)
                         323.132*cm};//beampipe ends at 323.132 per G3 geometry file.

  /*GPARVOL32  'BPI2' 209  'HALL'    0.    0.  311.566 0  'TUBE'  3   0.     5.08  11.566*/

  G4VSolid * BPUpstreamAlum[9];
  G4VSolid * BPUpstreamVac[9];
  G4LogicalVolume * BPAlLogVol[9];
  G4LogicalVolume * BPVacLogVol[9];
  for(G4int i = 0; i < 9; i++){
    G4String solidName = "BPalum_" + std::to_string(i);
    G4String vacName   = "BPvacm_" + std::to_string(i);
    BPUpstreamAlum[i]  = new G4Tubs( solidName, 0.*cm, pBPRout, (pBPpos[i+1] - pBPpos[i]) / 2., 0.0 * deg, 360.0 * deg );
    BPUpstreamVac[i]   = new G4Tubs( vacName,   0.*cm, pBPRin,  (pBPpos[i+1] - pBPpos[i]) / 2., 0.0 * deg, 360.0 * deg );
    G4String solidLogName = "BPalum_logical_" + std::to_string(i);
    G4String vacLogName   = "BPvacm_logical_" + std::to_string(i);
    BPAlLogVol[i]         = new G4LogicalVolume( BPUpstreamAlum[i], MolPol_Aluminum, solidLogName , 0, 0, 0);
    BPAlLogVol[i]->SetVisAttributes(BPVisAtt);
    BPVacLogVol[i]        = new G4LogicalVolume( BPUpstreamVac[i] , MolPol_Vacuum,   vacLogName, 0, 0, 0);
    BPVacLogVol[i]->SetVisAttributes(BPVisAtt);
    G4String solidPVName  = "BeamPipeAlmnum_" + std::to_string(i);
    G4String vacPVName    = "BeamPipeVacuum_" + std::to_string(i);
    if( i % 2 == 0 ) {
        new G4PVPlacement(0, G4ThreeVector(0,0,((pBPpos[i+1] + pBPpos[i]) / 2.)), BPAlLogVol[i],  solidPVName, world_log,     0, 0, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0,0,0),                                BPVacLogVol[i], vacPVName,   BPAlLogVol[i], 0, 0, fCheckOverlaps);
    } else {
        if( i == 1){
            new G4PVPlacement(0, G4ThreeVector(0,0,0), BPAlLogVol[i],  solidPVName, Q1MagLogical,  0, 0, fCheckOverlaps);
            new G4PVPlacement(0, G4ThreeVector(0,0,0), BPVacLogVol[i], vacPVName,   BPAlLogVol[i], 0, 0, fCheckOverlaps);
        }
        if( i == 3){
            new G4PVPlacement(0, G4ThreeVector(0,0,0), BPAlLogVol[i],  solidPVName, Q2MagLogical,  0, 0, fCheckOverlaps);
            new G4PVPlacement(0, G4ThreeVector(0,0,0), BPVacLogVol[i], vacPVName,   BPAlLogVol[i], 0, 0, fCheckOverlaps);
        }
        if( i == 5){
            new G4PVPlacement(0, G4ThreeVector(0,0,0), BPAlLogVol[i],  solidPVName, Q3MagLogical,  0, 0, fCheckOverlaps);
            new G4PVPlacement(0, G4ThreeVector(0,0,0), BPVacLogVol[i], vacPVName,   BPAlLogVol[i], 0, 0, fCheckOverlaps);
        }
        if( i == 7){
            new G4PVPlacement(0, G4ThreeVector(0,0,0), BPAlLogVol[i],  solidPVName, Q4MagLogical,  0, 0, fCheckOverlaps);
            new G4PVPlacement(0, G4ThreeVector(0,0,0), BPVacLogVol[i], vacPVName,   BPAlLogVol[i], 0, 0, fCheckOverlaps);
        }

     }
  }


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Planes for Virtual Detectors
  G4double         VBSolidHLZ = 0.00001*mm;
  G4double         VBSolidRad = 4.78*cm;
  G4VSolid*        VBSolid   = new G4Tubs("VBSolid",0,VBSolidRad, VBSolidHLZ ,0.0,360.0*deg);
  Q1ENLogical = new G4LogicalVolume(VBSolid, MolPol_Vacuum, "Q1ENLogical",0,0,0);
  Q1EXLogical = new G4LogicalVolume(VBSolid, MolPol_Vacuum, "Q1EXLogical",0,0,0);
  Q2ENLogical = new G4LogicalVolume(VBSolid, MolPol_Vacuum, "Q2ENLogical",0,0,0);
  Q2EXLogical = new G4LogicalVolume(VBSolid, MolPol_Vacuum, "Q2EXLogical",0,0,0);
  Q3ENLogical = new G4LogicalVolume(VBSolid, MolPol_Vacuum, "Q3ENLogical",0,0,0);
  Q3EXLogical = new G4LogicalVolume(VBSolid, MolPol_Vacuum, "Q3EXLogical",0,0,0);
  Q4ENLogical = new G4LogicalVolume(VBSolid, MolPol_Vacuum, "Q4ENLogical",0,0,0);
  Q4EXLogical = new G4LogicalVolume(VBSolid, MolPol_Vacuum, "Q4EXLogical",0,0,0);
  DipELogical = new G4LogicalVolume(VBSolid, MolPol_Vacuum, "DIPELogical",0,0,0);

  Q1ENLogical->SetVisAttributes(VacVisAtt);
  Q1EXLogical->SetVisAttributes(VacVisAtt);
  Q2ENLogical->SetVisAttributes(VacVisAtt);
  Q2EXLogical->SetVisAttributes(VacVisAtt);
  Q3ENLogical->SetVisAttributes(VacVisAtt);
  Q3EXLogical->SetVisAttributes(VacVisAtt);
  Q4ENLogical->SetVisAttributes(VacVisAtt);
  Q4EXLogical->SetVisAttributes(VacVisAtt);
  DipELogical->SetVisAttributes(VacVisAtt);

  MolPolDetector* Q1ENSD = new MolPolDetector("q1en", 1);
  MolPolDetector* Q1EXSD = new MolPolDetector("q1ex", 2);
  MolPolDetector* Q2ENSD = new MolPolDetector("q2en", 3);
  MolPolDetector* Q2EXSD = new MolPolDetector("q2ex", 4);
  MolPolDetector* Q3ENSD = new MolPolDetector("q3en", 5);
  MolPolDetector* Q3EXSD = new MolPolDetector("q3ex", 6);
  MolPolDetector* Q4ENSD = new MolPolDetector("q4en", 7);
  MolPolDetector* Q4EXSD = new MolPolDetector("q4ex", 8);

  MolPolDetector* DETSD  = new MolPolDetector("det",  9);
  //MolPolDetector* DETSD2 = new MolPolDetector("det2", 10);//No longer used
  MolPolDetector* APPSD1 = new MolPolDetector("app1", 11);
  MolPolDetector* APPSD2 = new MolPolDetector("app2", 12);
  MolPolDetector* DETVP  = new MolPolDetector("vp",   13);
  MolPolDetector* DPIN   = new MolPolDetector("dpin", 14);
  MolPolDetector* DPOUT  = new MolPolDetector("dpout",15);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  SDman->AddNewDetector(Q1ENSD);
  SDman->AddNewDetector(Q1EXSD);
  SDman->AddNewDetector(Q2ENSD);
  SDman->AddNewDetector(Q2EXSD);
  SDman->AddNewDetector(Q3ENSD);
  SDman->AddNewDetector(Q3EXSD);
  SDman->AddNewDetector(Q4ENSD);
  SDman->AddNewDetector(Q4EXSD);
  SDman->AddNewDetector(DETSD );
  //SDman->AddNewDetector(DETSD2);//No longer used
  SDman->AddNewDetector(APPSD1);
  SDman->AddNewDetector(APPSD2);
  SDman->AddNewDetector(DETVP );
  SDman->AddNewDetector(DPIN  );
  SDman->AddNewDetector(DPOUT );
  // Sensitive detector assignment for detectors 1-8 and 13-15 controlled by fluxVPs macro command
  // Detectors 9, 11, 12 are always active

  new G4PVPlacement(0, G4ThreeVector( 0, 0, -1.*pQ1HL + VBSolidHLZ ),       Q1ENLogical,   "VP.Q1.Entr", BPVacLogVol[1], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0,     pQ1HL - VBSolidHLZ ),       Q1EXLogical,   "VP.Q1.Exit", BPVacLogVol[1], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0, -1.*pQ2HL + VBSolidHLZ ),       Q2ENLogical,   "VP.Q2.Entr", BPVacLogVol[3], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0,     pQ2HL - VBSolidHLZ ),       Q2EXLogical,   "VP.Q2.Exit", BPVacLogVol[3], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0, -1.*pQ3HL + VBSolidHLZ ),       Q3ENLogical,   "VP.Q3.Entr", BPVacLogVol[5], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0,     pQ3HL - VBSolidHLZ ),       Q3EXLogical,   "VP.Q3.Exit", BPVacLogVol[5], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0, -1.*pQ4HL + VBSolidHLZ ),       Q4ENLogical,   "VP.Q4.Entr", BPVacLogVol[7], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0,     pQ4HL - VBSolidHLZ ),       Q4EXLogical,   "VP.Q4.Exit", BPVacLogVol[7], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0, ((pBPpos[9]-pBPpos[8])/2.)-VBSolidHLZ),DipELogical,"VP.Dp.Entr",BPVacLogVol[8],0,0,fCheckOverlaps);

  G4cout << std::setprecision(5) << "Z-Position of Q1ENLogical: " << -1.*pQ1HL + VBSolidHLZ + pQ1Pos_Z << G4endl;
  G4cout << std::setprecision(5) << "Z-Position of Q1EXLogical: " <<     pQ1HL - VBSolidHLZ + pQ1Pos_Z << G4endl;
  G4cout << std::setprecision(5) << "Z-Position of Q2ENLogical: " << -1.*pQ2HL + VBSolidHLZ + pQ2Pos_Z << G4endl;
  G4cout << std::setprecision(5) << "Z-Position of Q2EXLogical: " <<     pQ2HL - VBSolidHLZ + pQ2Pos_Z << G4endl;
  G4cout << std::setprecision(5) << "Z-Position of Q3ENLogical: " << -1.*pQ3HL + VBSolidHLZ + pQ3Pos_Z << G4endl;
  G4cout << std::setprecision(5) << "Z-Position of Q3EXLogical: " <<     pQ3HL - VBSolidHLZ + pQ3Pos_Z << G4endl;
  G4cout << std::setprecision(5) << "Z-Position of Q4ENLogical: " << -1.*pQ4HL + VBSolidHLZ + pQ4Pos_Z << G4endl;
  G4cout << std::setprecision(5) << "Z-Position of Q4EXLogical: " <<     pQ4HL - VBSolidHLZ + pQ4Pos_Z << G4endl;
  G4cout << std::setprecision(5) << "Z-Position of DipELogical: " <<  ((pBPpos[9] + pBPpos[8]) / 2.) + ((pBPpos[9]-pBPpos[8]) / 2.) - VBSolidHLZ << G4endl;


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DIPOLE Virtual Planes - NOTE: Dipole Entrance VP integrated into beampipe construction above.
  G4double pVP2HLX    = 6.00 * cm;   G4double pVP2HLY   = 20.00 * cm;  G4double pVP2HLZ   = 0.001   * cm;
  G4double pVP3Pos_X  = 0.00  * cm;  G4double pVP3Pos_Y = -9.00  * cm;  G4double pVP3Pos_Z = (537.0*cm - 14.0*cm - pVP2HLZ);

  G4VSolid* VP3Solid  = new G4Box( "VP3BOX",  pVP2HLX, pVP2HLY, pVP2HLZ );
  VP3Logical = new G4LogicalVolume(VP3Solid, MolPol_Vacuum, "VP3Logical", 0,0,0);
  // Sensitive detector assignment controlled by fluxVPs macro command
  VP3Logical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0,G4ThreeVector(pVP3Pos_X, pVP3Pos_Y, pVP3Pos_Z), VP3Logical, "VP.Dp.Exit", world_log, 0,0, fCheckOverlaps);

  G4cout << std::setprecision(5) << "Z-Position of DipXLogical: " << pVP3Pos_Z << G4endl;


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Virtual Planes inside of dipole box :: this block can safely be removed when no longer needed. -Eric King
  G4double DPIVPHLX = 2.295 * cm / 2;
  G4double DPIVPHLY = 15 * cm;
  G4double DPIVPHLZ = 0.00001 * cm;
  G4double DPIxloc = (3 + 2.295 * 0.5) * cm;// will be negative for LHS and pos for RHS
  G4double DPIzloc[10] = {-82.*cm,-63.*cm,-45.*cm,-27.*cm,-9.*cm,9.*cm,27.*cm,45.*cm,63.*cm,82.*cm};
  G4double DPIyloc = 0.0 * cm;
  G4VSolid * DPIVPbox = new G4Box("DipoleInternalVPlane",DPIVPHLX,DPIVPHLY,DPIVPHLZ);
  DP0L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP0L",0,0,0);
  DP0R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP0R",0,0,0);
  DP1L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP1L",0,0,0);
  DP1R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP1R",0,0,0);
  DP2L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP2L",0,0,0);
  DP2R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP2R",0,0,0);
  DP3L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP3L",0,0,0);
  DP3R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP3R",0,0,0);
  DP4L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP4L",0,0,0);
  DP4R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP4R",0,0,0);
  DP5L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP5L",0,0,0);
  DP5R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP5R",0,0,0);
  DP6L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP6L",0,0,0);
  DP6R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP6R",0,0,0);
  DP7L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP7L",0,0,0);
  DP7R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP7R",0,0,0);
  DP8L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP8L",0,0,0);
  DP8R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP8R",0,0,0);
  DP9L = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP9L",0,0,0);
  DP9R = new G4LogicalVolume(DPIVPbox, MolPol_Vacuum, "DP9R",0,0,0);
  DP0L->SetVisAttributes(VacVisAtt);
  DP0R->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[0]),DP0L,"DP0L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[0]),DP0R,"DP0R",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[1]),DP1L,"DP1L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[1]),DP1R,"DP1R",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[2]),DP2L,"DP2L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[2]),DP2R,"DP2R",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[3]),DP3L,"DP3L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[3]),DP3R,"DP3R",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[4]),DP4L,"DP4L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[4]),DP4R,"DP4R",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[5]),DP5L,"DP5L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[5]),DP5R,"DP5R",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[6]),DP6L,"DP6L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[6]),DP6R,"DP6R",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[7]),DP7L,"DP7L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[7]),DP7R,"DP7R",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[8]),DP8L,"DP8L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[8]),DP8R,"DP8R",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-DPIxloc,DPIyloc,DPIzloc[9]),DP9L,"DP9L",DLogical,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(+DPIxloc,DPIyloc,DPIzloc[9]),DP9R,"DP9R",DLogical,0,0,fCheckOverlaps);
  MolPolDetector* DP0Lsd = new MolPolDetector("dp0L",100);
  MolPolDetector* DP0Rsd = new MolPolDetector("dp0R",101);
  MolPolDetector* DP1Lsd = new MolPolDetector("dp1L",110);
  MolPolDetector* DP1Rsd = new MolPolDetector("dp1R",111);
  MolPolDetector* DP2Lsd = new MolPolDetector("dp2L",120);
  MolPolDetector* DP2Rsd = new MolPolDetector("dp2R",121);
  MolPolDetector* DP3Lsd = new MolPolDetector("dp3L",130);
  MolPolDetector* DP3Rsd = new MolPolDetector("dp3R",131);
  MolPolDetector* DP4Lsd = new MolPolDetector("dp4L",140);
  MolPolDetector* DP4Rsd = new MolPolDetector("dp4R",141);
  MolPolDetector* DP5Lsd = new MolPolDetector("dp5L",150);
  MolPolDetector* DP5Rsd = new MolPolDetector("dp5R",151);
  MolPolDetector* DP6Lsd = new MolPolDetector("dp6L",160);
  MolPolDetector* DP6Rsd = new MolPolDetector("dp6R",161);
  MolPolDetector* DP7Lsd = new MolPolDetector("dp7L",170);
  MolPolDetector* DP7Rsd = new MolPolDetector("dp7R",171);
  MolPolDetector* DP8Lsd = new MolPolDetector("dp8L",180);
  MolPolDetector* DP8Rsd = new MolPolDetector("dp8R",181);
  MolPolDetector* DP9Lsd = new MolPolDetector("dp9L",190);
  MolPolDetector* DP9Rsd = new MolPolDetector("dp9R",191);
  SDman->AddNewDetector(DP0Lsd);
  SDman->AddNewDetector(DP0Rsd);
  SDman->AddNewDetector(DP1Lsd);
  SDman->AddNewDetector(DP1Rsd);
  SDman->AddNewDetector(DP2Lsd);
  SDman->AddNewDetector(DP2Rsd);
  SDman->AddNewDetector(DP3Lsd);
  SDman->AddNewDetector(DP3Rsd);
  SDman->AddNewDetector(DP4Lsd);
  SDman->AddNewDetector(DP4Rsd);
  SDman->AddNewDetector(DP5Lsd);
  SDman->AddNewDetector(DP5Rsd);
  SDman->AddNewDetector(DP6Lsd);
  SDman->AddNewDetector(DP6Rsd);
  SDman->AddNewDetector(DP7Lsd);
  SDman->AddNewDetector(DP7Rsd);
  SDman->AddNewDetector(DP8Lsd);
  SDman->AddNewDetector(DP8Rsd);
  SDman->AddNewDetector(DP9Lsd);
  SDman->AddNewDetector(DP9Rsd);
  //////////////////////////////////////////////////////////////////////////////




  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Rotations information to be used for subsequent sections

  /* C RMATR07  90.   0.  80.  90.  10.  270.     Moller detector
   * RMATR07  90.   0.  82.7  90.  7.0  270.     Moller detector
   * RMATR09  90. 270.   0.   0.  90.  180.      I is oppos Y,II along Z,III oppos X */

  G4RotationMatrix* pRot9 = new G4RotationMatrix();
  pRot9->rotateY(90.*deg);
  pRot9->rotateZ(90.*deg);

  G4RotationMatrix* pRot7 = new G4RotationMatrix();
  pRot7->rotateX(-7.05*deg);




  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Shielding | 15: AIR; 9: Al; 13: Pb

  /* GPARVOL120 'S1LD'  13  'HALL'    0.  -29.  537.   0  'BOX '  3  30.  25. 14.0
     GPARVOL121 'S1H1'  15  'S1LD'   -4.1   9.5   0.   9  'PARA'  6  10.5 14. 1.6 10. 0. 0.
     GPARVOL122 'S1H2'  15  'S1LD'    4.1   9.5   0.   9  'PARA'  6  10.5 14. 1.6 10. 0. 0.
     GPARVOL123 'S2LD'  13  'HALL'    0.   -7.5 610.   0  'BOX '  3  17.   3.5 50.0 !!!***
     GPARVOL124 'S21B'   9  'HALL'    0.  -13.3 620.   0  'BOX '  3  31.   1.9  3.8 !!!***
     GPARVOL125 'S21H'  15  'S21B'    0.   -0.5   0.   0  'BOX '  3  31.   1.4  3.3
     GPARVOL126 'S3LD'  13  'HALL'    0.   -7.5 723.   0  'BOX '  3  20.   2.0 57.0  !!!***  */

  G4double pS1LDHLX   = 30.0 * cm;  G4double pS1LDHLY   = 25.0 * cm;  G4double pS1LDHLZ = 14.0 * cm;
  G4double pS1H1HLX   = 10.5 * cm;  G4double pS1H1HLY   = 14.0 * cm;  G4double pS1H1HLZ = 1.6  * cm;
  G4double pS1H1alpha = 10.0 * deg; G4double pS1H1theta = 0.0  * deg; G4double pS1H1phi = 0.0  * deg;
  G4double pS1H2HLX   = 10.5 * cm;  G4double pS1H2HLY   = 14.0 * cm;  G4double pS1H2HLZ = 1.6  * cm;
  G4double pS1H2alpha = 10.0 * deg; G4double pS1H2theta = 0.0  * deg; G4double pS1H2phi = 0.0  * deg;
  G4double pS2LDHLX   = 17.0 * cm;  G4double pS2LDHLY   = 3.5  * cm;  G4double pS2LDHLZ = 50.0 * cm;
  G4double pS21BHLX   = 31.0 * cm;  G4double pS21BHLY   = 1.9  * cm;  G4double pS21BHLZ = 3.8  * cm;
  G4double pS21HHLX   = 31.0 * cm;  G4double pS21HHLY   = 1.4  * cm;  G4double pS21HHLZ = 3.3  * cm;
  //G4double pS3LDHLX   = 20.0 * cm;  G4double pS3LDHLY   = 2.0  * cm;  G4double pS3LDHLZ = 57.0 * cm;//Unneeded as Shield4 (S3LD in G3) is not used

  G4double pS1LDPos_X   = 0.0  * cm;  G4double pS1LDPos_Y   = -29.0 * cm;  G4double pS1LDPos_Z = 537.0 * cm;
  G4double pS1H1Pos_X   = -4.1 * cm;  G4double pS1H1Pos_Y   = 9.5   * cm;  G4double pS1H1Pos_Z = 0     * cm;
  G4double pS1H2Pos_X   = 4.1  * cm;  G4double pS1H2Pos_Y   = 9.5   * cm;  G4double pS1H2Pos_Z = 0     * cm;
  G4double pS2LDPos_X   = 0.0  * cm;  G4double pS2LDPos_Y   = -7.5  * cm;  G4double pS2LDPos_Z = 610.0 * cm;
  G4double pS21BPos_X   = 0.0  * cm;  G4double pS21BPos_Y   = -13.3 * cm;  G4double pS21BPos_Z = 620.0 * cm;
  G4double pS21HPos_X   = 0.0  * cm;  G4double pS21HPos_Y   = -0.5  * cm;  G4double pS21HPos_Z = 0.0   * cm;
  //G4double pS3LDPos_X   = 0.0  * cm;  G4double pS3LDPos_Y   = -7.5  * cm;  G4double pS3LDPos_Z = 723.0 * cm;//Unneeded as Shield4 (S3LD in G3) is not used

  G4VSolid* S1LDSolid  = new G4Box ( "S1LDBox"  , pS1LDHLX, pS1LDHLY,  pS1LDHLZ );
  G4VSolid* S1H1Solid  = new G4Para( "S1H1Para" , pS1H1HLX, pS1H1HLY,  pS1H1HLZ, pS1H1alpha, pS1H1theta, pS1H1phi);
  G4VSolid* S1H2Solid  = new G4Para( "S1H2Para" , pS1H2HLX, pS1H2HLY,  pS1H2HLZ, pS1H2alpha, pS1H2theta, pS1H2phi);
  G4VSolid* S2LDSolid  = new G4Box ( "S2LDBox"  , pS2LDHLX, pS2LDHLY,  pS2LDHLZ );
  G4VSolid* S21BSolid  = new G4Box ( "S21BBox"  , pS21BHLX, pS21BHLY,  pS21BHLZ );
  G4VSolid* S21HSolid  = new G4Box ( "S21HBox"  , pS21HHLX, pS21HHLY,  pS21HHLZ );
  //G4VSolid* S3LDSolid  = new G4Box ( "S3LDBox"  , pS3LDHLX, pS3LDHLY,  pS3LDHLZ );//Unneeded as Shield4 (S3LD in G3) is not used

  G4SubtractionSolid* subSHL1 = new G4SubtractionSolid("subSHL1", S1LDSolid, S1H1Solid, pRot9, G4ThreeVector(pS1H1Pos_X, pS1H1Pos_Y, pS1H1Pos_Z) );
  G4SubtractionSolid* subSHL2 = new G4SubtractionSolid("subSHL2", subSHL1, S1H2Solid, pRot9, G4ThreeVector(pS1H2Pos_X, pS1H2Pos_Y, pS1H2Pos_Z) );
  G4SubtractionSolid* subSHL3 = new G4SubtractionSolid("subSHL3", S21BSolid, S21HSolid, 0, G4ThreeVector(pS21HPos_X, pS21HPos_Y, pS21HPos_Z) );

  G4LogicalVolume* SHL1Logical = new G4LogicalVolume(subSHL2,   MolPol_Lead, "Shield1", 0,0,0);
  G4LogicalVolume* SHL2Logical = new G4LogicalVolume(S2LDSolid, MolPol_Lead, "Shield2", 0,0,0);
  G4LogicalVolume* SHL3Logical = new G4LogicalVolume(subSHL3,   MolPol_Aluminum, "Shield3", 0,0,0);
  SHL1Logical->SetVisAttributes(LeadVisAtt);
  SHL2Logical->SetVisAttributes(LeadVisAtt);
  SHL3Logical->SetVisAttributes(AlumVisAtt);

  new G4PVPlacement(0, G4ThreeVector(pS1LDPos_X, pS1LDPos_Y, pS1LDPos_Z), SHL1Logical, "Shield1", world_log, 0,0,fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(pS2LDPos_X, pS2LDPos_Y, pS2LDPos_Z), SHL2Logical, "Shield2", world_log, 0,0,fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(pS21BPos_X, pS21BPos_Y, pS21BPos_Z), SHL3Logical, "Shield3", world_log, 0,0,fCheckOverlaps);

  //G4LogicalVolume* SHL4Logical = new G4LogicalVolume(S3LDSolid, MolPol_Lead, "Shield4", 0,0,0);
  //SHL4Logical->SetVisAttributes(LeadVisAtt);
  //new G4PVPlacement(0, G4ThreeVector(pS3LDPos_X, pS3LDPos_Y, pS3LDPos_Z), SHL4Logical, "Shield4", world_log, 0,0,fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DETECTOR AND BOX
  /*
   * C --     Detector box
   * GPARVOL80  'MDBX'  13  'HALL'    0.  -50.1 724.   0  'BOX '  3  37.   47.1   61.2  ***7cm up
   * GPARVOL81  'MDBA'  15  'MDBX'    0.    1.1   1.8  0  'BOX '  3  18.5  29.8   40.0
   * GPARVOL82  'MDBW'  15  'MDBX'    0.    8.5 -49.7  9  'PARA'  6  13.7  11.5   6.45   6.0  0.  0.
   * GPARVOL83  'MDBL'  15  'MDBX'    0.    7.5 -49.7  9  'PARA'  6  13.7  11.5   6.45  10.0  0.  0.
   * C --     Detector
   * GPARVOL90  'MDET'  15  'MDBA'    0.    0.    0.   7  'BOX '  3   9.2  17.0    34.0
   * GPARVOL91  'DLGB'   9  'MDET'    0.    0.  -11.   0  'BOX '  3   8.2  16.2    20.1
   * C
   * C --  Syracuse hodoscope, SPACAL, auxill. artificial plane
   * C
   * GPARVOL101 'HOD1'  15  'MDET'    0.    0.  -32.5  0  'BOX '  3   9.2  15.5     0.75
   * GPARVOL102 'HOD2'   9  'MDET'    0.    0.  -16.0  0  'BOX '  3   9.2  16.5    15.2
   * GPARVOL103 'HOD3'  15  'MDET'    0.    0.  -33.5  0  'BOX '  3   9.2  17.0     0.1
   */

  G4double pMDBXHLX   = 37.00 * cm;  G4double pMDBXHLY   = 47.10 * cm;  G4double pMDBXHLZ   = 61.20 * cm;
  G4double pMDBAHLX   = 18.50 * cm;  G4double pMDBAHLY   = 29.80 * cm;  G4double pMDBAHLZ   = 40.00 * cm;
  G4double pMDBWHLX   = 14.60 * cm;  G4double pMDBWHLY   = 11.50 * cm;  G4double pMDBWHLZ   = 6.305 * cm;
  G4double pMDBWalpha =  6.0  *deg;  G4double pMDBWtheta = 0.0   *deg;  G4double pMDBWphi   =  0.0  *deg;
  G4double pMDBLHLX   = 14.60 * cm;  G4double pMDBLHLY   = 11.50 * cm;  G4double pMDBLHLZ   = 6.305 * cm;
  G4double pMDBLalpha = 10.0  *deg;  G4double pMDBLtheta = 0.0   *deg;  G4double pMDBLphi   =  0.0  *deg;
  G4double pMDETHLX   =  9.20 * cm;  G4double pMDETHLY   = 17.00 * cm;  G4double pMDETHLZ   = 34.00 * cm;
  G4double pDLGBHLX   =  9.20 * cm;  G4double pDLGBHLY   = 16.20 * cm;  G4double pDLGBHLZ   = 20.10 * cm;
  G4double pSPCLHLX   =  9.00 * cm;  G4double pSPCLHLY   = 15.00 * cm;  G4double pSPCLHLZ   = 15.00 * cm;

  G4double pMDBXPos_X   =  0.055 * cm; G4double pMDBXPos_Y  = -46.70 * cm;  G4double pMDBXPos_Z   =723.20 * cm;//Adjusted for 2019 survey
  G4double pMDBAPos_X   =  0.00 * cm;  G4double pMDBAPos_Y   =  1.10 * cm;  G4double pMDBAPos_Z   =  1.80 * cm;
  G4double pMDBWPos_X   =  0.00 * cm;  G4double pMDBWPos_Y   =  8.49 * cm;  G4double pMDBWPos_Z   =-49.70 * cm;
  G4double pMDBLPos_X   =  0.00 * cm;  G4double pMDBLPos_Y   =  7.67 * cm;  G4double pMDBLPos_Z   =-49.70 * cm;
  G4double pMDETPos_X   =  0.405* cm;  G4double pMDETPos_Y   =  2.46 * cm;  G4double pMDETPos_Z   = -4.08 * cm;
  G4double pDLGBPos_X   =  0.00 * cm;  G4double pDLGBPos_Y   =  0.00 * cm;  G4double pDLGBPos_Z   =-11.00 * cm;//Calorimeter Case
  G4double pSPCLPos_X   =  0.00 * cm;  G4double pSPCLPos_Y   =  0.00 * cm;  G4double pSPCLPos_Z   =  -5.1 * cm;//Calorimeter Block

  G4VSolid* MDBXSolid  = new G4Box ( "MDBXBox"  , pMDBXHLX, pMDBXHLY, pMDBXHLZ );
  G4VSolid* MDBASolid  = new G4Box ( "MDBABox"  , pMDBAHLX, pMDBAHLY, pMDBAHLZ );
  G4VSolid* MDBWSolid  = new G4Para( "MDBWPara" , pMDBWHLX, pMDBWHLY, pMDBWHLZ, pMDBWalpha, pMDBWtheta, pMDBWphi);
  G4VSolid* MDBLSolid  = new G4Para( "MDBLPara" , pMDBLHLX, pMDBLHLY, pMDBLHLZ, pMDBLalpha, pMDBLtheta, pMDBLphi);
  G4VSolid* MDETSolid  = new G4Box ( "MDETBox"  , pMDETHLX, pMDETHLY, pMDETHLZ );
  //DLGB is the case which holds the spaghetti calorimeter; this is made of aluminum
  //Spaghetti calorimeter [SPCL] is of the dimensions specified in the Manual[2019].
  G4VSolid* DLGBSolid  = new G4Box ( "DLGBBox"  , pDLGBHLX, pDLGBHLY, pDLGBHLZ );
  G4VSolid* SPCLSolid  = new G4Box ( "SPCLBox"  , pSPCLHLX, pSPCLHLY, pSPCLHLZ );

  G4SubtractionSolid* sub5 = new G4SubtractionSolid("sub5", MDBXSolid, MDBASolid, 0    , G4ThreeVector(pMDBAPos_X, pMDBAPos_Y, pMDBAPos_Z) );
  G4SubtractionSolid* sub6 = new G4SubtractionSolid("sub6", sub5     , MDBWSolid, pRot9, G4ThreeVector(pMDBWPos_X, pMDBWPos_Y, pMDBWPos_Z) );
  G4SubtractionSolid* sub7 = new G4SubtractionSolid("sub7", sub6     , MDBLSolid, pRot9, G4ThreeVector(pMDBLPos_X, pMDBLPos_Y, pMDBLPos_Z) );
  G4LogicalVolume* MDBXLogical = new G4LogicalVolume ( sub7, MolPol_Lead, "Detector_MDBX", 0, 0, 0);
  G4VisAttributes *MDBXVisAtt(PbVisAtt);
  MDBXVisAtt->SetForceWireframe(true);
  MDBXLogical->SetVisAttributes(MDBXVisAtt);
  new G4PVPlacement(0 , G4ThreeVector(pMDBXPos_X, pMDBXPos_Y, pMDBXPos_Z) , MDBXLogical , "Detector_MDBX" , world_log , 0 , 0 , fCheckOverlaps);

  //This is the 'MDET' volume in G3. It's the mother volume for all the detector pieces which will be placed and then rotated by the 7.3*deg as a whole.
  //This must be positioned at the MDBA positions (In G3 this was just a daughter volume, in G4 it's part of a boolean subtraction).
  G4LogicalVolume * MDETLogical = new G4LogicalVolume( MDETSolid , MolPol_Air , "DetectorContainer" , 0 , 0 , 0);
  MDETLogical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(pRot7 , G4ThreeVector(pMDBAPos_X+pMDETPos_X+pMDBXPos_X,
                                          pMDBAPos_Y+pMDETPos_Y+pMDBXPos_Y,
                                          pMDBAPos_Z+pMDETPos_Z+pMDBXPos_Z), MDETLogical , "Detector_MDET" , world_log , 0 , 0 , fCheckOverlaps);

  G4LogicalVolume* DLGBLogical = new G4LogicalVolume(DLGBSolid, MolPol_Aluminum, "DLGBLogical",0,0,0);
  new G4PVPlacement(0 , G4ThreeVector(pDLGBPos_X, pDLGBPos_Y, pDLGBPos_Z) , DLGBLogical , "Detector_CAL_CASE" , MDETLogical , 0 , 0 , fCheckOverlaps);

  G4LogicalVolume* SPCLLogical = new G4LogicalVolume(SPCLSolid, MolPol_Vacuum, "SPCLLogical",0,0,0);
  new G4PVPlacement(0 , G4ThreeVector(pSPCLPos_X, pSPCLPos_Y, pSPCLPos_Z ) , SPCLLogical , "Detector_SpagCalor" , DLGBLogical , 0 , 0 , fCheckOverlaps);
  //Place virtual plane at the very front of calorimeter as a daughter volume for now and let that be sensitive detector #9
  G4double CalVPHLZ = 0.00001 * cm;
  G4VSolid* SpagCalVPSolid  = new G4Box ( "SpagCalVPSolid"  , pSPCLHLX, pSPCLHLY, CalVPHLZ );
  G4LogicalVolume* SpagCalVpLogical = new G4LogicalVolume(SpagCalVPSolid, MolPol_Vacuum, "SpagCalVpLogical",0,0,0);
  //This sets the flux plane that acts as our detector as active 
  SpagCalVpLogical->SetSensitiveDetector(DETSD);
  new G4PVPlacement(0 , G4ThreeVector(0. , 0. , -1.*pSPCLHLZ + CalVPHLZ ) , SpagCalVpLogical , "SpagCalor_VP" , SPCLLogical , 0 , 0 , fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Hodoscope stuffs
  /* AS WOULD BE CONSTRUCTED IN G3
   * G4double pHOD1HLX   = 9.2 * cm;   G4double pHOD1HLY   = 15.5 * cm;   G4double pHOD1HLZ   = 0.75 * cm;
   * G4double pHOD2HLX   = 9.2 * cm;   G4double pHOD2HLY   = 16.5 * cm;   G4double pHOD2HLZ   = 15.2 * cm;
   * G4double pHOD3HLX   = 9.2 * cm;   G4double pHOD3HLY   = 17.5 * cm;   G4double pHOD3HLZ   = 0.1  * cm;
   * G4double pHOD1Pos_X = 0.0 * cm;   G4double pHOD1Pos_Y = 0.0  * cm;   G4double pHOD1Pos_Z = -32.5* cm;
   * G4double pHOD2Pos_X = 0.0 * cm;   G4double pHOD2Pos_Y = 0.0  * cm;   G4double pHOD2Pos_Z = -16.0* cm;
   * G4double pHOD3Pos_X = 0.0 * cm;   G4double pHOD3Pos_Y = 0.0  * cm;   G4double pHOD3Pos_Z = -33.5* cm;
   * G4VSolid* HOD1Solid = new G4Box( "HOD1Box", pHOD1HLX, pHOD1HLY, pHOD1HLZ );
   * G4VSolid* HOD2Solid = new G4Box( "HOD2Box", pHOD2HLX, pHOD2HLY, pHOD2HLZ );
   * G4VSolid* HOD3Solid = new G4Box( "HOD3Box", pHOD3HLX, pHOD3HLY, pHOD3HLZ );
   * G4LogicalVolume* HOD1Logical = new G4LogicalVolume(HOD1Solid, MolPol_Air,      "DETLogical", 0, 0, 0);
   * G4LogicalVolume* HOD2Logical = new G4LogicalVolume(HOD2Solid, MolPol_Aluminum, "DETLogical", 0, 0, 0);
   * G4LogicalVolume* HOD3Logical = new G4LogicalVolume(HOD3Solid, MolPol_Air,      "DETLogical", 0, 0, 0);
   * new G4PVPlacement(0, G4ThreeVector(pHOD1Pos_X,pHOD1Pos_Y,pHOD1Pos_Z), HOD1Logical, "Detector_HOD1", MDETLogical, 0,0, fCheckOverlaps);
   * new G4PVPlacement(0, G4ThreeVector(pHOD2Pos_X,pHOD2Pos_Y,pHOD2Pos_Z), HOD2Logical, "Detector_HOD2", MDETLogical, 0,0, fCheckOverlaps);
   * new G4PVPlacement(0, G4ThreeVector(pHOD3Pos_X,pHOD3Pos_Y,pHOD3Pos_Z), HOD3Logical, "Detector_HOD3", MDETLogical, 0,0, fCheckOverlaps);
   */

  G4double pAPP1HLX   = 2.0 * cm;   G4double pAPP1HLY   = 15. * cm;   G4double pAPP1HLZ   =  0.65 * cm;// NOTE: 0.65 IS THE CORRECT VALUE FOR THE Z HALF LENGTH
  G4double pAPP1Pos_X = 4.4 * cm;   G4double pAPP1Pos_Y = 0.0 * cm;    G4double pAPP1Pos_Z = -32.5 * cm;
  G4VSolid* APP1LSolid = new G4Box( "APP1LBOX", pAPP1HLX, pAPP1HLY, pAPP1HLZ );
  G4VSolid* APP1RSolid = new G4Box( "APP1RBOX", pAPP1HLX, pAPP1HLY, pAPP1HLZ );
  G4LogicalVolume* APP1LLogical = new G4LogicalVolume(APP1LSolid, MolPol_Scint, "APP1LLogical", 0,0,0);
  // APP1LLogical->SetSensitiveDetector(APPSD1);
  APP1LLogical->SetVisAttributes(ScintVisAtt);
  G4LogicalVolume* APP1RLogical = new G4LogicalVolume(APP1RSolid, MolPol_Scint, "APP1RLogical", 0,0,0);
  // APP1RLogical->SetSensitiveDetector(APPSD2);
  APP1RLogical->SetVisAttributes(ScintVisAtt);
  new G4PVPlacement(0, G4ThreeVector(    pAPP1Pos_X,pAPP1Pos_Y,pAPP1Pos_Z), APP1LLogical, "Detector_APP1L", MDETLogical, 0,0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(-1.*pAPP1Pos_X,pAPP1Pos_Y,pAPP1Pos_Z), APP1RLogical, "Detector_APP1R", MDETLogical, 0,0, fCheckOverlaps);

  G4double  pVPHODHLZ = 0.0001 * cm;
  G4VSolid* VPHOD1Solid = new G4Box( "VPHOD1", pAPP1HLX, pAPP1HLY, pVPHODHLZ);
  G4VSolid* VPHOD2Solid = new G4Box( "VPHOD2", pAPP1HLX, pAPP1HLY, pVPHODHLZ);
  VPHOD1Logical = new G4LogicalVolume(VPHOD1Solid, MolPol_Vacuum, "VPHOD1Logical", 0,0,0);
  VPHOD1Logical->SetSensitiveDetector(APPSD1);
  VPHOD1Logical->SetVisAttributes(VacVisAtt);
  VPHOD2Logical = new G4LogicalVolume(VPHOD2Solid, MolPol_Vacuum, "VPHOD2Logical", 0,0,0);
  VPHOD2Logical->SetSensitiveDetector(APPSD2);
  VPHOD2Logical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0, G4ThreeVector(    pAPP1Pos_X,pAPP1Pos_Y,pAPP1Pos_Z-pAPP1HLZ-pVPHODHLZ), VPHOD1Logical, "APP1L-VP", MDETLogical, 0,0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(-1.*pAPP1Pos_X,pAPP1Pos_Y,pAPP1Pos_Z-pAPP1HLZ-pVPHODHLZ), VPHOD2Logical, "APP1R-VP", MDETLogical, 0,0, fCheckOverlaps);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DETECTOR Box Virtual Plane
  G4double pVP1HLX    = 37.00 * cm;  G4double pVP1HLY   = 50.00 * cm;  G4double pVP1HLZ   = 0.0001   * cm;
  G4double pVP1Pos_X  = 0.00  * cm;  G4double pVP1Pos_Y = -46.86* cm;  //G4double pVP1Pos_Z = 660.0 * cm;--Unneeded as it is 'pasted' onto the front of the Detector Box
  G4VSolid* VP1Solid  = new G4Box( "VP1BOX",  pVP1HLX, pVP1HLY, pVP1HLZ );
  VP1Logical = new G4LogicalVolume(VP1Solid, MolPol_Vacuum, "VPDetector", 0,0,0);
  // Sensitive detector assignment controlled by fluxVPs macro command
  VP1Logical->SetVisAttributes( VacVisAtt );
  new G4PVPlacement(0,G4ThreeVector(pVP1Pos_X, pVP1Pos_Y, pMDBXPos_Z - pMDBXHLZ - pVP1HLZ ), VP1Logical, "VP.Detector.Entr", world_log, 0,0, fCheckOverlaps);

  G4cout << "Z-position of DetectorBox Virtual Plane " << pMDBXPos_Z - pMDBXHLZ - pVP1HLZ << G4endl;

  return world_phys;
}

void MolPolDetectorConstruction::ConstructMaterials(){
  G4double a, z, density, pressure, temperature;
  G4int nelements, natoms;

  G4Element* N  = new G4Element("Nitrogen"  , "N" , z=7 , a=14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"    , "O" , z=8 , a=16.00*g/mole);
  //G4Element* H  = new G4Element("Hydrogen"  , "H" , z=1 , a=1.01 *g/mole); // Unneeded at this time
  G4Element* C  = new G4Element("Carbon"    , "C" , z=6 , a=12.01*g/mole);

  // USED VALUES FROM WOLFRAMALPHA PERIODIC TABLE DATA | ERIC KING
  G4Element* P  = new G4Element("Phosphorus", "P" , z=15, a=30.974*g/mole);
  G4Element* S  = new G4Element("Sulfur"    , "S" , z=16, a=32.06 *g/mole);

  G4Element* Al = new G4Element("Aluminum"  , "Al", z=13, a=26.98 *g/mole);
  G4Element* Fe = new G4Element("Iron"      , "Fe", z=26, a=55.845*g/mole);
  G4Element* Si = new G4Element("Silicon"   , "Si", z=14, a=28.09 *g/mole);
  //G4Element* Pb = new G4Element("Lead"      , "Pb", z=82, a=207.19*g/mole); // Unneeded at this time

  // USED VALUES FROM WOLFRAMALPHA PERIODIC TABLE DATA | ERIC KING
  G4Element* Mn = new G4Element("Manganese" , "Mn", z=25, a=54.938*g/mole);
  G4Element* Cr = new G4Element("Chromium"  , "Cr", z=24, a=51.966*g/mole);
  G4Element* Ni = new G4Element("Nickel"    , "Ni", z=28, a=58.693*g/mole);
  G4Element* Mo = new G4Element("Molybdenum", "Mo", z=42, a=95.95 *g/mole);

  // INFORMATION FROM SANGHWA
  density = 7.93 *g/cm3;
  G4Material* ss304 = new G4Material("MP_StainlessSteel304", density, 9);
  ss304->AddElement(Fe, 0.65);
  ss304->AddElement(Cr, 0.19);
  ss304->AddElement(Ni, 0.10);
  ss304->AddElement(Mn, 0.02);
  ss304->AddElement(Si, 0.01);
  ss304->AddElement(Mo, 0.0284);
  ss304->AddElement(C,  0.0008);
  ss304->AddElement(P,  0.0005);
  ss304->AddElement(S,  0.0003);

  density = 7.87 * g/cm3;
  a = 55.847 * g /mole;
  //G4Material* iron = new G4Material("MP_Iron", z=26, a, density);
  new G4Material("MP_Iron", z=26, a, density);

  density = 7.65 *g/cm3;
  G4Material* siliconsteel = new G4Material("MP_SiliconSteel", density, nelements=2);
  siliconsteel->AddElement(Fe, natoms=11);
  siliconsteel->AddElement(Si, natoms=1);

  density = 2.70 *g/cm3;
  G4Material* aluminum = new G4Material("MP_Aluminum", density, nelements=1);
  aluminum->AddElement(Al, natoms=1);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  //G4Material* lead = new G4Material("MP_Lead", z=82, a, density);
  new G4Material("MP_Lead", z=82, a, density);

  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  //G4Material* Cu = new G4Material("MP_Copper" , z=29., a, density);
  new G4Material("MP_Copper" , z=29., a, density);

  density = 1.032*g/cm3;
  a = 12.01*g/mole;
  //G4Material* scint = new G4Material("MP_Scint", z=6., a, density);
  new G4Material("MP_Scint", z=6., a, density);

  density = 1.e-6/760.0 * 1.29*mg/cm3; //0.001 of air density
  pressure = 1.e-6/760.0 *atmosphere;
  temperature = 293.15 *kelvin;
  a = 28.97 *g/mole;

  //G4Material* titanium = new G4Material("MP_Titanium", 22, 47.867*g/mole, 4.54*g/cm3);
  new G4Material("MP_Titanium", 22, 47.867*g/mole, 4.54*g/cm3);

  //G4Material* Vacuum = new G4Material("MP_Vacuum",z=1,a,density,kStateGas,temperature,pressure);
  new G4Material("MP_Vacuum",z=1,a,density,kStateGas,temperature,pressure);

  G4Material* Air    = new G4Material("MP_Air",    density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 79.0*perCent);
  Air->AddElement(O, 21.0*perCent);

  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


void MolPolDetectorConstruction::SetTargetThickness(G4double val){

  fTargetFullLength = val;
  fTargetVSolidTubs->SetZHalfLength( val / 2.0 );

  if (!( fTargetPhysVolume )) {
      G4cerr << "Target has not yet been constructed." << G4endl;
      return;
  }

  // tell G4RunManager that we change the geometry
  //G4RunManager::GetRunManager()->ReinitializeGeometry();  // THIS RESETS EVERYTHING, USE IF NEXT IS NOT SUFFICIENT
  G4RunManager::GetRunManager()->GeometryHasBeenModified(); // IS THIS SUFFICIENT FOR SOLID CHANGE...???
}


void MolPolDetectorConstruction::SetTargetZPosition(G4double val){
  fTargetBeamlinePz = val;

  if (!( fTargetPhysVolume )) { //Check one of the jaws.
      G4cerr << "Target has not yet been constructed." << G4endl;
      return;
  }

  fTargetPhysVolume->SetTranslation( G4ThreeVector(fTargetBeamlinePx,fTargetBeamlinePy,fTargetBeamlinePz - fHelmCoilBeamPosZ) );
  G4cout << "fTargetBeamlinePz: " << fTargetBeamlinePz << ", fHelmCoilBeamPosZ: " << fHelmCoilBeamPosZ << " [should be 6.9*cm]" << G4endl;

  // tell G4RunManager that we changed the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


void MolPolDetectorConstruction::SetDipolePbJawsGap(G4double val){
  fLeadJawGapWidth = val;

  if (!( fLeadJawsPhysicalLT && fLeadJawsPhysicalLB && fLeadJawsPhysicalRT && fLeadJawsPhysicalRB )) {
      G4cerr << "Jaws have not yet been constructed." << G4endl;
      return;
  }

  fLeadJawsPhysicalLT->SetTranslation( G4ThreeVector(-fLeadJawsXOrigin, (fLeadJawsYOrigin + (fLeadJawsHLength + (fLeadJawGapWidth/2.0))), fLeadJawsZOrigin) );
  fLeadJawsPhysicalLB->SetTranslation( G4ThreeVector(-fLeadJawsXOrigin, (fLeadJawsYOrigin - (fLeadJawsHLength + (fLeadJawGapWidth/2.0))), fLeadJawsZOrigin) );
  fLeadJawsPhysicalRT->SetTranslation( G4ThreeVector( fLeadJawsXOrigin, (fLeadJawsYOrigin + (fLeadJawsHLength + (fLeadJawGapWidth/2.0))), fLeadJawsZOrigin) );
  fLeadJawsPhysicalRB->SetTranslation( G4ThreeVector( fLeadJawsXOrigin, (fLeadJawsYOrigin - (fLeadJawsHLength + (fLeadJawGapWidth/2.0))), fLeadJawsZOrigin) );

  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void MolPolDetectorConstruction::SetDPVPSensitiveDetectors(){
  if(!DP0L){
    G4cerr << "Dipole internal virtual planes have not yet been constructed." << G4endl;
    return;
  }

  if(!fEnableDipoleInternalVPlanes) {
    DP0L->SetSensitiveDetector(nullptr);
    DP0R->SetSensitiveDetector(nullptr);
    DP1L->SetSensitiveDetector(nullptr);
    DP1R->SetSensitiveDetector(nullptr);
    DP2L->SetSensitiveDetector(nullptr);
    DP2R->SetSensitiveDetector(nullptr);
    DP3L->SetSensitiveDetector(nullptr);
    DP3R->SetSensitiveDetector(nullptr);
    DP4L->SetSensitiveDetector(nullptr);
    DP4R->SetSensitiveDetector(nullptr);
    DP5L->SetSensitiveDetector(nullptr);
    DP5R->SetSensitiveDetector(nullptr);
    DP6L->SetSensitiveDetector(nullptr);
    DP6R->SetSensitiveDetector(nullptr);
    DP7L->SetSensitiveDetector(nullptr);
    DP7R->SetSensitiveDetector(nullptr);
    DP8L->SetSensitiveDetector(nullptr);
    DP8R->SetSensitiveDetector(nullptr);
    DP9L->SetSensitiveDetector(nullptr);
    DP9R->SetSensitiveDetector(nullptr);
    return;
  }

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // Safety check to ensure the SDs are registered before applying
  if(!SDman->FindSensitiveDetector("dp0L", false) || !SDman->FindSensitiveDetector("dp0R", false)) {
    G4cerr << "Dipole internal SD(s) not registered yet" << G4endl;
    return;
  }

  DP0L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp0L") );
  DP0R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp0R") );
  DP1L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp1L") );
  DP1R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp1R") );
  DP2L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp2L") );
  DP2R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp2R") );
  DP3L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp3L") );
  DP3R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp3R") );
  DP4L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp4L") );
  DP4R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp4R") );
  DP5L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp5L") );
  DP5R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp5R") );
  DP6L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp6L") );
  DP6R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp6R") );
  DP7L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp7L") );
  DP7R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp7R") );
  DP8L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp8L") );
  DP8R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp8R") );
  DP9L->SetSensitiveDetector( SDman->FindSensitiveDetector("dp9L") );
  DP9R->SetSensitiveDetector( SDman->FindSensitiveDetector("dp9R") );

  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void MolPolDetectorConstruction::SetDipoleInternalVPlanes(G4String val){
  G4cout << "internalDipoleVPs macro value receive: " << val << G4endl;
  if(val == "true"){
    fEnableDipoleInternalVPlanes = true;
  } 
  else if(val == "false") {
    fEnableDipoleInternalVPlanes = false;
  } else {
    G4cerr << "Invalid value for dipole internal virtual planes: " << val << G4endl;
    return;
  }
  SetDPVPSensitiveDetectors();  // Update detector state (enables or disables based on flag)
  if(fEnableDipoleInternalVPlanes) {
    G4cout << "Dipole internal virtual planes are enabled." << G4endl;
  } else {
    G4cout << "Dipole internal virtual planes are disabled." << G4endl;
  }
}

void MolPolDetectorConstruction::SetFluxVPSensitiveDetectors(){
  if(!Q1ENLogical){
    G4cerr << "Flux virtual planes have not yet been constructed." << G4endl;
    return;
  }

  if(!fEnableFluxVPlanes) {
    Q1ENLogical->SetSensitiveDetector(nullptr);
    Q1EXLogical->SetSensitiveDetector(nullptr);
    Q2ENLogical->SetSensitiveDetector(nullptr);
    Q2EXLogical->SetSensitiveDetector(nullptr);
    Q3ENLogical->SetSensitiveDetector(nullptr);
    Q3EXLogical->SetSensitiveDetector(nullptr);
    Q4ENLogical->SetSensitiveDetector(nullptr);
    Q4EXLogical->SetSensitiveDetector(nullptr);
    VP1Logical->SetSensitiveDetector(nullptr);
    DipELogical->SetSensitiveDetector(nullptr);
    VP3Logical->SetSensitiveDetector(nullptr);
    return;
  }

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // Safety check to ensure the SDs are registered before applying
  if(!SDman->FindSensitiveDetector("q1en", false) || !SDman->FindSensitiveDetector("q1ex", false)) {
    G4cerr << "Flux VP SD(s) not registered yet" << G4endl;
    return;
  }

  Q1ENLogical->SetSensitiveDetector( SDman->FindSensitiveDetector("q1en") );
  Q1EXLogical->SetSensitiveDetector( SDman->FindSensitiveDetector("q1ex") );
  Q2ENLogical->SetSensitiveDetector( SDman->FindSensitiveDetector("q2en") );
  Q2EXLogical->SetSensitiveDetector( SDman->FindSensitiveDetector("q2ex") );
  Q3ENLogical->SetSensitiveDetector( SDman->FindSensitiveDetector("q3en") );
  Q3EXLogical->SetSensitiveDetector( SDman->FindSensitiveDetector("q3ex") );
  Q4ENLogical->SetSensitiveDetector( SDman->FindSensitiveDetector("q4en") );
  Q4EXLogical->SetSensitiveDetector( SDman->FindSensitiveDetector("q4ex") );
  VP1Logical->SetSensitiveDetector( SDman->FindSensitiveDetector("vp") );
  DipELogical->SetSensitiveDetector( SDman->FindSensitiveDetector("dpin") );
  VP3Logical->SetSensitiveDetector( SDman->FindSensitiveDetector("dpout") );

  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void MolPolDetectorConstruction::SetFluxVPlanes(G4String val){
  G4cout << "fluxVPs macro value receive: " << val << G4endl;
  if(val == "true"){
    fEnableFluxVPlanes = true;
  } 
  else if(val == "false") {
    fEnableFluxVPlanes = false;
  } else {
    G4cerr << "Invalid value for flux virtual planes: " << val << G4endl;
    return;
  }
  SetFluxVPSensitiveDetectors();  // Update detector state (enables or disables based on flag)
  if(fEnableFluxVPlanes) {
    G4cout << "Flux virtual planes are enabled." << G4endl;
  } else {
    G4cout << "Flux virtual planes are disabled." << G4endl;
  }
}


void MolPolDetectorConstruction::DefineGeometryCommands(){
  fMessenger = new G4GenericMessenger(this,"/MolPol/Geo/","Geometry control");

  // jaw width command
  auto& dipolePbJawWidthCmd = fMessenger->DeclareMethodWithUnit("jawWidth","mm", &MolPolDetectorConstruction::SetDipolePbJawsGap, "Set Pb Jaw Width at Dipole Entrance in cm.");
  dipolePbJawWidthCmd.SetParameterName("jawWidth", true);
  dipolePbJawWidthCmd.SetRange("jawWidth >= 0 && jawWidth <=38 ");//BETWEEN 0 and 36mm

  // target position command
  auto& targetPositionZCmd  = fMessenger->DeclareMethodWithUnit("targetPosition","mm", &MolPolDetectorConstruction::SetTargetZPosition, "Sets target position according to the /MolPol/fz");
  targetPositionZCmd.SetParameterName("targetPosition", true);
  targetPositionZCmd.SetRange("targetPosition>=50 && targetPosition<=80"); //BETWEEN 66 and 70 mm

  // target thickness command
  auto& targetThicknessCmd  = fMessenger->DeclareMethodWithUnit("targetThickness","mm", &MolPolDetectorConstruction::SetTargetThickness, "Set Pb Jaw Width at Dipole Entrance in cm.");
  targetThicknessCmd.SetParameterName("targetThickness", true);
  targetThicknessCmd.SetRange("targetThickness >= 0. && targetThickness < 100");//TODO: SET LIMIT BETWEEN 1 and 100 microns AFTER TESTING

  // dipole internal virtual planes command
  auto& dipoleInternalVPlanesCmd =
  fMessenger->DeclareMethod("internalDipoleVPs",
                            &MolPolDetectorConstruction::SetDipoleInternalVPlanes,
                            "Enable dipole internal VP sensitive detectors with 'true'");

  // flux virtual planes command
  auto& fluxVPlanesCmd =
  fMessenger->DeclareMethod("fluxVPs",
                            &MolPolDetectorConstruction::SetFluxVPlanes,
                            "Enable flux VP sensitive detectors with 'true'");

}
