#include "MolPolDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "MolPolDetector.hh"
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
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "MolPolEMFieldSetup.hh"

MolPolDetectorConstruction::MolPolDetectorConstruction():
  mEMFieldSetup(0), fCheckOverlaps(true),
  fTargetMaterial(0), fTargetFullLength(0.)
{
}

MolPolDetectorConstruction::~MolPolDetectorConstruction(){
  if(mEMFieldSetup) delete mEMFieldSetup;
}

G4VPhysicalVolume* MolPolDetectorConstruction::Construct() {

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  //// Define Visual Attributes
  G4double alphaVacuum = 0.15;
  G4double alphaMatStd = 0.50;
  G4double alphaTarget = 0.85;
  G4VisAttributes* IronVisAtt = new G4VisAttributes( G4Colour( 10./255., 10./255.,10./255.,alphaTarget) );
  G4VisAttributes* LeadVisAtt = new G4VisAttributes( G4Colour(149./255.,149./255.,100./255.,alphaMatStd) );
  G4VisAttributes* SteelVisAtt= new G4VisAttributes( G4Colour(  0./255., 80./255.,225./255.,alphaMatStd) );
  G4VisAttributes* AlumVisAtt = new G4VisAttributes( G4Colour(  0./255.,237./255.,  0./255.,alphaMatStd) );
  G4VisAttributes* BPVisAtt   = new G4VisAttributes( G4Colour(175./255.,175./255.,175./255.,alphaMatStd) );
  G4VisAttributes* VacVisAtt  = new G4VisAttributes( G4Colour(255./255.,255./255.,255./255.,alphaVacuum) );
  G4VisAttributes* CuVisAtt   = new G4VisAttributes( G4Colour(178./255.,102./255., 26./255.,alphaMatStd) );
  G4VisAttributes* ScintVisAtt= new G4VisAttributes( G4Colour(  0./255.,100./255.,100./255.,alphaMatStd) );
  G4VisAttributes* DipVisAtt  = new G4VisAttributes( G4Colour(  0./255., 80./255.,225./255.,alphaVacuum) );


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  //// Materials Setup
  G4double a, z, density, pressure, temperature;
  G4int nelements, natoms;

  G4Element* N  = new G4Element("Nitrogen"  , "N" , z=7 , a=14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"    , "O" , z=8 , a=16.00*g/mole);
  G4Element* H  = new G4Element("Hydrogen"  , "H" , z=1 , a=1.01 *g/mole);
  G4Element* C  = new G4Element("Carbon"    , "C" , z=6 , a=12.01*g/mole);

  // USED VALUES FROM WOLFRAMALPHA PERIODIC TABLE DATA | ERIC KING
  G4Element* P  = new G4Element("Phosphorus", "P" , z=15, a=30.974*g/mole);
  G4Element* S  = new G4Element("Sulfur"    , "S" , z=16, a=32.06 *g/mole);

  G4Element* Al = new G4Element("Aluminum"  , "Al", z=13, a=26.98 *g/mole);
  G4Element* Fe = new G4Element("Iron"      , "Fe", z=26, a=55.845*g/mole);
  G4Element* Si = new G4Element("Silicon"   , "Si", z=14, a=28.09 *g/mole);
  G4Element* Pb = new G4Element("Lead"      , "Pb", z=82, a=207.19*g/mole);

  // USED VALUES FROM WOLFRAMALPHA PERIODIC TABLE DATA | ERIC KING
  G4Element* Mn = new G4Element("Manganese" , "Mn", z=25, a=54.938*g/mole);
  G4Element* Cr = new G4Element("Chromium"  , "Cr", z=24, a=51.966*g/mole);
  G4Element* Ni = new G4Element("Nickel"    , "Ni", z=28, a=58.693*g/mole);
  G4Element* Mo = new G4Element("Molybdenum", "Mo", z=42, a=95.95 *g/mole);

  // INFORMATION FROM SANGHWA
  density = 7.93 *g/cm3;
  G4Material* stainlesssteel304 = new G4Material("StainlessSteel304", density, 9);
  // AddElement(elementPointer,weightedRatio)
  stainlesssteel304->AddElement(Fe, 0.65);
  stainlesssteel304->AddElement(Cr, 0.19);
  stainlesssteel304->AddElement(Ni, 0.10);
  stainlesssteel304->AddElement(Mn, 0.02);
  stainlesssteel304->AddElement(Si, 0.01);
  stainlesssteel304->AddElement(Mo, 0.0284);
  stainlesssteel304->AddElement(C,  0.0008);
  stainlesssteel304->AddElement(P,  0.0005);
  stainlesssteel304->AddElement(S,  0.0003);

  density = 7.87 * g/cm3;
  a = 55.847 * g /mole;
  G4Material* iron = new G4Material("iron", z=26, a, density);

  density = 7.65 *g/cm3;
  G4Material* siliconsteel = new G4Material("SiliconSteel", density, nelements=2);
  siliconsteel->AddElement(Fe, natoms=11);
  siliconsteel->AddElement(Si, natoms=1);

  density = 2.70 *g/cm3;
  G4Material* aluminum = new G4Material("aluminum", density, nelements=1);
  aluminum->AddElement(Al, natoms=1);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  G4Material* lead = new G4Material("lead", z=82, a, density);

  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  G4Material* Cu = new G4Material("Copper" , z=29., a, density);

  density = 1.032*g/cm3;
  a = 12.01*g/mole;
  G4Material* scint = new G4Material("scint", z=6., a, density);

  density = 1.e-6/760.0 * 1.29*mg/cm3; //0.001 of air density
  pressure = 1.e-6/760.0 *atmosphere;
  temperature = 293.15 *kelvin;  //room temperature
  a = 28.97 *g/mole;

  G4Material* Vacuum = new G4Material("Vacuum",z=1,a,density,kStateGas,temperature,pressure);
  G4Material* Air    = new G4Material("Air",    density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 79.0*perCent);
  Air->AddElement(O, 21.0*perCent);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Build world
  G4double world_x = 10*m;  G4double world_y = 10*m;  G4double world_z = 10*m;
  G4Box* world_box = new G4Box("World",world_x,world_y,world_z);
  G4LogicalVolume* world_log = new G4LogicalVolume(world_box,Vacuum,"World",0,0,0);
  world_log->SetVisAttributes(G4VisAttributes::Invisible);
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
  G4double pQ6Rin  =  0    * cm;  G4double pQ6Rout =  25.4 * cm;  G4double pQ6HL   = 38.1  * cm;  G4double pQ6Pos_z=   6.9 * cm;
  G4VSolid* Q6MagSolid = new G4Tubs( "Q6MagTubs", pQ6Rin, pQ6Rout, pQ6HL, 0.0, 360.0 * deg);
  G4LogicalVolume* Q6MagLogical = new G4LogicalVolume(Q6MagSolid, Vacuum, "Q6Mag", 0,0,0);
  Q6MagLogical->SetFieldManager(Q6FieldManager, allLocal);
  Q6MagLogical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0, G4ThreeVector(0, 0, pQ6Pos_z), Q6MagLogical, "SolenoidMag", world_log, 0, 0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Target Beam pipe - Alternative setup inside solenoid magnetic volume only.
  G4double pBPRin  = 4.78 * cm;
  G4double pBPRout = 5.08 * cm;
  G4VSolid * BPITAlum = new G4Tubs( "BPITAlum", 0., pBPRout, pQ6HL, 0.0 * deg, 360.0 * deg );
  G4VSolid * BPITVac  = new G4Tubs( "BPITVacm", 0., pBPRin,  pQ6HL, 0.0 * deg, 360.0 * deg );
  G4LogicalVolume * BPITLogical = new G4LogicalVolume(BPITAlum, aluminum, "BPalum_Targ", 0, 0, 0);
  BPITLogical->SetVisAttributes(BPVisAtt);
  G4LogicalVolume * BPVTLogical = new G4LogicalVolume(BPITVac,  Vacuum,   "BPvacm_Targ",  0, 0, 0);
  BPVTLogical->SetVisAttributes(BPVisAtt);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),BPITLogical, "BeamPipeAlmnum_Targ", Q6MagLogical, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),BPVTLogical, "BeamPipeVacuum_Targ", BPITLogical, 0, 0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Target
  G4double pMTATRin   = 0.0 * cm; G4double pMTATRout  = 1.5 * cm;   G4double pMTATHLZ = 0.0062 * cm;
  G4double pMTATPos_X = 0.0 * cm; G4double pMTATPos_Y = 0.0 * cm; G4double pMTATPos_Z = 6.9 * cm;
  G4VSolid* MTATSolid = new G4Tubs( "MTATTube", pMTATRin, pMTATRout, pMTATHLZ, 0.0, 360.0 * deg );

  G4LogicalVolume* TargetLogical = new G4LogicalVolume(MTATSolid, iron, "Target", 0, 0, 0);
  TargetLogical->SetVisAttributes(IronVisAtt);
  new G4PVPlacement(0, G4ThreeVector(0,0,pMTATPos_Z - pQ6Pos_z), TargetLogical, "Target", BPVTLogical, 0, 0, fCheckOverlaps);

  fTargetFullLength = pMTATHLZ * 2;
  fTargetMaterial = TargetLogical->GetMaterial();

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Helmholtz Coil 'Physical' Volume
  G4double pHLMZRin = 5.10 * cm;   G4double pHLMZRout = 15.0 * cm;   G4double pHLMZHLZ = 5.0 * cm;
  G4double pHLMZ1Pos_X = 0.0 * cm;   G4double pHLMZ1Pos_Y = 0.0 * cm;   G4double pHLMZ1Pos_Z = -5.7 * cm;
  G4double pHLMZ2Pos_X = 0.0 * cm;   G4double pHLMZ2Pos_Y = 0.0 * cm;   G4double pHLMZ2Pos_Z = 19.5 * cm;

  G4VSolid* HLMZ1Solid = new G4Tubs( "HLMZ1Tube", pHLMZRin, pHLMZRout, pHLMZHLZ, 0.0, 360.0 * deg );
  G4VSolid* HLMZ2Solid = new G4Tubs( "HLMZ2Tube", pHLMZRin, pHLMZRout, pHLMZHLZ, 0.0, 360.0 * deg );

  G4LogicalVolume* HLMZ1Logical = new G4LogicalVolume(HLMZ1Solid, Cu, "Helmholtz1", 0, 0, 0);
  G4LogicalVolume* HLMZ2Logical = new G4LogicalVolume(HLMZ2Solid, Cu, "Helmholtz2", 0, 0, 0);
  HLMZ1Logical->SetVisAttributes(CuVisAtt);
  HLMZ2Logical->SetVisAttributes(CuVisAtt);

  new G4PVPlacement(0, G4ThreeVector(0, 0, pHLMZ1Pos_Z - pQ6Pos_z), HLMZ1Logical, "Helmholtz1", Q6MagLogical, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0, 0, pHLMZ2Pos_Z - pQ6Pos_z), HLMZ2Logical, "Helmholtz2", Q6MagLogical, 0, 0, fCheckOverlaps);


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

  G4LogicalVolume* Q1Logical = new G4LogicalVolume(Q1Solid,siliconsteel,"Q1Logical",0,0,0);
  G4LogicalVolume* Q2Logical = new G4LogicalVolume(Q2Solid,siliconsteel,"Q2Logical",0,0,0);
  G4LogicalVolume* Q3Logical = new G4LogicalVolume(Q3Solid,siliconsteel,"Q3Logical",0,0,0);
  G4LogicalVolume* Q4Logical = new G4LogicalVolume(Q4Solid,siliconsteel,"Q4Logical",0,0,0);
  G4LogicalVolume* Q1MagLogical = new G4LogicalVolume(Q1MagSolid,Vacuum,"Q1MagLogical",0,0,0);
  G4LogicalVolume* Q2MagLogical = new G4LogicalVolume(Q2MagSolid,Vacuum,"Q2MagLogical",0,0,0);
  G4LogicalVolume* Q3MagLogical = new G4LogicalVolume(Q3MagSolid,Vacuum,"Q3MagLogical",0,0,0);
  G4LogicalVolume* Q4MagLogical = new G4LogicalVolume(Q4MagSolid,Vacuum,"Q4MagLogical",0,0,0);

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
  G4LogicalVolume* DLogical = new G4LogicalVolume ( DMagSolid, Vacuum, "DipoleMag", 0, 0, 0);
  DLogical->SetFieldManager(DFieldManager,allLocal);
  DLogical->SetVisAttributes(VacVisAtt);

  new G4PVPlacement(0,G4ThreeVector(pDMagPos_X, pDMagPos_Y - 9*cm, pDMagPos_Z), DLogical,"DipoleMag",world_log,0,0,fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DIPOLE BOX
  G4double pDBI1HLX =  6.0   * cm;  G4double pDBI1HLY = 16.5 * cm;  G4double pDBI1HLZ = 98.5  * cm;
  G4double pDBV1HLX =  5.295 * cm;  G4double pDBV1HLY = 16.0 * cm;  G4double pDBV1HLZ = 98.5  * cm;
  G4double pDBW1HLX = 11.66  * cm;  G4double pDBW1HLY = 21.4 * cm;  G4double pDBW1HLZ =  0.64 * cm;
  G4double pDBW2Rin =  0.0   * cm;  G4double pDBW2Rout=  3.0 * cm;  G4double pDBW2HL  =  0.64 * cm;
  G4double pDBW3HLX =  1.18  * cm;  G4double pDBW3HLY =  8.0 * cm;  G4double pDBW3HLZ =  0.64 * cm;
  G4double pDBW4HLX =  1.18  * cm;  G4double pDBW4HLY =  8.0 * cm;  G4double pDBW4HLZ =  0.64 * cm;
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
  G4double pDBT3Pos_X =  0.0   * cm;  G4double pDBT3Pos_Y =  0.0 * cm;  G4double pDBT3Pos_Z =  0.0  * cm;
  G4double pDBT4Pos_X =  0.0   * cm;  G4double pDBT4Pos_Y =  0.0 * cm;  G4double pDBT4Pos_Z =  0.0  * cm;
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

  G4LogicalVolume* sub4Logical = new G4LogicalVolume ( sub4, siliconsteel, "DipoleVacuumBox", 0, 0, 0);
  sub4Logical->SetVisAttributes(SteelVisAtt);

  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),sub4Logical,"DipoleVacuumBox",DLogical,0,0,fCheckOverlaps);

  G4SubtractionSolid* sub9 = new G4SubtractionSolid("sub9"  , DBW1Solid, DBW2Solid, 0, G4ThreeVector(pDBW2Pos_X, pDBW2Pos_Y, pDBW2Pos_Z) );
  G4SubtractionSolid* sub10= new G4SubtractionSolid("sub10" , sub9     , DBW3Solid, 0, G4ThreeVector(pDBW3Pos_X, pDBW3Pos_Y, pDBW3Pos_Z) );
  G4SubtractionSolid* sub11= new G4SubtractionSolid("sub11" , sub10    , DBW4Solid, 0, G4ThreeVector(pDBW4Pos_X, pDBW4Pos_Y, pDBW4Pos_Z) );

  G4LogicalVolume* sub11Logical = new G4LogicalVolume ( sub11, siliconsteel, "sub11Logical", 0, 0, 0);
  sub11Logical->SetVisAttributes(LeadVisAtt);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // COLLIMATORS
  G4double pDCOLHLX = 1.00   * cm;  G4double pDCOLHLY = 4.00 * cm;  G4double pDCOLHLZ = 3.00  * cm;
  G4double pDSLOHLX = 1.00   * cm;  G4double pDSLOHLY = 1.50 * cm;  G4double pDSLOHLZ = 3.00  * cm;
  G4double pDCOLPos_X =-4.00   * cm;  G4double pDCOLPos_Y = 9.00 * cm;  G4double pDCOLPos_Z =-92.00  * cm;
  G4double pDSLOPos_X = 0.00   * cm;  G4double pDSLOPos_Y = 0.00 * cm;  G4double pDSLOPos_Z =  0.00  * cm;
  G4VSolid* DCOLSolid = new G4Box ( "DCOLBox" , pDCOLHLX , pDCOLHLY  , pDCOLHLZ );
  G4VSolid* DSLOSolid = new G4Box ( "DSLOBox" , pDSLOHLX , pDSLOHLY  , pDSLOHLZ );
  G4SubtractionSolid* subcol = new G4SubtractionSolid("subcol", DCOLSolid, DSLOSolid, 0, G4ThreeVector(pDSLOPos_X,pDSLOPos_Y,pDSLOPos_Z) );
  G4LogicalVolume* subcolLogical = new G4LogicalVolume ( subcol, lead, "Collimator", 0, 0, 0);
  subcolLogical ->SetVisAttributes(LeadVisAtt);
  new G4PVPlacement(0,G4ThreeVector( pDCOLPos_X + pDBV1Pos_X + pDBI1Pos_X,
             pDCOLPos_Y + pDBV1Pos_Y + pDBI1Pos_Y,
             pDCOLPos_Z + pDBV1Pos_Z + pDBI1Pos_Z),subcolLogical,"Collimator1",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(-pDCOLPos_X + pDBV1Pos_X + pDBI1Pos_X,
            pDCOLPos_Y + pDBV1Pos_Y + pDBI1Pos_Y,
            pDCOLPos_Z + pDBV1Pos_Z + pDBI1Pos_Z),subcolLogical,"Collimator2",world_log,0,0,fCheckOverlaps);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Upstream Flange attached to dipole box
  G4double pDUPFLRin = 0.5 * 10.16 * cm;  G4double pDUPFLRout = 0.5 * 17.145 * cm;  G4double pDUPFLHLZ = 0.5 * 2.1336 * cm;
  G4VSolid * DUpstreamFlange = new G4Tubs( "DUpstreamFlange", pDUPFLRin , pDUPFLRout , pDUPFLHLZ  , 0.0, 360.0 * deg );
  G4LogicalVolume * DUpstreamFlangeLogical = new G4LogicalVolume(DUpstreamFlange,stainlesssteel304,"DUpstreamFlangeLogical",0,0,0);
  DUpstreamFlangeLogical->SetVisAttributes(LeadVisAtt);
  new G4PVPlacement(0,G4ThreeVector( 0 , 0 , pDMagPos_Z - 98.5*cm - pDUPFLHLZ ),DUpstreamFlangeLogical,"DipoleUpstreamFlange",world_log,0,0,0);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Dipole Magnets Physical
  G4double pDHLX   = 10.00 * cm;  G4double pDHLY   = 40.00 * cm;  G4double pDHLZ   = 76.00 * cm;
  G4double pD1Pos_X=-16.00 * cm;  G4double pD2Pos_X= 16.00 * cm;  G4double pDPos_Y = -9.00 * cm;  G4double pDPos_Z =422.80 * cm;
  G4VSolid* DSolid  = new G4Box ( "DBox"  , pDHLX , pDHLY  , pDHLZ );
  G4LogicalVolume* DMagLogical  = new G4LogicalVolume(DSolid ,siliconsteel, "DMagLogical" ,0,0,0);
  DMagLogical->SetVisAttributes(DipVisAtt);

  new G4PVPlacement(0,G4ThreeVector(pD1Pos_X,pDPos_Y,pDPos_Z),DMagLogical,"Dipole1",world_log,0,0,0);
  new G4PVPlacement(0,G4ThreeVector(pD2Pos_X,pDPos_Y,pDPos_Z),DMagLogical,"Dipole2",world_log,0,0,0);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Beam pipe - Alternative setup from solenoid to dipole box
  G4double pBPpos[10] = {pQ6Pos_z + pQ6HL,      //end of solenoid (world coord's)
                         pQ1Pos_Z - pQ1HL,      //begin Q1 (world coord's)
                         pQ1Pos_Z + pQ1HL,      //end Q1   (world coord's)
                         pQ2Pos_Z - pQ2HL,      //begin Q2 (world coord's)
                         pQ2Pos_Z + pQ2HL,      //end Q2   (world coord's)
                         pQ3Pos_Z - pQ3HL,      //begin Q3 (world coord's)
                         pQ3Pos_Z + pQ3HL,      //end Q3   (world coord's)
                         pQ4Pos_Z - pQ4HL,      //begin Q4 (world coord's)
                         pQ4Pos_Z + pQ4HL,      //end Q4   (world coord's)
                         pDMagPos_Z - pDBI1HLZ};//position of dipole mag volume minus half length of initial dipole box geometry
  G4VSolid * BPUpstreamAlum[9];
  G4VSolid * BPUpstreamVac[9];
  G4LogicalVolume * BPAlLogVol[9];
  G4LogicalVolume * BPVacLogVol[9];
  for(G4int i = 0; i < 9; i++){
    G4double startAt   = pBPpos[i];
    G4double endAt     = pBPpos[i+1];
    G4String solidName = "BPalum_" + std::to_string(i);
    G4String vacName   = "BPvacm_" + std::to_string(i);
    BPUpstreamAlum[i]  = new G4Tubs( solidName, 0.*cm, pBPRout, (pBPpos[i+1] - pBPpos[i]) / 2., 0.0 * deg, 360.0 * deg );
    BPUpstreamVac[i]   = new G4Tubs( vacName,   0.*cm, pBPRin,  (pBPpos[i+1] - pBPpos[i]) / 2., 0.0 * deg, 360.0 * deg );
    G4String solidLogName = "BPalum_logical_" + std::to_string(i);
    G4String vacLogName   = "BPvacm_logical_" + std::to_string(i);
    BPAlLogVol[i]         = new G4LogicalVolume( BPUpstreamAlum[i], aluminum, solidLogName , 0, 0, 0);
    BPAlLogVol[i]->SetVisAttributes(BPVisAtt);
    BPVacLogVol[i]        = new G4LogicalVolume( BPUpstreamVac[i] , Vacuum,   vacLogName, 0, 0, 0);
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
  G4LogicalVolume* Q1ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q1ENLogical",0,0,0);
  G4LogicalVolume* Q1EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q1EXLogical",0,0,0);
  G4LogicalVolume* Q2ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q2ENLogical",0,0,0);
  G4LogicalVolume* Q2EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q2EXLogical",0,0,0);
  G4LogicalVolume* Q3ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q3ENLogical",0,0,0);
  G4LogicalVolume* Q3EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q3EXLogical",0,0,0);
  G4LogicalVolume* Q4ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q4ENLogical",0,0,0);
  G4LogicalVolume* Q4EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q4EXLogical",0,0,0);
  G4LogicalVolume* DipELogical = new G4LogicalVolume(VBSolid, Vacuum, "DIPELogical",0,0,0);

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
  MolPolDetector* DETSD2 = new MolPolDetector("det2", 10);
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
  SDman->AddNewDetector(DETSD2);
  SDman->AddNewDetector(APPSD1);
  SDman->AddNewDetector(APPSD2);
  SDman->AddNewDetector(DETVP );
  SDman->AddNewDetector(DPIN  );
  SDman->AddNewDetector(DPOUT );

  Q1ENLogical->SetSensitiveDetector(Q1ENSD);
  Q1EXLogical->SetSensitiveDetector(Q1EXSD);
  Q2ENLogical->SetSensitiveDetector(Q2ENSD);
  Q2EXLogical->SetSensitiveDetector(Q2EXSD);
  Q3ENLogical->SetSensitiveDetector(Q3ENSD);
  Q3EXLogical->SetSensitiveDetector(Q3EXSD);
  Q4ENLogical->SetSensitiveDetector(Q4ENSD);
  Q4EXLogical->SetSensitiveDetector(Q4EXSD);
  DipELogical->SetSensitiveDetector( DPIN );

  new G4PVPlacement(0, G4ThreeVector( 0, 0, -1.*pQ1HL + VBSolidHLZ ),       Q1ENLogical,   "VP.Q1.Entr", BPVacLogVol[1], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0,     pQ1HL - VBSolidHLZ ),       Q1EXLogical,   "VP.Q1.Exit", BPVacLogVol[1], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0, -1.*pQ2HL + VBSolidHLZ ),       Q2ENLogical,   "VP.Q2.Entr", BPVacLogVol[3], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0,     pQ2HL - VBSolidHLZ ),       Q2EXLogical,   "VP.Q2.Exit", BPVacLogVol[3], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0, -1.*pQ3HL + VBSolidHLZ ),       Q3ENLogical,   "VP.Q3.Entr", BPVacLogVol[5], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0,     pQ3HL - VBSolidHLZ ),       Q3EXLogical,   "VP.Q3.Exit", BPVacLogVol[5], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0, -1.*pQ4HL + VBSolidHLZ ),       Q4ENLogical,   "VP.Q4.Entr", BPVacLogVol[7], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0,     pQ4HL - VBSolidHLZ ),       Q4EXLogical,   "VP.Q4.Exit", BPVacLogVol[7], 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( 0, 0, ((pBPpos[9]-pBPpos[8])/2.)-VBSolidHLZ),DipELogical,"VP.Dp.Entr",BPVacLogVol[8],0,0,fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DIPOLE Virtual Planes
  G4double pVP2HLX    = 6.00 * cm;   G4double pVP2HLY   = 20.00 * cm;  G4double pVP2HLZ   = 0.001   * cm;
  G4double pVP2Pos_X  = 0.00  * cm;  G4double pVP2Pos_Y = -9.00  * cm;  G4double pVP2Pos_Z = (422.8 - 98.5 - 2) * cm;
  G4double pVP3Pos_X  = 0.00  * cm;  G4double pVP3Pos_Y = -9.00  * cm;  G4double pVP3Pos_Z = (422.8 + 98.5 + 1) * cm;

  G4VSolid* VP3Solid  = new G4Box( "VP3BOX",  pVP2HLX, pVP2HLY, pVP2HLZ );
  G4LogicalVolume* VP3Logical = new G4LogicalVolume(VP3Solid, Vacuum, "VP3Logical", 0,0,0);
  VP3Logical->SetSensitiveDetector( DPOUT );
  VP3Logical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0,G4ThreeVector(pVP3Pos_X, pVP3Pos_Y, pVP3Pos_Z), VP3Logical, "VP.Dp.Exit", world_log, 0,0, fCheckOverlaps);

  // ** MOVED THIS VIRTUAL PLANE INTO THE LAST SECTION OF PRE-DIPOLE BEAMPIPE.
  //G4VSolid* VP2Solid  = new G4Box( "VP2BOX",  pVP2HLX, pVP2HLY, pVP2HLZ );
  //G4LogicalVolume* VP2Logical = new G4LogicalVolume(VP2Solid, Vacuum, "VP2Logical", 0,0,0);
  //VP2Logical->SetSensitiveDetector( DPIN );
  //VP2Logical->SetVisAttributes(VacVisAtt);
  //new G4PVPlacement(0,G4ThreeVector(pVP2Pos_X, pVP2Pos_Y, pVP2Pos_Z), VP2Logical, "VP.Dipole.Entr", world_log, 0,0, fCheckOverlaps);

    //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Virtual Planes inside of dipole box :: this block can safely be removed when no longer needed. -Eric King
  G4double DPIVPHLX = 2.295 * cm / 2;
  G4double DPIVPHLY = 15 * cm;
  G4double DPIVPHLZ = 0.00001 * cm;
  G4double DPIxloc = (3 + 2.295 * 0.5) * cm; // will be negative for LHS and pos for RHS
  G4double DPIzloc[10] = {-82.*cm,-63.*cm,-45.*cm,-27.*cm,-9.*cm,9.*cm,27.*cm,45.*cm,63.*cm,82.*cm};
  G4double DPIyloc = 0.0 * cm;
  G4VSolid * DPIVPbox = new G4Box("DipoleInternalVPlane",DPIVPHLX,DPIVPHLY,DPIVPHLZ);
  G4LogicalVolume* DP0L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP0L",0,0,0);
  G4LogicalVolume* DP0R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP0R",0,0,0);
  G4LogicalVolume* DP1L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP1L",0,0,0);
  G4LogicalVolume* DP1R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP1R",0,0,0);
  G4LogicalVolume* DP2L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP2L",0,0,0);
  G4LogicalVolume* DP2R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP2R",0,0,0);
  G4LogicalVolume* DP3L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP3L",0,0,0);
  G4LogicalVolume* DP3R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP3R",0,0,0);
  G4LogicalVolume* DP4L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP4L",0,0,0);
  G4LogicalVolume* DP4R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP4R",0,0,0);
  G4LogicalVolume* DP5L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP5L",0,0,0);
  G4LogicalVolume* DP5R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP5R",0,0,0);
  G4LogicalVolume* DP6L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP6L",0,0,0);
  G4LogicalVolume* DP6R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP6R",0,0,0);
  G4LogicalVolume* DP7L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP7L",0,0,0);
  G4LogicalVolume* DP7R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP7R",0,0,0);
  G4LogicalVolume* DP8L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP8L",0,0,0);
  G4LogicalVolume* DP8R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP8R",0,0,0);
  G4LogicalVolume* DP9L = new G4LogicalVolume(DPIVPbox, Vacuum, "DP9L",0,0,0);
  G4LogicalVolume* DP9R = new G4LogicalVolume(DPIVPbox, Vacuum, "DP9R",0,0,0);
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
  DP0L->SetSensitiveDetector(DP0Lsd);
  DP0R->SetSensitiveDetector(DP0Rsd);
  DP1L->SetSensitiveDetector(DP1Lsd);
  DP1R->SetSensitiveDetector(DP1Rsd);
  DP2L->SetSensitiveDetector(DP2Lsd);
  DP2R->SetSensitiveDetector(DP2Rsd);
  DP3L->SetSensitiveDetector(DP3Lsd);
  DP3R->SetSensitiveDetector(DP3Rsd);
  DP4L->SetSensitiveDetector(DP4Lsd);
  DP4R->SetSensitiveDetector(DP4Rsd);
  DP5L->SetSensitiveDetector(DP5Lsd);
  DP5R->SetSensitiveDetector(DP5Rsd);
  DP6L->SetSensitiveDetector(DP6Lsd);
  DP6R->SetSensitiveDetector(DP6Rsd);
  DP7L->SetSensitiveDetector(DP7Lsd);
  DP7R->SetSensitiveDetector(DP7Rsd);
  DP8L->SetSensitiveDetector(DP8Lsd);
  DP8R->SetSensitiveDetector(DP8Rsd);
  DP9L->SetSensitiveDetector(DP9Lsd);
  DP9R->SetSensitiveDetector(DP9Rsd);
  //////////////////////////////////////////////////////////////////////////////




  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Rotations information to be used for subsequent sections

  /* RMATR07  90.   0.  80.  90.  10.  270.     Moller detector
   * RMATR09  90. 270.   0.   0.  90.  180.     I is oppos Y,II along Z,III oppos X */

  G4RotationMatrix* pRot9 = new G4RotationMatrix();
  pRot9->rotateY(90.*deg);
  pRot9->rotateZ(90.*deg);

  G4RotationMatrix* pRot7 = new G4RotationMatrix();
  pRot7->rotateX(-7.3*deg);




  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Shielding

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
  G4double pS3LDHLX   = 20.0 * cm;  G4double pS3LDHLY   = 2.0  * cm;  G4double pS3LDHLZ = 57.0 * cm;

  G4double pS1LDPos_X   = 0.0  * cm;  G4double pS1LDPos_Y   = -29.0 * cm;  G4double pS1LDPos_Z = 537.0 * cm;
  G4double pS1H1Pos_X   = -4.1 * cm;  G4double pS1H1Pos_Y   = 9.5   * cm;  G4double pS1H1Pos_Z = 0     * cm;
  G4double pS1H2Pos_X   = 4.1  * cm;  G4double pS1H2Pos_Y   = 9.5   * cm;  G4double pS1H2Pos_Z = 0     * cm;
  G4double pS2LDPos_X   = 0.0  * cm;  G4double pS2LDPos_Y   = -7.5  * cm;  G4double pS2LDPos_Z = 610.0 * cm;
  G4double pS21BPos_X   = 0.0  * cm;  G4double pS21BPos_Y   = -13.3 * cm;  G4double pS21BPos_Z = 620.0 * cm;
  G4double pS21HPos_X   = 0.0  * cm;  G4double pS21HPos_Y   = -0.5  * cm;  G4double pS21HPos_Z = 0.0   * cm;
  G4double pS3LDPos_X   = 0.0  * cm;  G4double pS3LDPos_Y   = -7.5  * cm;  G4double pS3LDPos_Z = 723.0 * cm;

  G4VSolid* S1LDSolid  = new G4Box ( "S1LDBox"  , pS1LDHLX, pS1LDHLY,  pS1LDHLZ );
  G4VSolid* S1H1Solid  = new G4Para( "S1H1Para" , pS1H1HLX, pS1H1HLY,  pS1H1HLZ, pS1H1alpha, pS1H1theta, pS1H1phi);
  G4VSolid* S1H2Solid  = new G4Para( "S1H2Para" , pS1H2HLX, pS1H2HLY,  pS1H2HLZ, pS1H2alpha, pS1H2theta, pS1H2phi);
  G4VSolid* S2LDSolid  = new G4Box ( "S2LDBox"  , pS2LDHLX, pS2LDHLY,  pS2LDHLZ );
  G4VSolid* S21BSolid  = new G4Box ( "S21BBox"  , pS21BHLX, pS21BHLY,  pS21BHLZ );
  G4VSolid* S21HSolid  = new G4Box ( "S21HBox"  , pS21HHLX, pS21HHLY,  pS21HHLZ );
  G4VSolid* S3LDSolid  = new G4Box ( "S3LDBox"  , pS3LDHLX, pS3LDHLY,  pS3LDHLZ );

  G4SubtractionSolid* subSHL1 = new G4SubtractionSolid("subSHL1", S1LDSolid, S1H1Solid, pRot9, G4ThreeVector(pS1H1Pos_X, pS1H1Pos_Y, pS1H1Pos_Z) );
  G4SubtractionSolid* subSHL2 = new G4SubtractionSolid("subSHL2", subSHL1, S1H2Solid, pRot9, G4ThreeVector(pS1H2Pos_X, pS1H2Pos_Y, pS1H2Pos_Z) );
  G4SubtractionSolid* subSHL3 = new G4SubtractionSolid("subSHL3", S21BSolid, S21HSolid, 0, G4ThreeVector(pS21HPos_X, pS21HPos_Y, pS21HPos_Z) );

  // FIXME!!!!!  Something here may require attention. Are materials correct? Does not match VisAttributes
  G4LogicalVolume* SHL1Logical = new G4LogicalVolume(subSHL2,   lead, "Shield1", 0,0,0);
  G4LogicalVolume* SHL2Logical = new G4LogicalVolume(S2LDSolid, lead, "Shield2", 0,0,0);
  G4LogicalVolume* SHL3Logical = new G4LogicalVolume(subSHL3,   aluminum, "Shield3", 0,0,0);
  G4LogicalVolume* SHL4Logical = new G4LogicalVolume(S3LDSolid, lead, "Shield4", 0,0,0);
  // FIXME!!!!!  Something here may require attention. Are VisAttributes correct? Does not match materials.
  SHL1Logical->SetVisAttributes(LeadVisAtt);
  SHL2Logical->SetVisAttributes(AlumVisAtt);
  SHL3Logical->SetVisAttributes(LeadVisAtt);
  SHL4Logical->SetVisAttributes(LeadVisAtt);

  new G4PVPlacement(0, G4ThreeVector(pS1LDPos_X, pS1LDPos_Y, pS1LDPos_Z), SHL1Logical, "Shield1", world_log, 0,0,fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(pS2LDPos_X, pS2LDPos_Y, pS2LDPos_Z), SHL2Logical, "Shield2", world_log, 0,0,fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(pS21BPos_X, pS21BPos_Y, pS21BPos_Z), SHL3Logical, "Shield3", world_log, 0,0,fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(pS3LDPos_X, pS3LDPos_Y, pS3LDPos_Z), SHL4Logical, "Shield4", world_log, 0,0,fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DETECTOR AND BOX
  G4double pMDBXHLX   = 37.00 * cm;  G4double pMDBXHLY   = 47.10 * cm;  G4double pMDBXHLZ   = 61.20 * cm;
  G4double pMDBAHLX   = 18.50 * cm;  G4double pMDBAHLY   = 29.80 * cm;  G4double pMDBAHLZ   = 40.00 * cm;
  G4double pMDBWHLX   = 13.70 * cm;  G4double pMDBWHLY   = 11.50 * cm;  G4double pMDBWHLZ   =  6.45 * cm;
  G4double pMDBWalpha =  6.0  *deg;  G4double pMDBWtheta = 0.0  *deg;  G4double pMDBWphi   =  0.0  *deg;
  G4double pMDBLHLX   = 13.70 * cm;  G4double pMDBLHLY   = 11.50 * cm;  G4double pMDBLHLZ   =  6.45 * cm;
  G4double pMDBLalpha = 10.0  *deg;  G4double pMDBLtheta = 0.0  *deg;  G4double pMDBLphi   =  0.0  *deg;
  G4double pMDETHLX   =  9.20 * cm;  G4double pMDETHLY   = 17.00 * cm;  G4double pMDETHLZ   = 34.00 * cm;
  G4double pDLGBHLX   =  8.20 * cm;  G4double pDLGBHLY   = 16.20 * cm;  G4double pDLGBHLZ   = 20.10 * cm;

  G4double pMDBXPos_X   =  0.00 * cm;  G4double pMDBXPos_Y   =-46.86 * cm;  G4double pMDBXPos_Z   =724.00 * cm; //10cm up
  G4double pMDBAPos_X   =  0.00 * cm;  G4double pMDBAPos_Y   =  1.10 * cm;  G4double pMDBAPos_Z   =  1.80 * cm;
  G4double pMDBWPos_X   =  0.00 * cm;  G4double pMDBWPos_Y   =  8.50 * cm;  G4double pMDBWPos_Z   =-49.70 * cm;
  G4double pMDBLPos_X   =  0.00 * cm;  G4double pMDBLPos_Y   =  7.50 * cm;  G4double pMDBLPos_Z   =-49.70 * cm;
  G4double pMDETPos_X   =  0.00 * cm;  G4double pMDETPos_Y   =  0.00 * cm;  G4double pMDETPos_Z   =  0.00 * cm;
  G4double pDLGBPos_X   =  0.00 * cm;  G4double pDLGBPos_Y   =  0.00 * cm;  G4double pDLGBPos_Z   =-11.00 * cm;

  G4VSolid* MDBXSolid  = new G4Box ( "MDBXBox"  , pMDBXHLX, pMDBXHLY, pMDBXHLZ );
  G4VSolid* MDBASolid  = new G4Box ( "MDBABox"  , pMDBAHLX, pMDBAHLY, pMDBAHLZ );
  G4VSolid* MDBWSolid  = new G4Para( "MDBWPara" , pMDBWHLX, pMDBWHLY, pMDBWHLZ, pMDBWalpha, pMDBWtheta, pMDBWphi);
  G4VSolid* MDBLSolid  = new G4Para( "MDBLPara" , pMDBLHLX, pMDBLHLY, pMDBLHLZ, pMDBLalpha, pMDBLtheta, pMDBLphi);
  G4VSolid* MDETSolid  = new G4Box ( "MDETBox"  , pMDETHLX, pMDETHLY, pMDETHLZ );
  G4VSolid* DLGBSolid  = new G4Box ( "DLGBBox"  , pDLGBHLX, pDLGBHLY, pDLGBHLZ );

  G4SubtractionSolid* sub5 = new G4SubtractionSolid("sub5", MDBXSolid, MDBASolid, 0    , G4ThreeVector(pMDBAPos_X, pMDBAPos_Y, pMDBAPos_Z) );
  G4SubtractionSolid* sub6 = new G4SubtractionSolid("sub6", sub5     , MDBWSolid, pRot9, G4ThreeVector(pMDBWPos_X, pMDBWPos_Y, pMDBWPos_Z) );
  G4SubtractionSolid* sub7 = new G4SubtractionSolid("sub7", sub6     , MDBLSolid, pRot9, G4ThreeVector(pMDBLPos_X, pMDBLPos_Y, pMDBLPos_Z) );
  G4SubtractionSolid* sub8 = new G4SubtractionSolid("sub8", sub7     , MDETSolid, pRot7, G4ThreeVector(pMDETPos_X, pMDETPos_Y, pMDETPos_Z) );

  //  G4LogicalVolume* DETLogical = new G4LogicalVolume(DLGBSolid, aluminum, "DETLogical",0,0,0);
  //  DETLogical->SetVisAttributes(AlumVisAtt);
  G4LogicalVolume* DETLogical = new G4LogicalVolume(DLGBSolid, Vacuum, "DETLogical",0,0,0);
  DETLogical->SetSensitiveDetector(DETSD);

  G4LogicalVolume* sub8Logical = new G4LogicalVolume ( sub8, lead, "DetectorBox", 0, 0, 0);
  G4VisAttributes *sub8VisAtt(AlumVisAtt);
  sub8VisAtt->SetForceWireframe(true);
  sub8Logical->SetVisAttributes(sub8VisAtt);
  sub8Logical->SetSensitiveDetector(DETSD2 );

  new G4PVPlacement(0,G4ThreeVector(pMDBXPos_X, pMDBXPos_Y, pMDBXPos_Z),sub8Logical,"sub8",world_log,0,0,fCheckOverlaps);

  // Note: The name of this PV is used in SteppingAction. Please do not change.
  new G4PVPlacement(pRot7,G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pDLGBPos_X,
            pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pDLGBPos_Y,
            pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pDLGBPos_Z),
            DETLogical,"virtualBoundaryPhys_det",world_log,0,0,fCheckOverlaps);
  /*
    C --     Detector box
    GPARVOL80  'MDBX'  13  'HALL'    0.  -50.1 724.   0  'BOX '  3  37.   47.1   61.2  ***7cm up
    GPARVOL81  'MDBA'  15  'MDBX'    0.    1.1   1.8  0  'BOX '  3  18.5  29.8   40.0
    GPARVOL82  'MDBW'  15  'MDBX'    0.    8.5 -49.7  9  'PARA'  6  13.7  11.5   6.45 6.0 0. 0.
    GPARVOL83  'MDBL'  15  'MDBX'    0.    7.5 -49.7  9  'PARA'  6  13.7  11.5   6.45 10.0 0. 0.
    C --     Detector
    GPARVOL90  'MDET'  15  'MDBA'    0.    0.    0.   7  'BOX '  3   9.2  17.0    34.0
    GPARVOL91  'DLGB'   9  'MDET'    0.    0.  -11.   0  'BOX '  3   8.2  16.2    20.1

    C --  Syracuse hodoscope, SPACAL, auxill. artificial plane
    C
    GPARVOL101 'HOD1'  15  'MDET'    0.    0.  -32.5  0  'BOX '  3   9.2  15.5     0.75
    GPARVOL102 'HOD2'   9  'MDET'    0.    0.  -16.0  0  'BOX '  3   9.2  16.5    15.2
    GPARVOL103 'HOD3'  15  'MDET'    0.    0.  -33.5  0  'BOX '  3   9.2  17.0     0.1
  */


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DETECTOR Virtual Plane
  G4double pVP1HLX    = 37.00 * cm;  G4double pVP1HLY   = 50.00 * cm;  G4double pVP1HLZ   = 0.0001   * cm;
  G4double pVP1Pos_X  = 0.00  * cm;  G4double pVP1Pos_Y = -46.86* cm;  G4double pVP1Pos_Z = 660.0 * cm;
  G4VSolid* VP1Solid  = new G4Box( "VP1BOX",  pVP1HLX, pVP1HLY, pVP1HLZ );
  G4LogicalVolume* VP1Logical = new G4LogicalVolume(VP1Solid, Vacuum, "VPDetector", 0,0,0);
  VP1Logical->SetSensitiveDetector( DETVP );
  VP1Logical->SetVisAttributes(VacVisAtt);

  //new G4PVPlacement(0,G4ThreeVector(pVP1Pos_X, pVP1Pos_Y, pVP1Pos_Z), VP1Logical, "VP.Detector.Entr", world_log, 0,0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(pVP1Pos_X, pVP1Pos_Y, pMDBXPos_Z - pMDBXHLZ - pVP1HLZ ), VP1Logical, "VP.Detector.Entr", world_log, 0,0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // FIXME!!!!! This section of the geometry seems incomplete.
  // Hodoscope
  G4double pHOD1HLX   = 9.2 * cm;   G4double pHOD1HLY   = 15.5 * cm;   G4double pHOD1HLZ   = 0.75 * cm;
  G4double pHOD2HLX   = 9.2 * cm;   G4double pHOD2HLY   = 16.5 * cm;   G4double pHOD2HLZ   = 15.2 * cm;
  G4double pHOD3HLX   = 9.2 * cm;   G4double pHOD3HLY   = 17.5 * cm;   G4double pHOD3HLZ   = 0.1  * cm;

  G4double pHOD1Pos_X = 0.0 * cm;   G4double pHOD1Pos_Y = 0.0  * cm;   G4double pHOD1Pos_Z = -32.5* cm;
  G4double pHOD2Pos_X = 0.0 * cm;   G4double pHOD2Pos_Y = 0.0  * cm;   G4double pHOD2Pos_Z = -16.0* cm;
  G4double pHOD3Pos_X = 0.0 * cm;   G4double pHOD3Pos_Y = 0.0  * cm;   G4double pHOD3Pos_Z = -33.5* cm;

  G4VSolid* HOD1Solid = new G4Box( "HOD1Box", pHOD1HLX, pHOD1HLY, pHOD1HLZ );
  G4VSolid* HOD2Solid = new G4Box( "HOD2Box", pHOD2HLX, pHOD2HLY, pHOD2HLZ );
  G4VSolid* HOD3Solid = new G4Box( "HOD3Box", pHOD3HLX, pHOD3HLY, pHOD3HLZ );

  G4double offset = 3.5 * cm;
  // G4double offset = 0 * cm;

  G4double pAPP1HLX   = 2.0 * cm;   G4double pAPP1HLY   = 15.5 * cm;          G4double pAPP1HLZ   = 0.75 * cm;
  G4double pAPP1Pos_X = 4.4 * cm;   G4double pAPP1Pos_Y = 0.0 * cm + offset;  G4double pAPP1Pos_Z = 0.0 * cm;

  G4VSolid* APP1LSolid = new G4Box( "APP1LBOX", pAPP1HLX, pAPP1HLY, pAPP1HLZ );
  G4VSolid* APP1RSolid = new G4Box( "APP1RBOX", pAPP1HLX, pAPP1HLY, pAPP1HLZ );

  // G4LogicalVolume* APP1LLogical = new G4LogicalVolume(APP1LSolid, scint, "APP1LLogical", 0,0,0);
  // APP1LLogical->SetSensitiveDetector(APPSD1);
  G4LogicalVolume* APP1LLogical = new G4LogicalVolume(APP1LSolid, Vacuum, "APP1LLogical", 0,0,0);
  APP1LLogical->SetVisAttributes(ScintVisAtt);

  // G4LogicalVolume* APP1RLogical = new G4LogicalVolume(APP1RSolid, scint, "APP1RLogical", 0,0,0);
  // APP1RLogical->SetSensitiveDetector(APPSD2);
  G4LogicalVolume* APP1RLogical = new G4LogicalVolume(APP1RSolid, Vacuum, "APP1RLogical", 0,0,0);
  APP1RLogical->SetVisAttributes(ScintVisAtt);

  new G4PVPlacement(pRot7, G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pHOD1Pos_X + pAPP1Pos_X,
             pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pHOD1Pos_Y + pAPP1Pos_Y,
             pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pHOD1Pos_Z + pAPP1Pos_Z),
             APP1LLogical, "APP1L", world_log, 0,0, fCheckOverlaps);

  new G4PVPlacement(pRot7, G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pHOD1Pos_X - pAPP1Pos_X,
             pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pHOD1Pos_Y + pAPP1Pos_Y,
             pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pHOD1Pos_Z + pAPP1Pos_Z),
             APP1RLogical, "APP1R", world_log, 0,0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // HODOSCOPE Virtual Planes
  G4double pVPHODHLZ = 0.0001 * cm;
  G4double pVPHODPos_Z = -1 * pAPP1HLZ;
  // Rotated coordinates to keep VP's on front of HODs
  G4double VPHODPos_Zr = (pVPHODPos_Z - pVPHODHLZ) * cos(7.3*deg);
  G4double VPHodPos_Yr = -1 * (pVPHODPos_Z - pVPHODHLZ) * sin(7.3*deg);

  G4VSolid* VPHOD1Solid = new G4Box( "VPHOD1", pAPP1HLX, pAPP1HLY, pVPHODHLZ);
  G4VSolid* VPHOD2Solid = new G4Box( "VPHOD2", pAPP1HLX, pAPP1HLY, pVPHODHLZ);

  G4LogicalVolume* VPHOD1Logical = new G4LogicalVolume(VPHOD1Solid, Vacuum, "VPHOD1Logical", 0,0,0);
  VPHOD1Logical->SetSensitiveDetector(APPSD1);
  VPHOD1Logical->SetVisAttributes(VacVisAtt);
  G4LogicalVolume* VPHOD2Logical = new G4LogicalVolume(VPHOD2Solid, Vacuum, "VPHOD2Logical", 0,0,0);
  VPHOD2Logical->SetSensitiveDetector(APPSD2);
  VPHOD2Logical->SetVisAttributes(VacVisAtt);

  new G4PVPlacement(pRot7, G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pHOD1Pos_X + pAPP1Pos_X,
             pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pHOD1Pos_Y + pAPP1Pos_Y + VPHodPos_Yr,
             pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pHOD1Pos_Z + pAPP1Pos_Z + VPHODPos_Zr),
             VPHOD1Logical, "VP.Hodo.Left", world_log, 0,0, fCheckOverlaps);

  new G4PVPlacement(pRot7, G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pHOD1Pos_X - pAPP1Pos_X,
             pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pHOD1Pos_Y + pAPP1Pos_Y + VPHodPos_Yr,
             pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pHOD1Pos_Z + pAPP1Pos_Z + VPHODPos_Zr),
             VPHOD2Logical, "VP.Hodo.Right", world_log, 0,0, fCheckOverlaps);


  return world_phys;
}
