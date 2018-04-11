#ifndef MolPolIO_HH
#define MolPolIO_HH

#include "TROOT.h"
#include "TObject.h"
#include "G4Run.hh"
#include "MolPoltypes.hh"

#include "G4String.hh"

#include "MolPolEvent.hh"

class TFile;
class TTree;

class MolPolDetectorHit;
class MolPolScintDetectorHit;

#define __IO_MAXHIT 10000
#define __FILENAMELEN 255

// Units for output
#define __E_UNIT GeV
#define __L_UNIT m
#define __T_UNIT ns
#define __ANG_UNIT rad
#define __ASYMM_SCALE 1e-9 // ppb

class MolPolIO {
public:
  MolPolIO();
  ~MolPolIO();

  void SetFilename(G4String  fn);
  G4String GetFilename(){return fFilename;}

  void FillTree();
  void Flush();
  void WriteTree();

  void WriteRun();

  void InitializeTree();

private:
  TFile *fFile;
  TTree *fTree;

  char fFilename[__FILENAMELEN];

  //  Interfaces and buffers to the tree
  //  This is potentially going to get very long...

  // Event data
public:
  void SetEventData(MolPolEvent *);
private:

  G4int fNEvPart;
  G4double fEvEffXs;
  G4double fEvAsym;
  G4double fEvThCoM;
  G4double fEvPhCoM;
  G4double fUnpolWght;
  G4double fpolPlusWghtX;
  G4double fpolPlusWghtY;
  G4double fpolPlusWghtZ;
  G4double fpolMinusWghtX;
  G4double fpolMinusWghtY;
  G4double fpolMinusWghtZ;

  G4int fEvPart_PID[__IO_MAXHIT];

  G4double fEvPart_X[__IO_MAXHIT];
  G4double fEvPart_Y[__IO_MAXHIT];
  G4double fEvPart_Z[__IO_MAXHIT];
  G4double fEvPart_lX[__IO_MAXHIT];
  G4double fEvPart_lY[__IO_MAXHIT];
  G4double fEvPart_lZ[__IO_MAXHIT];
  G4double fEvPart_P[__IO_MAXHIT];
  G4double fEvPart_Px[__IO_MAXHIT];
  G4double fEvPart_Py[__IO_MAXHIT];
  G4double fEvPart_Pz[__IO_MAXHIT];
  G4double fEvPart_Th[__IO_MAXHIT];
  G4double fEvPart_Ph[__IO_MAXHIT];



  //  DetectorHit
public:
  void AddDetectorHit(MolPolDetectorHit *);
private:
  G4int fNDetHit;
  G4int fDetHit_det[__IO_MAXHIT];
  G4int fDetHit_id[__IO_MAXHIT];

  G4int fDetHit_trid[__IO_MAXHIT];
  G4int fDetHit_pid[__IO_MAXHIT];
  G4int fDetHit_gen[__IO_MAXHIT];
  G4int fDetHit_mtrid[__IO_MAXHIT];

  G4double fDetHit_X[__IO_MAXHIT];
  G4double fDetHit_Y[__IO_MAXHIT];
  G4double fDetHit_Z[__IO_MAXHIT];

  G4double fDetHit_lX[__IO_MAXHIT];
  G4double fDetHit_lY[__IO_MAXHIT];
  G4double fDetHit_lZ[__IO_MAXHIT];

  G4double fDetHit_Px[__IO_MAXHIT];
  G4double fDetHit_Py[__IO_MAXHIT];
  G4double fDetHit_Pz[__IO_MAXHIT];
  G4double fDetHit_P[__IO_MAXHIT];
  G4double fDetHit_E[__IO_MAXHIT];
  G4double fDetHit_M[__IO_MAXHIT];

  G4double fDetHit_Vx[__IO_MAXHIT];
  G4double fDetHit_Vy[__IO_MAXHIT];
  G4double fDetHit_Vz[__IO_MAXHIT];
  G4double fDetHit_Vdx[__IO_MAXHIT];
  G4double fDetHit_Vdy[__IO_MAXHIT];
  G4double fDetHit_Vdz[__IO_MAXHIT];

  //  ScintDetectorHit
public:
  void AddScintDetectorHit(MolPolScintDetectorHit *);
private:
  G4int fNScintDetHit;
  G4int fScintDetHit_det[__IO_MAXHIT];
  G4int fScintDetHit_id[__IO_MAXHIT];

  G4double fScintDetHit_edep[__IO_MAXHIT];

};

#endif//MolPolIO_HH
