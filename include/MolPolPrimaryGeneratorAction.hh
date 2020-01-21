#ifndef MolPolPrimaryGeneratorAction_h
#define MolPolPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "MolPolIO.hh"
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>

class MolPolDetectorConstruction;
class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;
class Simulation;
class remollMultScatt;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MolPolPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  MolPolPrimaryGeneratorAction();
  virtual ~MolPolPrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);
  void SetAngle(double deg){
    double rad = deg*M_PI/180.;
    angle = deg;
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,sin(rad),cos(rad)));
  }
  void SetIO( MolPolIO *io ){ fIO = io; }
  void rand();
  double xPos(){return particleGun->GetParticlePosition().x();}
  double yPos(){return particleGun->GetParticlePosition().y();}
  double GetAngle(){return angle;}
  G4double fXmin, fXmax, fYmin, fYmax;
  G4double fXsmear, fYsmear; 
  G4double fX, fY, fZ;
  G4double fTargLen;
  G4double fTargPol;//[0,1]
  G4double fBeamE;
  G4double fEmin, fEmax;
  G4double fthetaMin, fthetaMax;
  G4double fthetaComMin, fthetaComMax;
  G4double fphiMin, fphiMax;
  G4String fBeamPol;
  void SourceModeSet(G4int );
  void SetGenerator(G4String genname){ gentype = genname; }
  G4bool fLevchukFlag;
  G4bool fRadCorrFlag;
  G4bool fRemollMSFlag;
  G4double x1,x2,x3,x4,u1,u2,u3,u4,s;

  G4double fBeamRotZX = 0.00 * rad;
  G4double fBeamRotZY = 0.00 * rad;

private:
  MolPolEvent* fDefaultEvent;
  G4ParticleGun* particleGun; //pointer a to G4  class
  MolPolIO *fIO;
  G4String rndmFlag;     //flag for a rndm impact point
  double angle;//in deg
  G4String gentype;

  //Levchuk effect
  static const G4int eMomDistN = 150;
  G4double eMomDist[2][eMomDistN];
  void LevchukEffect();
  G4double fLEcorFac, fLEtgtPol;
  G4double SampleTargetMomentum(G4bool);
  void InitTargetMomentum();
  G4double GetTmpUnpolDist(const G4double p[8],const G4double refMom[8],const G4int k);
  G4double GetElectronStructFct(G4double&, const G4double);
  G4bool CheckLUNDFile(G4String LUNDfile_name);

  G4double fTargetA;
  G4double fTargetZ;
  G4double fTargetDensity;
  G4int    fNLUNDLines;
  std::ifstream LUNDfile;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
