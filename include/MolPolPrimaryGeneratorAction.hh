#ifndef MolPolPrimaryGeneratorAction_h
#define MolPolPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "MolPolIO.hh"
#define _USE_MATH_DEFINES
#include <cmath>
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
  G4double fZ;
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

private:
  MolPolEvent* fDefaultEvent;
  G4ParticleGun* particleGun; //pointer a to G4  class
  MolPolIO *fIO;
  PrimaryGeneratorMessenger* gunMessenger;   //messenger of this class
  G4String rndmFlag;     //flag for a rndm impact point
  double angle;//in deg
  G4String gentype;

  //Initial target multiple scattering
  remollMultScatt *fMS;

  //Levchuk effect
  G4bool fLevchukFlag;
  static const G4int eMomDistN = 150;
  G4double eMomDist[2][eMomDistN];
  void LevchukEffect();
  G4double fLEcorFac, fLEtgtPol;
  G4double SampleTargetMomentum(G4bool);
  void InitTargetMomentum();
  G4double GetTmpUnpolDist(const G4double p[8],const G4double refMom[8],const G4int k);
  G4double GetElectronStructFct(const G4double, const G4double);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
