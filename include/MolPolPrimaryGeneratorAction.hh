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

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
