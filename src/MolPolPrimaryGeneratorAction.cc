#include "MolPolPrimaryGeneratorAction.hh"

//MolPol headers
#include "MolPolMessenger.hh"

//G4 headers
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

//To use CLHEP variables (recent version of G4)
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include "MolPolIO.hh"

MolPolPrimaryGeneratorAction::MolPolPrimaryGeneratorAction()
  :rndmFlag("off")
{

  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  angle = 0;

  fBeamE = 1.063*GeV;
  particleGun->SetParticleEnergy(fBeamE*GeV);

  //Requirement (BPM values): <0.1mm

  fBeamPol = "unpol";

  fXmin = 0.0*mm;
  fXmax = 0.0*mm;

  fYmin = 0.0*mm;
  fYmax = 0.0*mm;

  fthetaMin = 1.0*deg;   //unit: rad (x* (pi/180))
  fthetaMax = 3.0*deg;

  fthetaComMin = 70*deg;
  fthetaComMax = 110*deg;

  fphiMin = -8.0*deg;
  fphiMax = 8.0*deg;

  fZ = 0.0;

  gentype = "flat";

  fDefaultEvent = new MolPolEvent();

}

MolPolPrimaryGeneratorAction::~MolPolPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
  delete fDefaultEvent;
}

void MolPolPrimaryGeneratorAction::rand(){
}

void MolPolPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fDefaultEvent->Reset();

  double xpos, ypos, zpos, thetaPos, phiPos;
  double pX, pY, pZ, eff_sigma, Azz;

  if( gentype == "moller" )
    {
      double beamE = fBeamE;
      double me = electron_mass_c2; //electron mass in MeV (0.51099906)

      double beta_com = sqrt( (beamE - me)/(beamE + me) );
      double gamma_com = 1.0/sqrt(1.0 - beta_com*beta_com);
      double e_com = me*gamma_com;
      double thcom = acos(G4RandFlat::shoot(cos(fthetaComMax), cos(fthetaComMin)));
      double phcom = G4RandFlat::shoot(fphiMin, fphiMax); //deg

      double sigma = fine_structure_const*fine_structure_const*pow(3.0+cos(thcom)*cos(thcom),2.0)*hbarc*hbarc/pow(sin(thcom),4.0)/(2.0*me*beamE); // units of area (mm^2)

      // double V = 2.0*pi*(cos(fthetaMin) - cos(fthetaMax));
      // double Z = 26.; //Z for Fe

      // Multiply by Z because we have Z electrons
      // here we must also divide by two because we are double covering
      // phasespace because of identical particles
      // eff_sigma = sigma*V*Z/2.0;

      // leave it for now
      eff_sigma = sigma;

      Azz = -((7+pow(cos(thcom),2))*pow(sin(thcom),2))/pow(3+cos(thcom)*cos(thcom),2);

      fDefaultEvent->SetEffCrossSection(eff_sigma);
      fDefaultEvent->SetAsymmetry(Azz);
      fDefaultEvent->SetThCoM(thcom);
      fDefaultEvent->SetPhCoM(phcom);

      //perpendicular and parallel momentum components
      double pperp = e_com * sin(thcom);
      double ppar  = e_com * cos(thcom);

      xpos = 0.0;
      ypos = 0.0;
      zpos = fZ;

      //generate first electron
      pX = pperp*cos(phcom);
      pY = pperp*sin(phcom);
      pZ = gamma_com*(ppar + e_com*beta_com);

      double kin = sqrt(pX*pX + pY*pY + pZ*pZ + me*me) - me;

      particleGun->SetParticlePosition( G4ThreeVector(xpos, ypos, zpos) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( pX, pY, pZ ).unit() );

      fDefaultEvent->ProduceNewParticle(G4ThreeVector(xpos, ypos, zpos),
                                        G4ThreeVector(pX, pY, pZ ),
                                        particleGun->GetParticleDefinition()->GetParticleName() );

      if(fBeamPol == "long")
        particleGun->SetParticlePolarization((G4ThreeVector(pX, pY, pZ).unit()));

      particleGun->SetParticleEnergy(kin);
      particleGun->GeneratePrimaryVertex(anEvent);

      //generate second electron
      pX = -pperp*cos(phcom);
      pY = -pperp*sin(phcom);
      pZ = gamma_com*(-ppar + e_com*beta_com);
      kin = sqrt(pX*pX + pY*pY + pZ*pZ + me*me) - me;

      particleGun->SetParticlePosition( G4ThreeVector(xpos, ypos, zpos) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( pX, pY, pZ ).unit() );

      fDefaultEvent->ProduceNewParticle(G4ThreeVector(xpos, ypos, zpos),
                                        G4ThreeVector(pX, pY, pZ ),
                                        particleGun->GetParticleDefinition()->GetParticleName() );

      fIO->SetEventData(fDefaultEvent);
      particleGun->SetParticleEnergy(kin);
      particleGun->GeneratePrimaryVertex(anEvent);


    }
  else
    {
      xpos     = G4RandFlat::shoot( fXmin, fXmax );
      ypos     = G4RandFlat::shoot( fYmin, fYmax );
      zpos     = fZ;

      thetaPos = G4RandFlat::shoot( fthetaMin, fthetaMax );
      phiPos   = G4RandFlat::shoot( fphiMin, fphiMax );

      double me = electron_mass_c2;
      double beamE = fBeamE;
      double p = sqrt( beamE*beamE - me*me); //unit:MeV

      pX = cos( phiPos ) * sin( thetaPos ) * p;
      pY = sin( phiPos ) * sin( thetaPos ) * p;
      pZ = cos( thetaPos ) * p;
      eff_sigma = 1;
      Azz = 1;

      fDefaultEvent->SetEffCrossSection(eff_sigma);
      fDefaultEvent->SetAsymmetry(Azz);

      double kinE = sqrt(pX*pX + pY*pY + pZ*pZ + me*me) - me;
      particleGun->SetParticleEnergy( kinE );

      particleGun->SetParticlePosition( G4ThreeVector(xpos, ypos, zpos) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( pX, pY, pZ ).unit() );


      fDefaultEvent->ProduceNewParticle(G4ThreeVector(xpos, ypos, zpos),
                                        G4ThreeVector(pX, pY, pZ ),
                                        particleGun->GetParticleDefinition()->GetParticleName() );
      fIO->SetEventData(fDefaultEvent);
      particleGun->GeneratePrimaryVertex(anEvent);

    }

}

void MolPolPrimaryGeneratorAction::SourceModeSet(G4int i) {
  //      SourceModeSet(0); // point to the one below with default settings = 0. // should I just use default parameters?
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MolPolPrimaryGeneratorAction::LevchukEffect(G4double targetPol) {
  G4double rand  = G4RandFlat::shoot( 0., 1. );

  G4double corFac(0), tgtPol(0);
  if (rand <= targetPol){
    G4double aaa = G4RandFlat::shoot(-1.,1.);
    G4double tgtMom = SampleTargetMomentum(1);
    //G4double pPol = 1e6 * tgtMom * aaa; //Used for printout in fortran
    corFac = 1 + tgtMom * aaa / electron_mass_c2;
    tgtPol = 1;
  }else{
    G4double aaa = G4RandFlat::shoot(-1.,1.);
    G4double tgtMom = SampleTargetMomentum(0);
    //G4double pUnpol = 1e6 * tgtMom * aaa;//Used for printout in fortran
    corFac = 1 + tgtMom * aaa / electron_mass_c2;
    tgtPol = 0;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MolPolPrimaryGeneratorAction::SampleTargetMomentum(G4bool polarizationFlag) {
  G4double rand  = G4RandFlat::shoot( 0., 1. );
  G4double eMom(0),p0(0),p1(0);
  for(int i = 0; i < eMomDistN; i++){
    if(eMomDist[polarizationFlag][i] >= rand){
      eMom = 2*i;
      p0 = eMomDist[polarizationFlag][i-1];
      p1 = eMomDist[polarizationFlag][i];
      break;
    }
  }

  //divide by 1e6 to get to GeV/c
  return (eMom + 2 * (rand - p0) / (p1 - p0)) / 1e6;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MolPolPrimaryGeneratorAction::InitTargetMomentum() {
  //array of 150 elements to hold momentum distribution
  G4double refMom[] = { 95.1, 76.5, 35.4, 5.60, 98.8, 80.2, 37.3, 5.60};

  for(int i = 0; i < eMomDistN; i++){
    eMomDist[0][i]=0;
    eMomDist[1][i]=0;
  }

  for(int i = 0; i < eMomDistN; i++){
    G4double maxMom = 2*(i+1);
    G4int nStep(0);
    if(maxMom < 5) nStep = 20;
    else nStep = 4*maxMom;

    G4double xStep = maxMom / nStep;
    G4double xSum[2] = {0, 0};
    for(int j = 0; j <= nStep; j++){
      G4double x = j * xStep;
      G4double p[8];
      for(int k = 0; k < 8 ; k++)
	p[k] = pow(x/refMom[k],2);
      G4double xTmp1 = GetTmpUnpolDist(p,refMom,0);
      G4double xTmp2 = GetTmpUnpolDist(p,refMom,4);
      G4double xInt1 = (xTmp1 + xTmp2) / 48.495;
      G4double xInt2 = ( 570281 * pow(p[2],3) / pow(1 + 9 * p[2],8) / refMom[2] +
			 570281 * pow(p[6],3) / pow(1 + 9 * p[6],8) / refMom[6])/2;
      G4double coeff = 2 + j%2 * 2;
      if(j==0 || j==nStep) coeff=1;

      xSum[0] += coeff*xInt1;
      xSum[1] += coeff*xInt2;
    }

    eMomDist[0][i] = xStep/3 * xSum[0];
    eMomDist[1][i] = xStep/3 * xSum[1];
  }

  eMomDist[0][0] = 0;
  eMomDist[1][0] = 0;
  for(int i = 1; i < eMomDistN; i++){
    eMomDist[0][i] /= eMomDist[0][eMomDistN - 1];
    eMomDist[1][i] /= eMomDist[1][eMomDistN - 1];
  }
  eMomDist[0][eMomDistN - 1] = 1;
  eMomDist[1][eMomDistN - 1] = 1;
}

G4double MolPolPrimaryGeneratorAction::GetTmpUnpolDist(const G4double p[8], const G4double refMom[8], const G4int k){
  G4double oddValue = 3.79;
  if(k == 4) oddValue = 4.705;
  return
    (2 * 10.1859 * p[0+k]) / (pow(1 + p[0+k],4) * refMom[0+k]) +
    (2 * 325.949 * p[1+k] * pow(-1 + 4 * p[1+k],2) + 6 * 1738.4 * pow(p[1+k],2)) /(pow(1 + 4 * p[1+k],6) * refMom[1+k]) +
    (2 * 275.02  * p[2+k] * pow(-1 + 4 * pow(-1 + 9 * p[2+k],2) / pow(1 + 9 * p[2+k],2),2)
     / pow(1 + 9 * p[2+k],4) + (6 * 79205.7 * pow(p[2+k],2) * pow(-1 + 9 * p[2+k],2) + oddValue * 570281 * pow(p[2+k],3) )
     /pow(1 + 9 * p[2+k],8) )/ refMom[2+k] +
    (2 * 651.899 * p[3+k] * pow(-1 + 16 * p[3+k],2) * pow(-4 + 8 * pow(-1 + 16 * p[3+k],2) / pow(1 + 16 * p[3+k],2) ,2) )
    / ( pow(1 + 16 * p[3+k],6) * refMom[3+k] );
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
