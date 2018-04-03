#include "MolPolPrimaryGeneratorAction.hh"

//MolPol headers
#include "MolPolMessenger.hh"

//G4 headers
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "remollMultScatt.hh"

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

  fBeamPol = "long";

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

  gentype = "moller";

  fDefaultEvent = new MolPolEvent();

  fTargLen = 1 * mm;
  //initialize target momentum distribution for Levchuck effect
  fLevchukFlag = false;
  fTargPol = 0.99;
  InitTargetMomentum();
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

      //direction at face of Target
      double thcom = acos(G4RandFlat::shoot(cos(fthetaComMax), cos(fthetaComMin)));
      double phcom = G4RandFlat::shoot(fphiMin, fphiMax); //deg

      double zpos = - fTargLen/2 + G4UniformRand() * fTargLen;

      //Multiple Scattering until the vertex position
      const G4int nTgtMat = 1;
      G4double msA[nTgtMat], msZ[nTgtMat], msThick[nTgtMat];
      msA[0] = 57.14;//FIXME
      msZ[0] = 26.43;//FIXME

      //~~ Iron density in g/cm3; divide by g/cm2 at the end to get back into G4units
      G4double ironDensity = 8.15; //FIXME density of Iron -- taken from DG fortran code
      msThick[0] = ((zpos*cm + fTargLen*cm/2) * ironDensity )/ (g/cm2);

      //~~ Sample multiple scattering + angles
      G4double msth(0), msph(0);
      G4double bmth(0), bmph(0);
      fMS->Init( fBeamE, nTgtMat, msThick, msA, msZ );
      msth = fMS->GenerateMSPlane();
      msph = fMS->GenerateMSPlane();

      assert( !std::isnan(msth) && !std::isnan(msph) );

      //~~FIXME do we need to take care of a raster?!?!
      G4ThreeVector direction = G4ThreeVector(0,0, 1);
      direction.rotateY( bmth); // Positive th pushes to positive X (around Y-axis)
      direction.rotateX(-bmph); // Positive ph pushes to positive Y (around X-axis)

      //~~apply MS
      direction.rotateY(msth);
      direction.rotateX(msph);

      //Energy Loss until the vertex position
      //~~ Sample beam energy based on radiation
      //~~ We do this so it doesn't affect the weighting

      G4double Ekin = fBeamE - electron_mass_c2;

      //~~ 13.84 g/cm2 is the radiation length of Fe
      G4double bt   = (msThick[0]*(g/cm2))/(13.84) * 4.0/3.0;
      G4double prob_sample, eloss, sample, env, value, ref;
      G4double fEcut = 1e-6*MeV; // copied over from remollBeamTarget.cc

      G4double Euler = 0.5772157;
      G4double prob = 1.- pow(fEcut/Ekin,bt) - bt/(bt+1.)*(1.- pow(fEcut/Ekin,bt+1.))
        + 0.75*bt/(2.+bt)*(1.- pow(fEcut/Ekin,bt+2.));
      prob = prob/(1.- bt*Euler + bt*bt/2.*(Euler*Euler+pi*pi/6.)); /* Gamma function */

      prob_sample = G4UniformRand();

      if (prob_sample <= prob) {
        do {
          sample = G4UniformRand();
          eloss = fEcut*pow(Ekin/fEcut,sample);
          env = 1./eloss;
          value = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,bt);

          sample = G4UniformRand();
          ref = value/env;
        } while (sample > ref);

        beamE = fBeamE - eloss;
        assert( fSampE > electron_mass_c2 );
      } else {
        beamE = fBeamE;
      }

      assert( beamE >= electron_mass_c2 );

      //Levchuck effect
      if(fLevchukFlag)
        LevchukEffect();

      G4double pBeam = beamE - electron_mass_c2;

      //Internal initial state radiation
      G4double s0 = 2 * electron_mass_c2 * pBeam * fLEcorFac;

      //~~The correct scale for the bremsstrahlung is the minimum of T and U
      G4double TUmin = 0.5 * s0 * ( 1 - abs(cos(direction.getTheta())) );

      //~~The constant HBETA is 1/2 of the bremsstrahlung constant beta
      G4double hBeta = fine_structure_const/pi * (log(TUmin / pow(electron_mass_c2,2))-1);

      //~~Define the minimum value of the photon energy fraction
      G4double uMin = 1e-30;

      //~~Calculate the corresponding minimum random number (for U1,U2 gen)
      G4double rMin = pow(uMin , hBeta);

      //~~Choose X1 and X2 to account for initial state bremsstrahlung
      //~~First, choose the photon energy fractions according to roughly the
      //~~correct distribution (U1,U2 are MUCH more important that X1,X2)

      G4double s(0), u1(0), u2(0);
      G4double x1(0), x2(0);
      do{
        G4double rand = G4UniformRand();
        if( rand < rMin )
          u1 = uMin;
        else
          u1 = pow(rand, 1/hBeta);

        rand = G4UniformRand();
        if( rand < rMin )
          u2 = uMin;
        else
          u2 = pow(rand, 1/hBeta);

        //~~Now convert them to X1,X2, and S
        x1 = 1 - u1;
        x2 = 1 - u2;

        s = s0 * x1 * x2;
      }while(s < 1e-6);//~~The cross section formally diverges at s=0, protect against

      G4double p1 = pBeam/2 * x1 * ( 1 + cos(direction.getTheta()) );
      G4double p2 = pBeam/2 * x2 * ( 1 - cos(direction.getTheta()) );
      G4double theta1 = sqrt(fLEcorFac * 2 * electron_mass_c2 * x2 * (1/p1 - 1/(x1*pBeam)) );
      G4double theta2 = sqrt(fLEcorFac * 2 * electron_mass_c2 * x2 * (1/p2 - 1/(x1*pBeam)) );

      //Internal final state radiation
      G4double u3(0), u4(0);
      G4double x3(0), x4(0);
      G4double rand = G4UniformRand();
      if( rand < rMin )
	u3 = uMin;
      else
	u3 = pow(rand, 1/hBeta);

      rand = G4UniformRand();
      if( rand < rMin )
	u4 = uMin;
      else
	u4 = pow(rand, 1/hBeta);

      p1 *= x3;
      p2 *= x4;

      //CALCULATE THE M0LLER CROSS SECTION
      G4double cos2t = pow(cos(direction.getTheta()),2);
      G4double sin2t = 1 - cos2t;

      G4double sigma = pow(fine_structure_const,2) / s * pow(hbarc,2) *
	pow(3 + cos2t,2)/pow(sin2t,2) / (2 * electron_mass_c2 * beamE);

      G4double strFct =
	pow(GetElectronStructFct(u1, TUmin) * u1, (1 - hBeta)/hBeta) *
	pow(GetElectronStructFct(u2, TUmin) * u2, (1 - hBeta)/hBeta) *
	pow(GetElectronStructFct(u3, TUmin) * u3, (1 - hBeta)/hBeta) *
	pow(GetElectronStructFct(u4, TUmin) * u4, (1 - hBeta)/hBeta);

      G4double dPhaseSpace = 1. * (cos(fthetaMax) - cos(fthetaMin));
      G4double zLum = msZ[0] * ironDensity * (zpos + fTargLen/2) * Avogadro / msA[0];
      G4double weight = zLum * dPhaseSpace * sigma * strFct;

      //CALCULATE THE EVENT WEIGHT FOR THE VARIOUS POLARIZATION DIRECTIONS
      G4double wTT = fLEtgtPol * sin2t / pow(3 + cos2t,2);//FIXME do we want to record these

      G4double wZZ = wTT * (7 + cos2t);
      G4double polPlusWeightZ = weight * (1 + wZZ);//FIXME do we want to record these
      G4double polMinusWeightZ = weight * (1 - wZZ);//FIXME do we want to record these

      //THE TGT SPIN IS ASSUMED TO BE ALONG THE X AXIS
      G4double wXX = wTT * sin2t * cos(2 * direction.getPhi());
      G4double polMinusWeightX = weight * (1 - wXX);//FIXME do we want to record these
      G4double polPlusWeightX = weight * (1 + wXX);//FIXME do we want to record these
      G4double wYY = wTT * sin2t * sin(2 * direction.getPhi());
      G4double polMinusWeightY = weight * (1 - wYY);//FIXME do we want to record these
      G4double polPlusWeightY = weight * (1 + wYY);//FIXME do we want to record these

      //STORE THE RESULTS OF THE EVENT GENERATION
      G4double tX = direction.getX()/pBeam;
      G4double tY = direction.getY()/pBeam;
      G4double tZ = direction.getZ()/pBeam;

      //ADD THE ELECTRON DIRECTION
      G4double tX1 = tX + theta1 * cos(direction.getPhi());
      G4double tY1 = tY + theta1 * sin(direction.getPhi());
      G4double arg = 1 - tX1*tX1 - tY1*tY1;
      G4double tZ1(0);
      assert(arg >= 0 and "your tX1^2 + tY1^2 > 1");
      if(arg>0) tZ1 = sqrt(arg);

      particleGun->SetParticlePosition( G4ThreeVector(0, 0, zpos) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( tX1, tY1, tZ1 ).unit() );
      fDefaultEvent->ProduceNewParticle(G4ThreeVector(0, 0, zpos),
                                        G4ThreeVector(tX1, tY1, tZ1 ) * p1,
                                        particleGun->GetParticleDefinition()->GetParticleName() );


      G4double tX2 = tX + theta2 * cos(direction.getPhi() + pi);
      G4double tY2 = tY + theta2 * sin(direction.getPhi() + pi);
      arg = 1 - tX2*tX2 - tY2*tY2;
      G4double tZ2(0);
      assert(arg >= 0 and "your tX2^2 + tY2^2 > 1");
      if(arg>0) tZ2 = sqrt(arg);

      particleGun->SetParticlePosition( G4ThreeVector(0, 0, zpos) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( tX2, tY2, tZ2 ).unit() );
      fDefaultEvent->ProduceNewParticle(G4ThreeVector(0, 0, zpos),
                                        G4ThreeVector(tX2, tY2, tZ2 ) * p2,
                                        particleGun->GetParticleDefinition()->GetParticleName() );

      eff_sigma = sigma;
      Azz = -((7+pow(cos(thcom),2))*pow(sin(thcom),2))/pow(3+cos(thcom)*cos(thcom),2);
      fDefaultEvent->SetEffCrossSection(eff_sigma);
      fDefaultEvent->SetAsymmetry(Azz);
      fDefaultEvent->SetThCoM(thcom);
      fDefaultEvent->SetPhCoM(phcom);

      fIO->SetEventData(fDefaultEvent);
      particleGun->SetParticleEnergy(pBeam);
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

void MolPolPrimaryGeneratorAction::LevchukEffect() {
  G4double rand  = G4RandFlat::shoot( 0., 1. );

  if (rand <= fTargPol){
    G4double aaa = G4RandFlat::shoot(-1.,1.);
    G4double tgtMom = SampleTargetMomentum(1);
    //G4double pPol = 1e6 * tgtMom * aaa; //Used for printout in fortran
    fLEcorFac = 1 + tgtMom * aaa / electron_mass_c2;
    fLEtgtPol = 1;
  }else{
    G4double aaa = G4RandFlat::shoot(-1.,1.);
    G4double tgtMom = SampleTargetMomentum(0);
    //G4double pUnpol = 1e6 * tgtMom * aaa;//Used for printout in fortran
    fLEcorFac = 1 + tgtMom * aaa / electron_mass_c2;
    fLEtgtPol = 0;
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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

///**********************************************************
///This routine calculates the electron structure function
///as given in Alexander et al., Phys.Rev.D37, 56 (1988)
///U is 1-X...for precision reasons.
///**********************************************************
G4double MolPolPrimaryGeneratorAction::GetElectronStructFct(G4double &u, const G4double s){

  //If X is too close to 1, reset it to a more sensible value
  if( u < 1e-30) u = 1e-30;

  G4double l = log(s/pow(electron_mass_c2,2));
  G4double l1 = l + 2 * log(u);
  G4double beta = 2 * fine_structure_const/pi * (l - 1);
  G4double eBeam = sqrt(s)/2;

  // Calculate e-photon structure function
  G4double dePhot = beta/2 * pow(u, beta/2 - 1)
    * ( 1 + 3 * beta/8 - pow(beta,2) / 48 * (l/3 + pow(pi,2) - 47/8.))
    - beta/4 * (2 - u) + pow(beta,2)/32
    * (4 * (2 - u) * log(1/u) - (1 + 3 * pow(1 - u , 2) ) / u * log(1. - u) - 6 + u);

  //Calculate e-e pair structure function
  G4double dee(0);
  if(u > 2*electron_mass_c2/eBeam)
    dee = pow(fine_structure_const/pi,2) *
      ( pow(u - 2 * electron_mass_c2/eBeam, beta/2) / u
	* pow(l1 - 5/3.,2) / 12 * (1 + pow(1 - u,2) + beta * (l1 - 5/3.) / 6 )
	+ pow(l,2) / 4 * (2/3. * (1 - pow(1 - u,3)) / (1 - u) + u/2 + (2 - u) * log(1 - u)));

  //Sum both contributions to the electron structure function
  return dePhot + dee;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
