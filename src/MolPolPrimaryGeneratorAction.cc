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

#include <cassert>

#include <cmath> /* Used for isnan, isinf, debugging... */

MolPolPrimaryGeneratorAction::MolPolPrimaryGeneratorAction()
  :rndmFlag("off")
{
  particleGun  = new G4ParticleGun();
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  angle = 0;

  //SET SOME DEFAULT VALUES IF NONE ARE SET IN THE MACRO
  fBeamE = 1.063*GeV;
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
  fTargLen = 1.0*mm;
  gentype = "moller";

  //radiative correction factors calculations, default to 1.0 for no usage.
  fBeamRadCorrFlag = false;
  x1 = 1.0;
  x2 = 1.0;
  fElectronsRadCorrFlag = false;
  x3 = 1.0;
  x4 = 1.0;
  u1 = 0.;
  u2 = 0.;
  u3 = 0.;
  u4 = 0.;
  s = 0.;

  //CREATE NEW MOLPOLEVENT OBJECT TO STORE EVENT INFORMATION FOR RECORDING
  fDefaultEvent = new MolPolEvent();

  //initialize target momentum distribution for Levchuck effect
  fLevchukFlag = false;
  fTargPol = 0.;
  fLEcorFac = 1.;

  //intialTargetMomentaCalculated = false;
  InitTargetMomentum();
}

MolPolPrimaryGeneratorAction::~MolPolPrimaryGeneratorAction()
{
  delete particleGun;
  delete fDefaultEvent;
}

void MolPolPrimaryGeneratorAction::rand(){
}

void MolPolPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  /*
  if(intialTargetMomentaCalculated == false){
    InitTargetMomentum();
  intialTargetMomentaCalculated = true;
  }*/

  fDefaultEvent->Reset();

  double xpos, ypos, zpos, thetaPos, phiPos;
  double pX, pY, pZ, eff_sigma, Azz;

  if( gentype == "moller" )
    {
      G4double beamE = fBeamE;

      //direction at face of Target
      G4double thcom = acos(G4RandFlat::shoot(cos(fthetaComMax), cos(fthetaComMin)));
      G4double phcom = G4RandFlat::shoot(fphiMin, fphiMax); //deg

      ///zpos: random position within the target
      G4double zpos = ( G4UniformRand() - 0.5 ) * fTargLen;

      //Multiple Scattering until the vertex position
      const G4int nTgtMat = 1;
      G4double msA[nTgtMat], msZ[nTgtMat], msThick[nTgtMat];
      msA[0] = 57.14;//
      msZ[0] = 26.43;//

      //~~ Iron density in g/cm3; divide by g/cm2 at the end to get back into G4units
      G4double ironDensity = 8.15*g/cm3; //FIXME
      msThick[0] = ((zpos + fTargLen/2) * ironDensity )/(g/cm2);
      //~~ Sample multiple scattering + angles
      G4double msth(0), msph(0);
      fMS = new remollMultScatt();
      fMS->Init(fBeamE,nTgtMat,msThick,msA,msZ);

      msth = fMS->GenerateMSPlane();
      msph = fMS->GenerateMSPlane();

      assert( !std::isnan(msth) && !std::isnan(msph) );

      //~~FIXME do we need to take care of a raster?!?!
      G4ThreeVector direction = G4ThreeVector(0,0,1.);

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
      if(fLevchukFlag) LevchukEffect();

      G4double pBeam = beamE - electron_mass_c2;

      //Internal initial state radiation
      G4double s0 = 2 * electron_mass_c2 * pBeam * fLEcorFac;
      //~~The correct scale for the bremsstrahlung is the minimum of T and U
      G4double TUmin = 0.5 * s0 * ( 1 - abs(cos(thcom)) );
      //~~The constant HBETA is 1/2 of the bremsstrahlung constant beta
      G4double hBeta = fine_structure_const/pi * (log(TUmin / pow(electron_mass_c2,2))-1);
      //~~Define the minimum value of the photon energy fraction
      G4double uMin = 1e-30;

      //~~Calculate the corresponding minimum random number (for U1,U2 gen)
      G4double rMin = pow(uMin , hBeta);
      //if(std::isinf(rMin)) G4Exception("rMin->inf","inf",FatalException,"");


      if( fBeamRadCorrFlag ){
        //~~Choose X1 and X2 to account for initial state bremsstrahlung
        //~~First, choose the photon energy fractions according to roughly the
        //~~correct distribution (U1,U2 are MUCH more important that X1,X2)
        do{
          G4double rand = G4UniformRand();
          if( rand < rMin ) u1 = uMin;
          else u1 = pow(rand, 1/hBeta);
          rand = G4UniformRand();
          if( rand < rMin ) u2 = uMin;
          else u2 = pow(rand, 1/hBeta);
          //~~Now convert them to X1,X2, and S
          x1 = 1 - u1;
          x2 = 1 - u2;
          s = s0 * x1 * x2;
        } while(s < 1e-6);//~~The cross section formally diverges at s=0, protect against
      } else {
        s = s0;
      }

      //CALCULATE LAB FRAME MOMENTA AND SCATTERING ANGLES
      G4double p1 = pBeam/2 * x1 * ( 1 + cos(thcom) );
      G4double p2 = pBeam/2 * x1 * ( 1 - cos(thcom) );
      G4double theta1 = sqrt(fLEcorFac * 2 * electron_mass_c2 * x2 * (1/p1 - 1/(x1*pBeam)) );
      G4double theta2 = sqrt(fLEcorFac * 2 * electron_mass_c2 * x2 * (1/p2 - 1/(x1*pBeam)) );

      //INTERNAL FINAL STATE RADIATION
      if( fElectronsRadCorrFlag ){
        G4double u3(0), u4(0);
        G4double rand = G4UniformRand();
        if( rand < rMin ) u3 = uMin;
        else u3 = pow(rand, 1/hBeta);
        rand = G4UniformRand();
        if( rand < rMin ) u4 = uMin;
        else u4 = pow(rand, 1/hBeta);
        x3 = 1 - u3;
        x4 = 1 - u4;
        p1 *= x3;
        p2 *= x4;
      }

      //CALCULATE THE M0LLER CROSS SECTION
      G4double cos2t = pow(cos(thcom),2);
      G4double sin2t = 1 - cos2t;

      G4double sigma = pow(fine_structure_const,2) / s * pow(hbarc,2) *
	pow(3 + cos2t,2)/pow(sin2t,2) / (2 * electron_mass_c2 * beamE);

      //CALCULATE INDV ELECTRON STRUCTURE FACTORS AND COMPLETE strFct
      G4double strFct1,strFct2,strFct3,strFct4;
      //in the case of no radiative corrections the u's initialize to zero.
      //GetElectronStructFct will set that to a nominal value to calculate
      //strFct product. strFct necessarily computes to zero unless these
      //individual factors for each are calculated first and the u's are set
      //to a non-zero value. Final values are still zero... at least for current
      //hBeta...
      strFct1 = GetElectronStructFct(u1, TUmin);
      strFct2 = GetElectronStructFct(u2, TUmin);
      strFct3 = GetElectronStructFct(u3, TUmin);
      strFct4 = GetElectronStructFct(u4, TUmin);
      G4double strFct =
        pow( strFct1*u1 , (1 - hBeta)/hBeta ) *
	pow( strFct2*u2 , (1 - hBeta)/hBeta ) *
	pow( strFct3*u3 , (1 - hBeta)/hBeta ) *
	pow( strFct4*u4 , (1 - hBeta)/hBeta );

      G4double dPhaseSpace = 1. * (cos(fthetaMax) - cos(fthetaMin));
      G4double zLum = msZ[0] * ironDensity * (zpos + fTargLen/2) * Avogadro / msA[0];
      G4double weight = zLum * dPhaseSpace * sigma * strFct;
      fDefaultEvent->fUnpolWght = weight;

      //CALCULATE THE EVENT WEIGHT FOR THE VARIOUS POLARIZATION DIRECTIONS
      G4double wTT = fLEtgtPol * sin2t / pow(3 + cos2t,2);
      G4double wZZ = wTT * (7 + cos2t);
      fDefaultEvent->fpolPlusWghtZ = weight * (1 + wZZ);
      fDefaultEvent->fpolMinusWghtZ  = weight * (1 - wZZ);

      //THE TGT SPIN IS ASSUMED TO BE ALONG THE X AXIS
      G4double wXX = wTT * sin2t * cos(2 * phcom );
      fDefaultEvent->fpolPlusWghtX = weight * (1 + wXX);
      fDefaultEvent->fpolMinusWghtX  = weight * (1 - wXX);
      G4double wYY = wTT * sin2t * sin(2 * phcom );
      fDefaultEvent->fpolPlusWghtY = weight * (1 + wYY);
      fDefaultEvent->fpolMinusWghtY  = weight * (1 - wYY);

      //STORE THE RESULTS OF THE EVENT GENERATION
      G4double tX = direction.getX();
      G4double tY = direction.getY();

      //CALCULATE EACH ELECTRONS MOMENTUM DIRECTION COMPONENTS
      //electron #1
      G4double tX1 = tX + theta1 * cos(phcom);
      G4double tY1 = tY + theta1 * sin(phcom);
      G4double arg1 = 1 - tX1*tX1 - tY1*tY1;
      G4double tZ1(0);
      assert(arg1 >= 0 and "your tX1^2 + tY1^2 > 1");
      if(arg1 > 0) tZ1 = sqrt(arg1);
      //electron #2
      G4double tX2 = tX + theta2 * cos(phcom + pi);
      G4double tY2 = tY + theta2 * sin(phcom + pi);
      G4double arg2 = 1 - tX2*tX2 - tY2*tY2;
      G4double tZ2(0);
      assert(arg2 >= 0 and "your tX2^2 + tY2^2 > 1");
      if(arg2 > 0) tZ2 = sqrt(arg2);

      //STORE EVENT QUANTITIES OF IMPORTANCE
      eff_sigma = sigma;
      Azz = -((7+pow(cos(thcom),2))*pow(sin(thcom),2))/pow(3+cos(thcom)*cos(thcom),2);

      fDefaultEvent->SetEffCrossSection(eff_sigma);
      fDefaultEvent->SetAsymmetry(Azz);
      fDefaultEvent->SetThCoM(thcom);
      fDefaultEvent->SetPhCoM(phcom);

      //GENERATE PARTICLES
      //electron #1
      G4double KE1 = electron_mass_c2 * sqrt( 1 + pow(p1/electron_mass_c2,2)) - electron_mass_c2;
      particleGun->SetParticleEnergy( KE1 );
      particleGun->SetParticlePosition( G4ThreeVector(0, 0, zpos) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( tX1, tY1, tZ1 ).unit() );
      fDefaultEvent->ProduceNewParticle(G4ThreeVector(0, 0, zpos),
                                        G4ThreeVector(tX1, tY1, tZ1 ) * p1,
                                        particleGun->GetParticleDefinition()->GetParticleName() );
      particleGun->GeneratePrimaryVertex(anEvent);
      //electron#2
      G4double KE2 = electron_mass_c2 * sqrt( 1 + pow(p2/electron_mass_c2,2)) - electron_mass_c2;
      particleGun->SetParticleEnergy( KE2 );
      particleGun->SetParticlePosition( G4ThreeVector(0, 0, zpos) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( tX2, tY2, tZ2 ).unit() );
      fDefaultEvent->ProduceNewParticle(G4ThreeVector(0, 0, zpos),
                                        G4ThreeVector(tX2, tY2, tZ2 ) * p2,
                                        particleGun->GetParticleDefinition()->GetParticleName() );
      particleGun->GeneratePrimaryVertex(anEvent);

      //FINALIZE EVENT DATA FOR RECORDING
      fIO->SetEventData(fDefaultEvent);
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
      particleGun->GeneratePrimaryVertex(anEvent);

      //FINALIZE EVENT DATA FOR RECORDING
      fIO->SetEventData(fDefaultEvent);
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
    fDefaultEvent->SetTargetMomentum( tgtMom );
    //G4double pPol = 1e6 * tgtMom * aaa; //Used for printout in fortran
    fLEcorFac = 1 + tgtMom * aaa / electron_mass_c2;
    fLEtgtPol = 1;
  }else{
    G4double aaa = G4RandFlat::shoot(-1.,1.);
    G4double tgtMom = SampleTargetMomentum(0);
    fDefaultEvent->SetTargetMomentum( tgtMom );
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

  //This is a (sort of) linear interpolation to get the momentum from the
  //dice roll; divide by 1e6 to get to GeV/c
  G4double sample = (eMom + 2 * (rand - p0) / (p1 - p0)) / 1e6;
  return sample;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MolPolPrimaryGeneratorAction::InitTargetMomentum() {
  //refMom[] is the conversion factor from dimensionless momentum to keV for each
  //successive shell (assuming screening). Iron atom specific.
  G4double refMom[] = { 95.1, 76.5, 35.4, 5.60, 98.8, 80.2, 37.3, 5.60};

  //initialize the eMomDist array, there are 150 bins (i) to be stored for polarized
  //and unpolarized electrons (hence the [2][150] array).
  for(int i = 0; i < eMomDistN; i++){
    eMomDist[0][i]=0;
    eMomDist[1][i]=0;
  }

  //Calculates and stores the momentum distribution, this was SUMMOM in DG's original code.
  for(int i = 0; i < eMomDistN; i++){
    G4double maxMom = 2*(i+1);
    G4int nStep(0);
    //Define the number of steps to calculate the integrals of the electron momentum
    //distribution for the FE atom. This must always be an even number of steps.
    if(maxMom < 5){
      nStep = 20;
    } else {
      nStep = 4*maxMom;
    }
    //xStep is the step size for the integration.
    G4double xStep = maxMom / nStep;
    //xSum is the integral
    G4double xSum[2] = {0, 0};
    //NOTE: THE INTEGRATION IS DONE ACCORDING TO SIMPSON'S COMPOSITE RULE...
    for(int j = 0; j <= nStep; j++){
      //self-explanator
      G4double x = j * xStep;
      //not sully sure that these 'momenta' are exactly... some sort of unitless momenta... :/
      G4double p[8];
      for(int k = 0; k < 8 ; k++){
        p[k] = pow(x/refMom[k],2);
      }
      //xTmp1/2 are partial calculation of the integrand (hence temporary)
      G4double xTmp1 = GetTmpUnpolDist(p,refMom,0);
      G4double xTmp2 = GetTmpUnpolDist(p,refMom,4);
      //xInt1 is the unpolarized distribution.  Although not noted in DG's original FORTRAN code, I assume that 
      //xInt2 is the polarized distribtuion.
      G4double xInt1 = (xTmp1 + xTmp2) / 48.495;
      G4double xInt2 = ( 570281 * pow(p[2],3) / pow(1 + 9 * p[2],8) / refMom[2] +
                         570281 * pow(p[6],3) / pow(1 + 9 * p[6],8) / refMom[6])/2;
      //Simpson's composite fules for coefficient selection: 2 for all (n/2 -1); 4 for all (n/2); and
      //1 for the first and last terms.
      G4double coeff = 2 + j%2 * 2;
      if(j==0 || j==nStep) coeff=1;
      //Sum up the integrands.
      xSum[0] += coeff*xInt1;
      xSum[1] += coeff*xInt2;
    }//END THE FOR(j) LOOP
    //Multiply by 1/3 of the step size to finish up Simpson's rule, and record the calculation.
    eMomDist[0][i] = xStep/3 * xSum[0];
    eMomDist[1][i] = xStep/3 * xSum[1];
  }//END THE FOR(i) LOOP

  //These remaining lines appear to be a normalizaton to get the continuous
  //momentum distribution intialization arrays on a scale of [0,1] so that a momentum can
  //be generated by rolling the die.
  eMomDist[0][0] = 0;
  eMomDist[1][0] = 0;
  for(int i = 1; i < eMomDistN; i++){
    eMomDist[0][i] /= eMomDist[0][eMomDistN - 1];
    eMomDist[1][i] /= eMomDist[1][eMomDistN - 1];
  }
  eMomDist[0][eMomDistN - 1] = 1;
  eMomDist[1][eMomDistN - 1] = 1;

  /*
  G4cout << "Election Momenta Table" << G4endl;
  G4cout << "Bin, Unpolarized,  Polarized" << G4endl;
  G4cout << "==================================" << G4endl;
  for(int i = 0; i < eMomDistN; i++){
     G4cout << i << "," << std::setprecision(9) << eMomDist[0][i] << "," << std::setprecision(9) << eMomDist[1][i] << G4endl;
  }
  */

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
