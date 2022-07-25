#include "MolPolPrimaryGeneratorAction.hh"

//MolPol headers
#include "MolPolMessenger.hh"
#include "MolPolDetectorConstruction.hh"
#include "MolPolIO.hh"

//G4 headers
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "remollMultScatt.hh"
#include "G4Material.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

//To use CLHEP variables (recent version of G4)
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

#include "G4ios.hh"

#include <cassert>

#include <iostream>
#include <fstream>
#include <string>
#include <cmath> /* Used for isnan, isinf, debugging... */

#define _unused(x) ((void)(x)) /* Avoid unused variable warning for variables actually used in assert() */

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

  //Flat generation range
  fXmin = 0.0*mm;
  fXmax = 0.0*mm;
  fYmin = 0.0*mm;
  fYmax = 0.0*mm;

  //Gaussian XY beam profile for the Moller generator
  fX = 0.0*mm;
  fY = 0.0*mm;
  fXsmear = 0.0*mm;
  fYsmear = 0.0*mm;
  fZ = 6.9*cm;

  fthetaMin = 1.0*deg;   //unit: rad (x* (pi/180))
  fthetaMax = 3.0*deg;
  fthetaComMin = 70*deg;
  fthetaComMax = 110*deg;
  fphiMin = -8.0*deg;
  fphiMax = 8.0*deg;

  gentype = "moller";

  //TARGET VARIABLES
  fTargLen = .0124*mm;
  fTargetA = 55.847;
  fTargetZ = 26;
  fTargetDensity = 7.87 * g/cm3;

  //CREATE NEW MOLPOLEVENT OBJECT TO STORE EVENT INFORMATION FOR RECORDING
  fDefaultEvent = new MolPolEvent();

  //LEVCHUCK FLAG AND VALUES - INITIALIZE FLAG TO 'true' BY DEFAULT AND POLARIZATION TO SOME VALUE
  //ADDITIONALLY, fLEcorFac MUST INITIALIZE TO 1.
  fLevchukFlag = true;
  fTargPol = 0.;
  fLEcorFac = 1.;
  //RADIATIVE CORRECTION CALCULATIONS - DEFAULT TO TRUE UNLESS OVERRIDDEN IN MACRO
  fRadCorrFlag = true;
  //MOLPOL MSC ON/OFF FLAG -- DEFAULTS TO TRUE/ON
  fRemollMSFlag = true;

  //INITIALIZE TARGET MOMENTUM DISTRIBUTION FOR LEVCHUK EFFECT
  InitTargetMomentum();

  fNLUNDLines = 0;
}

MolPolPrimaryGeneratorAction::~MolPolPrimaryGeneratorAction()
{
  delete particleGun;
  delete fDefaultEvent;

  if( gentype == "moller") LUNDfile.close();
}

void MolPolPrimaryGeneratorAction::rand(){
}

void MolPolPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

 fDefaultEvent->Reset();

  if( gentype == "moller" )
    {
      G4double beamE = fBeamE;

      //direction at face of Target
      G4double thcom = acos(G4RandFlat::shoot(cos(fthetaComMax), cos(fthetaComMin)));
      G4double phcom = G4RandFlat::shoot(fphiMin, fphiMax); //deg

      //zpos: random position within the target of length fTargLen (local coordinate)
      // vertex z position is calculated as vtx_z = zpos + fZ (target center) later
      G4double xpos     = G4RandGauss::shoot( fX, fXsmear );
      G4double ypos     = G4RandGauss::shoot( fY, fYsmear );
      G4double zpos = ( G4UniformRand() - 0.5 ) * fTargLen;

      //Multiple Scattering until the vertex position
      const G4int nTgtMat = 1;
      G4double msA[nTgtMat], msZ[nTgtMat], msThick[nTgtMat];
      msA[0] = fTargetA;
      msZ[0] = fTargetZ;

      //~~ Iron density in g/cm3; divide by g/cm2 at the end to get back into G4units
      G4double ironDensity = fTargetDensity;
      msThick[0] = ((zpos + fTargLen/2) * ironDensity ); // in (g/cm2);
      //~~ Sample multiple scattering + angles
      G4double msth(0), msph(0);

  	  //G4cout << "Remoll MS next..." << G4endl;
  	  //if(!fRemollMSFlag) G4cout << "Remoll MSc off..." << G4endl;
      if(fRemollMSFlag){
                  remollMultScatt fMS;
  		  fMS.Init(fBeamE,nTgtMat,msThick,msA,msZ);
  		  msth = fMS.GenerateMSPlane();
  		  msph = fMS.GenerateMSPlane();
  		  //G4cout << "Remoll MSc on..." << G4endl;
      }

      assert( !std::isnan(msth) && !std::isnan(msph) );

      //No Raster (actual data taking done mostly without raster)
      G4ThreeVector direction = G4ThreeVector(0,0,1.);

      //Rotate the directional vector to match that specified for the beam
      //Remember we're rotating the axis, so an angle from ZtoX is made by turning the Y-axis, good for very small angles only
      direction.rotateX(-fBeamRotZY);
      direction.rotateY( fBeamRotZX);

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
      G4double prob = 1.- pow(fEcut/Ekin,bt) - bt/(bt+1.)*(1.- pow(fEcut/Ekin,bt+1.)) + 0.75*bt/(2.+bt)*(1.- pow(fEcut/Ekin,bt+2.));
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
        assert( beamE > electron_mass_c2 );
      } else {
        beamE = fBeamE;
      }

      assert( beamE >= electron_mass_c2 );

      //Levchuck effect
      if(fLevchukFlag) LevchukEffect();

      G4double pBeam = sqrt(beamE * beamE - electron_mass_c2 * electron_mass_c2);

      //Internal initial state radiation
      G4double s0 = 2 * electron_mass_c2 * pBeam * fLEcorFac;
      //~~The correct scale for the bremsstrahlung is the minimum of T and U
      G4double TUmin = 0.5 * s0 * ( 1 - fabs(cos(thcom)) );
      //~~The constant HBETA is 1/2 of the bremsstrahlung constant beta
      G4double hBeta = fine_structure_const/pi * (log(TUmin / pow(electron_mass_c2,2))-1);
      //~~Define the minimum value of the photon energy fraction
      G4double uMin = 1e-30;

      //~~Calculate the corresponding minimum random number (for U1,U2 gen)
      G4double rMin = pow(uMin , hBeta);
      //if(std::isinf(rMin)) G4Exception("rMin->inf","inf",FatalException,"");

      G4double s(0), u1(0), u2(0);
      G4double x1(1), x2(1);
      if( fRadCorrFlag ){
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
      G4double u3(0), u4(0);
      G4double x3(1), x4(1);
      if( fRadCorrFlag ){
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
      G4double strFct1,strFct2,strFct3,strFct4,strFct;
      if(fRadCorrFlag){
        strFct1 = GetElectronStructFct(u1, TUmin);
        strFct2 = GetElectronStructFct(u2, TUmin);
        strFct3 = GetElectronStructFct(u3, TUmin);
        strFct4 = GetElectronStructFct(u4, TUmin);
        strFct =
          strFct1 * pow(u1 , (1. - hBeta) ) / hBeta *
          strFct2 * pow(u2 , (1. - hBeta) ) / hBeta *
          strFct3 * pow(u3 , (1. - hBeta) ) / hBeta *
          strFct4 * pow(u4 , (1. - hBeta) ) / hBeta;
      } else {
        strFct = 1.;
      }

      G4double dPhaseSpace = 1. * fabs(cos(fthetaComMax) - cos(fthetaComMin));
      G4double zLum = msZ[0] * ironDensity * (cm3/g) * (zpos + fTargLen/2) / cm * Avogadro / msA[0];
      G4double weight = 1. * zLum * dPhaseSpace * sigma * strFct;
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
      G4double eff_sigma = sigma;
      G4double Azz = -((7+pow(cos(thcom),2))*pow(sin(thcom),2))/pow(3+cos(thcom)*cos(thcom),2);

      fDefaultEvent->SetEffCrossSection(eff_sigma);
      fDefaultEvent->SetAsymmetry(Azz);
      fDefaultEvent->SetThCoM(thcom);
      fDefaultEvent->SetPhCoM(phcom);

      //GENERATE PARTICLES
      //electron #1
      G4double vtx_z = zpos + fZ;
      G4double KE1 = electron_mass_c2 * sqrt( 1 + pow(p1/electron_mass_c2,2)) - electron_mass_c2;
      particleGun->SetParticleEnergy( KE1 );
      particleGun->SetParticlePosition( G4ThreeVector(xpos, ypos, vtx_z) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( tX1, tY1, tZ1 ).unit() );
      fDefaultEvent->ProduceNewParticle(G4ThreeVector(xpos, ypos, vtx_z),
                                        G4ThreeVector(tX1, tY1, tZ1 ) * p1,
                                        particleGun->GetParticleDefinition()->GetParticleName() );
      particleGun->GeneratePrimaryVertex(anEvent);
      //electron#2
      G4double KE2 = electron_mass_c2 * sqrt( 1 + pow(p2/electron_mass_c2,2)) - electron_mass_c2;
      particleGun->SetParticleEnergy( KE2 );
      particleGun->SetParticlePosition( G4ThreeVector(xpos, ypos, vtx_z) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( tX2, tY2, tZ2 ).unit() );
      fDefaultEvent->ProduceNewParticle(G4ThreeVector(xpos, ypos, vtx_z),
                                        G4ThreeVector(tX2, tY2, tZ2 ) * p2,
                                        particleGun->GetParticleDefinition()->GetParticleName() );
      particleGun->GeneratePrimaryVertex(anEvent);

      //FINALIZE EVENT DATA FOR RECORDING
      fIO->SetEventData(fDefaultEvent);
    }
      else if( gentype == "LUND")
    {
      if( fNLUNDLines < 1 )
  	  {
  	    G4bool is_good_input = CheckLUNDFile("../g3input.dat");
        _unused( is_good_input ); //Avoids unused variable warning during compiling, var is used in assert()
  	    assert(is_good_input and "ERROR: bad LUND input file");

  	    const G4Run* aRun = G4RunManager::GetRunManager()->GetCurrentRun();
  	    G4int NTOTAL_TO_PROCESS = aRun->GetNumberOfEventToBeProcessed();
        _unused( NTOTAL_TO_PROCESS ); //Avoids unused variable warning during compiling, var is used in assert()
  	    assert(fNLUNDLines >= NTOTAL_TO_PROCESS and "Requested number of events exceeds LUND input limit");
  	  }

      // get kinematic variables
      G4double beamE;
      G4double xpos, ypos, zpos;
      G4double thcom, phcom;
      G4double p1[3];
      G4double p2[3];
      G4double Azz;

      LUNDfile >> beamE >> thcom >> phcom >> xpos >> ypos >> zpos >> p1[0] >> p1[1] >> p1[2] >> p2[0] >> p2[1] >> p2[2] >> Azz;

      // Attach units
      beamE = beamE * GeV;
      xpos *= cm;
      ypos *= cm;
      zpos *= cm;
      thcom *= radian;
      phcom *= radian;
      for(int i=0; i<3; i++)
  	  {
  	    p1[i] = p1[i] * GeV;
  	    p2[i] = p2[i] * GeV;
  	  }

      G4double me = electron_mass_c2;
      G4double eff_sigma = 1;

      G4double kinE1 = sqrt(p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] + me*me) - me;
      G4double kinE2 = sqrt(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2] + me*me) - me;

      fDefaultEvent->SetEffCrossSection(eff_sigma);
      fDefaultEvent->SetAsymmetry(Azz);
      fDefaultEvent->SetThCoM(thcom);
      fDefaultEvent->SetPhCoM(phcom);

      particleGun->SetParticleEnergy( kinE1 );
      particleGun->SetParticlePosition( G4ThreeVector(xpos, ypos, zpos) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( p1[0], p1[1], p1[2] ).unit() );
      fDefaultEvent->ProduceNewParticle(G4ThreeVector(xpos, ypos, zpos),
                                        G4ThreeVector(p1[0], p1[1], p1[2]),
                                        particleGun->GetParticleDefinition()->GetParticleName() );
      particleGun->GeneratePrimaryVertex(anEvent);

      particleGun->SetParticleEnergy( kinE2 );
      particleGun->SetParticlePosition( G4ThreeVector(xpos, ypos, zpos) );
      particleGun->SetParticleMomentumDirection( G4ThreeVector( p2[0], p2[1], p2[2] ).unit() );
      fDefaultEvent->ProduceNewParticle(G4ThreeVector(xpos, ypos, zpos),
                                        G4ThreeVector(p2[0], p2[1], p2[2]),
                                        particleGun->GetParticleDefinition()->GetParticleName() );
      particleGun->GeneratePrimaryVertex(anEvent);

      //FINALIZE EVENT DATA FOR RECORDING
      fIO->SetEventData(fDefaultEvent);
      }
        else if(gentype == "beam") // JUST A SINGLE BEAM
      {
      G4double xpos     = G4RandGauss::shoot( fX, fXsmear );
      G4double ypos     = G4RandGauss::shoot( fY, fYsmear );
      G4double zpos     = fZ;

      G4ThreeVector direction = G4ThreeVector(0,0,1.);

      direction.rotateX(-fBeamRotZY);
      direction.rotateY( fBeamRotZX);

      G4double me = electron_mass_c2;
      G4double beamE = fBeamE;
      G4double p = sqrt( beamE*beamE - me*me); //unit:MeV

      G4double pX = direction.x() * p;
      G4double pY = direction.y() * p;
      G4double pZ = direction.z() * p;

      fDefaultEvent->SetEffCrossSection(0);
      fDefaultEvent->SetAsymmetry(0);

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
  else
    {

      G4double xpos     = G4RandFlat::shoot( fXmin, fXmax );
      G4double ypos     = G4RandFlat::shoot( fYmin, fYmax );
      G4double zpos     = fZ;

      G4double thetaPos = G4RandFlat::shoot( fthetaMin, fthetaMax );
      G4double phiPos   = G4RandFlat::shoot( fphiMin, fphiMax );

      G4double me = electron_mass_c2;
      G4double beamE = fBeamE;
      G4double p = sqrt( beamE*beamE - me*me); //unit:MeV

      G4double pX = cos( phiPos ) * sin( thetaPos ) * p;
      G4double pY = sin( phiPos ) * sin( thetaPos ) * p;
      G4double pZ = cos( thetaPos ) * p;
      G4double eff_sigma = 1;
      G4double Azz = 1;

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
    fLEcorFac = 1 + tgtMom * aaa / (electron_mass_c2*0.001); //tgtMom returned in GeV ... need Me in GeV
    fLEtgtPol = 1;
  }else{
    G4double aaa = G4RandFlat::shoot(-1.,1.);
    G4double tgtMom = SampleTargetMomentum(0);
    fDefaultEvent->SetTargetMomentum( tgtMom );
    //G4double pUnpol = 1e6 * tgtMom * aaa;//Used for printout in fortran
    fLEcorFac = 1 + tgtMom * aaa / (electron_mass_c2*0.001); //tgtMom returned in GeV ... need Me in Gev
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

G4bool MolPolPrimaryGeneratorAction::CheckLUNDFile(G4String LUNDfile_name){

  // Get input kinematics variables
  // format
  // : BEAME, THETACOM, PHICOM, vtx[3], p1[3], p2[3], ANPOW

  G4bool is_good_lund = true;

  LUNDfile.open(LUNDfile_name);

  if( ! LUNDfile.good() )
    return false;

  G4int nlines = 0;
  G4String line;
  while(std::getline(LUNDfile, line))
    {
      nlines++;
    }

  LUNDfile.clear();
  LUNDfile.seekg(0, std::ios::beg);

  if( nlines < 1 )
    return false;

  fNLUNDLines = nlines;
  G4cout << "Total LUND generated: " << fNLUNDLines << G4endl;

  return is_good_lund;
}
