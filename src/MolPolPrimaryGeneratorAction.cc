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


  return;

}

void MolPolPrimaryGeneratorAction::SourceModeSet(G4int i) {
//      SourceModeSet(0); // point to the one below with default settings = 0. // should I just use default parameters?
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
