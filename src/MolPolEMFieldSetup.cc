// *************************************************************** (╯°□°）╯︵ ┻━┻
//
//	MolPolEMFieldSetup.cc
//
//  Sets up the integrated global field. Initalizes all 6 fields with zero strength
//  and boolean for each field usage set to false.
//
//  Takes values from MolPolFieldMessenger in UpdateConfiguration() method.  This
//  will turn on (or turn off) and update fields as specified.
//
//  Added G4cout lines for bug checking and testing.
//
//  Note: to turn off a field it is (believed to be) fully sufficient to simply
//        use a macro line for idealized field set to zero strength.
//
//
//	Eric King - 2018-11-19
//
// *****************************************************************************

#include "MolPolEMFieldSetup.hh"
#include "MolPolEMFieldMessenger.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ios.hh"

#include "MolPolQuad.hh"
#include "MolPolDipole.hh"
#include "MolPolSolenoid.hh"
#include "MolPolTOSCAField.hh"

#include <iostream>
using namespace std;


MolPolEMFieldSetup::MolPolEMFieldSetup()
  : fFieldManager(0),
    fChordFinder(0),
    fStepper(0),
    fIntgrDriver(0),
    fFieldMessenger(0),
    fStepperType(0),
    fMinStep(0)
{
  InitialseAll();
}

void MolPolEMFieldSetup::InitialseAll()
{
  // Initializes field messenger.
  fFieldMessenger = new MolPolEMFieldMessenger(this);

  // Constructs the global field object and assigns necessary parameters.
  fEMfield = new MolPolEMField();
  fEquation = new G4EqMagElectricField(fEMfield);
  fMinStep  = 0.01*mm ; // minimal step of 1 miron, default is 0.01 mm: Doesn't seem to make much difference here
  fStepperType = 4 ;    // ClassicalRK4 -- the default stepper
  fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fChordFinder = 0;
  UpdateField(); // This is the updater for the GLOBAL FIELD.

  // Creates zero strength solenoid object set to status off.  Will be updated with update field.
  G4cout << "Intializing Solenoid EMField Setup:" << G4endl;
  fEMfield->setSolenoidStatus(0); //Solenoid Off
  MolPolSolenoid* sol = new MolPolSolenoid( 0.0 );
  fEMfield->setSolenoidObject( sol ); // Assigns object. At zero strength. Will update after information gained from field messenger to proper strength.

  // Creates ideal quad objects for fEMfield with status set to off.  Updates with update field command to whatever is wanted.
  G4cout << "Intializing Quads in EMField:" << G4endl;
  G4cout << "  Defaulting Quad1 to ideal..." << G4endl;
  fEMfield->setQuad1Status(0);
  MolPolQuad* q1 = new MolPolQuad( 0.0 , G4ThreeVector(0,0,ORIGINQ1) , NOROT , BORERADIUS , 2*18.29*cm );
  fEMfield->setQuad1Object(1,q1);
  G4cout << "  Defaulting Quad2 to ideal..." << G4endl;
  fEMfield->setQuad2Status(0);
  MolPolQuad* q2 = new MolPolQuad( 0.0 , G4ThreeVector(0,0,ORIGINQ2) , NOROT , BORERADIUS , 2*22.30*cm );
  fEMfield->setQuad2Object(1,q2);
  G4cout << "  Defaulting Quad3 to ideal..." << G4endl;
  fEMfield->setQuad3Status(0);
  MolPolQuad* q3 = new MolPolQuad( 0.0 , G4ThreeVector(0,0,ORIGINQ3) , NOROT , BORERADIUS , 2*18.37*cm );
  fEMfield->setQuad3Object(1,q3);
  G4cout << "  Defaulting Quad4 to ideal..." << G4endl;
  fEMfield->setQuad4Status(0);
  MolPolQuad* q4 = new MolPolQuad( 0.0 , G4ThreeVector(0,0,ORIGINQ4) , NOROT , BORERADIUS , 2*18.37*cm );
  fEMfield->setQuad4Object(1,q4);

  // Creates ideal dipole objects for fEMfield with status set to off.  Updates with update field command to desired choice.
  G4cout << "Intializing Dipole in EMField:" << G4endl;
  G4cout << "  Dipole defaults to ideal for setup..." << G4endl;
  fEMfield->setDipoleStatus(0);
  MolPolDipole* dip = new MolPolDipole(0.0 , 82.25*cm); //(Tesla strength , Zeff) both from Geometry setup.
  fEMfield->setDipoleObject(1,dip); // Initializes as type ideal. Since no values are yet passed from the field macro this is fine. Will Update()
}

/////////////////////////////////////////////////////////////////////////////////
//

MolPolEMFieldSetup::~MolPolEMFieldSetup()
{
  if(fChordFinder)    delete fChordFinder;
  if(fStepper)        delete fStepper;
  if(fEquation)       delete fEquation;
  if(fEMfield)        delete fEMfield;
  if(fFieldMessenger) delete fFieldMessenger;
}

void MolPolEMFieldSetup::UpdateConfiguration(){
  // Tesla are the default unit. Currents, which are given to the same variable,
  // will be passed as unitless (when divided out by teslas). Additionally, it's
  // perfectly fine to assign the tesla unit to TOSCA map strengths and
  dDipRelevantStr   *= tesla;
  dDipToscaMapStr   *= tesla;
  dQuad4RelevantStr *= tesla;
  dQuad4ToscaMapStr *= tesla;
  dQuad3RelevantStr *= tesla;
  dQuad3ToscaMapStr *= tesla;
  dQuad2RelevantStr *= tesla;
  dQuad2ToscaMapStr *= tesla;
  dQuad1RelevantStr *= tesla;
  dQuad1ToscaMapStr *= tesla;
  dSolRelevantStr   *= tesla;

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Field information recieved from macro.  Work on this.
  G4cout << "\n>>>>>>>>>> Received values from macro <<<<<<<<<<" << G4endl;
  G4cout << "   SCS Strength: " << dSolRelevantStr / tesla << " tesla" << G4endl;
  G4cout << "------------------------------------------------" << G4endl;
  if(iQuad1Type  == 0){
    G4cout << "   Quad1 Type: " << iQuad1Type << G4endl;
    G4cout << "   Quad1 Amps: " << dQuad1RelevantStr / tesla << " amps" << G4endl;
  } else if(iQuad1Type  == 1){
    G4cout << "   Quad1 Type: " << iQuad1Type << G4endl;
    G4cout << "  Quad1 Field: " << dQuad1RelevantStr / tesla << " tesla" << G4endl;
  } else if(iQuad1Type  == 2){
    G4cout << "   Quad1 Type: " << iQuad1Type << G4endl;
    G4cout << "    Quad1 Map: " << strQuad1MapLoc << G4endl;
    G4cout << "    Map Scale: " << dQuad1ToscaMapStr / tesla << G4endl;
    G4cout << "Desired Scale: " << dQuad1RelevantStr / tesla << G4endl;
  }
  G4cout << "------------------------------------------------" << G4endl;
  if(iQuad2Type  == 0){
    G4cout << "   Quad2 Type: " << iQuad2Type << G4endl;
    G4cout << "   Quad2 Amps: " << dQuad2RelevantStr / tesla << " amps" << G4endl;
  } else if(iQuad2Type  == 1){
    G4cout << "   Quad2 Type: " << iQuad2Type << G4endl;
    G4cout << "  Quad2 Field: " << dQuad2RelevantStr / tesla << " tesla" << G4endl;
  } else if(iQuad2Type  == 2){
    G4cout << "   Quad2 Type: " << iQuad2Type << G4endl;
    G4cout << "    Quad2 Map: " << strQuad2MapLoc << G4endl;
    G4cout << "    Map Scale: " << dQuad2ToscaMapStr / tesla << G4endl;
    G4cout << "    Des Scale: " << dQuad2RelevantStr / tesla << G4endl;
  }
  G4cout << "------------------------------------------------" << G4endl;
  if(iQuad3Type  == 0){
    G4cout << "   Quad3 Type: " << iQuad3Type << G4endl;
    G4cout << "   Quad3 Amps: " << dQuad3RelevantStr / tesla << " amps" << G4endl;
  } else if(iQuad3Type  == 1){
    G4cout << "   Quad3 Type: " << iQuad3Type << G4endl;
    G4cout << "  Quad3 Field: " << dQuad3RelevantStr / tesla << " tesla" << G4endl;
  } else if(iQuad3Type  == 2){
    G4cout << "   Quad3 Type: " << iQuad3Type << G4endl;
    G4cout << "    Quad3 Map: " << strQuad3MapLoc << G4endl;
    G4cout << "    Map Scale: " << dQuad3ToscaMapStr / tesla << G4endl;
    G4cout << "Desired Scale: " << dQuad3RelevantStr / tesla << G4endl;
  }
  G4cout << "------------------------------------------------" << G4endl;
  if(iQuad4Type  == 0){
    G4cout << "   Quad4 Type: " << iQuad4Type << G4endl;
    G4cout << "   Quad4 Amps: " << dQuad4RelevantStr / tesla << " amps" << G4endl;
  } else if(iQuad4Type  == 1){
    G4cout << "   Quad4 Type: " << iQuad4Type << G4endl;
    G4cout << "  Quad4 Field: " << dQuad4RelevantStr / tesla << " tesla" << G4endl;
  } else if(iQuad4Type  == 2){
    G4cout << "   Quad4 Type: " << iQuad4Type << G4endl;
    G4cout << "    Quad4 Map: " << strQuad4MapLoc << G4endl;
    G4cout << "    Map Scale: " << dQuad4ToscaMapStr / tesla << G4endl;
    G4cout << "Desired Scale: " << dQuad4RelevantStr / tesla << G4endl;
  }
  G4cout << "------------------------------------------------" << G4endl;
  if(iDipoleType == 0){
    G4cout << "  Dipole Type: " << iDipoleType << G4endl;
    G4cout << "  Dipole Amps: " << dDipRelevantStr / tesla << " amps" << G4endl;
  } else if(iDipoleType == 1){
    G4cout << "  Dipole Type: " << iDipoleType << G4endl;
    G4cout << " Dipole Field: " << dDipRelevantStr / tesla << " tesla" << G4endl;
  } else if(iDipoleType == 2){
    G4cout << "  Dipole Type: " << iDipoleType << G4endl;
    G4cout << "   Dipole Map: " << strDipoleMapLoc << G4endl;
    G4cout << "    Map Scale: " << dDipToscaMapStr / tesla << G4endl;
    G4cout << "Desired Scale: " << dDipRelevantStr / tesla << G4endl;
  }

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Solenoid ... Only need the desired strength for solenoid field. 0.0 to turn off.
  if (dSolRelevantStr == 0.){
    G4cout << "\nUpdating Solenoid Configuration Status: OFF" << G4endl;
    fEMfield->setSolenoidStatus(0); //Solenoid Off
  } else {
    G4cout << "\nUpdating Solenoid Configuration..." << G4endl
           << "  Status: ON;" << G4endl
           << "  Desired Strength: " << dSolRelevantStr / tesla << " tesla" << G4endl;
    fEMfield->setSolenoidStatus(1); //Solenoid On
    fEMfield->setSolenoidStrength(dSolRelevantStr / tesla);
  }

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Quad Fields Updated ... What would be desirable here is if the specs of the
  // Quad1
  if(dQuad1RelevantStr != 0.){
    if(iQuad1Type < 2){
      G4double gradient = 0.0;
      if(iQuad1Type == 0){
        gradient = CalA2T(dQuad1RelevantStr/tesla,1) / BORERADIUS;
        G4cout << "\nUpdating Quad1 Configuration... " << G4endl
               << "    Status: ON; " << G4endl
               << "      Type: Ideal by Current;" << G4endl
               << "   Current: " << dQuad1RelevantStr / tesla << "tesla;" << G4endl // this is amps
               << "  Gradient: " << gradient / (tesla / m) << " tesla/meter" << G4endl;
      }
      if(iQuad1Type == 1){
        gradient = dQuad1RelevantStr / BORERADIUS;
        G4cout << "\nUpdating Quad1 Configuration... " << G4endl
               << "    Status: ON; " << G4endl
               << "      Type: Ideal by Field;" << G4endl
               << "     Field: " << dQuad1RelevantStr / tesla << "tesla;" << G4endl
               << "  Gradient: " << gradient / (tesla / m) << " tesla/meter" << G4endl;
      }
      fEMfield->setQuad1Status(1); //Quad1 On
      MolPolQuad* q1 = new MolPolQuad( gradient , G4ThreeVector(0,0,ORIGINQ1) , NOROT , BORERADIUS , 2*18.29*cm );
      fEMfield->setQuad1Object(1,q1);
    }
    if(iQuad1Type == 2){
      fEMfield->setQuad1Status(1); //Quad1 On
      G4double scale = dQuad1RelevantStr / dQuad1ToscaMapStr;
      G4cout << "Got here.." << G4endl;
      G4cout << "strQuad1MapLoc: " << strQuad1MapLoc << G4endl;
      MolPolTOSCAField* q1 = new MolPolTOSCAField( strQuad1MapLoc, scale, ORIGINQ1 );
      G4cout << "But not here..." << G4endl;
      fEMfield->setQuad1Object(2,q1);
      G4cout << "\nUpdating Quad1 Configuration... " << G4endl
             << "       Status: ON; " << G4endl
             << "         Type: TOSCA Map;" << G4endl
             << "    Map Scale: " << dQuad1ToscaMapStr << ";" << G4endl
             << "Desired Scale: " << dQuad1RelevantStr << G4endl;
    }
  } else {
    G4cout << "\nUpdating Quad1 Configuration Status: OFF" << G4endl;
    fEMfield->setQuad1Status(0);
  }
  //Quad2
  if(dQuad2RelevantStr != 0.){
    if(iQuad2Type < 2){
      G4double gradient = 0.0;
      if(iQuad2Type == 0){
        gradient = CalA2T(dQuad2RelevantStr/tesla,2) / BORERADIUS;
        G4cout << "\nUpdating Quad2 Configuration... " << G4endl
               << "    Status: ON; " << G4endl
               << "      Type: Ideal by Current;" << G4endl
               << "   Current: " << dQuad2RelevantStr / tesla << " amps;" << G4endl
               << "  Gradient: " << gradient / (tesla / m) << " tesla/meter" << G4endl;
      }
      if(iQuad2Type == 1){
        gradient = dQuad2RelevantStr / BORERADIUS;
        G4cout << "\nUpdating Quad2 Configuration... " << G4endl
               << "    Status: ON; " << G4endl
               << "      Type: Ideal by Field;" << G4endl
               << "     Field: " << dQuad2RelevantStr / tesla << " tesla;" << G4endl
               << "  Gradient: " << gradient / (tesla / m) << " tesla/meter" << G4endl;
      }
      fEMfield->setQuad2Status(1); //Quad2 On
      MolPolQuad* q1 = new MolPolQuad( gradient , G4ThreeVector(0,0,ORIGINQ2) , NOROT , BORERADIUS , 2*22.30*cm );
      fEMfield->setQuad2Object(1,q1);
    }
    if(iQuad2Type == 2){
      fEMfield->setQuad2Status(1); //Quad2 On
      G4double scale = dQuad2RelevantStr / dQuad2ToscaMapStr;
      MolPolTOSCAField* q1 = new MolPolTOSCAField( strQuad2MapLoc, scale, ORIGINQ2 );
      fEMfield->setQuad2Object(2,q1);
      G4cout << "\nUpdating Quad2 Configuration... " << G4endl
             << "       Status: ON; " << G4endl
             << "         Type: TOSCA Map;" << G4endl
             << "    Map Scale: " << dQuad2ToscaMapStr << ";" << G4endl
             << "Desired Scale: " << dQuad2RelevantStr << G4endl;
    }
  } else {
    G4cout << "\nUpdating Quad2 Configuration Status: OFF" << G4endl;
    fEMfield->setQuad2Status(0);
  }
  //Quad3
  if(dQuad3RelevantStr != 0.){
    if(iQuad3Type < 2){
      G4double gradient = 0.0;
      if(iQuad3Type == 0){
        gradient = CalA2T(dQuad3RelevantStr/tesla,3) / BORERADIUS;
        G4cout << "\nUpdating Quad3 Configuration... " << G4endl
               << "    Status: ON; " << G4endl
               << "      Type: Ideal by Current;" << G4endl
               << "   Current: " << dQuad3RelevantStr / tesla << " amps;" << G4endl
               << "  Gradient: " << gradient / (tesla / m) << " tesla/meter" << G4endl;
      }
      if(iQuad3Type == 1){
        gradient = dQuad3RelevantStr / BORERADIUS;
        G4cout << "\nUpdating Quad3 Configuration... " << G4endl
               << "    Status: ON; " << G4endl
               << "      Type: Ideal by Field;" << G4endl
               << "     Field: " << dQuad3RelevantStr / tesla << " tesla;" << G4endl
               << "  Gradient: " << gradient / (tesla / m) << " tesla/meter" << G4endl;
      }
      fEMfield->setQuad3Status(1); //Quad3 On
      MolPolQuad* q1 = new MolPolQuad( gradient , G4ThreeVector(0,0,ORIGINQ3) , NOROT , BORERADIUS , 2*18.37*cm );
      fEMfield->setQuad3Object(1,q1);
    }
    if(iQuad3Type == 2){
      fEMfield->setQuad3Status(1); //Quad3 On
      G4double scale = dQuad3RelevantStr / dQuad3ToscaMapStr;
      MolPolTOSCAField* q1 = new MolPolTOSCAField( strQuad3MapLoc, scale, ORIGINQ3 );
      fEMfield->setQuad3Object(2,q1);
      G4cout << "\nUpdating Quad3 Configuration... " << G4endl
             << "       Status: ON; " << G4endl
             << "         Type: TOSCA Map;" << G4endl
             << "    Map Scale: " << dQuad3ToscaMapStr << ";" << G4endl
             << "Desired Scale: " << dQuad3RelevantStr << G4endl;
    }
  } else {
    G4cout << "\nUpdating Quad3 Configuration Status: OFF" << G4endl;
    fEMfield->setQuad3Status(0);
  }
  //Quad4
  if(dQuad4RelevantStr != 0.){
    if(iQuad4Type < 2){
      G4double gradient = 0.0;
      if(iQuad4Type == 0){
        gradient = CalA2T(dQuad4RelevantStr/tesla,4) / BORERADIUS;
        G4cout << "\nUpdating Quad4 Configuration... " << G4endl
               << "    Status: ON; " << G4endl
               << "      Type: Ideal by Current;" << G4endl
               << "   Current: " << dQuad4RelevantStr / tesla << " amps;" << G4endl
               << "  Gradient: " << gradient / (tesla / m) << " tesla/meter" << G4endl;
      }
      if(iQuad4Type == 1){
        gradient = dQuad4RelevantStr / BORERADIUS;
        G4cout << "\nUpdating Quad4 Configuration... " << G4endl
               << "    Status: ON; " << G4endl
               << "      Type: Ideal by Field;" << G4endl
               << "     Field: " << dQuad4RelevantStr / tesla << " tesla;" << G4endl
               << "  Gradient: " << gradient / (tesla / m) << " tesla/meter" << G4endl;
      }
      fEMfield->setQuad4Status(1); //Quad4 On
      MolPolQuad* q1 = new MolPolQuad( gradient , G4ThreeVector(0,0,ORIGINQ4) , NOROT , BORERADIUS , 2*18.37*cm );
      fEMfield->setQuad4Object(1,q1);
    }
    if(iQuad4Type == 2){
      fEMfield->setQuad4Status(1); //Quad4 On
      G4double scale = dQuad4RelevantStr / dQuad4ToscaMapStr;
      MolPolTOSCAField* q1 = new MolPolTOSCAField( strQuad4MapLoc, scale, ORIGINQ4 );
      fEMfield->setQuad4Object(2,q1);
      G4cout << "\nUpdating Quad4 Configuration... " << G4endl
             << "       Status: ON; " << G4endl
             << "         Type: TOSCA Map;" << G4endl
             << "    Map Scale: " << dQuad4ToscaMapStr << ";" << G4endl
             << "Desired Scale: " << dQuad4RelevantStr << G4endl;
    }
  } else {
    G4cout << "\nUpdating Quad4 Configuration Status: OFF" << G4endl;
    fEMfield->setQuad4Status(0);
  }

  ///////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Dipole Field Updated ... Ideal if we don't have to reload map for TOSCA dipole.
  // Additional if statements will be needed at some point.
  // i.e. if was TOSCA and still is... just update.
  // i.e. if was TOSCA and now ideal...  create new object and set EMfield object to it.
  if(dDipRelevantStr != 0.){
    G4cout << "\nUpdating Dipole Configuration..." << G4endl;
    fEMfield->setDipoleStatus(1); //Dipole on
    if(iDipoleType == 1){
      MolPolDipole* dip = new MolPolDipole(dDipRelevantStr,G4ThreeVector(0,-9.0*cm,ORIGIND),NOROT,2*82.25*cm);
      fEMfield->setDipoleObject(1,(G4MagneticField*)dip );
      G4cout << "         Status: ON;" << G4endl
             << "           Type: Ideal;" << G4endl
             << "  Bend Strength: " << dDipRelevantStr/tesla << " tesla" << G4endl;
    }
    if(iDipoleType == 2){
      G4double scale = dDipRelevantStr / dDipToscaMapStr;
      MolPolTOSCAField* dip = new MolPolTOSCAField( strDipoleMapLoc, scale, ORIGIND );
      fEMfield->setDipoleObject(2,(G4MagneticField*)dip );
    }
  } else { //if for some odd reason the dipole is desired to be turned off
    fEMfield->setDipoleStatus(0); //Dipole off
  }

}

/////////////////////////////////////////////////////////////////////////////
// Register this field to 'global' Field Manager and
// Create Stepper and Chord Finder with predefined type, minstep (resp.)
void MolPolEMFieldSetup::UpdateField()
{
  fStepper = new G4ClassicalRK4( fEquation, 8 );
  fFieldManager->SetDetectorField(fEMfield);
  if(fChordFinder) delete fChordFinder;
  fIntgrDriver = new G4MagInt_Driver(fMinStep,fStepper,fStepper->GetNumberOfVariables());
  fChordFinder = new G4ChordFinder(fIntgrDriver);
  fFieldManager->SetChordFinder( fChordFinder );
}


/////////////////////////////////////////////////////////////////////////////
// Set stepper according to the stepper type
void MolPolEMFieldSetup::SetStepper()
{
  G4int nvar = 8;
  if(fStepper) delete fStepper;
  switch ( fStepperType ){
    case 0:
      fStepper = new G4ExplicitEuler( fEquation, nvar );
      break;
    case 1:
      fStepper = new G4ImplicitEuler( fEquation, nvar );
      break;
    case 2:
      fStepper = new G4SimpleRunge( fEquation, nvar );
      break;
    case 3:
      fStepper = new G4SimpleHeum( fEquation, nvar );
      break;
    case 4:
      fStepper = new G4ClassicalRK4( fEquation, nvar );
      break;
    case 5:
      fStepper = new G4CashKarpRKF45( fEquation, nvar );
      break;
    default: fStepper = 0;
  }
}


///////////////////////////////////////////////////////////////////////////////
// Current to Field Calculation (paprameters from Sasha)
// Return field at the pole tip (Quadrupole)
G4double MolPolEMFieldSetup::CalA2T(G4double current, G4int magnet)
{

  G4double gl1 = 0;
  G4double gl2 = 0;
  G4double gln = 0;
  G4double fld = 0;

  G4double cn = current / 300.;
  if(magnet == 1){
    // Moller Quad Q1/MQO1H01/LARGE/new/white
    gl1 = (.0110605+5.33237*cn-.0142794*pow(cn,2)+.259313*pow(cn,3));
    gln = (gl1+0.0058174*pow(cn,4)-0.831887*pow(cn,5));
    fld = gln*10.*5.08/36.5723;
  } else if(magnet == 2){
    // Moller Quad Q2/PATSY/MQM1H02/SMALL/RED
    gl1=(0.0196438+5.35443*cn+0.0297273*pow(cn,2)+0.103505*pow(cn,3));
    gln=(gl1-0.0449275*pow(cn,4)-0.211868*pow(cn,5));
    fld=gln*10.*5.08/44.76;
  } else if(magnet == 3){
    // Moller Quad Q3/TESSA/MQO1H03/LARGE/BLUE
    gl1=(0.000632446+5.15178*cn-0.00262778*pow(cn,2));
    gl2=(-0.107635*pow(cn,3)+0.00209902*pow(cn,4));
    gln=(gl1+gl2-0.640635*pow(cn,5));
    fld=gln*10.*5.08/36.74 ;
  } else if(magnet == 4){
    // Moller Quad Q4/FELICIA/MQO1H03A/LARGE/BLUE
    gl1=(0.0001732+5.2119*cn-0.000732518*pow(cn,2));
    gl2=(-0.133423*pow(cn,3)+0.000618402*pow(cn,4));
    gln=(gl1+gl2-0.647082*pow(cn,5));
    fld=gln*10.*5.08/36.50;
  } else if(magnet == 5){
    // Moller Dipole LILLY/MMA1H01/Blue
    fld=(-0.39026E-04+0.027051*current-0.17799E-08*pow(current,2));
  } else{
    //wrong magnet setup
    fld = 0.0;
  }

  return fld * 0.1 * tesla;

}
