#include "Randomize.hh"

#include "MolPolRunAction.hh"
#include "MolPolPrimaryGeneratorAction.hh"
#include "MolPolEventAction.hh"
#include "MolPolSteppingAction.hh"
#include "MolPolDetectorConstruction.hh"
#include "MolPolIO.hh"
#include "MolPolMessenger.hh"

//Standard physics list
#include "G4PhysListFactory.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4RunManagerKernel.hh"

//GUI Things
#include <G4UImanager.hh>
#include <G4UIExecutive.hh>
#include <G4UIterminal.hh>
#include <G4VisExecutive.hh>

//Needed if GDML geomety ever implemented
//#include "G4GDMLParser.hh"

#include <sys/types.h>
#include <sys/stat.h>

int main(int argc, char** argv){

    // Initialize seed
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4int seconds =  time(NULL);
    G4Random::setTheSeed(seconds);

    MolPolIO *io = new MolPolIO();

    //-------------------------------
    // Initialization of Run manager
    //-------------------------------
    G4cout << "RunManager construction starting...." << G4endl;
    G4RunManager * runManager = new G4RunManager;

    MolPolMessenger *rmmess = new MolPolMessenger();
    rmmess->SetIO(io);

    // Detector geometry
    MolPolDetectorConstruction* detector = new MolPolDetectorConstruction();
    runManager->SetUserInitialization(detector);
    rmmess->SetDetCon( ((MolPolDetectorConstruction *) detector) );

    // Physics we want to use
    G4int verbose = 0;
    G4PhysListFactory factory;
    G4VModularPhysicsList* physlist = factory.GetReferencePhysList("FTFP_BERT");
    physlist->SetVerboseLevel(verbose);
    runManager->SetUserInitialization(physlist);
    //physlist->RegisterPhysics( new MolPolOpticalPhysics() );

    //-------------------------------
    // UserAction classes
    //-------------------------------
    G4UserRunAction* run_action = new MolPolRunAction;
    ((MolPolRunAction *) run_action)->SetIO(io);
    runManager->SetUserAction(run_action);

    G4VUserPrimaryGeneratorAction* gen_action = new MolPolPrimaryGeneratorAction();
    ((MolPolPrimaryGeneratorAction *) gen_action)->SetIO(io);
    rmmess->SetPriGen((MolPolPrimaryGeneratorAction *)gen_action);
    runManager->SetUserAction(gen_action);

    G4UserEventAction* event_action = new MolPolEventAction;
    ((MolPolEventAction *) event_action)->SetIO(io);

    runManager->SetUserAction(event_action);
    G4UserSteppingAction* stepping_action = new MolPolSteppingAction;
    runManager->SetUserAction(stepping_action);
    rmmess->SetStepAct((MolPolSteppingAction *) stepping_action);

    runManager->Initialize(); 

    /*
    // Export gdml file
    G4VPhysicalVolume* pWorld = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
    G4GDMLParser parser;
    parser.Write("g4test.gdml", pWorld);
    */

    //NEW VISUALIZATION CODE--ERIC 6/22/2021 ... ISSUES RUNNING OLD CODE ON 10.7
    // Visualization, if you choose to have it!
    //
    // Simple graded message scheme - give first letter or a digit:
    //  0) quiet,         // Nothing is printed.
    //  1) startup,       // Startup and endup messages are printed...
    //  2) errors,        // ...and errors...
    //  3) warnings,      // ...and warnings...
    //  4) confirmations, // ...and confirming messages...
    //  5) parameters,    // ...and parameters of scenes and views...
    //  6) all            // ...and everything available.

    G4UImanager *  uiManager  = G4UImanager::GetUIpointer();
    G4VisManager * visManager = new G4VisExecutive("quiet");
    visManager->Initialize();
    visManager->SetVerboseLevel("quiet");
	
    if(argc==1){
        //Visualization session :: qt, xm or win32 -- we could pass this as an option -- default to "qt" for now.
        G4UIExecutive * session = new G4UIExecutive(argc,argv,"qt");
        if( session->IsGUI() ) uiManager->ApplyCommand("/control/execute macros/gui.mac");
        session->SessionStart();
        delete session;
    } else {
	G4String command = "/control/execute ";
	G4String fileName = argv[1];
	uiManager->ApplyCommand(command+fileName);
    }

    //clean up pointers
    delete runManager;
    delete visManager;

    return 0;
}
