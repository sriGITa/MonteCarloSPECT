//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "ExN01DetectorConstruction.hh"
#include "ExN01PhysicsList.hh"
#include "ExN01PrimaryGeneratorAction.hh"
#include "ExN01DicomRegularDetectorConstruction.hh"
#include "ExN01DicomIntersectVolume.hh"

#include "ExN01TrackerHit.hh"
#include "ExN01TrackerSD.hh"
#include "ExN01EventAction.hh"
#include "ExN01RunAction.hh"
#include "ExN01DicomHandler.hh"
#include "ExN01TrackingAction.hh"
#include "ExN01SteppingAction.hh"
#include "ExN01SteppingVerbose.hh"

#include "Randomize.hh"
#include "time.h"
#include "G4Timer.hh"
#include "G4ScoringManager.hh"

#include "G4UImanager.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
int main()
//int main(int argc,char** argv)
{
  // Calculate the total time for this simulation
  G4Timer totaltime;
  totaltime.Start();

  /*
  // Choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  G4long seed=time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);
   */
  // Run manager
  G4RunManager* runManager = new G4RunManager;
  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();
  scoringManager->SetVerboseLevel(1);

  // G4VUserDetectorConstruction* detector = new ExN01DetectorConstruction();
  ExN01DicomHandler* dcmHandler = 0;

  // Treatment of DICOM images before creating the G4runManager
  dcmHandler = new ExN01DicomHandler;
  dcmHandler->CheckFileFormat();

  ExN01DetectorConstruction* detector=0;
  detector=new ExN01DicomRegularDetectorConstruction;
  runManager->SetUserInitialization(detector);

  G4VUserPhysicsList* physics = new ExN01PhysicsList;
  runManager->SetUserInitialization(physics);

  // set mandatory user action class
  ExN01PrimaryGeneratorAction* gen_action = new ExN01PrimaryGeneratorAction();

  ExN01RunAction *runaction = new ExN01RunAction(gen_action);
  runManager->SetUserAction(runaction);
  runManager->SetUserAction(gen_action);

  ExN01EventAction *eventaction = new ExN01EventAction;
  runManager->SetUserAction(eventaction);

  ExN01TrackingAction*   track = new  ExN01TrackingAction(runaction,eventaction);
  runManager->SetUserAction(track);

  ExN01SteppingAction *steppingaction = new ExN01SteppingAction;
  runManager->SetUserAction(steppingaction);

  // Initialize G4 kernel
  runManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/run/verbose 0");
  UImanager->ApplyCommand("/event/verbose 0");
  UImanager->ApplyCommand("/tracking/verbose 0");
  UImanager->ApplyCommand("/control/execute vis.mac");

  /*
  //Start a run
  G4int numberOfEvent = 1000000;
  runManager->BeamOn(numberOfEvent);
  */
  new ExN01DicomIntersectVolume();

  /*
  #ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
    #ifdef G4UI_USE
            G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        #ifdef G4VIS_USE
            UImanager->ApplyCommand("/control/execute vis.mac");
   //#else
   // UImanager->ApplyCommand("/control/execute init.mac");
        #endif
        if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");

        ui->SessionStart();
        delete ui;
    #endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
*/
  delete dcmHandler;
  delete runManager;

  //output the total time
  totaltime.Stop();
  std::cout<<"**************************************************"<<std::endl;
  std::cout<<"The total time for this simulation is: "<<totaltime<<std::endl;
  std::cout<<"**************************************************"<<std::endl;
  return 0;
 }
