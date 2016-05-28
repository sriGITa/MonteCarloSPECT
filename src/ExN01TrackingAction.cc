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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN01TrackingAction.hh"
#include "ExN01RunAction.hh"

#include "ExN01RunAction.hh"
#include "ExN01EventAction.hh"
#include "ExN01TrackingMessenger.hh"
#include "ExN01HistoManager.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01TrackingAction::ExN01TrackingAction(ExN01RunAction* RA,ExN01EventAction* EA)
:fRun(RA),fEvent(EA)
{
  fullChain = false;
  fTrackMessenger = new ExN01TrackingMessenger(this);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01TrackingAction::~ExN01TrackingAction()
{
  delete fTrackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  G4ParticleDefinition* particle = track->GetDefinition();
  G4String name   = particle->GetParticleName();
  fCharge = particle->GetPDGCharge();
  fMass   = particle->GetPDGMass();  

  G4double Ekin = track->GetKineticEnergy();
  G4int ID      = track->GetTrackID();
  
  G4bool condition = false;  
  //count particles
  //
  fRun->ParticleCount(name, Ekin);
  
  //energy spectrum
  //
  G4int ih = 0;
  if (particle == G4Electron::Electron()||
      particle == G4Positron::Positron())  ih = 1;
  else if (particle == G4NeutrinoE::NeutrinoE()||
           particle == G4AntiNeutrinoE::AntiNeutrinoE()) ih = 2;
  else if (particle == G4Gamma::Gamma()) ih = 3;
  else if (particle == G4Alpha::Alpha()) ih = 4;
  else if (fCharge > 2.) ih = 5;
  if (ih) G4AnalysisManager::Instance()->FillH1(ih, Ekin);
  
  //fullChain: stop ion and print decay chain
  //fCharge determines whether the output is gamma, electron or other particles
  //if (fCharge > 2.) {//uncomment to output the result Co57[0.0] ---> nu_e ---> Fe57[136.5] ---> Fe57[0.0] ---> gamma ---> gamma
  // G4Track* tr = (G4Track*) track;
  // if (fullChain) tr->SetTrackStatus(fStopButAlive);
  if (ID == 1) fEvent->AddDecayChain(name);
    else       fEvent->AddDecayChain(" ---> " + name);
  //
  
  //example of saving random number seed of this fEvent, under condition
  //
  ////condition = (ih == 3);

  if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  //keep only ions
  //

  if (fCharge < 3. ) return;

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  
  //get time
  //   
  G4double time = track->GetGlobalTime();
  G4int ID = track->GetTrackID();
  if (ID == 1) fRun->PrimaryTiming(time);        //time of life of primary ion  
      
  //energy and momentum balance (from secondaries)
  //
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  size_t nbtrk = (*secondaries).size();

  if (nbtrk) {
    //there are secondaries --> it is a decay
    //
    //force 'single' decay
    //1
   //if ((!fullChain)&&(ID > 1)) G4RunManager::GetRunManager()->AbortEvent();
    //

    //balance    
  G4double EkinTot = 0.;
  G4ThreeVector Pbalance = - track->GetMomentum();
  for (size_t itr=0; itr<nbtrk; itr++) {
     G4Track* trk = (*secondaries)[itr];
     EkinTot += trk->GetKineticEnergy();
     //exclude gamma deexcitation from momentum balance
     if (trk->GetDefinition() != G4Gamma::Gamma())
       Pbalance += trk->GetMomentum();
     }
    G4double Pbal = Pbalance.mag();  
    fRun->Balance(EkinTot,Pbal);  
    analysis->FillH1(6,EkinTot);
    analysis->FillH1(7,Pbal);
    
  }

  //no secondaries --> end of chain    
  //  
 if (nbtrk!=0) {
    fRun->EventTiming(time);                     //total time of life   
    analysis->FillH1(8,time);
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

