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
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN01RunAction_h
#define ExN01RunAction_h 1
#include <map>

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Timer.hh"
#include "G4ParticleDefinition.hh"

#include <iostream>
#include <fstream>
using namespace std;

extern int runorder;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class ExN01PrimaryGeneratorAction;
class ExN01HistoManager;

class ExN01RunAction : public G4UserRunAction
{
  public:
    ExN01RunAction(ExN01PrimaryGeneratorAction*);
   ~ExN01RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

    void ParticleCount(G4String, G4double);
    void Balance(G4double,G4double);
    void EventTiming(G4double);
    void PrimaryTiming(G4double);

  public:
    ExN01PrimaryGeneratorAction* fPrimary;

    std::map<G4String,G4int> fParticleCount;
    std::map<G4String,G4double> fEmean;
    std::map<G4String,G4double> fEmin;
    std::map<G4String,G4double> fEmax;
    G4int    fDecayCount, fTimeCount;
    G4double fEkinTot[3];
    G4double fPbalance[3];
    G4double fEventTime[3];
    G4double fPrimaryTime;

   private:
    G4Timer myTimer;
    ExN01HistoManager*           fHistoManager;
    G4ParticleDefinition* particle;
    G4String partName;
    G4double eprimary;

   private:
    char buf[80];
    char bufup[80];
    char filename[128];
    char pathandname[256];
    FILE* TEST;
    int gpsorgun;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif





