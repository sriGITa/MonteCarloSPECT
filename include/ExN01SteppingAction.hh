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

#ifndef ExN01SteppingAction_h
#define ExN01SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"

#include "iostream"
#include "vector"
#include <fstream>
using namespace std;

extern  int penetration;
extern  int status1;
extern  int countpenetration;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN01SteppingAction : public G4UserSteppingAction
{
  public:
    ExN01SteppingAction();
   ~ExN01SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    long int num,runID_temp;
    double energy_temp,energy_blur_temp;
    long int count,max,max_blurred,energy_spectrum[10000],energy_blurred[10000];
    vector <int> vi;
    vector <int> vj;
    int display;
    G4ThreeVector  direction3;
    G4ThreeVector  direction4;
    G4int line;
    fstream penetantionrate;
    fstream primary_energy;
    fstream counts_primary;

    G4Material* material2;
    G4int status2;
    char filename[128];
    char pathandname[256];
    char buf[80];


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
