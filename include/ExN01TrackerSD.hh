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

#ifndef ExN01TrackerSD_h
#define ExN01TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ExN01TrackerHit.hh"

#include "G4ThreeVector.hh"

#include "iostream"
#include "vector"
#include <fstream>
#include "map"
using namespace std;
extern int  totalnum[100];

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ExN01TrackerSD : public G4VSensitiveDetector
{

  public:
      ExN01TrackerSD(G4String name);
      ~ExN01TrackerSD();
  
      void   Initialize(G4HCofThisEvent* HCE);
      G4bool ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist);
      void   EndOfEvent(G4HCofThisEvent* HCE);

 //store the deposited energy and the total counts in each pixel
      G4double  store[256*128][3];
      G4double  storepenetration[256*128][2];
      G4int     count[100000][2];
      G4int     count2[100000][2];

//when using energy window to filter some unreasonable values
      vector<int> orderhit;
      vector<double> orderenergy;
      map<int,double> pixel;

  private:
      ExN01TrackerHitsCollection *trackerCollection;
      G4double  energy;
      G4int  numberHits;
//
      long int i,j;
      G4int ii,jj,kk,mm;
      int compare;

//write the results to some out prifles.
      fstream datafile;
      fstream vector1;
      fstream Energy;
      fstream pixelenergy;
      fstream countsperrun;

//to save 256*128 matrix for different runs
      fstream rotate;
      fstream countsforrun;

//to save 256*128 matrix for different runs with 10% energy window.
      fstream windowenergy;
      fstream windowcounts;

//to save 256*128 matrix for penetration.
      fstream penetrationstore;
      fstream windowpenetration;

//to save the deposited energy into 256*128 matrix with /without 10% energy window.
      G4double  diffrun[256*128][90];
      G4double  windowrun[256*128][90];
      G4double  diffrunofcounts[256*128][90];
      G4double  windowrunofcounts[256*128][90];
      G4int     countsdiffrun[180];
      G4int     effectivecounts[180];
      G4int     effectivecounts2[180];

      char filename[128];
      char pathandname[256];
      char buf[80];
      char bufup[80];
      G4int nuclide_type;
      FILE* fp;
      G4int count_effective;
};
#endif

