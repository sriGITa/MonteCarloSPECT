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
//


#ifndef ExN01PrimaryGeneratorAction_h
#define ExN01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4AdjointPrimaryGenerator.hh"
#include "ExN01DetectorConstruction.hh"
#include <fstream>
#include "iostream"
#include "vector"
using namespace std;

class  G4Event;
class  ExN01DetectorConstruction;
class  G4GeneralParticleSource;

class ExN01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExN01PrimaryGeneratorAction();
    ~ExN01PrimaryGeneratorAction();

  public:
 //! defines primary particles (mandatory)
    void GeneratePrimaries(G4Event* anEvent);
    G4ParticleGun* GetParticleGun(){return generalParticleGun;};
    G4GeneralParticleSource* GetParticleSource(){return generalParticleSource;};

  private:
    G4ParticleGun* generalParticleGun;
    G4GeneralParticleSource* generalParticleSource;
    ExN01DetectorConstruction* myDetector;
    ifstream ACTIVITY;
    ifstream numFile;
    G4int num,rows,columns,gpsorgun;
    G4String* fname;
    G4double*** index;
    G4int sum,subsum,b,x,y,z,ord,ii;
    G4double VoxelHalfDimX,VoxelHalfDimY,VoxelHalfDimZ,position;
    G4double ang;
    G4double activity_type;

    ifstream parameterofphantom;
    ifstream fnameofphantom;
    fstream  para;
    G4int phantomlength,phantomwidth,phantomheight,type,vh,II,JJ,KK,bb,N;
    G4float xsize,ysize,zsize,offset,posv,posh,aa,xx,yy,zz;
    char buf[80];
    char bufup[80];
    char filename[128];
    char pathandname[256];
    vector <float> X;
    vector <float> Y;
    vector <float> Z;

    FILE* slice;
    G4int count;
    FILE* test;

};

#endif


