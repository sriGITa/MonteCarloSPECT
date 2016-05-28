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
// * institutes,nor the agencies pr#include "DicomPhantomZSliceHeader.hh"oviding financial support for this *
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
#ifndef ExN01DetectorConstruction_H
#define ExN01DetectorConstruction_H 1

#include "ExN01DicomPhantomZSliceHeader.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4AssemblyVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"
#include "map"
#include <fstream>
#include "iostream"
using namespace std;

class G4Box;
class G4Tubs;
class G4Ellipsoid;

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSensitiveDetector;
class G4Material;
class G4VPVParameterisation;
class ExN01DetectorConstMessenger;
class G4UserLimits;
class ExN01DetectorMessenger;
extern G4double off;
extern G4int phantomtype;


class ExN01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    // Constructor
    ExN01DetectorConstruction();
    // Destructor
    ~ExN01DetectorConstruction();
    // Construct geometry of the setup
   G4VPhysicalVolume* Construct();

   void SetArmAngle(G4double val);
   inline G4double GetArmAngle() { return fArmAngle; }
   void SetCheckOverlaps(G4bool );
   void ConstructPhantombyself();

  private:
     // Construct the world geometry
     G4Box*             solidWorld;    // pointer to the solid envelope 
     G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
     G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
  
     // Construct the collimator geometry
     G4Box*             solidCollimator;
     G4LogicalVolume*   logicCollimator;
     G4VPhysicalVolume* physiCollimator;
     G4VPhysicalVolume* physiCollimatortwo;

     G4Tubs*            solidCircle;
     G4LogicalVolume*   logicCircle;
     G4VPhysicalVolume* physiCircle;

     G4Tubs*            solidCircle2;
     G4LogicalVolume*   logicCircle2;
     G4VPhysicalVolume* physiCircle2;

     G4VPhysicalVolume* BOXphy;
     G4VPhysicalVolume* BOXphytwo;
     G4VPhysicalVolume* sidephy;
     G4VPhysicalVolume* sidephy2;
     G4VPhysicalVolume* sidephytwo;
     G4VPhysicalVolume* sidephy2two;
     G4VPhysicalVolume* leftphy;
     G4VPhysicalVolume* leftphy2;
     G4VPhysicalVolume* leftphytwo;
     G4VPhysicalVolume* leftphy2two;
//
    G4VPhysicalVolume* physiboxDetector1;
    G4VPhysicalVolume* physiboxDetector2;
    G4VPhysicalVolume* physiboxDetector3;
    G4VPhysicalVolume* physiboxDetector4;
    G4VPhysicalVolume* physibox11;
    G4VPhysicalVolume* physibox22;
    G4VPhysicalVolume* physibox33;
    G4VPhysicalVolume* physibox44;
    G4VPhysicalVolume* physiDetector;
    G4VPhysicalVolume* boxphy2;
    G4VPhysicalVolume* physiDetectortwo;
    G4VPhysicalVolume* polyhedraphy11;
    G4VPhysicalVolume* polyhedraphy12;
    G4VPhysicalVolume* polyhedraphy13;
    G4VPhysicalVolume* polyhedraphy14;
    G4VPhysicalVolume* polyhedraphy22;
    G4VPhysicalVolume* polyhedraphy31;
    G4VPhysicalVolume* polyhedraphy32;
    G4VPhysicalVolume* polyhedraphy41;
    G4VPhysicalVolume* polyhedraphy42;
  
//Ellipse

    G4Ellipsoid*       solidEllip;
    G4LogicalVolume*   logicEllip;
    G4VPhysicalVolume* physiEllip; 
    G4VPhysicalVolume* physicylinder;
    G4Ellipsoid* brain;
    G4LogicalVolume* logicBrain;
    G4VPhysicalVolume* physBrain;

    ExN01DetectorConstMessenger* fMessenger;
    G4double fArmAngle;
    G4RotationMatrix* fArmRotation;
    G4RotationMatrix* fArmRotation2;
    G4double x,y;
    G4ThreeVector* fArmRotation1;


protected:
    //create the original materials
    void InitialisationOfMaterial();

    // define materials
    G4Material *Air,*H2O,*Pb,*CdZnTe,*NaI,*Cu;
    G4Material *lunginhale,*lungexhale,*adiposeTissue,*breast,*soft;
    G4Material *water,*muscle,*liver,*trabecularBone,*denseBone;

    std::vector<G4Material*> fOriginalMaterials;
    std::vector<G4Material*> Materialsbysef;

    //read the DICOM files describing the phantom
    void ReadPhantomData();
    G4int fNoFiles;

    //Read one of the DICOM files describing the phantom (usually one per Z slice).
    //Build a DicomPhantomZSliceHeader for each file
    void ReadPhantomDataFile(const G4String& fname);

    // Density difference to distinguish material for each original material (by index)
    std::map<G4int,G4double> fDensityDiffs;
    size_t* fMateIDs; // index of material of each voxel

    // list of new materials created to distinguish different density voxels
    //that have the same original materials
    std::vector<G4Material*> fMaterials;

    //merge the slice headers of all the files
    void MergeZSliceHeaders();
    // z slice header resulted from merging all z slice headers
    ExN01DicomPhantomZSliceHeader* fZSliceHeaderMerged;

    // list of z slice header (one per DICOM files)
    std::vector<ExN01DicomPhantomZSliceHeader*> fZSliceHeaders;

    // build a new material if the density of the voxel is different to the other voxels
    G4Material* BuildMaterialWithChangingDensity( const G4Material* origMate, float density, G4String newMateName );

    // construct the phantom volumes. This method should be implemented for each of the derived classes

    G4Box* fContainer_solid;
    G4LogicalVolume* fContainer_logic;
    G4VPhysicalVolume* fContainer_phys;

    G4int fNVoxelX, fNVoxelY, fNVoxelZ;
    G4double fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ;
    ofstream saveactivity;

    void ConstructPhantomContainer();
    virtual void ConstructPhantom() = 0;
    // construct the phantom volumes.
    //This method should be implemented for each of the derived classes
private:
    G4double Leng,gap,length,thick;
    G4double pos,sourcetocollimator;
    G4double position,center;
    G4Material* WP;
    G4int materialofcollimator,materialofdetector,attenuation;
    G4int disorcon,ii,jj;
    G4LogicalVolume* BOXLog;

    // define the diameter for the hexgon hole
    G4double polyside1,polyside2;

    // define the diameter for the circle hole
    G4double circlerRadius1,circlerRadius2;

    // define the diameter for the rectangular hole
    G4double radius1,radius2;
    G4bool   fCheckOverlaps; // option to activate checking of volumes overlaps
    ifstream parameter;
    ifstream fname;
    fstream indexes;

    char buff[80];
    char buffup[80];
    char NameofFile[128];
    char NameofPath[256];
    FILE* Fp2;
    G4int number,gpsorgun;
    G4float radius;
    G4Material* detectormaterial;

};

#endif

