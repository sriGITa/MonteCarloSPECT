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

#include "ExN01DetectorConstruction.hh"
#include "ExN01DetectorConstMessenger.hh"
#include "ExN01TrackerSD.hh"
#include "ExN01TrackerHit.hh"
#include "ExN01RunAction.hh"

#include "G4Ellipsoid.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "globals.hh"
#include "G4PVReplica.hh"

#include "G4VSensitiveDetector.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "math.h"
#include "G4Polyhedra.hh"
#include "G4UIcommand.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "stdlib.h"
#include <stdio.h>
#include <unistd.h>
#include "G4AssemblyVolume.hh"
#include "ExN01DicomPhantomZSliceHeader.hh"

G4double off=0;
G4int phantomtype=2;

ExN01DetectorConstruction::ExN01DetectorConstruction()
// Initialize
 : solidWorld(0),  logicWorld(0),  physiWorld(0),
   solidCollimator(0), logicCollimator(0), physiCollimator(0),
   solidCircle(0), logicCircle(0), physiCircle(0),
   solidCircle2(0), logicCircle2(0), physiCircle2(0),
   BOXphy(0),sidephy(0),leftphy(0),physiboxDetector1(0),
    physiboxDetector2(0),physiboxDetector3(0),physiboxDetector4(0),
   physibox11(0), physibox22(0),physibox33(0),physibox44(0),
   physiDetector(0),boxphy2(0),physiDetectortwo(0),polyhedraphy11(0),polyhedraphy12(0),polyhedraphy13(0),
   polyhedraphy14(0),polyhedraphy22(0),polyhedraphy31(0),
   polyhedraphy32(0),polyhedraphy41(0),polyhedraphy42(0),
   solidEllip(0),logicEllip(0), physiEllip(0),physicylinder(0)
 {  
    sidephy2=0;leftphy2=0;
    fMessenger = new ExN01DetectorConstMessenger(this);
    x = 5.*cm * std::sin(fArmAngle);
    y = 5.*cm * std::cos(fArmAngle);
    fArmRotation = new G4RotationMatrix();
    fArmRotation2 = new G4RotationMatrix();
    //fArmRotation->rotateY(fArmAngle);
    fNoFiles=0;
    fMateIDs = 0;
    fCheckOverlaps=true;

    // get the current directory
    getcwd(buff, sizeof(buff));
    chdir("..");
    // get the parent directory
    getcwd(buffup,sizeof(buffup));
    chdir(buff);

    sprintf(NameofFile,"collimator_detector.txt");
    sprintf(NameofPath,"%s/run/%s",buffup,NameofFile);

    if (!(Fp2 = fopen(NameofPath,"r")))
    {
        printf("Error opening collimator_detector file %s!\n", NameofFile);
        fclose(Fp2);
        exit(1);
    }
    //fseek(Fp2,sizeof("gpsorgun=1"),0);
    fscanf(Fp2,"gpsorgun=%d\n",&gpsorgun);
    fscanf(Fp2,"sourcetocollimator=%lf\n",&sourcetocollimator);
    fscanf(Fp2,"halfLengofcollimator=%lf\n",&Leng);
    fscanf(Fp2,"materialofcollimator=%d\n",&materialofcollimator);
    fscanf(Fp2,"materialofdetector=%d\n",&materialofdetector);
    fscanf(Fp2,"attenuation=%d\n",&attenuation);
    fscanf(Fp2,"disorcon=%d\n",&disorcon);
    fscanf(Fp2,"number=%d\n",&number);
    fscanf(Fp2,"radius1=%lf\n",&radius1);
    fscanf(Fp2,"radius2=%lf\n",&radius2);
    fscanf(Fp2,"phantomtype=%d\n",&phantomtype);
    fscanf(Fp2,"radius=%f\n",&radius);
    fclose(Fp2);

    // units
    Leng=Leng*cm;         // unit is cm
    radius1 = radius1*mm; // unit is mm
    radius2 = radius2*mm; // unit is mm
    radius  = radius *cm; // unit is cm

    if(phantomtype>5)
        phantomtype=5;
    if(phantomtype<1)
        phantomtype=1;
    off=sourcetocollimator+Leng/cm;

    // store related information in indexes.txt
    indexes.open("indexes.txt",ios::out);//overlap
}

ExN01DetectorConstruction::~ExN01DetectorConstruction()
 {
   delete fMessenger;
   delete fArmRotation;
 }

G4VPhysicalVolume* ExN01DetectorConstruction::Construct()
 {

  //---------------------- materials--------------------------------
  InitialisationOfMaterial();

  //----------------------Volume--------------------------------
  //world

  solidWorld= new G4Box("Worldsolid",50*cm,50*cm,50*cm);
  logicWorld= new G4LogicalVolume( solidWorld, Air, "Worldlogic", 0, 0, 0);

  // must place the World Physical volume unrotated at (0,0,0).
  physiWorld = new G4PVPlacement(0,                 // no rotation
                                 G4ThreeVector(),   // at (0,0,0)
                                 logicWorld,        // its logical volume
                                 "Worldphysi",      // its name
                                 0,                 // its mother  volume
                                 fCheckOverlaps,    // no boolean operations
                                 0);                // copy number

  // attenuation=1, attenuation for hoffman phantom with dicom format
  if(attenuation==1)
   {
    ReadPhantomData();
    std::cout<<"The number of image files is "<<fNoFiles<<std::endl;
    ConstructPhantomContainer();
    ConstructPhantom();
  }

  G4int choose;
  choose=100;
  // attenuation=2, attenuation defined with text formats.
  if(attenuation==2)  //
   {
      ConstructPhantombyself();
   }

  //---------------------- collimator--------------------------------
  
  G4double collimatorlength;
  collimatorlength=10.24*cm;
  gap=0.025*cm;    // gap between the collimator and the detector.
  length=Leng+gap; // distance from the collimator center to the detector surface
  solidCollimator= new G4Box("Rectangular-solid-collimator",collimatorlength,collimatorlength,Leng);
  logicCollimator= new G4LogicalVolume( solidCollimator, WP, "Rectangle-logic-collimator", 0, 0, 0);
  physiCollimator = new G4PVPlacement(0,                                // no rotation
                                      G4ThreeVector(0*cm,0 *cm,0*cm),   // at (0,0,0)
                                      logicCollimator,                  // its logical volume
                                      "Rectangular-physi-collimator",   // its name
                                      logicWorld,                       // its mother  volume
                                      fCheckOverlaps,                   // no boolean operations
                                      0);                               // copy no.

  // the cycle holes for collimator
  circlerRadius1=radius1;
  circlerRadius2=radius2;

  solidCircle=new G4Tubs("Tubs",0,circlerRadius1,Leng,0,twopi);
  logicCircle=new G4LogicalVolume(solidCircle,Air,"tubs-log",0,0,0);

  solidCircle2=new G4Tubs("Tubs",0,circlerRadius2,Leng,0,twopi);
  logicCircle2=new G4LogicalVolume(solidCircle2,Air,"tubs-log",0,0,0);

  // the polygon hole
  polyside1=radius1;  // the cross-plat side of the polygon
  polyside2=radius2;

  G4double z_values[2] ={-Leng,Leng};
  G4double RMIN[2]={0.0*mm,0.0*mm};

 // G4BREPSolidPolyhedra* polyhedra = new G4BREPSolidPolyhedra("polyhedra",0,360 *deg,6,2,-10*mm,z_values,RMIN,RMAX1);
 // G4LogicalVolume* polyhedralog=new G4LogicalVolume(polyhedra,Air,"polyhedralog",0,0,0);
  double RMAX1[2]={polyside1,polyside1};
  G4Polyhedra* polyhedra1 = new G4Polyhedra("polyhedra1",0. *deg,360. *deg,6,2,z_values,RMIN,RMAX1);
  G4LogicalVolume* polyhedralog1=new G4LogicalVolume(polyhedra1,Air,"polyhedralog1",0,0,0);

  double RMAX2[2]={polyside2,polyside2};
  G4Polyhedra* polyhedra2 = new G4Polyhedra("polyhedra2",0,360 *deg,6,2,z_values,RMIN,RMAX2);
  G4LogicalVolume* polyhedralog2=new G4LogicalVolume(polyhedra2,Air,"polyhedralog2",0,0,0);

  double RMAX3[2]={polyside1,polyside2};
  G4Polyhedra* polyhedra3 = new G4Polyhedra("polyhedra3",0,360 *deg,6,2,z_values,RMIN,RMAX3);
  G4LogicalVolume* polyhedralog3=new G4LogicalVolume(polyhedra3,Air,"polyhedralog3",0,0,0);

  double RMAX4[2]={polyside2,polyside1};
  G4Polyhedra* polyhedra4 = new G4Polyhedra("polyhedra4",0,360 *deg,6,2,z_values,RMIN,RMAX4);
  G4LogicalVolume* polyhedralog4=new G4LogicalVolume(polyhedra4,Air,"polyhedralog4",0,0,0);

  // the rectangular holes for collimator
  G4VSolid* box11=new G4Box("Box11", radius1, radius1, Leng);
  G4LogicalVolume* box11log=new G4LogicalVolume(box11,Air,"box11log",0,0,0);

  G4VSolid* box22=new G4Box("Box22", radius2, radius2, Leng);
  G4LogicalVolume* box22log=new G4LogicalVolume(box22,Air,"box22log",0,0,0);

  G4VSolid* box33=new G4Box("Box33", radius1, radius2, Leng);
  G4LogicalVolume* box33log=new G4LogicalVolume(box33,Air,"box33log",0,0,0);

  G4VSolid* box44=new G4Box("Box44", radius2, radius1, Leng);
  G4LogicalVolume* box44log=new G4LogicalVolume(box44,Air,"box44log",0,0,0);

  //---------------------- detector--------------------------------

  thick=0.25 *cm;  // thickness of the detector
  G4double centerofdetector=thick+length;
  G4double detectorlength;
  detectorlength=10.24*cm;  // length and width of the detector

  G4VSolid* solidDetector=new G4Box("Box-Solid-detector",detectorlength,detectorlength,thick);
  G4LogicalVolume* logicDetector=new G4LogicalVolume(solidDetector,Cu,"Box-Log-detector",0,0,0);
  physiDetector=new G4PVPlacement(0,G4ThreeVector(0 *cm,0*cm,centerofdetector),logicDetector,"Box-phy-detector",logicWorld,fCheckOverlaps,0);

  // The width and length of the detector hole
  G4double holediameter1,holediameter2;

  // disorcon=1, discrete pixles placed in the detector geometry
  if(disorcon==1)
  {
  holediameter1=0.067*cm;  // the half length of the detector hole
  holediameter2=0.075*cm;  // the half width of the detector hole
  }
  // disorcon=2, there is no gaps which were made of copper.
  if(disorcon==2)
  {
  holediameter1=0.08*cm;  // continuous pixels
  holediameter2=0.08*cm;
  };

  G4VSolid* box1=new G4Box("box-1",holediameter1,holediameter1,thick);
  G4LogicalVolume* box1log=new G4LogicalVolume(box1,detectormaterial,"Box1",0,0,0);

  G4VSolid* box2=new G4Box("box-2",holediameter2,holediameter2 ,thick);
  G4LogicalVolume* box2log=new G4LogicalVolume(box2,detectormaterial,"BOX2",0,0,0);

  G4VSolid* box3=new G4Box("box-3",holediameter1,holediameter2 ,thick);
  G4LogicalVolume* box3log=new G4LogicalVolume(box3,detectormaterial,"Box3",0,0,0);

  G4VSolid* box4=new G4Box("box-4",holediameter2,holediameter1,thick);
  G4LogicalVolume* box4log=new G4LogicalVolume(box4,detectormaterial,"Box4",0,0,0);

  // shielding : back-compartment
  G4VSolid* BOX=new G4Box("BOX",14.24 *cm,14.24 *cm ,2*cm);
  G4LogicalVolume* BOXLog=new G4LogicalVolume(BOX,Pb,"BoxLog",0,0,0);
  BOXphy=new G4PVPlacement(0,G4ThreeVector(0 *cm,0*cm,(2+length/cm+2*thick/cm)*cm),BOXLog,"BOXphy",logicWorld,fCheckOverlaps,0);

  // shielding : up and down sides
  G4double uplength;
  uplength=(Leng/cm+thick/cm+(gap/2)/cm)*cm;
  G4VSolid* side=new G4Box("Side",10.24 *cm,2*cm,uplength);
  G4LogicalVolume* sidelog=new G4LogicalVolume(side,Pb,"Sidelog",0,0,0);
  sidephy=new G4PVPlacement(0,G4ThreeVector(0 *cm,12.24*cm,thick+gap/2),sidelog,"Sidephy",logicWorld,fCheckOverlaps,0);
  sidephy2=new G4PVPlacement(0,G4ThreeVector(0 *cm,-12.24*cm,thick+gap/2),sidelog,"Sidephy",logicWorld,fCheckOverlaps,1);

  // shielding : right and left sides
  G4VSolid* left=new G4Box("Left",2 *cm,14.24*cm,uplength);
  G4LogicalVolume* leftlog=new G4LogicalVolume(left,Pb,"Leftlog",0,0,0);
  leftphy=new G4PVPlacement(0,G4ThreeVector(12.24 *cm,0 *cm,thick+gap/2),leftlog,"Leftphy",logicWorld,fCheckOverlaps,0);
  leftphy2=new G4PVPlacement(0,G4ThreeVector(-12.24 *cm,0 *cm,thick+gap/2),leftlog,"Leftphy",logicWorld,fCheckOverlaps,1);

 // Place the detector holes in proper position
  if(disorcon==1)
  {
      for(int i=0;i<128;i++)
         for (int j=0;j<128;j++)
           {
            if((i%16==0)&&(j%16==0))
                physiboxDetector1 = new G4PVPlacement(0,
                                             G4ThreeVector((-10.152+(i/16)*2.560)*cm,(-10.152+(j/16)*2.560)*cm,0*cm),
                                             box1log,
                                             "BOX1",
                                             logicDetector,
                                             fCheckOverlaps,
                                             i+j*128);
            if((i%16==0)&&((j+1)%16==0))
                physiboxDetector1 = new G4PVPlacement(0,
                                             G4ThreeVector((-10.152+(i/16)*2.560)*cm,(-10.328+((j+1)/16)*2.560)*cm,0*cm),
                                             box1log,
                                             "BOX1",
                                             logicDetector,
                                             fCheckOverlaps,
                                             i+j*128);
            if(((i+1)%16==0)&&(j%16==0))
                physiboxDetector1 = new G4PVPlacement(0,
                                             G4ThreeVector((-10.328+((i+1)/16)*2.560)*cm,(-10.152+(j/16)*2.560)*cm,0*cm),
                                             box1log,
                                             "BOX1",
                                             logicDetector,
                                             fCheckOverlaps,
                                             i+j*128);
            if(((i+1)%16==0)&&((j+1)%16==0))
                physiboxDetector1 = new G4PVPlacement(0,
                                             G4ThreeVector((-10.328+((i+1)/16)*2.560)*cm,(-10.328+((j+1)/16)*2.560)*cm,0*cm),
                                             box1log,
                                             "BOX1",
                                             logicDetector,
                                             fCheckOverlaps,
                                             i+j*128);

            if(i%16==0&&(j+1)%16!=0&&j%16!=0)
                physiboxDetector3 = new G4PVPlacement(0,
                                             G4ThreeVector((-10.152+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                             box3log,
                                             "BOX3",
                                             logicDetector,
                                             fCheckOverlaps,
                                             i+j*128);

            if((i+1)%16==0&&(j+1)%16!=0&&j%16!=0)
                physiboxDetector3 = new G4PVPlacement(0,
                                             G4ThreeVector((-10.328+((i+1)/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                             box3log,
                                             "BOX3",
                                             logicDetector,
                                             fCheckOverlaps,
                                             i+j*128);


            if((i+1)%16!=0&&i%16!=0&&j%16==0)
                physiboxDetector4 = new G4PVPlacement(0,
                                             G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.152+(j/16)*2.560)*cm,0*cm),
                                             box4log,
                                             "BOX4",
                                             logicDetector,
                                             fCheckOverlaps,
                                             i+j*128);

            if((i+1)%16!=0&&i%16!=0&&(j+1)%16==0)
                physiboxDetector4 = new G4PVPlacement(0,
                                             G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.328+((j+1)/16)*2.560)*cm,0*cm),
                                             box4log,
                                             "BOX4",
                                             logicDetector,
                                             fCheckOverlaps,
                                             i+j*128);

           if((i+1)%16!=0&&i%16!=0&&(j+1)%16!=0&&j%16!=0)
            {    physiboxDetector2 = new G4PVPlacement(0,
                                             G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                             box2log,
                                             "BOX2",
                                             logicDetector,
                                             fCheckOverlaps,
                                             i+j*128);
             }
       }

  }

  if(disorcon==2)
  {
 for(int i=0;i<128;i++)
      for(int j=0;j<128;j++)
        {
              if((i%16==0)&&(j%16==0))
                  physiboxDetector1 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                                box1log,
                                               "BOX1",
                                               logicDetector,
                                               fCheckOverlaps,
                                               i+j*128);
              if((i%16==0)&&((j+1)%16==0))
                  physiboxDetector1 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                               box1log,
                                               "BOX1",
                                               logicDetector,
                                               fCheckOverlaps,
                                               i+j*128);
              if(((i+1)%16==0)&&(j%16==0))
                  physiboxDetector1 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                               box1log,
                                               "BOX1",
                                               logicDetector,
                                               fCheckOverlaps,
                                               i+j*128);
              if(((i+1)%16==0)&&((j+1)%16==0))
                  physiboxDetector1 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                               box1log,
                                               "BOX1",
                                               logicDetector,
                                               fCheckOverlaps,
                                               i+j*128);

              if(i%16==0&&(j+1)%16!=0&&j%16!=0)
                  physiboxDetector3 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                               box3log,
                                               "BOX3",
                                               logicDetector,
                                               fCheckOverlaps,
                                               i+j*128);

              if((i+1)%16==0&&(j+1)%16!=0&&j%16!=0)
                  physiboxDetector3 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                               box3log,
                                               "BOX3",
                                               logicDetector,
                                               fCheckOverlaps,
                                               i+j*128);


              if((i+1)%16!=0&&i%16!=0&&j%16==0)
                   physiboxDetector4 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                                box4log,
                                                "BOX4",
                                                logicDetector,
                                                fCheckOverlaps,
                                                i+j*128);

               if((i+1)%16!=0&&i%16!=0&&(j+1)%16==0)
                   physiboxDetector4 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                                box4log,
                                                "BOX4",
                                                logicDetector,
                                                fCheckOverlaps,
                                                i+j*128);

             if((i+1)%16!=0&&i%16!=0&&(j+1)%16!=0&&j%16!=0)
              {    physiboxDetector2 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                               box2log,
                                               "BOX2",
                                               logicDetector,
                                               fCheckOverlaps,
                                               i+j*128);
               }
       }
  }
 // to define different cuts in different sensitive regions
   G4Region* aRegion = new G4Region("BOX1");
   box1log->SetRegion(aRegion);
   aRegion->AddRootLogicalVolume(box1log);

   aRegion = new G4Region("BOX2");
   box2log->SetRegion(aRegion);
   aRegion->AddRootLogicalVolume(box2log);

   aRegion = new G4Region("BOX3");
   box3log->SetRegion(aRegion);
   aRegion->AddRootLogicalVolume(box3log);

   aRegion = new G4Region("BOX4");
   box4log->SetRegion(aRegion);
   aRegion->AddRootLogicalVolume(box4log);

  //choose the type of collimator hole
  // Pitch = 1.6mm for circle
  if(number==1)
   {
    for (ii=0;ii<128;ii++)
       for (jj=0;jj<128;jj++)
               {   physiCircle = new G4PVPlacement(0,
                                      G4ThreeVector((-10.16+ii*0.16)*cm,(-10.16+jj*0.16)*cm,0*cm),
                                      logicCircle,
                                      "rectang-circle",
                                      logicCollimator,
                                      fCheckOverlaps,
                                      0);
     }

   }

  // Pitch = 1.60mm for polygon
  if(number==2)
   {
      for(ii=0;ii<128;ii++)
          for(jj=0;jj<128;jj++)
          {
               polyhedraphy11 = new G4PVPlacement(0,
                                    G4ThreeVector((-10.16+ii*0.16)*cm,(-10.16+jj*0.160)*cm,0*cm),
                                    polyhedralog1,
                                    "rectang-polyhedra1",
                                    logicCollimator,
                                    fCheckOverlaps,
                                    0);
          }
  }

  // Pitch = 1.6mm for rectangular
  if(number==3)
   {
    //box
   for (ii=0;ii<128;ii++)
    for (jj=0;jj<128;jj++)

     {   physibox22 = new G4PVPlacement(0,
                          G4ThreeVector((-10.16+ii*0.16)*cm,(-10.16+jj*0.16)*cm,0*cm),
                          box22log,
                          "rectang-box",
                          logicCollimator,
                          fCheckOverlaps,
                          0);
     }
   }

  // match the detector pixels:circle
  if(number==4)
   {
       for(ii=0;ii<128;ii++)
           for (jj=0;jj<128;jj++)
              {
               if((ii%16==0)&&(jj%16==0))
                   physiCircle = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                logicCircle,
                                                "rectang-Circle",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if((ii%16==0)&&((jj+1)%16==0))
                   physiCircle = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                logicCircle,
                                                "rectang-Circle",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if(((ii+1)%16==0)&&(jj%16==0))
                   physiCircle = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                logicCircle,
                                                "rectang-Circle",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if(((ii+1)%16==0)&&((jj+1)%16==0))
                   physiCircle= new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                logicCircle,
                                                "rectang-Circle",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if(ii%16==0&&(jj+1)%16!=0&&jj%16!=0)
                   physiCircle = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                logicCircle,


                                                "rectang-Circle",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((ii+1)%16==0&&(jj+1)%16!=0&&jj%16!=0)
                   physiCircle = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                logicCircle,
                                                "rectang-Circle",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);


               if((ii+1)%16!=0&&ii%16!=0&&jj%16==0)
                   physiCircle = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                logicCircle,
                                                "rectang-Circle",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16==0)
                   physiCircle = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                logicCircle,
                                                "rectang-Circle",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

             if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16!=0&&jj%16!=0)
               {
                   physiCircle2 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                logicCircle2,
                                                "rectang-Circle2",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               }
           }
    }

  // match the detector pixels:polygon
  if(number==5)
  {
       for(ii=0;ii<128;ii++)
           for (jj=0;jj<128;jj++)
              {
               if((ii%16==0)&&(jj%16==0))
                  polyhedraphy11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog1,
                                                "rectang-polyhedra1",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if((ii%16==0)&&((jj+1)%16==0))
                   polyhedraphy12 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                polyhedralog1,
                                                "rectang-polyhedra1",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if(((ii+1)%16==0)&&(jj%16==0))
                   polyhedraphy13 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog1,
                                                "rectang-polyhedra1",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if(((ii+1)%16==0)&&((jj+1)%16==0))
                   polyhedraphy14 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                polyhedralog1,
                                                "rectang-polyhedra1",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if(ii%16==0&&(jj+1)%16!=0&&jj%16!=0)
                   polyhedraphy31 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog3,
                                               "rectang-polyhedra3",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((ii+1)%16==0&&(jj+1)%16!=0&&jj%16!=0)
                   polyhedraphy32 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog3,
                                                "rectang-polyhedra3",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);


               if((ii+1)%16!=0&&ii%16!=0&&jj%16==0)
                   polyhedraphy41 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog4,
                                                "rectang-polyhedra4",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16==0)
                   polyhedraphy42 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                polyhedralog4,
                                                "rectang-polyhedra4",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

             if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16!=0&&jj%16!=0)
               {
                   polyhedraphy22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog2,
                                                "rectang-polyhedra2",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               }
        }

 }

  // match the detector pixels:rectangular
  if(number==6)
   {

       for(ii=0;ii<128;ii++)
           for (jj=0;jj<128;jj++)
              {
               if((ii%16==0)&&(jj%16==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                box11log,
                                                "rectang-BOX111",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if((ii%16==0)&&((jj+1)%16==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                box11log,
                                                "rectang-BOX112",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if(((ii+1)%16==0)&&(jj%16==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                box11log,
                                                "rectang-BOX113",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if(((ii+1)%16==0)&&((jj+1)%16==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                box11log,
                                                "rectang-BOX114",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if(ii%16==0&&(jj+1)%16!=0&&jj%16!=0)
                   physibox33 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                box33log,
                                                "rectang-BOX33",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((ii+1)%16==0&&(jj+1)%16!=0&&jj%16!=0)
                   physibox33 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                box33log,
                                                "rectang-BOX33",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((ii+1)%16!=0&&ii%16!=0&&jj%16==0)
                   physibox44 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                box44log,
                                                "rectang-BOX44",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16==0)
                   physibox44 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                box44log,
                                                "rectang-BOX44",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

             if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16!=0&&jj%16!=0)
               {
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               }

       }
   }

  // honeycomb array of pixels:polyhedra
  if(number==7)
   {
     for(int kk=0;kk<157;kk++)
        for (int mm=0;mm<136;mm++) {
           if(kk%2==0)  {
               polyhedraphy11=new G4PVPlacement(0,
                                  G4ThreeVector((-9.925+kk*(0.165+0.05*sqrt(3))/2)*cm,(0.055*sqrt(3)/2-9.975+mm*(0.055*sqrt(3)+0.05))*cm,0*cm),
                                  polyhedralog1,
                                  "rectang-polyhedra",
                                  logicCollimator,
                                  fCheckOverlaps,
                                  0);
           }
          if(kk%2==1)   {
               polyhedraphy11=new G4PVPlacement(0,
                                  G4ThreeVector((-9.925+(0.165+0.05*sqrt(3))/2+(kk-1)*(0.165+0.05*sqrt(3))/2)*cm,(0.055*sqrt(3)-9.95+mm*(0.055*sqrt(3)+0.05))*cm,0*cm),
                                  polyhedralog1,
                                  "rectang-polyhedra",
                                  logicCollimator,
                                  fCheckOverlaps,
                                  0);
           }
        }
   }

  // four pixel-matching collimator
  if(number==8)
   {
       for(int i=0;i<64;i++)
           for(int j=0;j<64;j++)
           {
               if((i%8==0)&&(j%8==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.9470+(i/8)*2.5240)*cm,(-9.9470+(j/8)*2.5240)*cm,0*cm),
                                                box11log,
                                               "rectang-BOX111",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if((i%8==0)&&((j+1)%8==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.9470+(i/8)*2.5240)*cm,(-10.245+((j+1)/8)*2.5240)*cm,0*cm),
                                                 box11log,
                                                 "rectang-BOX112",
                                                 logicCollimator,
                                                 fCheckOverlaps,
                                                 0);
               if(((i+1)%8==0)&&(j%8==0))
                  physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.245+((i+1)/8)*2.5240)*cm,(-9.9470+(j/8)*2.5240)*cm,0*cm),
                                                box11log,
                                                "rectang-BOX113",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if(((i+1)%8==0)&&((j+1)%8==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.245+((i+1)/8)*2.5240)*cm,(-10.245+((j+1)/8)*2.5240)*cm,0*cm),
                                                box11log,
                                                "rectang-BOX114",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if(i%8==0&&(j+1)%8!=0&&j%8!=0)
                   physibox33 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.947+(i/8)*2.5240)*cm,(-9.947+(j%8)*0.3180+(j/8)*2.5240)*cm,0*cm),
                                                box33log,
                                                "rectang-BOX33",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((i+1)%8==0&&(j+1)%8!=0&&j%8!=0)
                   physibox33 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.245+((i+1)/8)*2.5240)*cm,(-9.947+(j%8)*0.3180+(j/8)*2.5240)*cm,0*cm),
                                                box33log,
                                                "rectang-BOX33",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);


               if((i+1)%8!=0&&i%8!=0&&j%8==0)
                   physibox44 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.947+(i%8)*0.3180+(i/8)*2.5240)*cm,(-9.947+(j/8)*2.5240)*cm,0*cm),
                                                box44log,
                                                "rectang-BOX44",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((i+1)%8!=0&&i%8!=0&&(j+1)%8==0)
                   physibox44 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.947+(i%8)*0.3180+(i/8)*2.5240)*cm,(-10.245+((j+1)/8)*2.5240)*cm,0*cm),
                                                box44log,
                                                "rectang-BOX44",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);

               if((i+1)%8!=0&&i%8!=0&&(j+1)%8!=0&&j%8!=0)
               {
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.947+(i%8)*0.3180+(i/8)*2.5240)*cm,(-9.947+(j%8)*0.3180+(j/8)*2.5240)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               }
           }
   }

  // the PGM pixels in the same dimension
  if(number==9)
  {
    for(int i=0;i<128;i++)
        for(int j=0;j<128;j++)
         {
               if((i%16==0)&&(j%16==0))
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                                 box22log,
                                                "rectang-BOX22",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if((i%16==0)&&((j+1)%16==0))
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                                  box22log,
                                                 "rectang-BOX22",
                                                 logicCollimator,
                                                 fCheckOverlaps,
                                                 0);
               if(((i+1)%16==0)&&(j%16==0))
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                0);
               if(((i+1)%16==0)&&((j+1)%16==0))
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                 logicCollimator,
                                                 fCheckOverlaps,
                                                 0);

               if(i%16==0&&(j+1)%16!=0&&j%16!=0)
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                 logicCollimator,
                                                 fCheckOverlaps,
                                                 0);

               if((i+1)%16==0&&(j+1)%16!=0&&j%16!=0)
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                 logicCollimator,
                                                 fCheckOverlaps,
                                                 0);


               if((i+1)%16!=0&&i%16!=0&&j%16==0)
                    physibox22 = new G4PVPlacement(0,
                                                 G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                                 box22log,
                                                 "rectang-BOX22",
                                                 logicCollimator,
                                                 fCheckOverlaps,
                                                 0);

                if((i+1)%16!=0&&i%16!=0&&(j+1)%16==0)
                    physibox22 = new G4PVPlacement(0,
                                                 G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                                 box22log,
                                                 "rectang-BOX22",
                                                 logicCollimator,
                                                 fCheckOverlaps,
                                                 0);

              if((i+1)%16!=0&&i%16!=0&&(j+1)%16!=0&&j%16!=0)
               {    physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                logicCollimator,
                                                fCheckOverlaps,
                                                 0);
                }
        }
}

  if(attenuation==3)
  choose=2;

  //choose has to be 1 or 2,1 for a water-filled ellipsoid and 2 for a water-filled cylinder
  if(choose==1)
  {
   solidEllip=new G4Ellipsoid("ellipse", 2*cm,8 *cm,10*cm);
   logicEllip=new G4LogicalVolume(solidEllip,H2O,"Ellipse", 0, 0, 0);
   physiEllip=new G4PVPlacement(fArmRotation,
                                G4ThreeVector(0,0,-11.25*cm),
                                logicEllip,
                                "Ellipse-phy",
                                logicWorld,
                                fCheckOverlaps,
                                0);
  }

  G4Tubs* cylinder=new G4Tubs("Tubs",0,radius,10*cm,0,twopi);
  G4LogicalVolume* logiccylinder=new G4LogicalVolume(cylinder,H2O,"tubs-log",0,0,0);
  //physicylinder=new G4PVPlacement(0,G4ThreeVector(0 *cm,0*cm,-off*cm),logiccylinder,"tubs-phy",logicWorld,fCheckOverlaps,0);
  if(choose==2)
  {
  G4RotationMatrix* rot;
  rot=new G4RotationMatrix();
  rot->rotateX(90*deg);
  fArmRotation2->rotateX(90*deg);
  physicylinder=new G4PVPlacement(fArmRotation2,G4ThreeVector(0 *cm,0*cm,-off*cm),logiccylinder,"tubs-phy",logicWorld,fCheckOverlaps,0);
  }

  if(choose==3)
  {

    // define the brain of the human phantom
     G4double ax = 6. * cm;
     G4double by= 9. * cm;
     G4double cz = 6.5 * cm;

     brain = new G4Ellipsoid("Brain", ax, by, cz);
     logicBrain =  new G4LogicalVolume(brain, soft,"logicalbrain", 0, 0, 0);
     // Define rotation and position here!
     physBrain = new G4PVPlacement(0,
                               G4ThreeVector(0.*cm, 0.*cm, -11.25* cm),
                               "physicalBrain",
                               logicBrain,
                               physiWorld,
                               fCheckOverlaps,
                               0, true);

     /*
        //define the head of the human phantom
        // MIRD male model
        // Ellipsoid
        G4double ax1 = 7.0 * cm;
        G4double by1 = 10.0 * cm;
        G4double cz1 = 8.50 * cm;
        G4double zcut1 = 0.0 * cm;
        G4double zcut2 = 8.5 * cm;

        G4Ellipsoid* head1 = new G4Ellipsoid("Head1", ax1, by1, cz1, zcut1, zcut2);

        G4double dx = 7.0 * cm;
        G4double dy = 10.0 * cm;
        G4double dz = 7.75 * cm;


        G4EllipticalTube* head2 = new G4EllipticalTube("Head2", dx, dy, dz);

        G4UnionSolid* head = new G4UnionSolid("Head",head2,head1,
                      0, // Rotation
                      G4ThreeVector(0.* cm, 0.*cm, 7.7500 * cm) );

        G4LogicalVolume* logicHead = new G4LogicalVolume(head, H2O,"logicalhead",
                                 0, 0,0);
        G4RotationMatrix* rm = new G4RotationMatrix();
        rm -> rotateX(90.* degree);

        // Define rotation and position here!
        G4VPhysicalVolume* physHead = new G4PVPlacement(rm,
                                G4ThreeVector(0.* cm,0 *cm, -11.25*cm),
                                logicHead,
                                "physicalHead",
                                logicWorld,
                                fCheckOverlaps,
                                0, true);

        // Testing Head Volume
        G4double HeadVol = logicHead->GetSolid()->GetCubicVolume();
        std::cout << "Volume of Head = " << HeadVol/cm3 << " cm^3" << std::endl;

        // Testing Head Material
        G4String HeadMat = logicHead->GetMaterial()->GetName();
        std::cout << "Material of Head = " << HeadMat << std::endl;

        // Testing Density
        G4double HeadDensity = logicHead->GetMaterial()->GetDensity();
        std::cout << "Density of Material = " << HeadDensity*cm3/g << " g/cm^3" << std::endl;

        // Testing Mass
        G4double HeadMass = (HeadVol)*HeadDensity;
        std::cout << "Mass of Head = " << HeadMass/gram << " g" << std::endl;

   */
  }

 //********************* Second detector *****************

  G4double position=-(sourcetocollimator*2+Leng*2/cm)*cm;
  G4double center=position-thick-length;
  G4VSolid* solidDetectortwo=new G4Box("Box-Solidtwo",detectorlength,detectorlength ,thick);
  G4LogicalVolume* logicDetectortwo=new G4LogicalVolume(solidDetectortwo,Cu,"Box-Logtwo",0,0,0);
  physiDetectortwo=new G4PVPlacement(0,G4ThreeVector(0 *cm,0*cm,center),logicDetectortwo,"Box-phytwo",logicWorld,fCheckOverlaps,0);

  // place the second detector
  if(disorcon==1)
  {
    for(int i=0;i<128;i++)
       for (int j=0;j<128;j++)
         {
          if((i%16==0)&&(j%16==0))
              physiboxDetector1 = new G4PVPlacement(0,
                                           G4ThreeVector((-10.152+(i/16)*2.560)*cm,(-10.152+(j/16)*2.560)*cm,0*cm),
                                           box1log,
                                           "BOX1",
                                           logicDetectortwo,
                                           fCheckOverlaps,
                                           i+j*128+16384);
          if((i%16==0)&&((j+1)%16==0))
              physiboxDetector1 = new G4PVPlacement(0,
                                           G4ThreeVector((-10.152+(i/16)*2.560)*cm,(-10.328+((j+1)/16)*2.560)*cm,0*cm),
                                           box1log,
                                           "BOX1",
                                           logicDetectortwo,
                                           fCheckOverlaps,
                                           i+j*128+16384);
          if(((i+1)%16==0)&&(j%16==0))
              physiboxDetector1 = new G4PVPlacement(0,
                                           G4ThreeVector((-10.328+((i+1)/16)*2.560)*cm,(-10.152+(j/16)*2.560)*cm,0*cm),
                                           box1log,
                                           "BOX1",
                                           logicDetectortwo,
                                           fCheckOverlaps,
                                           i+j*128+16384);
          if(((i+1)%16==0)&&((j+1)%16==0))
              physiboxDetector1 = new G4PVPlacement(0,
                                           G4ThreeVector((-10.328+((i+1)/16)*2.560)*cm,(-10.328+((j+1)/16)*2.560)*cm,0*cm),
                                           box1log,
                                           "BOX1",
                                           logicDetectortwo,
                                           fCheckOverlaps,
                                           i+j*128+16384);

          if(i%16==0&&(j+1)%16!=0&&j%16!=0)
              physiboxDetector3 = new G4PVPlacement(0,
                                           G4ThreeVector((-10.152+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                           box3log,
                                           "BOX3",
                                           logicDetectortwo,
                                           fCheckOverlaps,
                                           i+j*128+16384);

          if((i+1)%16==0&&(j+1)%16!=0&&j%16!=0)
              physiboxDetector3 = new G4PVPlacement(0,
                                           G4ThreeVector((-10.328+((i+1)/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                           box3log,
                                           "BOX3",
                                           logicDetectortwo,
                                           fCheckOverlaps,
                                           i+j*128+16384);


          if((i+1)%16!=0&&i%16!=0&&j%16==0)
              physiboxDetector4 = new G4PVPlacement(0,
                                           G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.152+(j/16)*2.560)*cm,0*cm),
                                           box4log,
                                           "BOX4",
                                           logicDetectortwo,
                                           fCheckOverlaps,
                                           i+j*128+16384);

          if((i+1)%16!=0&&i%16!=0&&(j+1)%16==0)
              physiboxDetector4 = new G4PVPlacement(0,
                                           G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.328+((j+1)/16)*2.560)*cm,0*cm),
                                           box4log,
                                           "BOX4",
                                           logicDetectortwo,
                                           fCheckOverlaps,
                                           i+j*128+16384);

         if((i+1)%16!=0&&i%16!=0&&(j+1)%16!=0&&j%16!=0)
          {    physiboxDetector2 = new G4PVPlacement(0,
                                           G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                           box2log,
                                           "BOX2",
                                           logicDetectortwo,
                                           fCheckOverlaps,
                                           i+j*128+16384);
           }
     }
  }

  if(disorcon==2)
  {
      for(int i=0;i<128;i++)
          for(int j=0;j<128;j++)
        {
              if((i%16==0)&&(j%16==0))
                  physiboxDetector1 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                                box1log,
                                               "BOX1",
                                               logicDetectortwo,
                                               fCheckOverlaps,
                                               i+j*128+16384);
              if((i%16==0)&&((j+1)%16==0))
                  physiboxDetector1 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                               box1log,
                                               "BOX1",
                                               logicDetectortwo,
                                               fCheckOverlaps,
                                               i+j*128+16384);
              if(((i+1)%16==0)&&(j%16==0))
                  physiboxDetector1 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                               box1log,
                                               "BOX1",
                                               logicDetectortwo,
                                               fCheckOverlaps,
                                               i+j*128+16384);
              if(((i+1)%16==0)&&((j+1)%16==0))
                  physiboxDetector1 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                               box1log,
                                               "BOX1",
                                               logicDetectortwo,
                                               fCheckOverlaps,
                                               i+j*128+16384);

              if(i%16==0&&(j+1)%16!=0&&j%16!=0)
                  physiboxDetector3 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                               box3log,
                                               "BOX3",
                                               logicDetectortwo,
                                               fCheckOverlaps,
                                               i+j*128+16384);

              if((i+1)%16==0&&(j+1)%16!=0&&j%16!=0)
                  physiboxDetector3 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                               box3log,
                                               "BOX3",
                                               logicDetectortwo,
                                               fCheckOverlaps,
                                               i+j*128+16384);


              if((i+1)%16!=0&&i%16!=0&&j%16==0)
                   physiboxDetector4 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                                box4log,
                                                "BOX4",
                                                logicDetectortwo,
                                                fCheckOverlaps,
                                                i+j*128+16384);

               if((i+1)%16!=0&&i%16!=0&&(j+1)%16==0)
                   physiboxDetector4 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                                box4log,
                                                "BOX4",
                                                logicDetectortwo,
                                                fCheckOverlaps,
                                                i+j*128+16384);

             if((i+1)%16!=0&&i%16!=0&&(j+1)%16!=0&&j%16!=0)
              {    physiboxDetector2 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                               box2log,
                                               "BOX2",
                                               logicDetectortwo,
                                               fCheckOverlaps,
                                               i+j*128+16384);
               }
       }

  }

 //********************* Second collimator *****************
  G4LogicalVolume* logicCollimatortwo= new G4LogicalVolume(solidCollimator, WP, "Rectangle-logictwo", 0, 0, 0);
  physiCollimatortwo = new G4PVPlacement(0,                              // no rotation
                                         G4ThreeVector(0*cm,0 *cm,position),
                                         logicCollimatortwo,               // its logical volume
                                         "Rectangular-physitwo",            // its name
                                         logicWorld,                     // its mother  volume
                                         fCheckOverlaps,                          // no boolean operations
                                         0);
 // Pitch=1.6mm for circle
  if(number==1)
  {
     for (ii=0;ii<128;ii++){
          for (jj=0;jj<128;jj++)
            {    physiCircle = new G4PVPlacement(0,
                                             G4ThreeVector((-10.16+ii*0.16)*cm,(-10.16+jj*0.16)*cm,0*cm),
                                             logicCircle,
                                             "rectang-circle2",
                                             logicCollimatortwo,
                                             fCheckOverlaps,
                                             0);
           }
       }
  }

 // Pitch=1.60mm for polygon
  if(number==2)
  {
      for(ii=0;ii<128;ii++){
          for(jj=0;jj<128;jj++)
          {
             polyhedraphy11 = new G4PVPlacement(0,
                                  G4ThreeVector((-10.16+ii*0.16)*cm,(-10.16+jj*0.160)*cm,0*cm),
                                  polyhedralog1,
                                  "rectang-polyhedra2",
                                  logicCollimatortwo,
                                  fCheckOverlaps,
                                  0);
          }
      }
  }

 // Pitch=1.6mm for rectangular
  if(number==3)
  {
      for(ii=0;ii<128;ii++)
          for(jj=0;jj<128;jj++)
           {
              physibox22 = new G4PVPlacement(0,
                               G4ThreeVector((-10.16+ii*0.16)*cm,(-10.16+jj*0.16)*cm,0*cm),
                               box22log,
                               "rectang-box2",
                                logicCollimatortwo,
                                fCheckOverlaps,
                                0);
      }
  }

 // match the detector pixels:circle
  if(number==4)
   {
      for(ii=0;ii<128;ii++)
          for (jj=0;jj<128;jj++)
             {
              if((ii%16==0)&&(jj%16==0))
                  physiCircle = new G4PVPlacement(0,
                                               G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                               logicCircle,
                                               "rectang-Circle",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              if((ii%16==0)&&((jj+1)%16==0))
                  physiCircle = new G4PVPlacement(0,
                                               G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                               logicCircle,
                                               "rectang-Circle",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              if(((ii+1)%16==0)&&(jj%16==0))
                  physiCircle = new G4PVPlacement(0,
                                               G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                               logicCircle,
                                               "rectang-Circle",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              if(((ii+1)%16==0)&&((jj+1)%16==0))
                  physiCircle= new G4PVPlacement(0,
                                               G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                               logicCircle,
                                               "rectang-Circle",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);

              if(ii%16==0&&(jj+1)%16!=0&&jj%16!=0)
                  physiCircle = new G4PVPlacement(0,
                                               G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                               logicCircle,
                                               "rectang-Circle",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);

              if((ii+1)%16==0&&(jj+1)%16!=0&&jj%16!=0)
                  physiCircle = new G4PVPlacement(0,
                                               G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                               logicCircle,
                                               "rectang-Circle",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);


              if((ii+1)%16!=0&&ii%16!=0&&jj%16==0)
                  physiCircle = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                               logicCircle,
                                               "rectang-Circle",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);

              if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16==0)
                  physiCircle = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                               logicCircle,
                                               "rectang-Circle",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);

            if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16!=0&&jj%16!=0)
              {
                  physiCircle2 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                               logicCircle2,
                                               "rectang-Circle",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              }
          }
  }

 // match the detector pixels:polygon
  if(number==5)
  {
       for(ii=0;ii<128;ii++)
           for (jj=0;jj<128;jj++)
              {
               if((ii%16==0)&&(jj%16==0))
                  polyhedraphy11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog1,
                                                "rectang-polyhedra1",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);
               if((ii%16==0)&&((jj+1)%16==0))
                   polyhedraphy12 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                polyhedralog1,
                                                "rectang-polyhedra1",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);
               if(((ii+1)%16==0)&&(jj%16==0))
                   polyhedraphy13 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog1,
                                                "rectang-polyhedra1",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);
               if(((ii+1)%16==0)&&((jj+1)%16==0))
                   polyhedraphy14 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                polyhedralog1,
                                                "rectang-polyhedra1",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

               if(ii%16==0&&(jj+1)%16!=0&&jj%16!=0)
                   polyhedraphy31 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog3,
                                               "rectang-polyhedra3",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

               if((ii+1)%16==0&&(jj+1)%16!=0&&jj%16!=0)
                   polyhedraphy32 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog3,
                                                "rectang-polyhedra3",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);


               if((ii+1)%16!=0&&ii%16!=0&&jj%16==0)
                   polyhedraphy41 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog4,
                                                "rectang-polyhedra4",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

               if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16==0)
                   polyhedraphy42 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                                polyhedralog4,
                                                "rectang-polyhedra4",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

             if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16!=0&&jj%16!=0)
               {
                   polyhedraphy22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                                polyhedralog2,
                                                "rectang-polyhedra2",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);
               }
        }

 }

 // match the detector pixels:rectangular
  if(number==6)
  {
      for(ii=0;ii<128;ii++)
          for (jj=0;jj<128;jj++)
             {
              if((ii%16==0)&&(jj%16==0))
                  physibox11 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                               box11log,
                                               "rectang-BOX111",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              if((ii%16==0)&&((jj+1)%16==0))
                  physibox11 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                               box11log,
                                               "rectang-BOX112",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              if(((ii+1)%16==0)&&(jj%16==0))
                  physibox11 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                               box11log,
                                               "rectang-BOX113",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              if(((ii+1)%16==0)&&((jj+1)%16==0))
                  physibox11 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                               box11log,
                                               "rectang-BOX114",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);

              if(ii%16==0&&(jj+1)%16!=0&&jj%16!=0)
                  physibox33 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.152+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                               box33log,
                                               "rectang-BOX33",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);

              if((ii+1)%16==0&&(jj+1)%16!=0&&jj%16!=0)
                  physibox33 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.328+((ii+1)/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                               box33log,
                                               "rectang-BOX33",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);

              if((ii+1)%16!=0&&ii%16!=0&&jj%16==0)
                  physibox44 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.152+(jj/16)*2.560)*cm,0*cm),
                                               box44log,
                                               "rectang-BOX44",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);

              if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16==0)
                  physibox44 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.328+((jj+1)/16)*2.560)*cm,0*cm),
                                               box44log,
                                               "rectang-BOX44",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);

            if((ii+1)%16!=0&&ii%16!=0&&(jj+1)%16!=0&&jj%16!=0)
              {
                  physibox22 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(ii%16)*0.160+(ii/16)*2.560)*cm,(-10.16+(jj%16)*0.160+(jj/16)*2.560)*cm,0*cm),
                                               box22log,
                                               "rectang-BOX22",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              }

      }
  }

 // honeycomb array of pixels:polyhedra
  if(number==7)
  {
    for(int kk=0;kk<157;kk++)
       for (int mm=0;mm<136;mm++) {
          if(kk%2==0)  {
              polyhedraphy11=new G4PVPlacement(0,
                                 G4ThreeVector((-9.925+kk*(0.165+0.05*sqrt(3))/2)*cm,(0.055*sqrt(3)/2-9.975+mm*(0.055*sqrt(3)+0.05))*cm,0*cm),
                                 polyhedralog1,
                                 "rectang-polyhedra",
                                 logicCollimatortwo,
                                 fCheckOverlaps,
                                 0);
          }
         if(kk%2==1)   {
              polyhedraphy11=new G4PVPlacement(0,
                                 G4ThreeVector((-9.925+(0.165+0.05*sqrt(3))/2+(kk-1)*(0.165+0.05*sqrt(3))/2)*cm,(0.055*sqrt(3)-9.95+mm*(0.055*sqrt(3)+0.05))*cm,0*cm),
                                 polyhedralog1,
                                 "rectang-polyhedra",
                                 logicCollimatortwo,
                                 fCheckOverlaps,
                                 0);
          }
       }
  }

 // four pixel-matching collimator
  if(number==8)
   {
       for(int i=0;i<64;i++)
           for(int j=0;j<64;j++)
           {
               if((i%8==0)&&(j%8==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.9470+(i/8)*2.5240)*cm,(-9.9470+(j/8)*2.5240)*cm,0*cm),
                                                box11log,
                                               "rectang-BOX111",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);
               if((i%8==0)&&((j+1)%8==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.9470+(i/8)*2.5240)*cm,(-10.245+((j+1)/8)*2.5240)*cm,0*cm),
                                                 box11log,
                                                 "rectang-BOX112",
                                                 logicCollimatortwo,
                                                 fCheckOverlaps,
                                                 0);
               if(((i+1)%8==0)&&(j%8==0))
                  physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.245+((i+1)/8)*2.5240)*cm,(-9.9470+(j/8)*2.5240)*cm,0*cm),
                                                box11log,
                                                "rectang-BOX113",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);
               if(((i+1)%8==0)&&((j+1)%8==0))
                   physibox11 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.245+((i+1)/8)*2.5240)*cm,(-10.245+((j+1)/8)*2.5240)*cm,0*cm),
                                                box11log,
                                                "rectang-BOX114",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

               if(i%8==0&&(j+1)%8!=0&&j%8!=0)
                   physibox33 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.947+(i/8)*2.5240)*cm,(-9.947+(j%8)*0.3180+(j/8)*2.5240)*cm,0*cm),
                                                box33log,
                                                "rectang-BOX33",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

               if((i+1)%8==0&&(j+1)%8!=0&&j%8!=0)
                   physibox33 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.245+((i+1)/8)*2.5240)*cm,(-9.947+(j%8)*0.3180+(j/8)*2.5240)*cm,0*cm),
                                                box33log,
                                                "rectang-BOX33",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);


               if((i+1)%8!=0&&i%8!=0&&j%8==0)
                   physibox44 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.947+(i%8)*0.3180+(i/8)*2.5240)*cm,(-9.947+(j/8)*2.5240)*cm,0*cm),
                                                box44log,
                                                "rectang-BOX44",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

               if((i+1)%8!=0&&i%8!=0&&(j+1)%8==0)
                   physibox44 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.947+(i%8)*0.3180+(i/8)*2.5240)*cm,(-10.245+((j+1)/8)*2.5240)*cm,0*cm),
                                                box44log,
                                                "rectang-BOX44",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

               if((i+1)%8!=0&&i%8!=0&&(j+1)%8!=0&&j%8!=0)
               {
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-9.947+(i%8)*0.3180+(i/8)*2.5240)*cm,(-9.947+(j%8)*0.3180+(j/8)*2.5240)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);
               }
           }
   }

 // the PGM pixels in the same dimension
  if(number==9)
  {
      for(int i=0;i<128;i++)
          for(int j=0;j<128;j++)
         {
              if((i%16==0)&&(j%16==0))
                  physibox22 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                                box22log,
                                               "rectang-BOX22",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              if((i%16==0)&&((j+1)%16==0))
                  physibox22 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                                 box22log,
                                                "rectang-BOX22",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);
              if(((i+1)%16==0)&&(j%16==0))
                  physibox22 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                               box22log,
                                               "rectang-BOX22",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                               0);
              if(((i+1)%16==0)&&((j+1)%16==0))
                  physibox22 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                               box22log,
                                               "rectang-BOX22",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

              if(i%16==0&&(j+1)%16!=0&&j%16!=0)
                  physibox22 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                               box22log,
                                               "rectang-BOX22",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

              if((i+1)%16==0&&(j+1)%16!=0&&j%16!=0)
                  physibox22 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.32+((i+1)/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                               box22log,
                                               "rectang-BOX22",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);


              if((i+1)%16!=0&&i%16!=0&&j%16==0)
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j/16)*2.560)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

               if((i+1)%16!=0&&i%16!=0&&(j+1)%16==0)
                   physibox22 = new G4PVPlacement(0,
                                                G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.32+((j+1)/16)*2.560)*cm,0*cm),
                                                box22log,
                                                "rectang-BOX22",
                                                logicCollimatortwo,
                                                fCheckOverlaps,
                                                0);

             if((i+1)%16!=0&&i%16!=0&&(j+1)%16!=0&&j%16!=0)
              {    physibox22 = new G4PVPlacement(0,
                                               G4ThreeVector((-10.16+(i%16)*0.160+(i/16)*2.560)*cm,(-10.16+(j%16)*0.160+(j/16)*2.560)*cm,0*cm),
                                               box22log,
                                               "rectang-BOX22",
                                               logicCollimatortwo,
                                               fCheckOverlaps,
                                                0);
               }
       }
}

 // shielding: rear
  G4LogicalVolume* BOXLogtwo=new G4LogicalVolume(BOX,Pb,"BoxLogtwo",0,0,0);
  BOXphytwo=new G4PVPlacement(0,G4ThreeVector(0 *cm,0*cm,(center/cm-thick/cm-2)*cm),BOXLogtwo,"BOXphytwo",logicWorld,fCheckOverlaps,0);

 // shielding: top and bottom
  sidephytwo=new G4PVPlacement(0,G4ThreeVector(0 *cm,12.24*cm,position-thick-gap/2),sidelog,"Sidephy",logicWorld,fCheckOverlaps,0);
  sidephy2two=new G4PVPlacement(0,G4ThreeVector(0 *cm,-12.24*cm,position-thick-gap/2),sidelog,"Sidephy",logicWorld,fCheckOverlaps,1);

 // shielding: left and right
  leftphytwo=new G4PVPlacement(0,G4ThreeVector(12.24 *cm,0 *cm,position-thick-gap/2),leftlog,"Leftphy",logicWorld,fCheckOverlaps,0);
  leftphy2two=new G4PVPlacement(0,G4ThreeVector(-12.24 *cm,0 *cm,position-thick-gap/2),leftlog,"Leftphy",logicWorld,fCheckOverlaps,1);

 // ******************* Visualization *******************

  G4Color
  yellow(1.0,1.0,0.0),
  green(0.0,1.0,0.0),
  blue(0.0,0.0,1.0),
  brown(0.4,0.4,0.1),
  red(1.0,0.0,0.0),
  gray(0.5,0.5,0.5),
  black(0.0,0.0,0.0),
  cyan(0.0,1.0,1.0),
  magenta(1.0,0.0,1.0),
  white(1.0,1.0,1.0);

  logicWorld -> SetVisAttributes(new G4VisAttributes(white)); //world

  G4VisAttributes* simple= new G4VisAttributes(yellow);
  simple->SetVisibility(false);
  simple->SetForceSolid(false);
  logicCollimator-> SetVisAttributes(simple);//collimator
  logicCollimatortwo->SetVisAttributes(simple);

  if(number==1)
   {
       G4VisAttributes* viscirle= new G4VisAttributes(black);//color: white
       viscirle->SetVisibility(false);
       viscirle->SetForceSolid(false);//object set to false when displayed as a wireframe
       logicCircle->SetVisAttributes(viscirle);
   }

  if(number==2)
  {
      G4VisAttributes* vispoly= new G4VisAttributes(black);
      vispoly->SetVisibility(false);
      vispoly->SetForceSolid(false);
      polyhedralog1->SetVisAttributes(vispoly);
  }

  if(number==3||number==9)
   {
       G4VisAttributes* visrect= new G4VisAttributes(black);
       visrect->SetVisibility(false);
       visrect->SetForceSolid(false);
       box22log->SetVisAttributes(visrect);

   }

  if(number==4)
  {
      G4VisAttributes* vis01= new G4VisAttributes(green);
      vis01->SetVisibility(false);
      vis01->SetForceSolid(false);
      logicCircle->SetVisAttributes(vis01);
      vis01= new G4VisAttributes(yellow);
      vis01->SetVisibility(false);
      vis01->SetForceSolid(false);
      logicCircle2->SetVisAttributes(vis01);
  }

  if(number==5)
  {
      G4VisAttributes* vis02= new G4VisAttributes(green);
      vis02->SetVisibility(false);
      vis02->SetForceSolid(false);
      polyhedralog1->SetVisAttributes(vis02);
      vis02= new G4VisAttributes(gray);
      vis02->SetVisibility(false);
      vis02->SetForceSolid(false);
      polyhedralog2->SetVisAttributes(vis02);
      vis02= new G4VisAttributes(blue);
      vis02->SetVisibility(false);
      vis02->SetForceSolid(false);
      polyhedralog3->SetVisAttributes(vis02);
      vis02= new G4VisAttributes(black);
      vis02->SetVisibility(false);
      vis02->SetForceSolid(false);
      polyhedralog4->SetVisAttributes(vis02);
  }

  if(number==7)
   {
       G4VisAttributes* vis22= new G4VisAttributes(brown);
       vis22->SetVisibility(true);
       vis22->SetForceSolid(false);
       polyhedralog1->SetVisAttributes(vis22);

   }

  if(number==6||number==8)
    {
      G4VisAttributes* vis03= new G4VisAttributes(white);
      vis03->SetVisibility(false);
      vis03->SetForceSolid(false);
      box11log->SetVisAttributes(vis03);

      vis03= new G4VisAttributes(green);
      vis03->SetVisibility(false);
      vis03->SetForceSolid(false);
      box22log->SetVisAttributes(vis03);

      vis03= new G4VisAttributes(magenta);
      vis03->SetVisibility(false);
      vis03->SetForceSolid(false);
      box33log->SetVisAttributes(vis03);

      vis03= new G4VisAttributes(yellow);
      vis03->SetVisibility(false);
      vis03->SetForceSolid(false);
      box44log->SetVisAttributes(vis03);
   }

   G4VisAttributes* simpleAlSVisAtt_Ellipse= new G4VisAttributes(red);
   simpleAlSVisAtt_Ellipse->SetVisibility(true);
   simpleAlSVisAtt_Ellipse->SetForceSolid(false);
   logicDetector->SetVisAttributes(simpleAlSVisAtt_Ellipse);
   logicDetectortwo->SetVisAttributes(simpleAlSVisAtt_Ellipse);

   G4VisAttributes* vis1= new G4VisAttributes(green);
   vis1->SetVisibility(false);
   vis1->SetForceSolid(false);
   box1log->SetVisAttributes(vis1);  // the first detector

   G4VisAttributes* vis2= new G4VisAttributes(red);
   vis2->SetVisibility(false);
   vis2->SetForceSolid(false);
   box2log->SetVisAttributes(vis2);

   G4VisAttributes* vis3= new G4VisAttributes(cyan);
   vis3->SetVisibility(false);
   vis3->SetForceSolid(false);
   box3log->SetVisAttributes(vis3);

   G4VisAttributes* vis4= new G4VisAttributes(yellow);
   vis4->SetVisibility(false);
   vis4->SetForceSolid(false);
   box4log->SetVisAttributes(vis4);

   G4VisAttributes* vis5= new G4VisAttributes(yellow);
   vis5->SetVisibility(true);
   vis5->SetForceSolid(false);

   if(choose==1)
    {  logicEllip-> SetVisAttributes(vis5);  }

   if(choose==2)
    {  logiccylinder->SetVisAttributes(vis5);}

   if(choose==3)
   {   logicBrain->SetVisAttributes(vis5);   }

   G4VisAttributes* vis6= new G4VisAttributes(cyan);
   vis6->SetVisibility(true);
   vis6->SetForceSolid(false);
   BOXLog-> SetVisAttributes(vis6);
   BOXLogtwo-> SetVisAttributes(vis6);

 // ************** Sensitivie Dectector *****************
   G4String trackerChamberSDname = "ExN01/TrackerChamberSD";
   G4VSensitiveDetector* sensitiveBox=new ExN01TrackerSD(trackerChamberSDname);
   G4SDManager* SDManager=G4SDManager::GetSDMpointer();
   SDManager->AddNewDetector(sensitiveBox);

// logicDetector->SetSensitiveDetector(sensitiveBox);
// logicDetectortwo->SetSensitiveDetector(sensitiveBox);
// BOXLog->SetSensitiveDetector(sensitiveBox);
// BOXLogtwo->SetSensitiveDetector(sensitiveBox);
   box1log->SetSensitiveDetector(sensitiveBox);
   box2log->SetSensitiveDetector(sensitiveBox);
   box3log->SetSensitiveDetector(sensitiveBox);
   box4log->SetSensitiveDetector(sensitiveBox);

// output some information
   indexes<<"sourcetocollimator="<<sourcetocollimator<<"cm"<<" // distance from source to the collimator surface (unit: cm)."<<std::endl;
   indexes<<"Lengofcollimator="<<Leng*2/mm<<"mm"<<" // length of the collimator (unit: cm)."<<std::endl;
   indexes<<"materialofcollimator="<<materialofcollimator<<" // 1 for W-19.34; 2 for W-17.0; 3 for lead "<<std::endl;
   indexes<<"materialofdetector="<<materialofdetector<<" // 1 for CZT; 2 for NaI "<<std::endl;
   indexes<<"attenuation="<<attenuation<<"  // 1 for Hoffman phantom attenuation (DICOM); 2 for user defined attenuation; 3 for water attenuation."<<std::endl;
   if(attenuation==3)  {
   indexes<<"choose="<<choose<<"  // 2 for water attenuation"<<std::endl;
   indexes<<"radius="<<radius/cm<<"cm"<<"  // radius of water phantom"<<std::endl;}

   indexes<<"disorcon="<<disorcon<<" // 1 for discrete detector; 2 for continuous detector"<<std::endl;
   indexes<<"number="<<number<<"  // 1 for uniform array of circle holes; 2 for uniform array of hexagonal holes ;3 for uniform array of rectangular holes;"<<std::endl;
   indexes<<"// 4 for PGM circle holes; 5 for PGM hexagonal holes; 6 for PGM rectangular holes."<<std::endl;
   if(number==1)
   indexes<<"radius of the circluar hole = "<<radius1/mm<<"mm"<<std::endl;
   if(number==2)
   indexes<<"radius of the hexagonal hole = "<<radius1/mm<<"mm"<<std::endl;
   if(number==3)
   indexes<<"radius of the rectangular hole = "<<radius2/mm<<"mm"<<std::endl;
   if(number==4){
   indexes<<"radius1 of circle holes="<<radius1/mm<<"mm"<<std::endl;
   indexes<<"radius2 of circle holes="<<radius2/mm<<"mm"<<std::endl;}
   if(number==5){
   indexes<<"radius1 of hexagonal holes="<<radius1/mm<<"mm"<<std::endl;
   indexes<<"radius2 of hexagonal holes="<<radius2/mm<<"mm"<<std::endl;}
   if(number==6){
   indexes<<"radius1 of rectangular holes="<<radius1/mm<<"mm"<<std::endl;
   indexes<<"radius2 of rectangular holes="<<radius2/mm<<"mm"<<std::endl;}
   if(gpsorgun==1){
   indexes<<"gun:"<<" "<<"phantomtype="<<phantomtype<<" // 1 for Hoffman phantom; 2 for cylinder phantom; 3 for one point; 4 for three points." <<std::endl;
   }
   if(gpsorgun==2){
   indexes<<"gps:"<<" "<<"phantomtype="<<phantomtype<<" // 1 for point source; 2 for cylinder source; 3 for Jaszczack source."<<std::endl;
   }

   indexes.close();

  return physiWorld;
}


void ExN01DetectorConstruction::InitialisationOfMaterial()
{
    // Creating elements :
    G4double z, a, density,fractionmass;
    G4String name, symbol;
    G4int nel,ncomponents,natoms;

    //Air
    G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
    G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
    G4Element* H = new G4Element("Hydrogen","H" , z= 1., a= 1.01*g/mole);

    Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
    Air->AddElement(N, 70*perCent);
    Air->AddElement(O, 30*perCent);

    //water
    H2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
    H2O->AddElement(H, natoms=2);
    H2O->AddElement(O, natoms=1);

    //tungsten
    // Change the density of the tungsten density here.
    // materialofcollimator=2;
     if(materialofcollimator==1)
    {  WP = new G4Material("Tungsten", z=74., a= 183.85*g/mole, density= 19.35*g/cm3);}
     if(materialofcollimator==2)
    {  WP = new G4Material("Tungsten", z=74., a= 183.85*g/mole, density= 17.0*g/cm3);}
     if(materialofcollimator==3)
    {//lead
       WP = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.34*g/cm3);
    }

    //Lead to shied
    Pb =  new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.34*g/cm3);
    //G4Material* Si =  new G4Material("Silicon", z=14., a= 28.09*g/mole, density= 2.33*g/cm3);

    //Cd,Zn,Te
    G4Element* Cd = new G4Element("Ge","Cd", z=48., a= 112.4*g/mole);
    G4Element* Zn = new G4Element("Zin","Zn", z=30., a= 65.39*g/mole);
    G4Element* Te = new G4Element("Tie","Te", z=52., a= 127.60*g/mole);

    //CdZnTe
    CdZnTe = new G4Material("CdZnTe", density= 5.84*g/cm3, ncomponents=3);
    CdZnTe->AddElement(Cd, fractionmass=45*perCent);
    CdZnTe->AddElement(Zn, fractionmass=5*perCent);
    CdZnTe->AddElement(Te, fractionmass=50*perCent);

    //Na,I
    G4Element* I = new G4Element("Indine", "I", z=53., a=126.91*g/mole);
    G4Element* Na = new G4Element("Sodium", "Na", z=11., a = 22.99*g/mole);

    //NaI
    NaI = new G4Material("NaI", density= 3.67*g/cm3, ncomponents=2);
    NaI->AddElement(Na, natoms=1);
    NaI->AddElement(I, natoms=1);

   // materialofdetector=1;
    if(materialofdetector==1)
        detectormaterial=CdZnTe;
    if(materialofdetector==2)
        detectormaterial=NaI;

    //Cu
    Cu=new G4Material("Cu", z=29., a=63.546*g/mole,density=8.92*g/cm3);

    // print the materical table
    G4cout<<"************The materical are listed:***********"<<G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    std::cout<<"The number of the materials is :"<<fOriginalMaterials.size()<<std::endl;

   // another way to define elements and materials
    G4Element* elC = new G4Element( name = "Carbon",
                                    symbol = "C",
                                    z = 6.0, a = 12.011 * g/mole );
    G4Element* elH = new G4Element( name = "Hydrogen",
                                    symbol = "H",
                                    z = 1.0, a = 1.008  * g/mole );
    G4Element* elN = new G4Element( name = "Nitrogen",
                                    symbol = "N",
                                    z = 7.0, a = 14.007 * g/mole );
    G4Element* elO = new G4Element( name = "Oxygen",
                                    symbol = "O",
                                    z = 8.0, a = 16.00  * g/mole );
    G4Element* elNa = new G4Element( name = "Sodium",
                                     symbol = "Na",
                                     z= 11.0, a = 22.98977* g/mole );
    G4Element* elS = new G4Element( name = "Sulfur",
                                    symbol = "S",
                                    z = 16.0,a = 32.065* g/mole );
    G4Element* elCl = new G4Element( name = "Chlorine",
                                     symbol = "P",
                                     z = 17.0, a = 35.453* g/mole );
    G4Element* elK = new G4Element( name = "Potassium",
                                    symbol = "P",
                                    z = 19.0, a = 30.0983* g/mole );
    G4Element* elP = new G4Element( name = "Phosphorus",
                                    symbol = "P",
                                    z = 30.0, a = 30.973976* g/mole );
    G4Element* elFe = new G4Element( name = "Iron",
                                     symbol = "Fe",
                                     z = 26, a = 56.845* g/mole );
    G4Element* elMg = new G4Element( name = "Magnesium",
                                     symbol = "Mg",
                                     z = 12.0, a = 24.3050* g/mole );
    G4Element* elCa = new G4Element( name="Calcium",
                                     symbol = "Ca",
                                     z = 20.0, a = 40.078* g/mole );
    // Creating Materials :
    G4int numberofElements;

    //  Lung Inhale
    lunginhale = new G4Material( "LungInhale",
                                 density = 0.217*g/cm3,
                                 numberofElements = 9);
    lunginhale->AddElement(elH,0.103);
    lunginhale->AddElement(elC,0.105);
    lunginhale->AddElement(elN,0.031);
    lunginhale->AddElement(elO,0.749);
    lunginhale->AddElement(elNa,0.002);
    lunginhale->AddElement(elP,0.002);
    lunginhale->AddElement(elS,0.003);
    lunginhale->AddElement(elCl,0.002);
    lunginhale->AddElement(elK,0.003);

    // Lung exhale
    lungexhale = new G4Material( "LungExhale",
                                 density = 0.508*g/cm3,
                                 numberofElements = 9 );
    lungexhale->AddElement(elH,0.103);
    lungexhale->AddElement(elC,0.105);
    lungexhale->AddElement(elN,0.031);
    lungexhale->AddElement(elO,0.749);
    lungexhale->AddElement(elNa,0.002);
    lungexhale->AddElement(elP,0.002);
    lungexhale->AddElement(elS,0.003);
    lungexhale->AddElement(elCl,0.002);
    lungexhale->AddElement(elK,0.003);

    // Adipose tissue
    adiposeTissue = new G4Material( "AdiposeTissue",
                                    density = 0.967*g/cm3,
                                    numberofElements = 7);
    adiposeTissue->AddElement(elH,0.114);
    adiposeTissue->AddElement(elC,0.598);
    adiposeTissue->AddElement(elN,0.007);
    adiposeTissue->AddElement(elO,0.278);
    adiposeTissue->AddElement(elNa,0.001);
    adiposeTissue->AddElement(elS,0.001);
    adiposeTissue->AddElement(elCl,0.001);

    // Breast
    breast = new G4Material( "Breast",
                             density = 0.990*g/cm3,
                             numberofElements = 8 );
    breast->AddElement(elH,0.109);
    breast->AddElement(elC,0.506);
    breast->AddElement(elN,0.023);
    breast->AddElement(elO,0.358);
    breast->AddElement(elNa,0.001);
    breast->AddElement(elP,0.001);
    breast->AddElement(elS,0.001);
    breast->AddElement(elCl,0.001);

     // Water
    water = new G4Material( "Water",
                              density = 1.0*g/cm3,
                              numberofElements = 2 );
    water->AddElement(elH,0.112);
    water->AddElement(elO,0.888);

    // Muscle
    muscle = new G4Material( "Muscle",
                             density = 1.061*g/cm3,
                             numberofElements = 9 );
    muscle->AddElement(elH,0.102);
    muscle->AddElement(elC,0.143);
    muscle->AddElement(elN,0.034);
    muscle->AddElement(elO,0.710);
    muscle->AddElement(elNa,0.001);
    muscle->AddElement(elP,0.002);
    muscle->AddElement(elS,0.003);
    muscle->AddElement(elCl,0.001);
    muscle->AddElement(elK,0.004);

    // Liver
    liver = new G4Material( "Liver",
                            density = 1.071*g/cm3,
                            numberofElements = 9);
    liver->AddElement(elH,0.102);
    liver->AddElement(elC,0.139);
    liver->AddElement(elN,0.030);
    liver->AddElement(elO,0.716);
    liver->AddElement(elNa,0.002);
    liver->AddElement(elP,0.003);
    liver->AddElement(elS,0.003);
    liver->AddElement(elCl,0.002);
    liver->AddElement(elK,0.003);

    // Trabecular Bone
    trabecularBone = new G4Material( "TrabecularBone",
                                     density = 1.159*g/cm3,
                                     numberofElements = 12 );
    trabecularBone->AddElement(elH,0.085);
    trabecularBone->AddElement(elC,0.404);
    trabecularBone->AddElement(elN,0.058);
    trabecularBone->AddElement(elO,0.367);
    trabecularBone->AddElement(elNa,0.001);
    trabecularBone->AddElement(elMg,0.001);
    trabecularBone->AddElement(elP,0.034);
    trabecularBone->AddElement(elS,0.002);
    trabecularBone->AddElement(elCl,0.002);
    trabecularBone->AddElement(elK,0.001);
    trabecularBone->AddElement(elCa,0.044);
    trabecularBone->AddElement(elFe,0.001);

    // Dense Bone
    denseBone = new G4Material( "DenseBone",
                                density = 1.575*g/cm3,
                                numberofElements = 11 );
    denseBone->AddElement(elH,0.056);
    denseBone->AddElement(elC,0.235);
    denseBone->AddElement(elN,0.050);
    denseBone->AddElement(elO,0.434);
    denseBone->AddElement(elNa,0.001);
    denseBone->AddElement(elMg,0.001);
    denseBone->AddElement(elP,0.072);
    denseBone->AddElement(elS,0.003);
    denseBone->AddElement(elCl,0.001);
    denseBone->AddElement(elK,0.001);
    denseBone->AddElement(elCa,0.146);

    //soft tissue
    G4double A;   G4int Z;
    A = 65.38*g/mole;
    G4Element* elZn = new G4Element("Zinc","Zn",Z = 30.,A);

    A = 85.47 *g/mole;
    G4Element* elRb = new G4Element("Rb","Rb",Z = 37.,A);

    A = 87.62 *g/mole;
    G4Element* elSr = new G4Element("Sr","Sr",Z = 38.,A);

    A = 91.22 *g/mole;
    G4Element* elZr = new G4Element("Zr","Zr",Z = 40.,A);

    A = 207.19 *g/mole;
    G4Element* elPb = new G4Element("Lead","Pb", Z = 82.,A);

    G4double d = 0.9869 *g/cm3;
    soft = new G4Material("soft_tissue",d,16);
    soft->AddElement(elH,0.1047);
    soft->AddElement(elC,0.2302);
    soft->AddElement(elN,0.0234);
    soft->AddElement(elO,0.6321);
    soft->AddElement(elNa,0.0013);
    soft->AddElement(elMg,0.00015);
    soft->AddElement(elP,0.0024);
    soft->AddElement(elS,0.0022);
    soft->AddElement(elCl,0.0014);
    soft->AddElement(elK,0.0021);
    soft->AddElement(elFe,0.000063);
    soft->AddElement(elZn,0.000032);
    soft->AddElement(elRb,0.0000057);
    soft->AddElement(elSr,0.00000034);
    soft->AddElement(elZr,0.000008);
    soft->AddElement(elPb,0.00000016);

    //----- Put the materials in a vector
    fOriginalMaterials.push_back(Air); // rho = 0.00129
//  fOriginalMaterials.push_back(lunginhale); // rho = 0.217
//  fOriginalMaterials.push_back(lungexhale); // rho = 0.508
//  fOriginalMaterials.push_back(adiposeTissue); // rho = 0.967
    fOriginalMaterials.push_back(soft);  //rho = 0.9869
//  fOriginalMaterials.push_back(breast ); // rho = 0.990
    fOriginalMaterials.push_back(water); // rho = 1.018
    fOriginalMaterials.push_back(muscle); // rho = 1.061
    fOriginalMaterials.push_back(liver); // rho = 1.071
    fOriginalMaterials.push_back(trabecularBone); // rho = 1.159
    fOriginalMaterials.push_back(denseBone); // rho = 1.575
}

void ExN01DetectorConstruction::ReadPhantomData()
{
    std::ifstream finDF("DICOM/Data.dat");
    saveactivity.open("DICOM/activity.txt",ios::out);
    G4String fname;
    if(finDF.good() != 1 ) {
      G4Exception(" DicomDetectorConstruction::ReadPhantomData",
                  "",
                  FatalException,
                  "Problem reading data file: Data.dat");
    }

    G4int compression;
    finDF >> compression; // not used here

    finDF >> fNoFiles;
    saveactivity<<fNoFiles<<endl;
    for(G4int i = 0; i < fNoFiles; i++ ) {
      finDF >> fname;
      saveactivity<<fname<<endl;
      //--- Read one data file
      fname += ".g4dcm";
      ReadPhantomDataFile(fname);
    }

    //----- Merge data headers
    MergeZSliceHeaders();

    finDF.close();

}
void ExN01DetectorConstruction::ReadPhantomDataFile(const G4String& fname)
{
#ifdef G4VERBOSE
  G4cout << " DicomDetectorConstruction::ReadPhantomDataFile opening file " << fname << G4endl;
#endif
  std::ifstream fin(fname.c_str(), std::ios_base::in);
  if( !fin.is_open() ) {
    G4Exception("DicomDetectorConstruction::ReadPhantomDataFile",
                "",
                FatalErrorInArgument,
                G4String("File not found " + fname ).c_str());
  }
  //----- Define density differences (maximum density difference to create a new material)
  char* part = getenv( "DICOM_CHANGE_MATERIAL_DENSITY" );

  G4double densityDiff = -1.;

  if( part ) densityDiff = G4UIcommand::ConvertToDouble(part);
  if( densityDiff != -1. ) {
    for( unsigned int ii = 0; ii < fOriginalMaterials.size(); ii++ ){
      fDensityDiffs[ii] = densityDiff; //currently all materials with same difference
    }
  }else {
    if( fMaterials.size() == 0 ) { // do it only for first slice
      for( unsigned int ii = 0; ii < fOriginalMaterials.size(); ii++ ){
        fMaterials.push_back( fOriginalMaterials[ii] );
      }
    }

  }
  //----- Read data header
  ExN01DicomPhantomZSliceHeader* sliceHeader = new  ExN01DicomPhantomZSliceHeader( fin );
  fZSliceHeaders.push_back( sliceHeader );

  //----- Read material indices
  G4int nVoxels = sliceHeader->GetNoVoxels(); //64*64=4096


  //--- If first slice, initialize fMateIDs
  if( fZSliceHeaders.size() == 1 ) {
    //fMateIDs = new unsigned int[fNoFiles*nVoxels];
    fMateIDs = new size_t[fNoFiles*nVoxels];
  }

  unsigned int mateID;
  G4int voxelCopyNo = (fZSliceHeaders.size()-1)*nVoxels; // number of voxels from previously read slices
  for( G4int ii = 0; ii < nVoxels; ii++, voxelCopyNo++ ){
    fin >> mateID;
    fMateIDs[voxelCopyNo] = mateID;
  }

  //----- Read material densities and build new materials if two voxels have same material but its density is in a
  //different density interval (size of density intervals defined by densityDiff)
  G4double density;
  voxelCopyNo = (fZSliceHeaders.size()-1)*nVoxels; // number of voxels from previously read slices
  for( G4int ii = 0; ii < nVoxels; ii++, voxelCopyNo++ ){
    fin >> density;

    //-- Get material from list of original materials
    mateID = fMateIDs[voxelCopyNo];

    G4Material* mateOrig  = fOriginalMaterials[mateID];

    //-- Get density bin: middle point of the bin in which the current density is included
    G4String newMateName = mateOrig->GetName();
    float densityBin = 0.;
   if( densityDiff != -1.) {
     densityBin = fDensityDiffs[mateID] * (G4int(density/fDensityDiffs[mateID])+0.5);
     //-- Build the new material name
      newMateName += G4UIcommand::ConvertToString(densityBin);
      std::cout<<newMateName<<"**"<<std::endl;
      sleep(2);
    }

    //-- Look if a material with this name is already created (because a previous voxel was already in this density bin)
    unsigned int im;
    for( im = 0; im < fMaterials.size(); im++ ){
      if( fMaterials[im]->GetName() == newMateName ) {
        break;
      }
    }
    //-- If material is already created use index of this material
    if( im != fMaterials.size() ) {
      fMateIDs[voxelCopyNo] = im;
    //-- else, create the material
    } else {
      if( densityDiff != -1.) {
        fMaterials.push_back( BuildMaterialWithChangingDensity( mateOrig, densityBin, newMateName ) );
        fMateIDs[voxelCopyNo] = fMaterials.size()-1;
      } else {
        G4cerr << " im " << im << " < " << fMaterials.size() << " name " << newMateName << G4endl;
        G4Exception("DicomDetectorConstruction::ReadPhantomDataFile",
                    "",
                    FatalErrorInArgument,
                    "Wrong index in material"); //it should never reach here
      }
    }
  }
}


void ExN01DetectorConstruction::MergeZSliceHeaders()
{

    //----- Images must have the same dimension ...
    fZSliceHeaderMerged = new ExN01DicomPhantomZSliceHeader( *fZSliceHeaders[0] );
    for( unsigned int ii = 1; ii < fZSliceHeaders.size(); ii++ ) {
      *fZSliceHeaderMerged += *fZSliceHeaders[ii];
    };
}

G4Material* ExN01DetectorConstruction::BuildMaterialWithChangingDensity( const G4Material* origMate, float density, G4String newMateName )
{
  //----- Copy original material, but with new density
  G4int nelem = origMate->GetNumberOfElements();
  G4Material* mate = new G4Material( newMateName, density*g/cm3, nelem, kStateUndefined, STP_Temperature );

  for( G4int ii = 0; ii < nelem; ii++ ){
    G4double frac = origMate->GetFractionVector()[ii];
    G4Element* elem = const_cast<G4Element*>(origMate->GetElement(ii));
    mate->AddElement( elem, frac );
  }
  std::cout<<density<<std::endl;
  return mate;
}

void ExN01DetectorConstruction::ConstructPhantomContainer()
{
  //---- Extract number of voxels and voxel dimensions
  fNVoxelX = fZSliceHeaderMerged->GetNoVoxelX();
  fNVoxelY = fZSliceHeaderMerged->GetNoVoxelY();
  fNVoxelZ = fZSliceHeaderMerged->GetNoVoxelZ();
  saveactivity<<fNVoxelX<<" "<<fNVoxelY<<endl;

  fVoxelHalfDimX = fZSliceHeaderMerged->GetVoxelHalfX()/mm;
  fVoxelHalfDimY = fZSliceHeaderMerged->GetVoxelHalfY()/mm;
  fVoxelHalfDimZ = fZSliceHeaderMerged->GetVoxelHalfZ()/mm;
  saveactivity<<fVoxelHalfDimX<<" "<<fVoxelHalfDimY<<" "<<fVoxelHalfDimZ<<endl;

#ifdef G4VERBOSE
  G4cout << " fNVoxelX " << fNVoxelX << " fVoxelHalfDimX " << fVoxelHalfDimX/cm <<G4endl;
  G4cout << " fNVoxelY " << fNVoxelY << " fVoxelHalfDimY " << fVoxelHalfDimY/cm <<G4endl;
  G4cout << " fNVoxelZ " << fNVoxelZ << " fVoxelHalfDimZ " << fVoxelHalfDimZ/cm <<G4endl;
  G4cout << " totalPixels " << fNVoxelX*fNVoxelY*fNVoxelZ <<  G4endl;
#endif

  //----- Define the volume that contains all the voxels
  fContainer_solid = new G4Box("phantomContainer",fNVoxelX*fVoxelHalfDimX,fNVoxelY*fVoxelHalfDimY,fNVoxelZ*fVoxelHalfDimZ);
  fContainer_logic =
    new G4LogicalVolume( fContainer_solid,
                         fMaterials[0],  //material is not important, it is completely filled by voxels
                         "phantomContainer",
                         0, 0, 0 );
  //--- Place it in the world
  G4double fOffsetX = (fZSliceHeaderMerged->GetMaxX() + fZSliceHeaderMerged->GetMinX() ) /2.;
  G4double fOffsetY = (fZSliceHeaderMerged->GetMaxY() + fZSliceHeaderMerged->GetMinY() ) /2.;
  G4double fOffsetZ = (fZSliceHeaderMerged->GetMaxZ() + fZSliceHeaderMerged->GetMinZ() ) /2.;
  fOffsetZ=-(sourcetocollimator*10+Leng/cm*10);  //112.5=100+12.5
  saveactivity<<fOffsetZ<<endl;
  fOffsetZ=fOffsetZ/2;
  G4ThreeVector posCentreVoxels(fOffsetX,fOffsetY,fOffsetZ);


#ifdef G4VERBOSE
  G4cout << " placing voxel container volume at " << posCentreVoxels << G4endl;
#endif
  G4RotationMatrix* rot=new G4RotationMatrix();
  rot->rotateX(90*deg);

  fContainer_phys =
          new G4PVPlacement(rot,           // rotation
                      posCentreVoxels*2,
                      fContainer_logic,     // The logic volume
                      "phantomContainer",   // Name
                      logicWorld,           // Mother
                      fCheckOverlaps,       // No op. bool.
                      1);                   // Copy number

 saveactivity.close();
}

void ExN01DetectorConstruction::SetArmAngle(G4double val)
{

  fArmAngle = val;
  *fArmRotation = G4RotationMatrix();  // make it unit vector
  int choose;
  choose=1;
  if(choose==1)
  {
  fArmRotation->rotateY(fArmAngle);
  }
  if(choose==2)
  {
  fArmRotation->rotateY(fArmAngle);
  fArmRotation->rotateX(-90*deg);
 // fArmRotation->rotateX(-90*deg);
  }
  G4double xx = (sourcetocollimator+Leng/cm)*cm*std::sin(-fArmAngle);
  G4double zz = (sourcetocollimator+Leng/cm)*cm*(std::cos(-fArmAngle)-1);
 // physiEllip->SetTranslation(G4ThreeVector(xx,yy,-14.5*cm));
 // fContainer_phys->SetRotation(fArmRotation);

 //  rotate the collimator
  physiCollimator->SetRotation(fArmRotation);
  physiCollimator->SetTranslation(G4ThreeVector(xx,0,zz));

 // rotate the detector
  xx=(sourcetocollimator+Leng/cm+length/cm+thick/cm)*cm*std::sin(-fArmAngle);
  zz=(sourcetocollimator+Leng/cm+length/cm+thick/cm)*cm*std::cos(-fArmAngle)-(sourcetocollimator+Leng/cm)*cm;
  physiDetector->SetRotation(fArmRotation);
  physiDetector->SetTranslation(G4ThreeVector(xx,0,zz));

  xx=(sourcetocollimator+Leng/cm+length/cm+2*thick/cm+2)*cm*std::sin(-fArmAngle);
  zz=(sourcetocollimator+Leng/cm+length/cm+2*thick/cm+2)*cm*std::cos(-fArmAngle)-(sourcetocollimator+Leng/cm)*cm;
  BOXphy->SetTranslation(G4ThreeVector(xx,0,zz));
  BOXphy->SetRotation(fArmRotation);

  xx=(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*cm*std::sin(-fArmAngle);
  zz=(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*cm*cos(-fArmAngle)-(sourcetocollimator+length/cm)*cm;
  sidephy->SetTranslation(G4ThreeVector(xx,12.24*cm,zz));
  sidephy->SetRotation(fArmRotation);
  sidephy2->SetTranslation(G4ThreeVector(xx,-12.24*cm,zz));
  sidephy2->SetRotation(fArmRotation);

  xx=std::sqrt(12.24*12.24+(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))*cm*std::sin(-fArmAngle-std::atan(12.24/(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)));
  zz=std::sqrt(12.24*12.24+(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))*cm*std::cos(-fArmAngle-std::atan(12.24/(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)))-(sourcetocollimator+Leng/cm)*cm;

  leftphy->SetTranslation(G4ThreeVector(xx,0*cm,zz));
  leftphy->SetRotation(fArmRotation);

  xx=std::sqrt(12.24*12.24+(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))*cm*std::sin(-fArmAngle+std::atan(12.24/(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)));
  zz=std::sqrt(12.24*12.24+(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))*cm*std::cos(-fArmAngle+std::atan(12.24/(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)))-(sourcetocollimator+Leng/cm)*cm;

  leftphy2->SetTranslation(G4ThreeVector(xx,0*cm,zz));
  leftphy2->SetRotation(fArmRotation);
//end of the first detector

  xx=(sourcetocollimator+Leng/cm+length/cm+thick/cm)*cm*std::sin(-fArmAngle-180*deg);
  zz=(sourcetocollimator+Leng/cm+length/cm+thick/cm)*cm*std::cos(-fArmAngle-180*deg)-(sourcetocollimator+Leng/cm)*cm;
  physiDetectortwo->SetRotation(fArmRotation);
  physiDetectortwo->SetTranslation(G4ThreeVector(xx,0,zz));

  xx=(sourcetocollimator+Leng/cm)*cm*std::sin(-fArmAngle-180*deg);
  zz=(sourcetocollimator+Leng/cm)*cm*std::cos(-fArmAngle-180*deg)-(sourcetocollimator+Leng/cm)*cm;
  physiCollimatortwo->SetRotation(fArmRotation);
  physiCollimatortwo->SetTranslation(G4ThreeVector(xx,0,zz));

  xx=(sourcetocollimator+Leng/cm+length/cm+2*thick/cm+2)*cm*std::sin(-fArmAngle-180*deg);
  zz=(sourcetocollimator+Leng/cm+length/cm+2*thick/cm+2)*cm*std::cos(-fArmAngle-180*deg)-(sourcetocollimator+Leng/cm)*cm;
  BOXphytwo->SetRotation(fArmRotation);
  BOXphytwo->SetTranslation(G4ThreeVector(xx,0,zz));

  xx=(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*cm*std::sin(-fArmAngle-180*deg);
  zz=(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*cm*cos(-fArmAngle-180*deg)-(sourcetocollimator+Leng/cm)*cm;
  sidephytwo->SetTranslation(G4ThreeVector(xx,12.24*cm,zz));
  sidephytwo->SetRotation(fArmRotation);
  sidephy2two->SetTranslation(G4ThreeVector(xx,-12.24*cm,zz));
  sidephy2two->SetRotation(fArmRotation);

  xx=std::sqrt(12.24*12.24+(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))*cm*std::sin(-fArmAngle-std::atan(12.24/(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))-180*deg);
  zz=std::sqrt(12.24*12.24+(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))*cm*std::cos(-fArmAngle-std::atan(12.24/(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))-180*deg)-(sourcetocollimator+Leng/cm)*cm;

  leftphytwo->SetTranslation(G4ThreeVector(xx,0*cm,zz));
  leftphytwo->SetRotation(fArmRotation);

  xx=std::sqrt(12.24*12.24+(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))*cm*std::sin(-fArmAngle+std::atan(12.24/(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))-180*deg);
  zz=std::sqrt(12.24*12.24+(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm)*(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))*cm*std::cos(-fArmAngle+std::atan(12.24/(sourcetocollimator+Leng/cm+thick/cm+(gap/2)/cm))-180*deg)-(sourcetocollimator+Leng/cm)*cm;

  leftphy2two->SetTranslation(G4ThreeVector(xx,0*cm,zz));
  leftphy2two->SetRotation(fArmRotation);

  // tell G4Manager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void ExN01DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{  fCheckOverlaps = checkOverlaps;}

void ExN01DetectorConstruction::ConstructPhantombyself()
{
    G4Box*             SolidPhantom;    // pointer to the solid envelope
    G4Box*             SolidPhantom2;
    G4LogicalVolume*   LogicPhantom0;    // pointer to the logical envelope
    G4LogicalVolume*   LogicPhantom1;
    G4LogicalVolume*   LogicPhantom3;
    G4LogicalVolume*   LogicPhantom4;
    G4VPhysicalVolume* PhysiPhantom;    // pointer to the physical envelope
    G4int phantomlength,phantomwidth,phantomheight;
    G4float xsize,ysize,zsize,offset,posv,posh;
    G4int num1,num2,num3,num4;
    char filename[128]={0};
    char pathandname[256]={0};
    char buf[80];char bufup[80];
    getcwd(buf, sizeof(buf));
    chdir("..");
    getcwd(bufup,sizeof(bufup));
    chdir(buf);
    printf("%s\n",bufup);

    G4double a, z; G4int i,j,k,type;
    G4double density;
    G4int nel,ncomponents,natoms;
    float aa=0.0; int bb=0;
    num1=0;num2=0;num3=0;num4=0;
    G4int vh; // vh=1,vertical profile;vh=2;horizonal profile;



    G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
    G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
    G4Element* H = new G4Element("Hydrogen","H" , z= 1., a= 1.01*g/mole);

    G4Material* H2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
    H2O->AddElement(H, natoms=2);
    H2O->AddElement(O, natoms=1);

    G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
    Air->AddElement(N, 70*perCent);
    Air->AddElement(O, 30*perCent);

    Materialsbysef.push_back(Air); // rho = 0.00129
    Materialsbysef.push_back(H2O);  //rho = 0.9869

    G4Color
    green(0.0,1.0,0.0),
    blue(0.0,0.0,1.0);

    sprintf(filename,"parameter.txt");
    sprintf(pathandname,"%s/phantom/%d/%s",bufup,phantomtype,filename);
    parameter.open(pathandname);
    parameter>>phantomlength>>phantomwidth>>phantomheight;
    parameter>>xsize>>ysize>>zsize;
    parameter>>type;
    parameter>>vh;
    parameter.close();

    SolidPhantom=new G4Box("SolidPhantom",xsize *cm,ysize *cm ,zsize *cm);
    SolidPhantom2=new G4Box("SolidPhantom",xsize *cm,zsize *cm ,ysize *cm);

    LogicPhantom0=new G4LogicalVolume(SolidPhantom,Air,"LogicPhantom",0,0,0);
    LogicPhantom1=new G4LogicalVolume(SolidPhantom,H2O,"LogicPhantom",0,0,0);
    LogicPhantom3=new G4LogicalVolume(SolidPhantom2,Air,"LogicPhantom",0,0,0);
    LogicPhantom4=new G4LogicalVolume(SolidPhantom2,H2O,"LogicPhantom",0,0,0);


    offset=-(Leng/cm+sourcetocollimator)*cm;// 10cm+1.25cm
    posv = (offset/cm-(phantomheight-1)*zsize)*cm;
    posh = -(phantomheight-1)*zsize*cm;

    for(i=0;i<phantomheight;i++)
     {
        if(i<10)
        //sprintf(filename,"attenuation0%d.txt",i);
          sprintf(filename,"slice0%d.txt",i);
        else
        //sprintf(filename,"attenuation%d.txt",i);
        sprintf(filename,"slice%d.txt",i);
        sprintf(pathandname,"%s/phantom/%d/attenuation/%s",bufup,phantomtype,filename);

        std::cout<<pathandname<<std::endl;
        sleep(2);

        parameter.open(pathandname);
        for(j=0;j<phantomwidth;j++)
       {
            for(k=0;k<phantomlength;k++)
            {
                if(type==1)
                parameter>>bb;
                if(type==2) {
                    parameter>>aa;
                    if(aa>0) bb=1;
                    else     bb=0;}

                if(bb>1)     bb=1;

                //if(bb>Materialsbysef.size())
                // bb=Materialsbysef.size();

         if(vh==1)
          {
           if(bb==0)
            { PhysiPhantom=new G4PVPlacement(0,
                              G4ThreeVector((k-phantomlength/2+0.5)*2*xsize*cm,(j-phantomwidth/2+0.5)*2*ysize*cm,(posv/cm+i*2*zsize)*cm),
                              LogicPhantom0,
                              "PhysiPhantom",
                              logicWorld,
                              fCheckOverlaps,
                              bb);
              num1++;}

           if(bb==1)
            { PhysiPhantom=new G4PVPlacement(0,
                                G4ThreeVector((k-phantomlength/2+0.5)*2*xsize*cm,(j-phantomwidth/2+0.5)*2*ysize*cm,(posv/cm+i*2*zsize)*cm),
                                LogicPhantom1,
                                "PhysiPhantom",
                                logicWorld,
                                fCheckOverlaps,
                                bb);
               num2++;
             }
           }
          if(vh==2)
         {
            if(bb==0)
            { PhysiPhantom=new G4PVPlacement(0,
                               G4ThreeVector((k-phantomlength/2+0.5)*2*xsize*cm,(posh/cm+i*2*zsize)*cm,((j-phantomwidth/2+0.5)*2*ysize+offset/cm)*cm),
                               LogicPhantom3,
                               "PhysiPhantom",
                               logicWorld,
                               fCheckOverlaps,
                               bb);
                num3++;
              }
             if(bb==1)
             { PhysiPhantom=new G4PVPlacement(0,
                               G4ThreeVector((k-phantomlength/2+0.5)*2*xsize*cm,(posh/cm+i*2*zsize)*cm,((j-phantomwidth/2+0.5)*2*ysize+offset/cm)*cm),
                               LogicPhantom4,
                               "PhysiPhantom",
                               logicWorld,
                               fCheckOverlaps,
                               bb);
                num4++;
              }

            }
        }
      }
       parameter.close();
    }
    G4VisAttributes* vis= new G4VisAttributes(green);//
    vis->SetVisibility(true);
    vis->SetForceSolid(false);//
    LogicPhantom0->SetVisAttributes(vis);  // the first detector
    LogicPhantom3->SetVisAttributes(vis);

    G4VisAttributes* vis2= new G4VisAttributes(blue);//
    vis2->SetVisibility(true);
    vis2->SetForceSolid(false);//
    LogicPhantom1->SetVisAttributes(vis2);  // the first detector
    LogicPhantom4->SetVisAttributes(vis2);

    std::cout<<num1<<" "<<num2<<" "<<num3<<" "<<num4<<std::endl;
    // Testing Head Volume
    G4double HeadVol = LogicPhantom4->GetSolid()->GetCubicVolume();
    std::cout << "Volume of Phantom4 = " << HeadVol/cm3*num4 << " cm^3" << std::endl;

    // Testing Head Material
    G4String HeadMat = LogicPhantom4->GetMaterial()->GetName();
    std::cout << "Material of Phantom4 = " << HeadMat << std::endl;

    // Testing Density
    G4double HeadDensity = LogicPhantom4->GetMaterial()->GetDensity();
    std::cout << "Density of Phantom4 = " << HeadDensity*cm3/g << " g/cm^3" << std::endl;

    // Testing Mass
    G4double HeadMass = (HeadVol)*HeadDensity;
    std::cout << "Mass of Phantom4 = " << HeadMass/gram*num4 << " g" << std::endl;
    sleep(2);
}


