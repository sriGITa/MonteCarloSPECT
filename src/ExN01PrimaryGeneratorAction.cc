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
#include "ExN01PrimaryGeneratorAction.hh"
#include "ExN01DetectorConstruction.hh"
#include "ExN01RunAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4GeneralParticleSource.hh"
#include "Randomize.hh"
#include "G4Geantino.hh"
#include "G4ChargedGeantino.hh"
#include "G4SPSPosDistribution.hh"
#include "time.h"
#include <unistd.h>  // get the current folder

//the class constructor PrimaryGeneratorAction::PrimaryGeneratorAction() and the method PrimaryGeneratorAction::GeneratePrimaries(G4Event* )

ExN01PrimaryGeneratorAction::ExN01PrimaryGeneratorAction()
{
  // get the current directory
  getcwd(buf, sizeof(buf));
  chdir("..");
  // get the parent directory
  getcwd(bufup,sizeof(bufup));
  chdir(buf);
  sprintf(filename,"collimator_detector.txt");
  sprintf(pathandname,"%s/run/%s",bufup,filename);
  test = fopen(pathandname,"r");
  fscanf(test,"gpsorgun=%d\n",&gpsorgun);
  fclose(test);
  //gpsorgun=1;//(1:generalParticleGun; 2:generalParticleSource)

  if(gpsorgun==1)
  {
    para.open("index.txt",ios::out);
    activity_type=3;   // activity_type=1: point source
                       // activity_type=2: import activity profiles with dicom format
                       // activity_type=3: phantom activity profiles without dicom format
   if(activity_type==2)
    {
    sum=0;
    ord=0;ii=0;
    x=0;y=0;z=0;

    numFile.open("DICOM/activity.txt");
    numFile>>num;
    fname=new G4String[num];
    for(int i=0;i<num;i++)
    {
        numFile>>fname[i];
    }
    numFile >> rows;
    numFile >> columns;

    index = new G4double**[rows];
    for ( G4int i = 0; i < rows; i ++ ) {
      index[i] = new G4double*[columns];
      for(G4int j=0;j<columns;j++)
      {index[i][j]=new G4double[num];}
    }
   chdir("DICOM/");
    for(G4int i = 0; i < num; i++ ) {
        fname[i] += ".activity.txt";
        ACTIVITY.open(fname[i]);
        for(int j=0;j<rows;j++)
          for(int k=0;k<columns;k++)
          {
              ACTIVITY>>index[j][k][i];
              if(index[j][k][i]>100)
             { std::cout<<(int)index[j][k][i]<<" "<<index[j][k][i]<<std::endl;
              sleep(2);
              }
              sum=sum+index[j][k][i];
          }
      ACTIVITY.close();
     }
    chdir("..");
    std::cout<<sum<<std::endl;
    sleep(2);
    subsum=index[0][0][0];
   }

   if(activity_type==3)
   {
       sprintf(filename,"parameter.txt");
       sprintf(pathandname,"%s/phantom/%d/%s",bufup,phantomtype,filename);
       parameterofphantom.open(pathandname);
       parameterofphantom>>phantomlength>>phantomwidth>>phantomheight;
       parameterofphantom>>xsize>>ysize>>zsize;
       parameterofphantom>>type;
       parameterofphantom>>vh;
       parameterofphantom.close();

       sum=0;
       offset=-off*cm;// 10cm+1.0cm
       posv = (offset/cm-(phantomheight-1)*zsize)*cm;
       posh = -(phantomheight-1)*zsize*cm;
       for(II=0;II<phantomheight;II++)
        {
           if(II<10)
           //sprintf(filename,"activity0%d.txt",II);
            sprintf(filename,"slice0%d.txt",II);
           else
           //sprintf(filename,"activity%d.txt",II);
           sprintf(filename,"slice%d.txt",II);
           sprintf(pathandname,"%s/phantom/%d/activity/%s",bufup,phantomtype,filename);
           parameterofphantom.open(pathandname);
            for(JJ=0;JJ<phantomwidth;JJ++)
          {
               for(KK=0;KK<phantomlength;KK++)
               {
                   if(type==1)
                   parameterofphantom>>bb;
                   if(type==2)
                   { parameterofphantom>>aa;
                       bb=(int)aa;}
                   if(bb==0)
                       bb=0;
                   else
                       bb=bb*1;
                    sum=sum+bb;

             if(vh==1)
              {
                 xx=(KK-phantomlength/2+0.5)*2*xsize*cm;
                 yy=(JJ-phantomwidth/2+0.5)*2*ysize*cm;
                 zz=(posv/cm+II*2*zsize)*cm;
              }
             if(vh==2)
              {
                 xx=(KK-phantomlength/2+0.5)*2*xsize*cm;
                 yy=(posh/cm+II*2*zsize)*cm;
                 zz=((JJ-phantomwidth/2+0.5)*2*ysize+offset/cm)*cm;

             }
             if(bb>0)
             {
             for(int ii=0;ii<bb;ii++)
             { X.push_back(xx);
               Y.push_back(yy);
               Z.push_back(zz);
             }

           }
          }
         }
          parameterofphantom.close();
       }
       std::cout<<sum<<std::endl;
       sleep(2);
   }

    G4int n_particle=1;
    generalParticleGun =new G4ParticleGun(n_particle);
    generalParticleGun->SetParticleEnergy(0*eV);

  //  generalParticleGun->SetParticlePosition(G4ThreeVector(0.0745*cm,0.0745*cm,-14 *cm));
    generalParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));

    if(activity_type==1)
        para<<"pointsource"<<std::endl;
    if(activity_type==2)
    {   para<<"Activities of Hoffman phantom (dicom)"<<std::endl;
        para<<"total counts: "<<sum<<std::endl;
        for(int i=0;i<num;i++)
        para<<"profile: "<<fname[i]<<std::endl;
        para<<"rowsize= "<<rows<<"\n"<<"columnsize= "<<columns<<std::endl;
     }
    if(activity_type==3)
    {   para<<"User defined source activities"<<std::endl;
        para<<"total counts: "<<sum<<std::endl;
        para<<"pixelNum_x="<<phantomlength<<"\n"<<"pixelNum_y="<<phantomwidth<<"\n"<<"pixelNum_z="<<phantomheight<<std::endl;
        para<<"pixel_x="<<xsize<<"cm"<<std::endl;
        para<<"pixel_y="<<ysize<<"cm"<<std::endl;
        para<<"pixel_z="<<zsize<<"cm"<<std::endl;
        para<<"data type: "<<type<<"  //type=1:int type; type=2:float type"<<std::endl;
        para<<"oriention: "<<vh<<"  // vh=1: horizontal profile; vh=2: vertical profile "<<std::endl;
     }
  para.close();
  }

  // particle source
  if(gpsorgun==2){
    generalParticleSource=new G4GeneralParticleSource();
   // generalParticleSource->GetCurrentSource()->GetEneDist()->SetMonoEnergy(140*keV);
  // generalParticleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0.0745*cm,0.0745*cm,15 *cm));
   // generalParticleSource->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
   // generalParticleSource->SetParticlePosition(G4ThreeVector(0.0745*cm,0.0745*cm,-14 *cm));
    if(activity_type==1)
        para<<"pointsource"<<std::endl;
    if(activity_type==2)
    {   para<<"Activities of cylinder phantom"<<std::endl;}
    if(activity_type==3)
    {   para<<"Activities of Jaszczak source"<<std::endl; }
  para.close();
 }

}

ExN01PrimaryGeneratorAction::~ExN01PrimaryGeneratorAction()
{
  if(gpsorgun==1)
    {
      delete generalParticleGun;
     }
  if(gpsorgun==2)
    {
      delete generalParticleSource;
    }
    if(gpsorgun==1)
  {
  if(activity_type==2)
  {
    for(int i=0;i<num;i++)
    {
        for(int j=0;j<rows;j++)
        {
            delete []index[i][j];
            index[i][j]=NULL;
        }
        delete []index[i];
        index[i]=NULL;
      }
    delete []index;
    index=NULL;
    delete []fname;
     }
   }
}

void ExN01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

   G4ThreeVector v(0,0.0745*cm,-off*cm);

   G4int i=anEvent->GetEventID();
 // generalParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,0.1));
  if(gpsorgun==1)
  {
   if (generalParticleGun->GetParticleDefinition() == G4ChargedGeantino::ChargedGeantino())
    {
      G4int Z = 27, A = 57;
      G4double ionCharge   = 0.*eplus;
      G4double excitEnergy = 0 *keV;
    
      G4ParticleDefinition* ion
       = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);

  //  ion->SetPDGStable( false );
  //  G4double lifetime=1*s;
  //  ion->SetPDGLifeTime( lifetime );

      generalParticleGun->SetParticleDefinition(ion);
      generalParticleGun->SetParticleCharge(ionCharge);
   }

   // choose the source, activity_type=1,then is the point source.

   if(activity_type==1)
   {
       switch(i%3)
       {
       case 0:
         v.setX(0.0*cm);
               break;
       case 1:
         v.setX(0.6*cm);
               break;
       case 2:
         v.setX(1.5*cm);
               break;
        }
       generalParticleGun->SetParticlePosition(v);
       generalParticleGun->GeneratePrimaryVertex(anEvent);

   }

   if(activity_type==2)
  {
    if(i==0&&runorder==0)
     {
      numFile>>VoxelHalfDimX>>VoxelHalfDimY>>VoxelHalfDimZ;
      numFile>>position;
      position=position-(num-1)*VoxelHalfDimZ;
    }
  if(i==0)
    {
      ord=0;
    }

   //create vertex
   // if(ii<(100-(int)index[x][y][z]))
   if(ord>=rows*columns*num)
    {ord=0;}

    z=ord/(rows*columns);
    x=(ord-z*rows*columns)/rows;
    y=ord-z*rows*columns-x*rows;

    v.setZ(position+(x-rows/2)*VoxelHalfDimX*2);
    v.setX((y-columns/2)*VoxelHalfDimY*2);
    v.setY((z-(num-1)/2)*2*VoxelHalfDimZ);

   if(ii<(int)index[x][y][z])
    {
       ii=ii+1;
    }
   else
    {
       ii=0;
       ord=ord+1;
    }

   generalParticleGun->SetParticlePosition(v);
   generalParticleGun->GeneratePrimaryVertex(anEvent);

  }

   if(activity_type==3)
   {
    v.setX(X[i%sum]);
    v.setY(Y[i%sum]);
    v.setZ(Z[i%sum]);
    generalParticleGun->SetParticlePosition(v);
    generalParticleGun->GeneratePrimaryVertex(anEvent);

   }

  }

  if(gpsorgun==2)
  {
      generalParticleSource->GeneratePrimaryVertex(anEvent);
  }

 }


