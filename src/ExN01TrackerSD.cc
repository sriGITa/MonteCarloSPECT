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

#include "ExN01TrackerSD.hh"
#include "ExN01TrackerHit.hh"
#include "ExN01PhysicsList.hh"
#include "ExN01RunAction.hh"
#include "ExN01SteppingAction.hh"
#include "ExN01EventAction.hh"

#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Region.hh"
#include "G4Material.hh"
#include "G4VProcess.hh"
#include "G4ParticleTypes.hh"
//#include "G4Electron.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "math.h"
#include "Randomize.hh"
#include <fstream>
#include "G4Track.hh"
#include  <iostream>
#include "vector"
#include <unistd.h>

using namespace std;
int totalnum[100]={0};

#ifdef HAVE_ROOT
    #include <TROOT.h> 
    #include <TSystem.h> 
    #include <TApplication.h>  
    #include <TStyle.h>     
    #include <TGraph.h> 
    #include <TGraphErrors.h> 
    #include <TH1.h> 
    #include <TH2.h> 
    #include <TCanvas.h> 
    #include <TNtuple.h>  
    #include <TFile.h>  
    #include <TPad.h>  
    #include <TF1.h>  
    #include <TLegend.h> 
    #include <TPaveText.h> 
    #include <TRandom.h>     
    #include <TStopwatch.h> 
    #include <TGaxis.h>
    #include <TLatex.h>
    #include <TPaveStats.h> 
    #include <TError.h>   
#endif


ExN01TrackerSD::ExN01TrackerSD(G4String name):G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");

  getcwd(buf, sizeof(buf));
  chdir("..");
  getcwd(bufup,sizeof(bufup));
  chdir(buf);
  sprintf(filename,"nuclide.txt");
  sprintf(pathandname,"%s/run/%s",bufup,filename);
  if (!(fp = fopen(pathandname,"r")))
  {
      printf("Error opening nuclide file %s!\n", filename);
      fclose(fp);
      exit(1);
  }
  fscanf(fp,"%d\n",&nuclide_type);
  fclose(fp);

  if(nuclide_type>6||nuclide_type<1)
  nuclide_type=1;

 //save each step of the gammas in each event
  sprintf(filename,"data.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  datafile.open(pathandname,ios::out);//overlap
  datafile<<"particles"<<" "<<"x"<<" "<<"y"<<" "<<"z"<<" "<<"deposited_energy"<<" "<<"Copy_no."<<" "<<"total_energy"<<" "<<"Steps"<<endl;
  datafile.close();

 //save the counts in each pixel
  sprintf(filename,"saveNum.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  vector1.open(pathandname,ios::out);

 //save the deposited energy in each pixel without energy window
  sprintf(filename,"saveEnergy.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  Energy.open(pathandname,ios::out);

 //output the deposited energy in each pixel with 10% energy window
  sprintf(filename,"savepixelenergy.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  pixelenergy.open(pathandname,ios::out);

 //save the deposited energy in each pixel for different runs
  sprintf(filename,"differentrun.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  rotate.open(pathandname,ios::out);
  rotate.close();

 //save the counts in each pixel for different runs
  sprintf(filename,"differentcounts.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  countsforrun.open(pathandname,ios::out);
  countsforrun.close();

 //save the deposited energy in each pixel for different runs with 10% energy window
  sprintf(filename,"windowenergy.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  windowenergy.open(pathandname,ios::out);
  windowenergy.close();

 //save the counts in each pixel for different runs with 10% energy window
  sprintf(filename,"windowcounts.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  windowcounts.open(pathandname,ios::out);
  windowcounts.close();

 //save the deposited energy in each pixel for penetration
  sprintf(filename,"storepenetration.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  penetrationstore.open(pathandname,ios::out);
  sprintf(filename,"penetration_window.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  windowpenetration.open(pathandname,ios::out);

 //save total counts for each run
  sprintf(filename,"counts_run.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  countsperrun.open(pathandname,ios::out);
  countsperrun<<"run_order"<<" "<<"Total_counts"<<" "<<"10%_windowed_counts"<<" "<<"in_one_pixel_windowed_counts"<<endl;
  countsperrun.close();

  ii=0;  jj=1;  kk=0;  mm=0;
}

ExN01TrackerSD::~ExN01TrackerSD( )
{

  windowpenetration.close();
  for(int kk=0;kk<256*128;kk++)
   {
    if(store[kk][2]==0)
    {
        if(kk==128*128)
        {
            vector1<<endl;
            vector1<<"SD2";
        }
        if (kk%128==0&&kk!=0)
          { vector1<<endl;  }
       vector1<<store[kk][1]<<" ";
    }

   }
   vector1<<endl;
   vector1<<endl;
   vector1<<endl;
   vector1.close();

  for(int mm=0;mm<256*128;mm++)
  {
      if(mm==128*128)
      {
          Energy<<endl;
          Energy<<"SD2";
      }
      if(store[mm][2]==0)
     {
       if (mm%128==0&&mm!=0)
       { Energy<<endl; }
       Energy<<store[mm][0]<<" ";
     }

  }
  Energy.close();

  sprintf(filename,"saveNum.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  vector1.open(pathandname,ios::out|ios::app);
  for(int ii=0;ii<100000;ii++)
   {
      if(count[ii][1]!=0)
      { vector1<<count[ii][0]<<" "<<count[ii][1]<<endl; }
   }
  vector1<<endl;
  vector1<<endl;
  vector1<<endl;
  vector1.close();

  for(int tt=0;tt<256*128;tt++)
  {
      if(tt==128*128)
      {   penetrationstore<<endl;
          penetrationstore<<"SD2";
      }
      if(tt%128==0&&tt!=0)
      {
          penetrationstore<<endl;
      }
          penetrationstore<<storepenetration[tt][0]<<" ";
  }
  penetrationstore.close();

  sprintf(filename,"saveNum.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  vector1.open(pathandname,ios::out|ios::app);
   for(int jj=0;jj<100000;jj++)
   {
      if(count2[jj][1]!=0)
      {vector1<<count2[jj][0]<<" "<<count2[jj][1]<<endl;}
   }
  vector1.close();

  map<int,double>::iterator iter;
  //for(iter=pixel.begin();iter!=pixel.end();iter++)
  //{std::cout<<iter->first<<"   "<<iter->second<<endl;}
  for(int b=0;b<256;b++)
     {
      for(int a=0;a<128;a++)
        {

          pixelenergy<<pixel[b*128+a]<<" ";
          if(b==128&&a==127)
          {
              pixelenergy<<endl;
              pixelenergy<<"SD2";
          }

        }
      pixelenergy<<endl;
     }

  sprintf(filename,"differentrun.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  rotate.open(pathandname,ios::out|ios::app);
  for(int j=0;j<=runorder;j++)
   {  rotate<<j<<"**"<<endl;
        for(int i=0;i<256*128;i++)
      {
         if(i==128*128)
            {
                rotate<<endl;
                rotate<<"SD2";
            }
         if (i%128==0&&i!=0)
           { rotate<<endl; }
           rotate<<diffrun[i][j]<<" ";

      }
        rotate<<endl;
  }
  rotate.close();

  sprintf(filename,"differentcounts.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  countsforrun.open(pathandname,ios::out|ios::app);
  for(int j=0;j<=runorder;j++)
   {  countsforrun<<j<<"**"<<endl;
        for(int i=0;i<256*128;i++)
      {
         if(i==128*128)
            {
                countsforrun<<endl;
                countsforrun<<"SD2";
            }
         if (i%128==0&&i!=0)
           { countsforrun<<endl; }
           countsforrun<<diffrunofcounts[i][j]<<" ";

      }
        countsforrun<<endl;
  }
  countsforrun.close();

  sprintf(filename,"windowenergy.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  windowenergy.open(pathandname,ios::out|ios::app);
  for(int j=0;j<=runorder;j++)
   {  windowenergy<<j<<"**"<<endl;
        for(int i=0;i<256*128;i++)
      {
            if (i==128*128)
            {
                windowenergy<<endl;
                windowenergy<<"SD2";
            }
            if (i%128==0&&i!=0)
           { windowenergy<<endl; }
           windowenergy<<windowrun[i][j]<<" ";
      }
        windowenergy<<endl;
  }
  windowenergy.close();

  sprintf(filename,"windowcounts.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  windowcounts.open(pathandname,ios::out|ios::app);
  for(int j=0;j<=runorder;j++)
   {  windowcounts<<j<<"**"<<endl;
        for(int i=0;i<256*128;i++)
      {
            if (i==128*128)
            {
                windowcounts<<endl;
                windowcounts<<"SD2";
            }
            if (i%128==0&&i!=0)
           { windowcounts<<endl; }
           windowcounts<<windowrunofcounts[i][j]<<" ";


      }
        windowcounts<<endl;
  }
  windowcounts.close();

  sprintf(filename,"counts_run.xls");
  sprintf(pathandname,"%s/result/%s",buf,filename);
  countsperrun.open(pathandname,ios::out|ios::app);
  for(int JJ=0;JJ<=runorder;JJ++)
  {
      countsperrun<<JJ<<" "<<countsdiffrun[JJ]<<" "<<effectivecounts[JJ]<<" "<<effectivecounts2[JJ]<<endl;

  }
  countsperrun.close();
}


void ExN01TrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  static int HCID = -1;
  trackerCollection = new ExN01TrackerHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,trackerCollection);
}


G4bool ExN01TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)//Step function is passed to record data inside SD
{    
    G4Track* theTrack=aStep->GetTrack();
    G4double edep = aStep->GetTotalEnergyDeposit()/keV;
        if(edep==0.)
          {  return false;  }

    G4double       N=theTrack->GetCurrentStepNumber();
    G4double       K=theTrack->GetVertexKineticEnergy()/keV;
    G4double       TL =theTrack->GetTrackLength()/um;
    G4double       TID=aStep->GetTrack()->GetTrackID();
    G4double       PID=aStep->GetTrack()->GetParentID();
    G4double       GlobalTime=aStep->GetTrack()->GetGlobalTime();
    G4double       LocalTime=aStep->GetTrack()->GetLocalTime();
    G4double       ProperTime=aStep->GetTrack()->GetProperTime();
    G4String       Name=aStep->GetTrack()->GetDefinition()->GetParticleName();
    G4double       SL = theTrack->GetStepLength();
    G4double       Velocity=aStep->GetTrack()->GetVelocity();
    G4double       time = aStep->GetPreStepPoint()->GetGlobalTime()/ns;
    G4ThreeVector  post=aStep->GetPostStepPoint()->GetPosition()/um;
    G4ThreeVector  pre=aStep->GetPreStepPoint()->GetPosition()/um;
    G4StepPoint*   preStepPoint = aStep->GetPreStepPoint();
    G4Material*    material = preStepPoint->GetMaterial();
    G4double       KineticE=theTrack->GetKineticEnergy()/keV;

    const G4VProcess* process 
                      = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    G4String  procName = process->GetProcessName();
   // std::cout<< procName<<std::endl;

    if(procName=="phot")
        kk=kk+1;
    if(procName=="compt")
        mm=mm+1;

  // output the percent of the compt effect
  // std::cout<<double(mm)/(mm+kk)<<std::endl;

    ExN01TrackerHit* newHit = new ExN01TrackerHit();
   // Example of coloring hits by particle type
   if(theTrack->GetParticleDefinition()->GetParticleName() == "gamma")
      newHit->SetHitColour(G4Colour(1.0, 0.0, 0.0, 1.0)); // RED
   // else if(theTrack->GetParticleDefinition()->GetParticleName() =="e-") // or use the follwing code.
      else if(theTrack->GetParticleDefinition() == G4Electron::ElectronDefinition())
      newHit->SetHitColour(G4Colour(0.0, 1.0, 0.0, 1.0)); // GREEN
      else
      newHit->SetHitColour(G4Colour(0.0, 0.0, 1.0, 1.0)); // BLUE

//     newHit->SetSL (SL);                                                  //Step length in every step
       newHit->SetEdep     (edep);                                         //Energy deposit in every step
       newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
       newHit->SetPos      (aStep->GetPreStepPoint()->GetPosition() );
//     newHit->SetPid      (aStep->GetTrack()->GetTrackID());
//     newHit->Setname     (aStep->GetTrack()->GetDefinition()->GetParticleName());
       G4double energyDeposit= aStep ->GetTotalEnergyDeposit();
       newHit->SetEnergyDeposit(energyDeposit);
       
    G4VPhysicalVolume* volume
         = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

    if(Name=="e-")
       { ii=1;}
    else
       { ii=0;}

 //make sure this event has gamma particles.
    jj=jj&&ii;
    j=volume->GetCopyNo();

 // 10% energy window to filter some energies
    orderhit.push_back(j);
    orderenergy.push_back(edep);

    if(jj==0)
    {
    store[j][0]=edep+store[j][0];
    store[j][1]= store[j][1]+1;
    store[j][2]=jj;
    diffrun[j][runorder]=edep+diffrun[j][runorder];
    }
    
    if(countpenetration==event_id)
    {
       storepenetration[j][0]=edep+storepenetration[j][0];
       storepenetration[j][1]=storepenetration[j][1]+1;
    }

    sprintf(filename,"data.xls");
    sprintf(pathandname,"%s/result/%s",buf,filename);
    datafile.open(pathandname,ios::out|ios::app);  //opens the file data.xls in append mode
    datafile<<Name<<" "<<pre.x()<<" "<<pre.y()<<" "<<pre.z()<<" "<<edep<<" "<<j<<" "<<store[j][0]<<" "<<store[j][1]<<endl;

    G4int print;
    print=0;
    if(print==1)
    {
       G4cout<<"  Name       " <<Name<<"   "<<G4endl;
       G4cout<<"  PreStep     " <<pre<<"    "<<G4endl;
       G4cout<<"  PostStep    " <<post<<"   "<<G4endl;
       G4cout<<"  GetTrackLength"<<TL<<"   " <<G4endl;
       G4cout<<"  StepCurrentN " <<N<<"   "<<G4endl;
       G4cout<<"  KineticE   " <<KineticE<<"   "<<G4endl;
       G4cout<<"  GetGlobalTime "<<time<<"   "<<G4endl;
       G4cout<<"  Vertecx    " <<K << "   " <<G4endl;
       G4cout<<"  edep       " <<edep << "   " <<G4endl;
       G4cout<<"  TrackID    " <<TID << "   " <<G4endl;
       G4cout<<"  StepLength " <<SL << "   " <<G4endl;
       G4cout<<"  ParentID   " <<PID<<"   "<<G4endl;
       G4cout<<"  GlobalTime " <<GlobalTime<<"   "<<G4endl;
       G4cout<<"  LocalTime  " <<LocalTime<<"   "<<G4endl;
       G4cout<<"  ProperTime " <<ProperTime<<"   "<<G4endl;
       G4cout<<"  Velocity   " <<Velocity<<"   "<<G4endl;
       std::cout<<"  Materical   " <<material->GetName()<<"   "<<std::endl;

       if (theTrack->GetTrackID() == 1)  //determine whether it is a primary particle（Primary Particle TrackID=1）
       std::cout << "Particle is a primary " << std::endl;   //if it’s G4cout，Qt is the output in cout
       if (theTrack->GetParentID() == 1)  //determine if the current particle is produced by a primary particle
       std::cout << "Parent was a primary " << std::endl;  //
If you use std :: cout, output goes to the terminal
   }

//  ExN01TrackerHit* newHit = new ExN01TrackerHit();
//  newHit->SetEdep( edep );
//  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
    datafile.close();
    trackerCollection->insert( newHit );
    return true;
}

G4int i=-1;
G4int c=0;
void ExN01TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  i=i+1;
  G4int out;
  out=0;

if(trackerCollection)
  {
   energy=0;
   numberHits= trackerCollection->entries();//number of acquired photons
   if(numberHits!=0)
   {

    for (int j= 0; j < numberHits; j++)
     {
      ExN01TrackerHit* hit= (*trackerCollection)[j];
      energy += hit ->GetEnergyDeposit()/keV;//accumulated energy deposition at each step

      if(out==1) //output the related message
       {
       std::cout<< numberHits<<std::endl;
       std::cout<< hit ->GetPos()<<std::endl;
       }
     }
    c=c+1;
    std::cout<<event_id<<" "<<energy<<std::endl;

    for(vector<int>::size_type i=0;i<orderhit.size();++i)
       {
        if(i==0)
       { diffrunofcounts[orderhit[i]][runorder]=1+diffrunofcounts[orderhit[i]][runorder];}
    }
    count[int(energy)][0]=int(energy);
    count[int(energy)][1]= count[int(energy)][1]+1;

    sprintf(filename,"data.xls");
    sprintf(pathandname,"%s/result/%s",buf,filename);
    datafile.open(pathandname,ios::out|ios::app);  //open the file in append mode
    datafile<<" "<<" "<<" "<<" "<<" "<<" "<<energy;
    G4double FWHM,sigma;
    FWHM=0.03*122;
    sigma=FWHM/2.3548;
    energy=energy+sigma*sqrt(-2.0*log(G4UniformRand()))*cos(2.0*pi*G4UniformRand());

    // the code for adding the energy window
    // Tc: 126-154   Co57: 110-150  I-123:143-175
    G4double energy_min1,energy_min2;
    G4double energy_max1,energy_max2;
    if(nuclide_type==1)
     {energy_min1=126;
      energy_min2=126;
      energy_max1=154;
      energy_max2=154;}
    if(nuclide_type==2)
     {energy_min1=143;
      energy_min2=143;
      energy_max1=175;
      energy_max2=175;}
    if(nuclide_type==3)
     {energy_min1=110;
      energy_min2=110;
      energy_max1=150;
      energy_max2=150;}
    if(nuclide_type==4)
    { energy_min1=154;
      energy_min2=188;
      energy_max1=220;
      energy_max2=270;}
    if(nuclide_type==5)
    { energy_min1=63;
      energy_min2=63;
      energy_max1=88;
      energy_max2=88;}
    if(nuclide_type==6)
     {energy_min1=328;
      energy_min2=328;
      energy_max1=400;
      energy_max2=400;
    }

    if((energy>=energy_min1&&energy<=energy_max1)||(energy>=energy_min2&&energy<=energy_max2))
    {
     for(vector<int>::size_type i=0;i<orderhit.size();++i)
        {
         pixel[orderhit[i]]=orderenergy[i]+pixel[orderhit[i]];
         if(i==orderhit.size()-1&&orderhit[i]==compare)
         {
          windowrun[orderhit[i]][runorder]=energy+windowrun[orderhit[i]][runorder];
          windowrunofcounts[orderhit[0]][runorder]=1+windowrunofcounts[orderhit[0]][runorder];
          effectivecounts2[runorder]=effectivecounts2[runorder]+1;
         }
         if(i==0)
         {compare=orderhit[i];}
        }
         effectivecounts[runorder]=effectivecounts[runorder]+1;
         count_effective=count_effective+1;
         if(countpenetration==event_id)
         {windowpenetration<<event_id<<" "<<energy<<endl;}
         //std::cout<<c<<" "<<count_effective<<" "<<energy<<std::endl;
     }
  orderhit.clear();
  orderenergy.clear();

  datafile<<" "<<energy<<endl;
  datafile.close();

  count2[int(energy)][0]=int(energy);
  count2[int(energy)][1]=count2[int(energy)][1]+1;

  totalnum[runorder]=totalnum[runorder]+1;
  countsdiffrun[runorder]=countsdiffrun[runorder]+1;
    }
  }
}

