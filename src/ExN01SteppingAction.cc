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

#include "ExN01SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Material.hh"
#include "G4VProcess.hh"
#include "ExN01EventAction.hh"

using namespace std;

int penetration=1;
int status1=0;
int countpenetration=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01SteppingAction::ExN01SteppingAction():num(0),count(0)
{
    display=0;
    direction3=G4ThreeVector(0,0,0);
    direction4=G4ThreeVector(0,0,0);
    line=1;
    getcwd(buf, sizeof(buf));
    sprintf(filename,"penetrationrate.xls");
    sprintf(pathandname,"%s/result/%s",buf,filename);
    penetantionrate.open(pathandname,ios::out);
    penetantionrate<<"The number of events penetrate the material"<<std::endl;
    G4int z,nel;
    G4double a,density;
    //Air
    G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
    G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
    material2 = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
    material2->AddElement(N, 70*perCent);
    material2->AddElement(O, 30*perCent);
    status1=0;

    sprintf(filename,"primary_energy.xls");
    sprintf(pathandname,"%s/result/%s",buf,filename);
    primary_energy.open(pathandname,ios::out);
    runID_temp=0;
    energy_temp=0;
    energy_blur_temp=0;

    sprintf(filename,"counts_primary.xls");
    sprintf(pathandname,"%s/result/%s",buf,filename);
    counts_primary.open(pathandname,ios::out);
    counts_primary<<"Counts"<<" "<<"Initial_energy"<<" "<<"Counts"<<" "<<"Blurred_energy"<<endl;
    max=-100;max_blurred=-100;
    for(int i=0;i<10000;i++){
    energy_spectrum[i]=0;
    energy_blurred[i]=0;
    }
 }

ExN01SteppingAction::~ExN01SteppingAction()
{
 penetantionrate.close();
 primary_energy.close();
 int max_min;
 if(max>max_blurred)
     max_min=max_blurred;
 else
     max_min=max;

 for(int i=0;i<=max_min;i++)
 counts_primary<<energy_spectrum[i]<<" "<<i<<" "<<energy_blurred[i]<<" "<<i<<std::endl;
 if(max_min==max){
 for(int j=max_min+1;j<=max_blurred;j++)
 counts_primary<<" "<<" "<<energy_blurred[j]<<" "<<j<<std::endl;
 }
 else{
     for(int j=max_min+1;j<=max;j++)
        counts_primary<<energy_spectrum[j]<<" "<<j<<std::endl;
 }


 counts_primary.close();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01SteppingAction::UserSteppingAction(const G4Step* Step)
{


   G4Track*       theTrack=Step->GetTrack();
   G4int          ID=theTrack->GetTrackID();
   G4double       depositedenergy = Step->GetTotalEnergyDeposit()/keV;
   G4String       Name=Step->GetTrack()->GetDefinition()->GetParticleName();
   G4Material*    material1 = Step->GetPreStepPoint()->GetMaterial();
   G4double       TID=Step->GetTrack()->GetTrackID();
   G4double       PID=Step->GetTrack()->GetParentID();
   G4double       K=theTrack->GetVertexKineticEnergy()/keV;
   G4ThreeVector  direction1=Step->GetPreStepPoint()->GetMomentumDirection();
   G4ThreeVector  direction2=Step->GetPostStepPoint()->GetMomentumDirection();
   G4ThreeVector  post11=Step->GetPostStepPoint()->GetPosition()/cm;
   G4ThreeVector  pre11=Step->GetPreStepPoint()->GetPosition()/cm;
   G4double       vecity=Step->GetPreStepPoint()->GetVelocity();

   if(display==1)
   {
       std::cout<<ID<<" "<<TID<<"  "<<PID<<"  "<<vecity<<"  "<<std::endl;
       std::cout<<direction1<<"  "<<direction2<<"  "<<direction3<<"  "<<std::endl;
       std::cout<<K<<std::endl;
   }

    status2=0;
   if((Name=="gamma")&&material2->GetName()=="Tungsten")
   {

     //theTrack->SetTrackStatus(fStopAndKill);
   if((direction2==direction1)&&(direction1==direction3))
      { line=1;}
   else
      { line=0;}
    penetration=penetration&&line;
    //std::cout<<penetration<<std::endl;
    status2=1;
   }
   status1=status1||status2;
   if(status1==1&&penetration==1&&material1->GetName()=="CdZnTe"&&material2->GetName()=="Tungsten")
    { num=num+1;
      vi.push_back(num);
      penetantionrate<<penetration<<"**"<<depositedenergy<<"**"<<num<<"**"<<material1->GetName()<<endl;
      countpenetration=event_id;
      std::cout<<event_id<<"**"<<Name<<std::endl;
      }
   direction3=direction1;

  if((Name=="e-")&&(material1->GetName()=="Tungsten"||material1->GetName()=="Air"||material1->GetName()=="Lead"))
  {
    // theTrack->SetTrackStatus(fStopAndKill);
  }

   if((Name=="nu_e")&&(material1->GetName()=="Tungsten"||material1->GetName()=="Lead"||material1->GetName()=="Air"))
   {
     count=count+1;
     vj.push_back(count);
     theTrack->SetTrackStatus(fStopAndKill);
    }
   material2=material1;
   if(runID_temp==event_id)
   {
      energy_temp+=depositedenergy;
    }
   else
   {
       primary_energy<<runID_temp<<" "<<energy_temp<<" ";
      if(runID_temp>0&&runID_temp%200==0)
          primary_energy<<endl;
      runID_temp=event_id;
      energy_spectrum[(int)(energy_temp+0.5)]+=1;

      G4double FWHM,sigma;
      FWHM=0.03*122;
      sigma=FWHM/2.3548;
      energy_blur_temp=energy_temp+sigma*sqrt(-2.0*log(G4UniformRand()))*cos(2.0*pi*G4UniformRand());
      energy_blurred[(int)(energy_blur_temp+0.5)]+=1;

      if((int)(energy_temp+0.5)>max)
          max=(int)(energy_temp+0.5);
      if((int)(energy_blur_temp+0.5)>max_blurred)
          max_blurred=(int)(energy_blur_temp+0.5);

      energy_temp=depositedenergy;
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
