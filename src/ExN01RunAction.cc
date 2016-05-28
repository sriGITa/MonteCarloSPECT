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
// *                                                                      *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN01RunAction.hh"
#include "ExN01PrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "ExN01HistoManager.hh"
#include "ExN01TrackerSD.hh"

#include "G4Run.hh"
#include <iomanip>

#include <fstream>
#include  <iostream>
using namespace std;

int runorder=0;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01RunAction::ExN01RunAction(ExN01PrimaryGeneratorAction * kin):fPrimary(kin)
{
fHistoManager = new ExN01HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01RunAction::~ExN01RunAction()
{
    delete fHistoManager;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::BeginOfRunAction(const G4Run* aRun)
{

    std::cout << "### Run " << aRun->GetRunID() << " start." <<" ###"<<std::endl;
    runorder=aRun->GetRunID();

    fDecayCount = fTimeCount = 1;
    for (G4int i=0; i<3; i++) fEkinTot[i] = fPbalance[i] = fEventTime[i] = 0. ;
    fPrimaryTime = 0.;
    myTimer.Start();

    //histograms
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if (analysisManager->IsActive() ) {
      analysisManager->OpenFile();}
    //inform the runManager to save random number seed
    //
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::Balance(G4double Ekin, G4double Pbal)
{
  fDecayCount++;
  fEkinTot[0] += Ekin;
  //update min max
  if (fDecayCount == 1) fEkinTot[1] = fEkinTot[2] = Ekin;
  if (Ekin < fEkinTot[1]) fEkinTot[1] = Ekin;
  if (Ekin > fEkinTot[2]) fEkinTot[2] = Ekin;
  fPbalance[0] += Pbal;
  //update min max
  if (fDecayCount == 1) fPbalance[1] = fPbalance[2] = Pbal;
  if (Pbal < fPbalance[1]) fPbalance[1] = Pbal;
  if (Pbal > fPbalance[2]) fPbalance[2] = Pbal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::EventTiming(G4double time)
{
  fTimeCount++;
  fEventTime[0] += time;
  if (fTimeCount == 1) fEventTime[1] = fEventTime[2] = time;
  if (time < fEventTime[1]) fEventTime[1] = time;
  if (time > fEventTime[2]) fEventTime[2] = time;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::PrimaryTiming(G4double ptime)
{
  fPrimaryTime += ptime;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::ParticleCount(G4String name, G4double Ekin)
{
  fParticleCount[name]++;
  fEmean[name] += Ekin;
  //update min max
  if (fParticleCount[name] == 1) fEmin[name] = fEmax[name] = Ekin;
  if (Ekin < fEmin[name]) fEmin[name] = Ekin;
  if (Ekin > fEmax[name]) fEmax[name] = Ekin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::EndOfRunAction(const G4Run* run)
{

    G4int nbEvents = run->GetNumberOfEvent();
    // get the current directory
    getcwd(buf, sizeof(buf));
    chdir("..");
    // get the parent directory
    getcwd(bufup,sizeof(bufup));
    chdir(buf);
    sprintf(filename,"collimator_detector.txt");
    sprintf(pathandname,"%s/run/%s",bufup,filename);
    TEST = fopen(pathandname,"r");
    fscanf(TEST,"gpsorgun=%d\n",&gpsorgun);
    fclose(TEST);

    if(gpsorgun==1)
    {
    particle = fPrimary->GetParticleGun()
                                            ->GetParticleDefinition();
    partName = particle->GetParticleName();
    eprimary = fPrimary->GetParticleGun()->GetParticleEnergy();
    }
    if(gpsorgun==2)
    {
    particle = fPrimary->GetParticleSource()
                                            ->GetParticleDefinition();
    partName = particle->GetParticleName();
    eprimary = fPrimary->GetParticleSource()->GetParticleEnergy();
    }

    std::cout << "\n ======================== run summary ======================";
    std::cout << "\n The run was " << nbEvents << " " << partName << " of "
           << G4BestUnit(eprimary,"Energy");
    std::cout << "\n ===========================================================\n";
    std::cout << std::endl;

    G4int prec = 4, wid = prec + 2;
    G4int dfprec = G4cout.precision(prec);

    //particle count
    //
    std::cout << " Nb of generated particles: \n" << std::endl;

    std::map<G4String,G4int>::iterator it;
    for (it = fParticleCount.begin(); it != fParticleCount.end(); it++) {
       G4String name = it->first;
       G4int count   = it->second;
       G4double eMean = fEmean[name]/count;
       G4double eMin = fEmin[name], eMax = fEmax[name];

       std::cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
              << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
              << "\t( "  << G4BestUnit(eMin, "Energy")
              << " --> " << G4BestUnit(eMax, "Energy")
              << ")" << std::endl;
    }

    //energy momentum balance
    //

    if (fDecayCount > 0) {
       G4double Ebmean = fEkinTot[0]/fDecayCount;
       G4double Pbmean = fPbalance[0]/fDecayCount;

       std::cout << "\n   Ekin Total (Q): mean = "
              << std::setw(wid) << G4BestUnit(Ebmean, "Energy")
              << "\t( "  << G4BestUnit(fEkinTot[1], "Energy")
              << " --> " << G4BestUnit(fEkinTot[2], "Energy")
              << ")" << std::endl;

       std::cout << "\n   Momentum balance (excluding gamma desexcitation): mean = "
              << std::setw(wid) << G4BestUnit(Pbmean, "Energy")
              << "\t( "  << G4BestUnit(fPbalance[1], "Energy")
              << " --> " << G4BestUnit(fPbalance[2], "Energy")
              << ")" << std::endl;
    }
    //total time of life
    //
    if (fTimeCount > 0) {
       G4double Tmean = fEventTime[0]/fTimeCount;
       G4double halfLife = Tmean*std::log(2.);

       std::cout << "\n   Total time of life : mean = "
              << std::setw(wid) << G4BestUnit(Tmean, "Time")
              << "  half-life = "
              << std::setw(wid) << G4BestUnit(halfLife, "Time")
              << "   ( "  << G4BestUnit(fEventTime[1], "Time")
              << " --> "  << G4BestUnit(fEventTime[2], "Time")
              << ")" << std::endl;
       std::cout<<"+++"<<halfLife/s<<"+++"<<std::endl;
      }

    //activity of primary ion
     //
     G4double pTimeMean = fPrimaryTime/nbEvents;


     G4double molMass = particle->GetAtomicMass()*g/mole;
     G4double molMass2= particle->GetAtomicMass()/mole;
     G4double nAtoms_perUnitOfMass = CLHEP::Avogadro/molMass;

     G4double Activity_perUnitOfMass = 0.0;
     if (pTimeMean > 0.0)
       { Activity_perUnitOfMass = nAtoms_perUnitOfMass/pTimeMean; }
     std::cout<<"("<<molMass2<<"-->"<<molMass*mole/g<<"-->"<<molMass*g/mole<<"-->"<<nAtoms_perUnitOfMass<<"-->"<<Activity_perUnitOfMass/becquerel<<")"<<std::endl;

     std::cout << "\n   Activity of " << partName << " = "
                << std::setw(wid)  << Activity_perUnitOfMass*g/becquerel
                << " Bq/g   ("     << Activity_perUnitOfMass*g/curie
                << " Ci/g) \n"
                << std::endl;
     std::cout << "\n   Activity of " << partName << " = "
                << std::setw(wid)  << Activity_perUnitOfMass*g/becquerel
                << " Bq/g   ("     << Activity_perUnitOfMass*g/curie
                << " Ci/g) \n"
                << std::endl;

     // restore default precision
     //
     G4cout.precision(dfprec);


     // remove all contents in fParticleCount
     //
     fParticleCount.clear();
     fEmean.clear();  fEmin.clear(); fEmax.clear();




      // restore default precision
      //
      G4cout.precision(dfprec);

      //normalize and save histograms
      //
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      G4double factor = 100./nbEvents;
      analysisManager->ScaleH1(1,factor);
      analysisManager->ScaleH1(2,factor);
      analysisManager->ScaleH1(3,factor);
      analysisManager->ScaleH1(4,factor);
      analysisManager->ScaleH1(5,factor);
      if ( analysisManager->IsActive() ) {
       analysisManager->Write();
       analysisManager->CloseFile();
      }

      //output the total time of this run
      myTimer.Stop();
      std::cout<<myTimer<<"  "<<"Events:"<<totalnum[runorder]<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



