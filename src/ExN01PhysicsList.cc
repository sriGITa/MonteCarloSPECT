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

#include "ExN01PhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "globals.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProductionCutsTable.hh"

#include "G4IonConstructor.hh"
#include "G4RadioactiveDecay.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01PhysicsList::ExN01PhysicsList():  G4VUserPhysicsList()
{
   G4LossTableManager::Instance()->SetVerbose(0);

   defaultCutValue =0.1*cm;

   cutForGamma     = defaultCutValue;
   cutForElectron  = defaultCutValue;
   cutForPositron  = defaultCutValue;

   //add new units for radioActive decays
   //
   const G4double minute = 60*second;
   const G4double hour   = 60*minute;
   const G4double day    = 24*hour;
   const G4double year   = 365*day;
   new G4UnitDefinition("minute", "min", "Time", minute);
   new G4UnitDefinition("hour",   "h",   "Time", hour);
   new G4UnitDefinition("day",    "d",   "Time", day);
   new G4UnitDefinition("year",   "y",   "Time", year);

   SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01PhysicsList::~ExN01PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();

  // ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01PhysicsList::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  // mu+/-
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  // nu_e
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  // nu_mu
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01PhysicsList::ConstructMesons()
{
  //  mesons
  //  light mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01PhysicsList::ConstructBaryons()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
  AddStepMax();

  G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
  radioactiveDecay->SetHLThreshold(-1.*s);
  radioactiveDecay->SetICM(false);                //Internal Conversion
  radioactiveDecay->SetARM(false);               //Atomic Rearrangement
  

 // G4ProcessManager* ph = G4GenericIon :: GenericIon()->GetProcessManager ();
 // ph->AddProcess(radioactiveDecay , 0, -1, 1);}
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
      
  // Deexcitation (in case of Atomic Rearrangement)
  //
  G4UAtomicDeexcitation* de = new G4UAtomicDeexcitation();
  de->SetFluo(false);
  de->SetAuger(false);
  de->SetPIXE(false);
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4EmLivermorePhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4Generator2BS.hh"

// e+
#include "G4eplusAnnihilation.hh"

// mu+-
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4hBremsstrahlungModel.hh"
#include "G4hPairProductionModel.hh"

// hadrons
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

// msc models
#include "G4UrbanMscModel95.hh"
#include "G4WentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01PhysicsList::ConstructEM()
{
    if (verboseLevel >0){
    std::cout << "PhysicsList::SetCuts:";
    std::cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << std::endl;
    }
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    G4double highEnergyLimit = 10*MeV;
    G4double LivermoreHighEnergyLimit = 10*MeV;

    if (particleName == "gamma") {
     // gamma
     // Photoelectric effect - define low-energy model
     G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
     G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel =
     new G4LivermorePhotoElectricModel();
     thePhotoElectricEffect->SetEmModel(theLivermorePhotoElectricModel);
     ph->RegisterProcess(thePhotoElectricEffect, particle);

     // Compton scattering - define low-energy model
     G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
     G4LivermoreComptonModel* theLivermoreComptonModel =
     new G4LivermoreComptonModel();
     theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
     theComptonScattering->SetEmModel(theLivermoreComptonModel, 1);
     ph->RegisterProcess(theComptonScattering, particle);

     // gamma conversion - define low-energy model
     G4GammaConversion* theGammaConversion = new G4GammaConversion();
     G4VEmModel* theLivermoreGammaConversionModel =
     new G4LivermoreGammaConversionModel();
     theGammaConversion->SetEmModel(theLivermoreGammaConversionModel, 1);
     ph->RegisterProcess(theGammaConversion, particle);

     // default Rayleigh scattering is Livermore
     G4RayleighScattering* theRayleigh = new G4RayleighScattering();
     ph->RegisterProcess(theRayleigh, particle);
      
    } else if (particleName == "e-") {
      // multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      G4UrbanMscModel95* msc1 = new G4UrbanMscModel95();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->AddEmModel(0, msc1);
      msc->AddEmModel(0, msc2);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1);
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);
      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(ss, particle);

      // Ionisation - Livermore should be used only for low energies
      G4eIonisation* eIoni = new G4eIonisation();
      G4LivermoreIonisationModel* theIoniLivermore = new
      G4LivermoreIonisationModel();
      theIoniLivermore->SetHighEnergyLimit(0.1*MeV);
      eIoni->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation() );
      eIoni->SetStepFunction(0.2, 100*um); //
      ph->RegisterProcess(eIoni, particle);

      // Bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      G4VEmModel* theBremLivermore = new G4LivermoreBremsstrahlungModel();
      theBremLivermore->SetHighEnergyLimit(1*GeV);
      //theBremLivermore->SetAngularDistribution(new G4Generator2BS());
      eBrem->SetEmModel(theBremLivermore,1);

      ph->RegisterProcess(eBrem, particle);

    } else if (particleName == "e+") {
        // Identical to G4EmStandardPhysics_option3

      // multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering;
      msc->SetStepLimitType(fUseDistanceToBoundary);
      G4UrbanMscModel95* msc1 = new G4UrbanMscModel95();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(highEnergyLimit);
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->AddEmModel(0, msc1);
      msc->AddEmModel(0, msc2);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1);
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);

      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);

      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
      ph->RegisterProcess(ss, particle);
    
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      ph->RegisterProcess(new G4MuMultipleScattering, particle);
      ph->RegisterProcess(new G4MuIonisation,         particle);
      ph->RegisterProcess(new G4MuBremsstrahlung,     particle);
      ph->RegisterProcess(new G4MuPairProduction,     particle);
             
    } else if( particleName == "proton" || 
               particleName == "pi-" ||
               particleName == "pi+"    ) {
      //proton  
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4hIonisation,         particle);
      ph->RegisterProcess(new G4hBremsstrahlung,     particle);
      ph->RegisterProcess(new G4hPairProduction,     particle);       
     
   } else if( particleName == "alpha" ||
	       particleName == "He3" )     {
      //alpha 
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4ionIonisation,       particle);
     
   } else if( particleName == "GenericIon" ) {
      //Ions 
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4ionIonisation,       particle);     
      
   } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4hIonisation,         particle);        
    }     
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

void ExN01PhysicsList::ConstructGeneral()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
  G4ParticleDefinition* particle = theParticleIterator->value();
  if (theDecayProcess->IsApplicable(*particle)) {
  ph->RegisterProcess(theDecayProcess, particle);
    }
  }
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

void ExN01PhysicsList::AddStepMax()
{
  // Step limitation seen as a process
  G4StepLimiter* stepLimiter = new G4StepLimiter();
  ////G4UserSpecialCuts* userCuts = new G4UserSpecialCuts();
  
  theParticleIterator->reset();
  while ((*theParticleIterator)()){
  G4ParticleDefinition* particle = theParticleIterator->value();
  G4ProcessManager* pmanager = particle->GetProcessManager();

  if (particle->GetPDGCharge() != 0.0)
   {
	  pmanager ->AddDiscreteProcess(stepLimiter);
	  ////pmanager ->AddDiscreteProcess(userCuts);
   }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

void ExN01PhysicsList::SetCuts()
{
  //G4VUserPhysicsList::SetCutsWithDefault method sets 
  //the default cut value for all particle types 
  //  SetCutsWithDefault();

    if (verboseLevel >1){
       std::cout << "PhysicsList::SetCuts:"<<std::endl;
       std::cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << std::endl;
     }

    // default production thresholds for the world volume
    G4double lowlimit=50*eV;
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit, 10.*MeV);

    SetCutsWithDefault();

    SetCutValue(cutForGamma, "gamma");
    SetCutValue(cutForElectron/10, "e-");
    SetCutValue(cutForPositron/10, "e+");

    // Production thresholds for detector regions
    G4Region* region;
    G4String regName;
    G4ProductionCuts* cuts;

    regName = "BOX1";
    region = G4RegionStore::GetInstance()->GetRegion(regName);
    cuts = new G4ProductionCuts;
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e-"));
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e+"));
    region->SetProductionCuts(cuts);

    regName = "BOX2";
    region = G4RegionStore::GetInstance()->GetRegion(regName);
    cuts = new G4ProductionCuts;
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e-"));
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e+"));
    region->SetProductionCuts(cuts);

    regName = "BOX3";
    region = G4RegionStore::GetInstance()->GetRegion(regName);
    cuts = new G4ProductionCuts;
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e-"));
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e+"));
    region->SetProductionCuts(cuts);

    regName = "BOX4";
    region = G4RegionStore::GetInstance()->GetRegion(regName);
    cuts = new G4ProductionCuts;
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e-"));
    cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e+"));
    region->SetProductionCuts(cuts);

    /*
     // set cut values for gamma at first and for e- second and next for e+,
     // because some processes for e+/e- need cut values for gamma
     SetCutValue(cutForGamma, "gamma");
     SetCutValue(cutForElectron, "e-");
     SetCutValue(cutForPositron, "e+");
     */
     G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(3*keV, 10*MeV);
     if (verboseLevel>0) DumpCutValuesTable();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01PhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

