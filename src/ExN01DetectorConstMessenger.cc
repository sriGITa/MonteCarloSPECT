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
// --------------------------------------------------------------
//
#include "ExN01DetectorConstMessenger.hh"
#include "ExN01DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"


ExN01DetectorConstMessenger::ExN01DetectorConstMessenger(ExN01DetectorConstruction* mpga)
:fTarget(mpga)
{
  // To rotate a fixed object by typing in commmonds.
  fMydetDirectory = new G4UIdirectory("/mydet/");
  fMydetDirectory->SetGuidance("ExN01 detector setup control commands.");

  fArmCmd = new G4UIcmdWithADoubleAndUnit("/mydet/armAngle",this);
  fArmCmd->SetGuidance("Rotation angle of the second arm.");
  fArmCmd->SetParameterName("angle",true);
  fArmCmd->SetRange("angle>=0. && angle<=360.");
  fArmCmd->SetDefaultValue(0);
  fArmCmd->SetDefaultUnit("deg");
}

ExN01DetectorConstMessenger::~ExN01DetectorConstMessenger()
{
  delete fArmCmd;
  delete fMydetDirectory;
}

void ExN01DetectorConstMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==fArmCmd )
  { fTarget->SetArmAngle(fArmCmd->GetNewDoubleValue(newValue)); }
}

G4String ExN01DetectorConstMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command==fArmCmd )
  { cv = fArmCmd->ConvertToString(fTarget->GetArmAngle(),"deg"); }

  return cv;
}
