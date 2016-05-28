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

#include "ExN01TrackerHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4AttValue.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"


G4Allocator<ExN01TrackerHit> ExN01TrackerHitAllocator;

ExN01TrackerHit::ExN01TrackerHit()
{;}

ExN01TrackerHit::~ExN01TrackerHit()
{;}

ExN01TrackerHit::ExN01TrackerHit(const ExN01TrackerHit &right)
  : G4VHit()
{
  edep = right.edep;
  pos = right.pos;

}

const ExN01TrackerHit& ExN01TrackerHit::operator=(const ExN01TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;

  return *this;
}

G4int ExN01TrackerHit::operator==(const ExN01TrackerHit &right) const
{
  return (this==&right) ? 1 : 0;
}

std::map<G4String,G4AttDef> ExN01TrackerHit::fAttDefs;

void ExN01TrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(3);
    circle.SetFillStyle(G4Circle::filled);
   // G4Colour colour(1.,0.,1.);
    G4VisAttributes attribs(HitColour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

const std::map<G4String,G4AttDef>* ExN01TrackerHit::GetAttDefs() const
{
  // G4AttDefs have to have long life.  Use static member...
  if (fAttDefs.empty()) {
      fAttDefs["HitType"] =
      G4AttDef("HitType","Type of hit","Physics","","G4String");
  }
  return &fAttDefs;
}

std::vector<G4AttValue>* ExN01TrackerHit::CreateAttValues() const
{
  // Create expendable G4AttsValues for picking...
  std::vector<G4AttValue>* attValues = new std::vector<G4AttValue>;
  attValues->push_back
    (G4AttValue("HitType","ExN01TrackerHit",""));
  //G4cout << "Checking...\n" << G4AttCheck(attValues, GetAttDefs());
  return attValues;
}

void ExN01TrackerHit::Print()
{;}


