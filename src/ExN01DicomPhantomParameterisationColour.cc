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


#include "ExN01DicomPhantomParameterisationColour.hh"

#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN01DicomPhantomParameterisationColour::ExN01DicomPhantomParameterisationColour()
{
  ReadColourData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN01DicomPhantomParameterisationColour::~ExN01DicomPhantomParameterisationColour()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01DicomPhantomParameterisationColour::ReadColourData()
{
  //----- Add a G4VisAttributes for materials not defined in file;
  G4VisAttributes* blankAtt = new G4VisAttributes;
  blankAtt->SetVisibility(true);
  fColours["Default"] = blankAtt;

  //----- Read file
  std::ifstream fin("DICOM/ColourMap.dat");
  G4int nMate;
  G4String mateName;
  G4double cred, cgreen, cblue, copacity;// opacity 不透明
  fin >> nMate;
  for( G4int ii = 0; ii < nMate; ii++ ){
    fin >> mateName >> cred >> cgreen >> cblue >> copacity;
    G4Colour colour( cred, cgreen, cblue, copacity );
    G4VisAttributes* visAtt = new G4VisAttributes( colour );
    // set the phantom to be invisibiled
    visAtt->SetVisibility(true);
    fColours[mateName] = visAtt;
  }

}
       
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Material* ExN01DicomPhantomParameterisationColour::
ComputeMaterial(const G4int copyNo, G4VPhysicalVolume * physVol, const G4VTouchable *) 
{ 
  G4Material* mate = G4PhantomParameterisation::ComputeMaterial( copyNo, physVol, 0 );
  if( physVol ) {
    G4String mateName = mate->GetName();
    std::string::size_type iuu = mateName.find("__");
    if( iuu != std::string::npos ) {
      mateName = mateName.substr( 0, iuu );
    }
    std::map<G4String,G4VisAttributes*>::const_iterator ite = fColours.find(mateName);
    if( ite != fColours.end() ){
      physVol->GetLogicalVolume()->SetVisAttributes( (*ite).second ); 
    } else {
      physVol->GetLogicalVolume()->SetVisAttributes( (*(fColours.begin()) ).second ); // set it as unseen
    }
  }  

  return mate;
}

