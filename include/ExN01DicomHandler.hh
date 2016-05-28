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
// The code was written by :
//      *Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Universit� Laval, Qu�bec (QC) Canada
//*******************************************************//

#ifndef ExN01DicomHandler_h
#define ExN01DicomHandler_h 1

#include <cstdio>
#include <map>
#include <fstream>
#include "iostream"
using namespace std;

#include "globals.hh"

//*******************************************************
/// Dicom Handler class
///        - Handling of DICM images
///        - Transforming *.dcm to *.g4 ( pixels->density )
///        - Reading headers and pixels
///        - Transforming pixel to density and creating *.g4
///          files
///        - Functions are in ExN01DicomHandler.cc
///
/// Base on previous code by :
///        Dragan Tubic <tdragan@gel.ulaval.ca>
//*******************************************************

class ExN01DicomHandler
{
public:

  ExN01DicomHandler();

    ~ExN01DicomHandler();

  G4int ReadFile(FILE *,char *);
  G4int ReadData(FILE *,char *); // note: always use readHeader 
                                    // before readData

  // use ImageMagick to display the image
  G4int displayImage(char[500]);

    void CheckFileFormat();
    G4double jj;

private:
    template <class Type> void GetValue(char *, Type &);

private:

    const int DATABUFFSIZE;
    const int LINEBUFFSIZE;
    const int FILENAMESIZE;

    void GetInformation(G4int &, char *);
    G4float Pixel2density(G4int pixel);
    void ReadMaterialIndices( std::ifstream& finData);
    unsigned int GetMaterialIndex( G4float density );
    void StoreData(std::ofstream& foutG4DCM);
    G4int read_defined_nested(FILE *, G4int);
    void read_undefined_nested(FILE *);
    void read_undefined_item(FILE *);

    short fCompression;
    G4int fNFiles;
    short fRows;
    short fColumns;
    short fBitAllocated;
    G4int fMaxPixelValue, fMinPixelValue;
    
    G4double fPixelSpacingX, fPixelSpacingY;
    G4double fSliceThickness;
    G4double fSliceLocation;
    
    G4int fRescaleIntercept, fRescaleSlope;
    
    G4bool fLittleEndian, fImplicitEndian;
    short fPixelRepresentation;

    G4int** fTab;
    std::map<G4float,G4String> fMaterialIndices;

    G4int numberoffiles;

};
#endif

