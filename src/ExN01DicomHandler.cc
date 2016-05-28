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
//*******************************************************
//
//*******************************************************
//
/// ExN01DicomHandler.cc :
///        - Handling of DICM images
///         - Reading headers and pixels
///        - Transforming pixel to density and creating *.g4dcm
///          files
//*******************************************************
//the normal size of the human head is 196mm*155mm*230mm
#include "ExN01DicomHandler.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cctype>
#include <cstring>
#include <iostream>
#include <sstream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01DicomHandler::ExN01DicomHandler()
    : DATABUFFSIZE(8192), LINEBUFFSIZE(5020), FILENAMESIZE(512),
      fCompression(0), fNFiles(0), fRows(0), fColumns(0),
      fBitAllocated(0), fMaxPixelValue(0), fMinPixelValue(0),
      fPixelSpacingX(0.), fPixelSpacingY(0.),
      fSliceThickness(0.), fSliceLocation(0.),
      fRescaleIntercept(0), fRescaleSlope(0),
      fLittleEndian(true), fImplicitEndian(false),
      fPixelRepresentation(0) {
    jj=0;
    numberoffiles=0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01DicomHandler::~ExN01DicomHandler()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int ExN01DicomHandler::ReadFile(FILE *dicom, char * filename2)
{
  G4cout << " ReadFile " << filename2 << G4endl;
    G4int returnvalue = 0; size_t rflag = 0;
    char * buffer = new char[LINEBUFFSIZE];

    fImplicitEndian = false;
    fLittleEndian = true;

    rflag = std::fread( buffer, 1, 128, dicom ); // The first 128 bytes 
                                         //are not important
    // Reads the "DICOM" letters
    rflag = std::fread( buffer, 1, 4, dicom );
    // if there is no preamble, the FILE pointer is rewinded.
    if(std::strncmp("DICM", buffer, 4) != 0) {
        std::fseek(dicom, 0, SEEK_SET);
        fImplicitEndian = true;
    }

    short readGroupId;    // identify the kind of input data 
    short readElementId;  // identify a particular type information
    short elementLength2; // deal with element length in 2 bytes
    //unsigned int elementLength4; // deal with element length in 4 bytes
    unsigned long elementLength4; // deal with element length in 4 bytes

    char * data = new char[DATABUFFSIZE];

 
    // Read information up to the pixel data
    while(true) {

        //Reading groups and elements :
        readGroupId = 0;
        readElementId = 0;
        // group ID
        rflag = std::fread(buffer, 2, 1, dicom);
        GetValue(buffer, readGroupId);
        // element ID
        rflag = std::fread(buffer, 2, 1, dicom);
        GetValue(buffer, readElementId);

        // Creating a tag to be identified afterward
        G4int tagDictionary = readGroupId*0x10000 + readElementId;

        // beginning of the pixels
            if(tagDictionary == 0x7FE00010) break;

      // VR or element length
        rflag = std::fread(buffer,2,1,dicom);
        GetValue(buffer, elementLength2);
          
         // If value representation (VR) is OB, OW, SQ, UN, added OF and UT 
        //the next length is 32 bits
        if((elementLength2 == 0x424f ||  // "OB"
            elementLength2 == 0x574f ||  // "OW"
            elementLength2 == 0x464f ||  // "OF"
            elementLength2 == 0x5455 ||  // "UT"
            elementLength2 == 0x5153 || //  "SQ"
            elementLength2 == 0x4e55) && // "UN"
           !fImplicitEndian ) {           // explicit VR

            rflag = std::fread(buffer, 2, 1, dicom); // Skip 2 reserved bytes

            // element length
            rflag = std::fread(buffer, 4, 1, dicom);
            GetValue(buffer, elementLength4);
            
            if(elementLength2 == 0x5153)
            {
             if(elementLength4 == 0xFFFFFFFF)           
             {
              read_undefined_nested( dicom );
              elementLength4=0;           
             }  else{
               if(read_defined_nested( dicom, elementLength4 )==0){
               G4cerr << "Function read_defined_nested() failed!" << G4endl;
               exit(-10);               }
              }
            } else  { 
            // Reading the information with data
            rflag = std::fread(data, elementLength4,1,dicom);
            }  

                
        }  else { 

                if(!fImplicitEndian || readGroupId == 2) {  //  explicit with VR different than previous ones
                    
                  //G4cout << "Reading  DICOM files with Explicit VR"<< G4endl;
                  // element length (2 bytes)
                  rflag = std::fread(buffer, 2, 1, dicom);
                  GetValue(buffer, elementLength2);
                  elementLength4 = elementLength2;
                  
                  rflag = std::fread(data, elementLength4, 1, dicom);
                
                } else {                                  // Implicit VR

                  //G4cout << "Reading  DICOM files with Implicit VR"<< G4endl;
   
                  // element length (4 bytes)
                  if(std::fseek(dicom, -2, SEEK_CUR) != 0) {
                      G4cerr << "[ExN01DicomHandler] fseek failed" << G4endl;
                      exit(-10);}

                  rflag = std::fread(buffer, 4, 1, dicom);
                  GetValue(buffer, elementLength4);

                  //G4cout <<  std::hex<< elementLength4 << G4endl;
              
                  if(elementLength4 == 0xFFFFFFFF) 
                      {
                      read_undefined_nested(dicom);
                     elementLength4=0;           
                  }  else{
                  rflag = std::fread(data, elementLength4, 1, dicom);
                 } 
                      
               } 
        }

        // NULL termination
        data[elementLength4] = '\0';

        // analyzing information
       GetInformation(tagDictionary, data);
    }

    // Creating files to store information
    std::ofstream foutG4DCM;
    G4String fnameG4DCM = G4String(filename2) + ".g4dcm";
    foutG4DCM.open(fnameG4DCM);
    G4cout << "### Writing of " << fnameG4DCM << " ### " << G4endl;

    foutG4DCM << fMaterialIndices.size() << G4endl;
    //--- Write materials
    unsigned int ii = 0;
    std::map<G4float,G4String>::const_iterator ite;
    for( ite = fMaterialIndices.begin(); ite != fMaterialIndices.end(); ite++, ii++ ){
      foutG4DCM << ii << " " << (*ite).second << G4endl;
    }
    //--- Write number of voxels (assume only one voxel in Z)
    foutG4DCM << fColumns/fCompression << " " << fRows/fCompression << " 1 " << G4endl;
    //--- Write minimum and maximum extensions
    foutG4DCM << -fPixelSpacingX*fColumns/4. << " " << fPixelSpacingX*fColumns/4. << G4endl;
    foutG4DCM << -fPixelSpacingX*fRows/4. << " " << fPixelSpacingX*fRows/4. <<
     G4endl;


 // foutG4DCM << fSliceLocation-fSliceThickness/2.<< " " << fSliceLocation+fSliceThickness/2 << G4endl;
    std::cout<<fSliceLocation-fSliceThickness/2.+jj-10<<"**"<<std::endl;
    foutG4DCM << fSliceLocation-fSliceThickness/2.+jj-10<< " " << fSliceLocation+fSliceThickness/2.+jj-10 << G4endl;
    jj=jj+2.0;   //2.0 equal to the thickness*2  10mm equal to the thickness

    // modify the thickness of the phantom, you should updata line 226,227,228,358
    // foutG4DCM << fCompression << G4endl;
    
    ReadData( dicom, filename2 );
    numberoffiles=numberoffiles+1;
    
    StoreData( foutG4DCM );

    foutG4DCM.close();

    //
    delete [] buffer;
    delete [] data;

    if (rflag) return returnvalue;
    return returnvalue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01DicomHandler::GetInformation(G4int & tagDictionary, char * data) {
    if(tagDictionary == 0x00280010 ) { // Number of Rows
        GetValue(data, fRows);
        std::printf("[0x00280010] Rows -> %i\n",fRows);

    } else if(tagDictionary == 0x00280011 ) { // Number of fColumns
        GetValue(data, fColumns);
        std::printf("[0x00280011] Columns -> %i\n",fColumns);

    } else if(tagDictionary == 0x00280102 ) { // High bits  ( not used )
        short highBits;
        GetValue(data, highBits);
        std::printf("[0x00280102] High bits -> %i\n",highBits);

    } else if(tagDictionary == 0x00280100 ) { // Bits allocated
        GetValue(data, fBitAllocated);
        std::printf("[0x00280100] Bits allocated -> %i\n", fBitAllocated);

    } else if(tagDictionary == 0x00280101 ) { //  Bits stored ( not used )
        short bitStored;
        GetValue(data, bitStored);
        std::printf("[0x00280101] Bits stored -> %i\n",bitStored);

    } else if(tagDictionary == 0x00280106 ) { //  Min. pixel value
        GetValue(data, fMinPixelValue);
        std::printf("[0x00280106] Min. pixel value -> %i\n", fMinPixelValue);

    } else if(tagDictionary == 0x00280107 ) { //  Max. pixel value
        GetValue(data, fMaxPixelValue);
        std::printf("[0x00280107] Max. pixel value -> %i\n", fMaxPixelValue);

    } else if(tagDictionary == 0x00281053) { //  Rescale slope
        fRescaleSlope = atoi(data);
        std::printf("[0x00281053] Rescale Slope -> %d\n", fRescaleSlope);

    } else if(tagDictionary == 0x00281052 ) { // Rescalse intercept
        fRescaleIntercept = atoi(data);
        std::printf("[0x00281052] Rescale Intercept -> %d\n", fRescaleIntercept );

    } else if(tagDictionary == 0x00280103 ) {
        //  Pixel representation ( functions not design to read signed bits )
        fPixelRepresentation = atoi(data); // 0: unsigned  1: signed 
        std::printf("[0x00280103] Pixel Representation -> %i\n", fPixelRepresentation);
        if(fPixelRepresentation == 1 ) {
            std::printf("### PIXEL REPRESENTATION = 1, BITS ARE SIGNED, ");
            std::printf("DICOM READING SCAN FOR UNSIGNED VALUE, POSSIBLE ");
            std::printf("ERROR !!!!!! -> \n");
        }

    } else if(tagDictionary == 0x00080060 ) { //  Modality
        std::printf("[0x00080006] Modality -> %s\n", data);

    } else if(tagDictionary == 0x00080070 ) { //  Manufacturer
        std::printf("[0x00080070] Manufacturer -> %s\n", data);

    } else if(tagDictionary == 0x00080080 ) { //  Institution Name
        std::printf("[0x00080080] Institution Name -> %s\n", data);

    } else if(tagDictionary == 0x00080081 ) { //  Institution Address
        std::printf("[0x00080081] Institution Address -> %s\n", data);

    } else if(tagDictionary == 0x00081040 ) { //  Institution Department Name
        std::printf("[0x00081040] Institution Department Name -> %s\n", data);

    } else if(tagDictionary == 0x00081090 ) { //  Manufacturer's Model Name
        std::printf("[0x00081090] Manufacturer's Model Name -> %s\n", data);

    } else if(tagDictionary == 0x00181000 ) { //  Device Serial Number
        std::printf("[0x00181000] Device Serial Number -> %s\n", data);

    } else if(tagDictionary == 0x00080008 ) { //  Image type ( not used )
        std::printf("[0x00080008] Image Types -> %s\n", data);
            
    } else if(tagDictionary == 0x00283000 ) { //  Modality LUT Sequence ( not used )
        std::printf("[0x00283000] Modality LUT Sequence SQ 1 -> %s\n", data);

    } else if(tagDictionary == 0x00283002 ) { // LUT Descriptor ( not used )
        std::printf("[0x00283002] LUT Descriptor US or SS 3 -> %s\n", data);

    } else if(tagDictionary == 0x00283003 ) { // LUT Explanation ( not used )
        std::printf("[0x00283003] LUT Explanation LO 1 -> %s\n", data);

    } else if(tagDictionary == 0x00283004 ) { // Modality LUT ( not used )
        std::printf("[0x00283004] Modality LUT Type LO 1 -> %s\n", data);

    } else if(tagDictionary == 0x00283006 ) { // LUT Data ( not used )
        std::printf("[0x00283006] LUT Data US or SS -> %s\n", data);

    } else if(tagDictionary == 0x00283010 ) { // VOI LUT ( not used )
        std::printf("[0x00283010] VOI LUT Sequence SQ 1 -> %s\n", data);

    } else if(tagDictionary == 0x00280120 ) { // Pixel Padding Value ( not used )
        std::printf("[0x00280120] Pixel Padding Value US or SS 1 -> %s\n", data);

    } else if(tagDictionary == 0x00280030 ) { // Pixel Spacing
      G4String datas(data);
      int iss = datas.find('\\');
      fPixelSpacingX = atof( datas.substr(0,iss).c_str() );
      fPixelSpacingX = fPixelSpacingX*0.7;
      fPixelSpacingY = atof( datas.substr(iss+2,datas.length()).c_str() );
      fPixelSpacingY = fPixelSpacingY*0.7;
    } else if(tagDictionary == 0x00200037 ) { // Image Orientation ( not used )
        std::printf("[0x00200037] Image Orientation (Phantom) -> %s\n", data);

    } else if(tagDictionary == 0x00200032 ) { // Image Position ( not used )
        std::printf("[0x00200032] Image Position (Phantom,mm) -> %s\n", data);

    } else if(tagDictionary == 0x00180050 ) { // Slice Thickness
        fSliceThickness = atof(data);
        fSliceThickness = fSliceThickness*2;  // change the z slice thickness

        std::printf("[0x00180050] Slice Thickness (mm) -> %f\n", fSliceThickness);

    } else if(tagDictionary == 0x00201041 ) { // Slice Location
        fSliceLocation = atof(data);
        std::printf("[0x00201041] Slice Location -> %f\n", fSliceLocation);

    } else if(tagDictionary == 0x00280004 ) { // Photometric Interpretation ( not used )
        std::printf("[0x00280004] Photometric Interpretation -> %s\n", data);

    } else if(tagDictionary == 0x00020010) { // Endian
        if(strcmp(data, "1.2.840.10008.1.2") == 0)
            fImplicitEndian = true;
        else if(strncmp(data, "1.2.840.10008.1.2.2", 19) == 0)
            fLittleEndian = false;
        //else 1.2.840..10008.1.2.1 (explicit little endian)
                   
        std::printf("[0x00020010] Endian -> %s\n", data);
    }

    // others
    else {
        //std::printf("[0x%x] -> %s\n", tagDictionary, data);
        ;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01DicomHandler::StoreData(std::ofstream& foutG4DCM) 
{
  G4int mean;
  G4double density;
  G4bool overflow = false;

  //---------------***------------

  //to change the attenuation from here ,change the density
  ifstream at;
  G4String name;
  stringstream o;
  o<<numberoffiles;
  o>>name;
  name="att"+name+".txt";
  at.open(name,ios::in);

  //----- Print indices of material 
  if(fCompression == 1) { // no fCompression: each pixel has a density value)
    for( G4int ww = 0; ww < fRows; ww++) {
      for( G4int xx = 0; xx < fColumns; xx++) {
        mean = fTab[ww][xx];
        density = Pixel2density(mean);
        at>>density;
        foutG4DCM << GetMaterialIndex( density ) << " ";
      }
      foutG4DCM << G4endl;
    }
    
  } else {
    // density value is the average of a square region of
    // fCompression*fCompression pixels
    for(G4int ww = 0; ww < fRows ;ww += fCompression ) {
      for(G4int xx = 0; xx < fColumns ;xx +=fCompression ) {
        overflow = false;
        mean = 0;
        for(int sumx = 0; sumx < fCompression; sumx++) {
          for(int sumy = 0; sumy < fCompression; sumy++) {
            if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
            mean += fTab[ww+sumy][xx+sumx];
          }
          if(overflow) break;
        }
        mean /= fCompression*fCompression;
        
        if(!overflow) {
          density = Pixel2density(mean);
          at>>density;
          foutG4DCM << GetMaterialIndex( density ) << " ";
        }
      }
      foutG4DCM << G4endl;
    }

  }
  at.close();

  //----- Print densities
  if(fCompression == 1) { // no fCompression: each pixel has a density value)
    for( G4int ww = 0; ww < fRows; ww++) {
      for( G4int xx = 0; xx < fColumns; xx++) {
        mean = fTab[ww][xx];
        density = Pixel2density(mean);
        foutG4DCM << density << " ";
        if( xx%8 == 3 ) foutG4DCM << G4endl; // just for nicer reading
      }
    }
    
  } else {
    // density value is the average of a square region of
    // fCompression*fCompression pixels
    for(G4int ww = 0; ww < fRows ;ww += fCompression ) {
      for(G4int xx = 0; xx < fColumns ;xx +=fCompression ) {
        overflow = false;
        mean = 0;
        for(int sumx = 0; sumx < fCompression; sumx++) {
          for(int sumy = 0; sumy < fCompression; sumy++) {
            if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
            mean += fTab[ww+sumy][xx+sumx];
          }
          if(overflow) break;
        }
        mean /= fCompression*fCompression;
        
        if(!overflow) {
          density = Pixel2density(mean);
          foutG4DCM << density  << " ";
          if( xx/fCompression%8 == 3 )
              foutG4DCM << G4endl; // just for nicer reading
        }
      }
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01DicomHandler::ReadMaterialIndices( std::ifstream& finData)
{
  unsigned int nMate;
  G4String mateName;
  G4float densityMax;
  finData >> nMate;
  if( finData.eof() ) return;

 // G4cout << " ReadMaterialIndices " << nMate << G4endl;
  for( unsigned int ii = 0; ii < nMate; ii++ ){
    finData >> mateName >> densityMax;
    fMaterialIndices[densityMax] = mateName;
   // G4cout << ii << " ReadMaterialIndices " << mateName << " " << densityMax << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
unsigned int ExN01DicomHandler::GetMaterialIndex( G4float density )
{
  std::map<G4float,G4String>::reverse_iterator ite;
  G4int ii = fMaterialIndices.size();

  for( ite = fMaterialIndices.rbegin(); ite != fMaterialIndices.rend(); ite++, ii-- ) {
    if( density >= (*ite).first ) {
      break;
    }
  }

  //-  G4cout << " GetMaterialIndex " << density << " = " << ii << G4endl;
  return  ii;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int ExN01DicomHandler::ReadData(FILE *dicom,char * filename2)
{
    G4int returnvalue = 0; size_t rflag = 0;


    //  READING THE PIXELS :
    G4int w = 0;
    
    fTab = new G4int*[fRows];
    for ( G4int i = 0; i < fRows; i ++ ) {
      fTab[i] = new G4int[fColumns];
    }

    if(fBitAllocated == 8) { // Case 8 bits :

        std::printf("@@@ Error! Picture != 16 bits...\n");
        std::printf("@@@ Error! Picture != 16 bits...\n"); 
        std::printf("@@@ Error! Picture != 16 bits...\n"); 

        unsigned char ch = 0;

        for(G4int j = 0; j < fRows; j++) {
            for(G4int i = 0; i < fColumns; i++) {
                w++;
                rflag = std::fread( &ch, 1, 1, dicom);
                fTab[j][i] = ch*fRescaleSlope + fRescaleIntercept;

            }
        }
        returnvalue = 1;

    } else { //  from 12 to 16 bits :
        char sbuff[2];
        short pixel;
        for( G4int j = 0; j < fRows; j++) {
            for( G4int i = 0; i < fColumns; i++) {
                w++;
                rflag = std::fread(sbuff, 2, 1, dicom);
                GetValue(sbuff, pixel);
                //The CT value
                fTab[j][i] = pixel*fRescaleSlope + fRescaleIntercept;
             }
        }

    }

    // Creation of .g4 files which contains averaged density data
    char * nameProcessed = new char[FILENAMESIZE];
    FILE* fileOut;

    std::sprintf(nameProcessed,"DICOM/%s.g4dcmb",filename2);
    fileOut = std::fopen(nameProcessed,"w+b");
    std::printf("### Writing of %s ###\n",nameProcessed);

    unsigned int nMate = fMaterialIndices.size();
    rflag = std::fwrite(&nMate, sizeof(unsigned int), 1, fileOut);
    //--- Write materials
    std::map<G4float,G4String>::const_iterator ite;
    for( ite = fMaterialIndices.begin(); ite != fMaterialIndices.end(); ite++ ){
      G4String mateName = (*ite).second;
      for( G4int ii = (*ite).second.length(); ii < 40; ii++ ) {
        mateName += " ";
      }         //mateName = const_cast<char*>(((*ite).second).c_str());

      const char* mateNameC = mateName.c_str();
      rflag = std::fwrite(mateNameC, sizeof(char),40, fileOut);
    }

    unsigned int fRowsC = fRows/fCompression;
    unsigned int fColumnsC = fColumns/fCompression;
    unsigned int planesC = 1;
    G4float pixelLocationXM = -fPixelSpacingX*fColumns/2.;
    G4float pixelLocationXP = fPixelSpacingX*fColumns/2.;
    G4float pixelLocationYM = -fPixelSpacingY*fRows/2.;
    G4float pixelLocationYP = fPixelSpacingY*fRows/2.;
    G4float fSliceLocationZM = fSliceLocation-fSliceThickness/2.;
    G4float fSliceLocationZP = fSliceLocation+fSliceThickness/2.;
    //--- Write number of voxels (assume only one voxel in Z)
    rflag = std::fwrite(&fRowsC, sizeof(unsigned int), 1, fileOut);
    rflag = std::fwrite(&fColumnsC, sizeof(unsigned int), 1, fileOut);
    rflag = std::fwrite(&planesC, sizeof(unsigned int), 1, fileOut);
    //--- Write minimum and maximum extensions
    rflag = std::fwrite(&pixelLocationXM, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&pixelLocationXP, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&pixelLocationYM, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&pixelLocationYP, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&fSliceLocationZM, sizeof(G4float), 1, fileOut);
    rflag = std::fwrite(&fSliceLocationZP, sizeof(G4float), 1, fileOut);
    // rflag = std::fwrite(&fCompression, sizeof(unsigned int), 1, fileOut);

    std::printf("%8i   %8i\n",fRows,fColumns);
    std::printf("%8f   %8f\n",fPixelSpacingX,fPixelSpacingY);
    std::printf("%8f\n", fSliceThickness);
    std::printf("%8f\n", fSliceLocation);
    std::printf("%8i\n", fCompression);

    G4int compSize = fCompression;
    G4int mean;
    G4float density;
    G4bool overflow = false;

    //----- Write index of material for each pixel
    if(compSize == 1) { // no fCompression: each pixel has a density value)
      for( G4int ww = 0; ww < fRows; ww++) {
        for( G4int xx = 0; xx < fColumns; xx++) {
          mean = fTab[ww][xx];
          density = Pixel2density(mean);
          unsigned int mateID = GetMaterialIndex( density );
          rflag = std::fwrite(&mateID, sizeof(unsigned int), 1, fileOut);
        }
      }


    } else {
      // density value is the average of a square region of
      // fCompression*fCompression pixels
      for(G4int ww = 0; ww < fRows ;ww += compSize ) {
        for(G4int xx = 0; xx < fColumns ;xx +=compSize ) {
          overflow = false;
          mean = 0;
          for(int sumx = 0; sumx < compSize; sumx++) {
            for(int sumy = 0; sumy < compSize; sumy++) {
              if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
              mean += fTab[ww+sumy][xx+sumx];
            }
            if(overflow) break;
          }
          mean /= compSize*compSize;
          
          if(!overflow) {
            density = Pixel2density(mean);
            unsigned int mateID = GetMaterialIndex( density );
            rflag = std::fwrite(&mateID, sizeof(unsigned int), 1, fileOut);
          }
        }
      }
    }

    //----- Write density for each pixel
    if(compSize == 1) { // no fCompression: each pixel has a density value)
      for( G4int ww = 0; ww < fRows; ww++) {
        for( G4int xx = 0; xx < fColumns; xx++) {
          mean = fTab[ww][xx];
          density = Pixel2density(mean);
          rflag = std::fwrite(&density, sizeof(G4float), 1, fileOut);
        }
      }
      
    } else {
      // density value is the average of a square region of
      // fCompression*fCompression pixels
      for(G4int ww = 0; ww < fRows ;ww += compSize ) {
        for(G4int xx = 0; xx < fColumns ;xx +=compSize ) {
          overflow = false;
          mean = 0;
          for(int sumx = 0; sumx < compSize; sumx++) {
            for(int sumy = 0; sumy < compSize; sumy++) {
              if(ww+sumy >= fRows || xx+sumx >= fColumns) overflow = true;
              mean += fTab[ww+sumy][xx+sumx];
            }
            if(overflow) break;
          }
          mean /= compSize*compSize;

          
          if(!overflow) {
            density = Pixel2density(mean);
            rflag = std::fwrite(&density, sizeof(G4float), 1, fileOut);
          }
        }
        
      }
      }
    
    rflag = std::fclose(fileOut);
    
    delete [] nameProcessed;

    /*    for ( G4int i = 0; i < fRows; i ++ ) {
      delete [] fTab[i];
    }
    delete [] fTab;
    */

    if (rflag) return returnvalue;
    return returnvalue;
}

/*
  G4int ExN01DicomHandler::displayImage(char command[300])
  {
  //   Display DICOM images using ImageMagick
  char commandName[500];
  std::sprintf(commandName,"display  %s",command);
  std::printf(commandName);
  G4int i = system(commandName);
  return (G4int )i;
  }

*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4float ExN01DicomHandler::Pixel2density(G4int pixel)
{
    G4float density = -1.;
    G4int nbrequali = 0;
    G4double deltaCT = 0;
    G4double deltaDensity = 0;

    // CT2Density.dat contains the calibration curve to convert CT (Hounsfield)
    // number to physical density
    std::ifstream calibration("DICOM/CT2Density.dat");
    calibration >> nbrequali;

    G4double * valuedensity = new G4double[nbrequali];
    G4double * valueCT = new G4double[nbrequali];

    if(!calibration) {
        G4cerr << "@@@ No value to transform pixels in density!" << G4endl;
        exit(1);

    } else { // calibration was successfully opened
        for(G4int i = 0; i < nbrequali; i++) { // Loop to store all the pts in CT2Density.dat
            calibration >> valueCT[i] >> valuedensity[i];
        }
    }
    calibration.close();

    for(G4int j = 1; j < nbrequali; j++) {
        if( pixel >= valueCT[j-1] && pixel < valueCT[j]) {

            deltaCT = valueCT[j] - valueCT[j-1];
            deltaDensity = valuedensity[j] - valuedensity[j-1];

            // interpolating linearly
            density = valuedensity[j] - ((valueCT[j] - pixel)*deltaDensity/deltaCT );

            break;
        }
    }

    if(density < 0.) {
        std::cout<<density<<std::endl;
        sleep(8);
        std::printf("@@@ Error density = %f && Pixel = %i (0x%x) && deltaDensity/deltaCT = %f\n",density,pixel,pixel, deltaDensity/deltaCT);
    }
    
    delete [] valuedensity;
    delete [] valueCT;

    return density;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01DicomHandler::CheckFileFormat()
{

    std::ifstream checkData("DICOM/Data.dat");
    char * oneLine = new char[128];


    if(!(checkData.is_open())) { //Check existance of Data.dat

        G4cout << "\nDicomG4 needs Data.dat :\n\tFirst line: number of image pixel for a "
               << "voxel (G4Box)\n\tSecond line: number of images (CT slices) to "
               << "read\n\tEach following line contains the name of a Dicom image except "
               << "for the .dcm extension\n";
        exit(0);
    }
    checkData >> fCompression;
    checkData >> fNFiles;
    G4String oneName;
    checkData.getline(oneLine,100);
    std::ifstream testExistence;
    G4bool existAlready = true;
    for(G4int rep = 0; rep < fNFiles; rep++) { 
      checkData.getline(oneLine,100);
      oneName = oneLine;
      oneName += ".g4dcm"; // create dicomFile.g4dcm
     // G4cout << fNFiles << " test file " << oneName << G4endl;
      chdir("DICOM/");
      testExistence.open(oneName.data());
      if(!(testExistence.is_open())) {
        existAlready = false;
        testExistence.clear();
        testExistence.close();
      }
      testExistence.clear();
      testExistence.close();
    }
    chdir("..");

    ReadMaterialIndices( checkData );

    checkData.close();
    delete [] oneLine;

    if( existAlready == false ) { // The files *.g4dcm have to be created

        G4cout << "\nAll the necessary images were not found in processed form, starting "
               << "with .dcm images\n";

        FILE * dicom;
        FILE * lecturePref;
        char * fCompressionc = new char[LINEBUFFSIZE];
        char * maxc = new char[LINEBUFFSIZE];
        //char name[300], inputFile[300];
        char * name = new char[FILENAMESIZE];
        char * inputFile = new char[FILENAMESIZE];
        G4int rflag;

        lecturePref = std::fopen("DICOM/Data.dat","r");
        rflag = std::fscanf(lecturePref,"%s",fCompressionc);
        fCompression = atoi(fCompressionc);
        rflag = std::fscanf(lecturePref,"%s",maxc);
        fNFiles = atoi(maxc);
        G4cout << " fNFiles " << fNFiles << G4endl;

        for( G4int i = 1; i <= fNFiles; i++ ) { // Begin loop on filenames

            rflag = std::fscanf(lecturePref,"%s",inputFile);
            std::sprintf(name,"DICOM/%s.dcm",inputFile);
            std::cout << "check 1: " << name << std::endl;
            //  Open input file and give it to gestion_dicom :
            std::printf("### Opening %s and reading :\n",name);
            dicom = std::fopen(name,"rb");
            // Reading the .dcm in two steps:
            //      1.  reading the header
            //        2. reading the pixel data and store the density in Moyenne.dat
            if( dicom != 0 ) {
                ReadFile(dicom,inputFile);
            } else {
                G4cout << "\nError opening file : " << name << G4endl;
            }
            rflag = std::fclose(dicom);
        }
        rflag = std::fclose(lecturePref);

        delete [] fCompressionc;
        delete [] maxc;
        delete [] name;
        delete [] inputFile;
        if (rflag) return;
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <class Type>
void ExN01DicomHandler::GetValue(char * _val, Type & _rval) {

#if BYTE_ORDER == BIG_ENDIAN
    if(fLittleEndian) {      // little endian
#else // BYTE_ORDER == LITTLE_ENDIAN
    if(!fLittleEndian) {     // big endian
#endif
        const int SIZE = sizeof(_rval);
        char ctemp;
        for(int i = 0; i < SIZE/2; i++) {
            ctemp = _val[i];
            _val[i] = _val[SIZE - 1 - i];
            _val[SIZE - 1 - i] = ctemp;
        }
    }
    _rval = *(Type *)_val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int ExN01DicomHandler::read_defined_nested(FILE * nested,G4int SQ_Length)
{ 
  //      VARIABLES
  unsigned short item_GroupNumber;
  unsigned short item_ElementNumber;
  G4int item_Length;
  G4int items_array_length=0;
  char * buffer= new char[LINEBUFFSIZE];
  size_t rflag = 0;

  while(items_array_length < SQ_Length)
  {
   rflag = std::fread(buffer, 2, 1, nested);
   GetValue(buffer, item_GroupNumber);
   
   rflag = std::fread(buffer, 2, 1, nested);
   GetValue(buffer, item_ElementNumber);
   
   rflag = std::fread(buffer, 4, 1, nested);
   GetValue(buffer, item_Length);
   
   rflag = std::fread(buffer, item_Length, 1, nested);
   
   items_array_length= items_array_length+8+item_Length;
  }
 
  delete [] buffer;
  
  if( SQ_Length>items_array_length )
   return 0;
  else
   return 1;
  if (rflag) return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01DicomHandler::read_undefined_nested(FILE * nested)
{
  //      VARIABLES
  unsigned short item_GroupNumber;
  unsigned short item_ElementNumber;
  unsigned int item_Length;
  char * buffer= new char[LINEBUFFSIZE];
  size_t rflag = 0;

  do
  {
   rflag = std::fread(buffer, 2, 1, nested);
   GetValue(buffer, item_GroupNumber);
   
   rflag = std::fread(buffer, 2, 1, nested);
   GetValue(buffer, item_ElementNumber);
   
   rflag = std::fread(buffer, 4, 1, nested);
   GetValue(buffer, item_Length);
   
   if(item_Length!=0xffffffff)
    rflag = std::fread(buffer, item_Length, 1, nested);
   else
    read_undefined_item(nested);
   
   
  } while(item_GroupNumber!=0xFFFE || item_ElementNumber!=0xE0DD || item_Length!=0);

  delete [] buffer;
  if (rflag) return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01DicomHandler::read_undefined_item(FILE * nested)
{
  //      VARIABLES
 unsigned short item_GroupNumber;
 unsigned short item_ElementNumber;
 G4int item_Length; size_t rflag = 0;
 char *buffer= new char[LINEBUFFSIZE];
 
 do
 {
  rflag = std::fread(buffer, 2, 1, nested);
  GetValue(buffer, item_GroupNumber);
   
  rflag = std::fread(buffer, 2, 1, nested);
  GetValue(buffer, item_ElementNumber);
   
  rflag = std::fread(buffer, 4, 1, nested);
  GetValue(buffer, item_Length);


  if(item_Length!=0)
   rflag = std::fread(buffer,item_Length,1,nested);

 }
 while(item_GroupNumber!=0xFFFE || item_ElementNumber!=0xE00D || item_Length!=0);
 
 delete [] buffer;
 if (rflag) return;
}
