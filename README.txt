=========================================================================================
Monte Carlo simulation toolkit to evaluate collimator design for energy-independent SPECT
=========================================================================================
 
Introduction: 
————————————
This application was developed to evaluate gamma camera performance for low and medium energy radionuclides with the goal of developing an energy-independent camera. The simulation calculates energy deposited and number of photon counts recorded in the detector. 
 
Input Variables: 
———————————————
- Distance from source to collimator 
- Length of collimator hole 
- Collimator material 
- Detector material 
- Phantom type 
- Collimator hole shape 
- Radius of collimator hole 
 
Phantoms: 
————————
- Single point source 
- Three point sources 
- Cylinder 
- Hoffman 
- Jaszczak 
 
How to Run: 
——————————
1. Compile after activating Geant4 environment variables (refer to Geant4 documentation: http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/InstallationGuide/html/ch03s02.html) 
2. Execute the script run_initial.sh to create 12 directories: 
	source run_initial.sh 
3. Then execute the script run_main.sh: 
	source run_main.sh 
4. Once the run is complete and the results folder is generated, execute the instructions in the scripts merge-all.sh and merge-total.sh: 
	source merge-all.sh and
	source merge-total.sh 
5. Execute run_clean.sh to initialize and run a new simulation: 
	source run_clean.sh 
 
Output files: 
————————————
- Projections of counts are stored in /Geant4txt_int/data and corresponding sinograms in /Geant4txt_int/sino 
- Projections of deposited energy are in /Geant4txt_float/data and the corresponding sinograms in /Geant4txt_float/sino 
- Total values for all projection angles for both counts and deposited energy are in data files names 'total' 
 
 
 
 
 