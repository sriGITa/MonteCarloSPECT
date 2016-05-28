cd Geant2txt_float/
cp ../Simulation_data/total-energy.txt .
cp ../Simulation_data/total-win-energy.txt .
make clean
make
./Geant2txt_float total-win-energy.txt
source merge-sino.sh
source merge-sino2.sh
cd ..
cd Geant2txt_int/

cp ../Simulation_data/total-counts.txt .
cp ../Simulation_data/total-win-counts.txt .
make clean
make
./Geant2txt_int total-win-counts.txt
source merge-sino.sh
source merge-sino2.sh
cd ..


