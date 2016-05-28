#include "G4VHit.hh"
class MyHit: public G4VHit{
         public:
           MyHit();
         virtual ~MyHit();
          inline void SetEnergyDeposit(double energy) 
           { energyDeposit= energy; }
          inline double GetEnergyDeposit() { return energyDeposit;}
         private:
	 G4double energyDeposit; // for recording energy
};
