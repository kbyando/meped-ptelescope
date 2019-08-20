// NOAA POES Monte Carlo Simulation v1.4, 08/04/2010e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2010
//
// "MEPED_EventAction.hh"

#ifndef MEPED_EventAction_h
#define MEPED_EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class MEPED_RunAction;
class MEPED_EventActionMessenger;

class MEPED_EventAction : public G4UserEventAction
{
  public:
    MEPED_EventAction(MEPED_RunAction*);
   ~MEPED_EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

    // method "StoreE()": works in conjunction with MEPED_SteppingAction.cc and
    // MEPED_EventAction.cc to store initial energy of particle (taken at step 0)
    void StoreE(G4double etot) {eTotal += etot;};

    // methods to receive initial position and momentum data (v1.2+)
    void StoreX(G4double dx, G4double dy, G4double dz) 
	{x = dx; y = dy; z = dz;};
    void StoreP(G4double dPx, G4double dPy, G4double dPz)
	{px = dPx; py = dPy; pz = dPz;};

    // Tally deposited energies
    void AddDete1(G4double de, G4double dl) {EnergyDete1 += de; TrackLDete1 +=dl;};
    void AddDete2(G4double de, G4double dl) {EnergyDete2 += de; TrackLDete2 +=dl;};

    void SetPrintModulo(G4int	val) {printModulo = val;};

  private:
    MEPED_RunAction* runAct;

    G4double eTotal;
    G4double x, y, z, px, py, pz;

    // Energy and Tracklength Collection for D1 & D2 (active areas only)
    G4double EnergyDete1, TrackLDete1;
    G4double EnergyDete2, TrackLDete2;

    G4int	printModulo;

    MEPED_EventActionMessenger* eventMessenger;
};

#endif
