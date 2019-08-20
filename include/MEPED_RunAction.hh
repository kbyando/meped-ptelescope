// NOAA POES Monte Carlo Simulation v1.3, 15/09/2008e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_RunAction.hh"

#ifndef MEPED_RunAction_h
#define MEPED_RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class MEPED_RunAction : public G4UserRunAction
{
  public:
    MEPED_RunAction();
   ~MEPED_RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

    void fillPerEvent(G4double, G4double, G4double, G4double);
    void incrHitCount(); 
    void incrNullCount();

  private:
    G4double sumEActiveD1, sum2EActiveD1;
    G4double sumEActiveD2, sum2EActiveD2;

    G4double sumLActiveD1, sum2LActiveD1;
    G4double sumLActiveD2, sum2LActiveD2;

    G4int hitCount, nullCount;
};

#endif
