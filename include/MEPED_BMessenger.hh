// NOAA POES Monte Carlo Simulation v1.3, 17/09/2008e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_BMessenger.hh"
// Created 01/29/2008 by Karl Yando, modified 04/14/2008 for v1.0

#ifndef MEPED_BMessenger_h
#define MEPED_BMessenger_h 1

class MEPED_BField;
class G4UIcmdWith3VectorAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class MEPED_BMessenger: public G4UImessenger
{
  public:
	MEPED_BMessenger(MEPED_BField* mpga);
	~MEPED_BMessenger();

  public:
	void SetNewValue(G4UIcommand * command,G4String newValues);
	G4String GetCurrentValue(G4UIcommand * command);

  private:
	MEPED_BField * target;

	G4UIcmdWith3VectorAndUnit*	fieldCmd;
};

#endif
