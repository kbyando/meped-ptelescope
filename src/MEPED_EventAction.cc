// NOAA POES Monte Carlo Simulation v1.4, 08/04/2008e
// MEPED - Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_EventAction.cc"

#include "MEPED_EventAction.hh"

#include "MEPED_RunAction.hh"
#include "MEPED_EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

MEPED_EventAction::MEPED_EventAction(MEPED_RunAction* run)
:runAct(run),printModulo(1),eventMessenger(0)
{
  eventMessenger = new MEPED_EventActionMessenger(this);
}

MEPED_EventAction::~MEPED_EventAction()
{
  delete eventMessenger;
}

void MEPED_EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) {
   // CLHEP::HepRandom::showEngineStatus();
 }

 // Event Initialization
  EnergyDete1 = EnergyDete2 = 0.;	// energy deposited in Dete1/Dete2
  TrackLDete1 = TrackLDete2 = 0.;	// length of particle track in ""

  eTotal = 0.;				// [incident] particle energy

  // Note: these values are updated from SteppingAction using the methods 
  //   StoreE(), AddDete1(), and AddDete2() defined in EventAction.hh
}

void MEPED_EventAction::EndOfEventAction(const G4Event* evt)
{
  // Accumulate Statistics
  runAct->fillPerEvent(EnergyDete1, EnergyDete2, TrackLDete1, TrackLDete2);

  // Number Per Event (modulo n)
  G4int evtNb = evt->GetEventID();

  // Output simulation data 
  if (evtNb%printModulo == 0) {

//    if (evtNb == 0) {				// if the first event, then..
//      G4RunManager::GetRunManager()->rndmSaveThisEvent(); // save seed
//    }

    if ((EnergyDete1 > 0.) || (EnergyDete2 > 0.)) {

      // hits were scored: report event information
      runAct->incrHitCount();  

      G4int prec = G4cout.precision(8);

      G4cout
       << evtNb << ", "		// eventID
       << x / mm << ", "  	// x-coordinate of position, in mm
       << y / mm << ", "  	// y-coordinate of position, in mm
       << z / mm << ", "  	// z-coordinate of position, in mm
       << px << ", " 		// x-direction of momentum
       << py << ", " 		// y-direction of momentum
       << pz << ", " 		// z-direction of momentum
       << eTotal / keV << ", "	// incident energy, in keV
       << EnergyDete1 / keV << ", "	// energy deposited in D1, in keV
       << EnergyDete2 / keV	// energy deposited in D2, in keV

       << G4endl;
    }
    else {
      // no hits were scored
      runAct->incrNullCount();
    }
  }
} 

