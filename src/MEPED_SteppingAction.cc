// "MEPED_SteppingAction.cc" 22/11/2007
// NOAA POES Monte Carlo Simulation v0.2.1
// MEPED (Medium Energy Proton and Electron Detector)
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2007
//

#include "MEPED_SteppingAction.hh"

#include "MEPED_DetectorConstruction.hh"
#include "MEPED_EventAction.hh"

#include "G4Step.hh"

MEPED_SteppingAction::MEPED_SteppingAction(MEPED_DetectorConstruction* det, MEPED_EventAction* evt)
:detector(det), eventaction(evt)
{ }

MEPED_SteppingAction::~MEPED_SteppingAction()
{ }

// Collect Relevant Data for Output
void MEPED_SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // Get Volume of the Current Step
  G4VPhysicalVolume* volume
  = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // Collect Energy and Track Length Step by Step
  G4double edep = aStep->GetTotalEnergyDeposit();

  // Data read from stepping dynamics and passed to EventAction.cc
  if (aStep->GetPreStepPoint()->GetGlobalTime() == 0.) {

   G4double eTotal = aStep->GetPreStepPoint()->GetKineticEnergy();
   eventaction->StoreE(eTotal);
   
   G4double x, y, z, px, py, pz;   // position and momentum
    x = aStep->GetPreStepPoint()->GetPosition().x();
    y = aStep->GetPreStepPoint()->GetPosition().y();
    z = aStep->GetPreStepPoint()->GetPosition().z();

    px = aStep->GetPreStepPoint()->GetMomentumDirection().x(); 
    py = aStep->GetPreStepPoint()->GetMomentumDirection().y();
    pz = aStep->GetPreStepPoint()->GetMomentumDirection().z();
    eventaction->StoreX(x, y, z); // Stores position data
    eventaction->StoreP(px, py, pz); // Stores momentum data
  } 

  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();

  // Declare "detector" volumes (recognized via reference to MEPED_DetectorC[...])
  if (volume == detector->GetDete1()) eventaction->AddDete1(edep,stepl);
  if (volume == detector->GetDete2()) eventaction->AddDete2(edep,stepl);
}
