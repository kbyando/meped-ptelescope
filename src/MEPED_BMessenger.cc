// NOAA POES Monte Carlo Simulation v1.3, 17/09/2008e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_BMessenger.cc"
// Created 01/29/2008 by Karl Yando, modified 04/14/2008 for v1.0
// Provides for dynamic manipulation of the module's magnetic field

#include "MEPED_BMessenger.hh"
#include "MEPED_BField.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

MEPED_BMessenger::MEPED_BMessenger(MEPED_BField * mpga)
:target(mpga)
{
  fieldCmd = new G4UIcmdWith3VectorAndUnit("/MEPED/Bfield/setFieldVector",this);
  fieldCmd->SetGuidance("Set Magnetic Field Strength, in Tesla.");
  fieldCmd->SetGuidance("3-Vector Field of the Sm-Co Sweep Magnet.");
  fieldCmd->SetParameterName("xcomp","ycomp","zcomp",true);
  fieldCmd->SetDefaultValue( G4ThreeVector(0.,0.,0.) ); // default 0.2 T
  fieldCmd->SetDefaultUnit("tesla"); // all parameters set elsewhere
  fieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

MEPED_BMessenger::~MEPED_BMessenger()
{
  delete fieldCmd;
}

void MEPED_BMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if( command==fieldCmd )
  { target->SetFieldValue(fieldCmd->GetNew3VectorValue(newValues)); }
}

G4String MEPED_BMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if(command==fieldCmd)
  { cv = fieldCmd->ConvertToString(target->GetFieldStrength(),"tesla"); }

  return cv;
}


