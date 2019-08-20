// NOAA POES Monte Carlo Simulation v1.3, 17/09/2008e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_BField.hh"
// Created 01/29/2008 by Karl Yando, modified 04/14/2008 to v1.0

#ifndef MEPED_BField_H
#define MEPED_BField_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"
class MEPED_BMessenger;

class MEPED_BField : public G4MagneticField
{
  public:
	//MEPED_BField();
	
	MEPED_BField(const G4ThreeVector& FieldVector);
	
	~MEPED_BField();
	virtual void SetFieldValue(const G4ThreeVector newFieldVector);
	virtual void GetFieldValue(const double Point[3],
					double *Bfield ) const;
		// Actually queries the field

  private:
	MEPED_BMessenger* messenger;
	G4double Bx;
	G4double By;
 	G4double Bz;
	G4double rmax_sq;
	G4double zmax;
	G4double xCenter;
	G4double yCenter;
	G4double zCenter;

  public:
	inline G4double GetFieldStrength() const { return Bz; }
		// Only here for us
};

#endif
