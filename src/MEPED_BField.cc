// NOAA POES Monte Carlo Simulation v1.3, 17/09/2008e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_BField.cc"
//
// Defines a Localized (Uniform) Magnetic Field.
// Created 01/29/2008 Karl Yando, modified 04/14/2008 to v1.0

#include "MEPED_BField.hh"
#include "MEPED_BMessenger.hh"

MEPED_BField::MEPED_BField(const G4ThreeVector& FieldVector)
{
  messenger = new MEPED_BMessenger(this); // for on-the-fly changes
  Bx = FieldVector.x();
  By = FieldVector.y();
  Bz = FieldVector.z();

// Magnetic Field Dimensions
  rmax_sq = sqr( 0.250*25.4*mm ); // Squares the field radius (default 0.250")
  zmax = (0.200)*25.4*mm; // Half-depth of field area (default 0.200")

// Placement of the Magnetic Field (in Absolute Coordinates)
  xCenter = -(1.975/2-0.062-0.250-0.522/2)*2.54*cm; // revised 4/15/2008
  yCenter = 0.;  
  zCenter = 0.;  
}  

void MEPED_BField::SetFieldValue(const G4ThreeVector newFieldVector)
{
  Bx = newFieldVector.x();
  By = newFieldVector.y();
  Bz = newFieldVector.z();
}

MEPED_BField::~MEPED_BField()
{ delete messenger; }

void MEPED_BField::GetFieldValue(const double Point[3],double *Bfield) const
{
//  Trig to determine whether our particle is affected
// 	(following the logic r2=x2+y2, such that our z component falls along the magnetic axis)
// 	Note: this orients the default WIRED view for viewing of deflection,
// 		*not* like the MEPED schematics
  if ( std::abs(Point[2] - zCenter) < zmax &&
	(sqr(Point[0] - xCenter) + sqr(Point[1] - yCenter)) < rmax_sq )
  { Bfield[0] = Bx; // x-comp
    Bfield[1] = By; // y-comp
    Bfield[2] = Bz; // z-comp
  } 
  else // alternative components (outside of field)
  { Bfield[0] = 0.;
    Bfield[1] = 0.;
    Bfield[2] = 0.;
  }  
} // a return value is evidently not needed with this
