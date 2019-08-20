// NOAA POES Monte Carlo Simulation v1.3, 22/10/2008e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_DetectorConstruction.cc"

// Our Geometry - MEPED Proton Telescope Assembly (ptAssembly)
// I. Physical Assembly (ShellStruct)
// II. Collimator & Sm-Co Magnets
//  1. Collimator 
//   A. Outer Collimator (OuterCollimator)
//   B. Inner Cone of Outer Collimator (OCInnerCone)
//  2. Magnetic Field Chamber (BChamber)
//   A. Sm-Co Magnet (1)
//   B. Sm-Co Magnet (2)
//   C. Fe Magnet Return Yoke
//   D. Baffles
//	i. Baffle (1)
//	ii. Baffle (2)
//	iii. Baffle (3)
// III. Misc Shielding
//   1. Front Cover
//   2. Partition
//   3. Back Cover
// IV. Detector Housing
//  1. Front Shield / Innermost Collimator
//  2. Inner Housing (InnerHousing)
//   A. Detector Stack
//      i. Belville Washer (1) (Belville1)
//      ii. Contact (1)
//      iii. Detector 1
//	  a. Ring Mount (front part)
// 	  b. Ring Mount (back part) 
//	  c. Al Film Deposit (SurfDeposit)       
//	  d. D1 Active Layer
//      iv. Contact (2)
//      v. Contact (3)
//      vi. Detector 2
//	  a. Ring Mount (front part)
//	  b. Ring Mount (back part)
//	  c. Au Film Deposit (SurfDeposit)
//	  d. D2 Active Layer
//      vii. Contact (4)
//      viii. Belville Washer (2) (Belville2)
//   B. Detector Insulator
//    n/a
//      i. Front
//      ii. Longitudinal
//      iii. Back
//   C. Housing Shield
//  3. Back Shield
// 

#include "MEPED_DetectorConstruction.hh"
#include "MEPED_DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "MEPED_BField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "math.h" // Added 01/13/2007 for direct calc of parameters

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MEPED_DetectorConstruction::MEPED_DetectorConstruction()
:DetectorMaterial(0),ShieldingMaterial(0),HousingMaterial(0),ContactMaterial(),BelvilleMaterial(),InsulatorMaterial(),defaultMaterial(),
 solidWorld(0),logicWorld(0),physiWorld(0),

 solid_ptAssembly(0),logic_ptAssembly(0),physi_ptAssembly(0),
// I. Physical Assembly (ShellStruct)
 solidShellStruct(0),logicShellStruct(0),physiShellStruct(0),
// II. Collimator & Sm-Co Magnets
 solidCollimatorChamber(0),logicCollimatorChamber(0),physiCollimatorChamber(0),
 solidCollimator(0),logicCollimator(0),physiCollimator(0),
 solidOCOuterCone(0),logicOCOuterCone(0),physiOCOuterCone(0),
 solidOCInnerCone(0),logicOCInnerCone(0),physiOCInnerCone(0),
 solidBChamber(0),logicBChamber(0),physiBChamber(0),
	solidSmCoMagnet1(0),logicSmCoMagnet1(0),physiSmCoMagnet1(0),
	solidSmCoMagnet2(0),logicSmCoMagnet2(0),physiSmCoMagnet2(0),
	solidReturnYokeA(0),logicReturnYokeA(0),physiReturnYokeA(0),
	solidReturnYokeB(0),logicReturnYokeB(0),physiReturnYokeB(0),
	solidReturnYokeCD(0),logicReturnYokeC(0),physiReturnYokeC(0),
			logicReturnYokeD(0),physiReturnYokeD(0),
	solidAlBaffle1(0),logicAlBaffle1(0),physiAlBaffle1(0),
	solidAlBaffle2(0),logicAlBaffle2(0),physiAlBaffle2(0),
	solidAlBaffle3(0),logicAlBaffle3(0),physiAlBaffle3(0),
// III. Misc Shielding
	solidFrontCover(0),logicFrontCover(0),physiFrontCover(0),
	solidPartition(0),logicPartition(0),physiPartition(0),
	solidRearCover(0),logicRearCover(0),physiRearCover(0),
// IV. Detector Housing
 solidDeteHousing(0),logicDeteHousing(0),physiDeteHousing(0),
 solidInnerCollimator(0),logicInnerCollimator(0),physiInnerCollimator(0),
 solidInnerHousing(0),logicInnerHousing(0),physiInnerHousing(0),
//   2A. Detector Stack
solidDStack(0),logicDStack(0),physiDStack(0),
 solidBelville1(0),logicBelville1(0),physiBelville1(0),
 solidContact1(0),logicContact1(0),physiContact1(0),
 solid_Detector1(0),logic_Detector1(0),physi_Detector1(0),
  solid_D1RingMntA(0),logic_D1RingMntA(0),physi_D1RingMntA(0),
  solid_D1RingMntB(0),logic_D1RingMntB(0),physi_D1RingMntB(0),
  solid_D1SurfDeposit(0),logic_D1SurfDeposit(0),physi_D1SurfDeposit(0),
  solid_D1ActiveLayer(0),logic_D1ActiveLayer(0),physi_D1ActiveLayer(0),
 solidContact2(0),logicContact2(0),physiContact2(0),
 solidContact3(0),logicContact3(0),physiContact3(0),
 solid_Detector2(0),logic_Detector2(0),physi_Detector2(0),
  solid_D2RingMntA(0),logic_D2RingMntA(0),physi_D2RingMntA(0),
  solid_D2RingMntB(0),logic_D2RingMntB(0),physi_D2RingMntB(0),
  solid_D2SurfDeposit(0),logic_D2SurfDeposit(0),physi_D2SurfDeposit(0),
  solid_D2ActiveLayer(0),logic_D2ActiveLayer(0),physi_D2ActiveLayer(0),
 solidContact4(0),logicContact4(0),physiContact4(0),
 solidBelville2(0),logicBelville2(0),physiBelville2(0),
//   2B. Detector Insulator
 solidDInsulator(0),logicDInsulator(0),physiDInsulator(0),
  solidFrontInsul(0),logicFrontInsul(0),physiFrontInsul(0),
  solidLongiInsul(0),logicLongiInsul(0),physiLongiInsul(0),
  solidBackInsul(0), logicBackInsul(0), physiBackInsul(0),
//   2C. Housing Shield
 solidHouseShielding(0),logicHouseShielding(0),physiHouseShielding(0),
//   3. Back Shield
 solidBackShield(0),logicBackShield(0),physiBackShield(0),
// Introducing Shielding and Detector Subunit 
 ptAssemblyAngle(90.*deg),magField(0)

{
  // Detector Characteristics - converted from inches
  DetectorRadius 	= (0.598/2)*2.54*cm; // physical detector
  DetectorThickness   	= 0.146*2.54*cm; // physical detector
  DActiveLayer		= 0.2*mm; // 200 microns
  D1Silicon_rad		= sqrt(25*mm2/3.14159); // A = pi*r^2,
  D2Silicon_rad		= sqrt(50*mm2/3.14159); // so r = sqrt(A/pi)

  // External Specifications for the Proton Telescope Assembly
  ptAssembly_rad 	= (37.8/2)*mm; // from 2004 rev
  ptAssembly_length	= 1.975*2.54*cm; // from orig. spec

  // Miscellaneous Components
  outerWall_thickness	= (0.125)*2.54*cm; // An estimate (1)
  partition_thickness	= (0.077)*2.54*cm; // revised spec (center partition)
  shield_thickness	= (0.250)*2.54*cm; // from orig. spec
  frontCover_thickness	= (0.062)*2.54*cm; // 	"
  rearCover_thickness	= (0.082)*2.54*cm; // 	"

  // Detector Housing
  DeteHousing_length	= (0.730 + 0.250)*2.54*cm; // revised spec
  DeteHousing_dia	= (2*ptAssembly_rad - 2*outerWall_thickness); // (1A)
  inner_DH_dia		= (DeteHousing_dia - 2*shield_thickness); // (1B)
  
  // Inner Housing
   // Element Thicknesses
  BelvilleThickness	= 0.012*2.54*cm; // orig. spec
  ContactThickness	= 0.015*2.54*cm; // "
  ContactDeposit	= 0.020*mg/cm2; // Pseudolength (20*µg/m2)
  AlDensity		= 2.700*g/cm3; // Contact Density (Aluminium)
   // Element Radii
  BelvilleInnerRadius	= (0.422/2)*2.54*cm; // revised spec, estimate (2)
  BelvilleOuterRadius	= DetectorRadius;
  ContactInnerRadius	= BelvilleInnerRadius; 	// estimate (2A)
  ContactOuterRadius	= DetectorRadius;
  // Misc Elements
  InsulatorThickness	= 0.050*2.54*cm;

  // Collimator Characteristics
  fullWidthAngle 	= 30.*deg; // from orig. spec
  fullWidthICAngle	= 13.*deg; // from orig. spec

  DIA1			= (0.562)*2.54*cm; // revised spec

  DIA2inner		= (0.462)*2.54*cm; // 	" " (OC_InnerConeF_rad)
  DIA2outer		= (0.874)*2.54*cm; //	" " (OC_OuterConeF_rad)

  DIA3inner		= (0.328)*2.54*cm; // 	" " (OC_InnerConeB_rad)
  DIA3outer		= (0.646)*2.54*cm; // 	" " (OC_OuterConeB_rad)

  DIA4			= (0.207)*2.54*cm; // 	" " (Structural Aperture)
  DIA5			= (0.188)*2.54*cm; // 	" " (IC / Front Shield)
  DIA6			= (0.130)*2.54*cm; // 	" " (IC / " ")

  AlColl1_dia		= (0.296)*2.54*cm; // 	"
  AlColl2_dia		= (0.265)*2.54*cm; // 	"
  AlColl3_dia		= (0.233)*2.54*cm; // 	"
  AlColl_odia		= (0.400)*2.54*cm; // Fe Ret. Yoke manu. spec

  AlColl_spacing	= (0.115)*2.54*cm; //	"
  AlColl_thickness	= (0.020)*2.54*cm; // 	"

  CollimatorSep		= (0.850)*2.54*cm; // from orig. spec
  CollimatorDia		= DeteHousing_dia; // (1C)
  CollimatorHeight	= (shield_thickness + CollimatorSep 
		- shield_thickness - partition_thickness); 

  // Magnetic Trap Elements
  SmCoMagnet_thickness	= (0.250)*2.54*cm; // guess (measured)
  SmCoMagnet_height	= (0.500)*2.54*cm; // guess (a diameter)

  FeReturnYoke_idia	= (1.040)*2.54*cm; // manufacturing spec
  FeReturnYoke_odia	= (1.220)*2.54*cm; // (trimmed from manufacturing spec of 0.1240")
  FeReturnYoke_length	= (0.520)*2.54*cm; // revised spec

  ComputeModuleParameters();
  
  // Defining Rotation Matrices 
  AssemblyRotation = new G4RotationMatrix();
  AssemblyRotation->rotateY(ptAssemblyAngle);

  MagTrapRotation = new G4RotationMatrix();
  MagTrapRotation->rotateY(90.*deg);
  // materials
  DefineMaterials();
  SetDetectorMaterial("Silicon");
  SetShieldingMaterial("Tungsten");
  SetHousingMaterial("Aluminium");
 
  SetContactMaterial("Brass"); // Detector Contacts
  SetD1CoatingMaterial("Aluminium");
  SetD2CoatingMaterial("Aluminium"); // Al in the original spec
  SetEpoxyMaterial("Kel-F"); // Not quite the case, but close

  SetBelvilleMaterial("Aluminium"); // Belville Washers
  SetInsulatorMaterial("Kel-F"); 

  SetMagnetMaterial("Sm-Co");
  SetReturnYokeMaterial("Iron");

  // create commands for interactive definition of the calorimeter
  detectorMessenger = new MEPED_DetectorMessenger(this);

  // Creating a Magnetic Field Object
  magneticField = new MEPED_BField(G4ThreeVector(0.,0.,0.2*tesla) ); // in absolute coordinates
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MEPED_DetectorConstruction::~MEPED_DetectorConstruction()
{ delete magneticField;
  delete detectorMessenger;
  delete AssemblyRotation;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* MEPED_DetectorConstruction::Construct()
{
  return ConstructModule(); // Create our Module's Geometry

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MEPED_DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String symbol;             //a=mass of a mole;
G4double a, z, density;      //z=mean number of protons;  
// G4int iz, n;                 //iz=number of protons  in an isotope; 
//                             // n=number of nucleons in an isotope;

G4int ncomponents, natoms;
G4double fractionmass;

// define Elements (metals described as "Materials")

G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);

G4Element* F  = new G4Element("Fluorine",symbol="F" , z= 9., a= 19.00*g/mole);
G4Element* Cl = new G4Element("Chlorine",symbol="Cl", z=17., a= 35.45*g/mole);

G4Element* Cu = new G4Element("Copper"  ,symbol="Cu", z=29., a=63.55*g/mole);
G4Element* Zn = new G4Element("Zinc"	,symbol="Zn", z=30., a=65.38*g/mole);

G4Element* Sm = new G4Element("Samarium",symbol="Sm", z=62., a=150.4*g/mole);
G4Element* Co = new G4Element("Cobalt"	,symbol="Co", z=27., a=58.93*g/mole);
	// References: 	http://periodic.lanl.gov/elements/9.html  (F)
	//		http://periodic.lanl.gov/elements/17.html (Cl)
	//		http://periodic.lanl.gov/elements/29.html (Cu)
	//		http://periodic.lanl.gov/elements/30.html (Zn)
	//		http://periodic.lanl.gov/elements/62.html (Sm)
	//		http://periodic.lanl.gov/elements/27.html (Co)

// define Simple Materials

new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.699*g/cm3);
new G4Material("Silicon", z= 14., a= 28.09*g/mole, density= 2.33*g/cm3);
new G4Material("Tungsten", z=74., a= 183.84*g/mole, density= 19.3*g/cm3);

new G4Material("Iron", z=26., a= 55.85*g/mole, density= 7.874*g/cm3);
new G4Material("Gold", z=79., a= 196.967*g/mole, density= 19.32*g/cm3);
	// References: 	http://periodic.lanl.gov/elements/26.html (Fe)
	//		http://periodic.lanl.gov/elements/79.html (Au)
	//	Densities from NIST STAR material composition database
	//	http://physics.nist.gov/PhysRefData/Star/Text/contents.html
// define a material from elements.   case 1: chemical molecule

G4Material* KelF =
new G4Material("Kel-F",density= 2.13*g/cm3, ncomponents=3);
KelF->AddElement(C, natoms=2);
KelF->AddElement(F, natoms=3);
KelF->AddElement(Cl, natoms=1);

G4Material* SmCo = 
new G4Material("Sm-Co", density= 9.3*g/cm3, ncomponents=2);
SmCo->AddElement(Sm, natoms=2);
SmCo->AddElement(Co, natoms=17);

// define a material from elements.   case 2: mixture by fractional mass

  G4Material* Air = 
  new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* Brass =
  new G4Material("Brass", density= 8.53*g/cm3, ncomponents=2);
  Brass->AddElement(Cu, fractionmass=0.7);
  Brass->AddElement(Zn, fractionmass=0.3);

// examples of vacuum

G4Material* Vacuum =
new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                           kStateGas, 2.73*kelvin, 3.e-18*pascal);

// Note: Following line commented out.  Will output a mess of geometry
// characteristics (mainly the properties of the defined materials).  Useful
// for easy reading.
//    
// G4cout << *(G4Material::GetMaterialTable()) << G4endl;

// A final material assignment
  defaultMaterial  = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* MEPED_DetectorConstruction::ConstructModule()

{
  // Clean old geometry, if any
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  ComputeModuleParameters();

  // Localized Magetic Field (BChamber)
  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized)
  {
	G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
	fieldMgr->SetDetectorField(magneticField);
	fieldMgr->CreateChordFinder(magneticField);
	fieldIsInitialized = true;
  }
   
  // World
  solidWorld = new G4Box("World",				//its name
                   1.*m,1.*m,1.*m);	//its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 logicWorld,		//its logical volume				 
                                 "World",		//its name
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

//
// Proton Telescope Assembly
//
 solid_ptAssembly=0;logic_ptAssembly=0; physi_ptAssembly=0;

if (ptAssembly_length > 0.)
  { solid_ptAssembly = new G4Tubs("Proton Telescope Assembly",0.*mm,
	ptAssembly_rad,ptAssembly_length/2,0.*deg,360.*deg);
    logic_ptAssembly = new G4LogicalVolume(solid_ptAssembly,defaultMaterial,
	"Proton Telescope Assembly");
    physi_ptAssembly = new G4PVPlacement(AssemblyRotation,
	G4ThreeVector(0.,0.,0.),logic_ptAssembly,"Proton Telescope Assembly",
    logicWorld,false,0);
  }

//
// Physical Assembly (Shell Structure)
//
 solidShellStruct=0; logicShellStruct=0; physiShellStruct=0;

if (ptAssembly_length > 0.)
  { solidShellStruct = new G4Tubs("(Shell) Assembly",DeteHousing_dia/2,
	ptAssembly_rad,(ptAssembly_length - (0.062 + 0.082)*25.4*mm)/2,
	0.*deg,360.*deg);
    logicShellStruct = new G4LogicalVolume(solidShellStruct,HousingMaterial,
	"(Shell) Assembly");
    physiShellStruct = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.01*25.4*mm),
	logicShellStruct,"(Shell) Assembly",logic_ptAssembly,false,0);
  }


//
// Collimator Chamber
// 
 solidCollimatorChamber=0;logicCollimatorChamber=0;physiCollimatorChamber=0;
if (ptAssembly_length > 0.)
  { G4double zdepth = (ptAssembly_length/2 - frontCover_thickness 
		- CollimatorHeight/2);
     solidCollimatorChamber = new G4Tubs("Collimator Chamber",0.*mm, CollimatorDia/2,
	(CollimatorHeight)/2,0.*deg,360.*deg);
    logicCollimatorChamber = new G4LogicalVolume(solidCollimatorChamber,
	defaultMaterial,"Collimator Chamber");
    physiCollimatorChamber = new G4PVPlacement(0,G4ThreeVector(0.,0.,zdepth),
	logicCollimatorChamber,"Collimator Chamber",logic_ptAssembly,false,0);
  }

//
// Collimator
//
 solidCollimator=0;logicCollimator=0;physiCollimator=0;

if(ptAssembly_length > 0.)
  { G4double zdepth = (CollimatorHeight/2 - shield_thickness/2);
    solidCollimator = new G4Cons("(Outer) Collimator",DIA3inner/2,
	CollimatorDia/2, DIA2inner/2, CollimatorDia/2, 
	shield_thickness/2,0.*deg,360.*deg);
    logicCollimator = new G4LogicalVolume(solidCollimator,
	defaultMaterial,"(Outer) Collimator");
    physiCollimator = new G4PVPlacement(0,G4ThreeVector(0.,0.,zdepth),
	logicCollimator,"(Outer) Collimator",logicCollimatorChamber,false,0);
  }

//
// Outer Cone of Outer Collimator
//
solidOCOuterCone=0;logicOCOuterCone=0;physiOCOuterCone=0;
if(ptAssembly_length > 0.)
  { solidOCOuterCone = new G4Cons("OCOuterCone",DIA3outer/2,CollimatorDia/2,
	DIA2outer/2,CollimatorDia/2,shield_thickness/2,0.*deg,360.*deg);
    logicOCOuterCone = new G4LogicalVolume(solidOCOuterCone,HousingMaterial,
	"OCOuterCone");
    physiOCOuterCone = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
	logicOCOuterCone,"OCOuterCone",logicCollimator,false,0);
  }

//
// Inner Cone of Outer Collimator
//
solidOCInnerCone=0;logicOCInnerCone=0;physiOCInnerCone=0;
if(ptAssembly_length > 0.)
  { solidOCInnerCone = new G4Cons("OCInnerCone",DIA3inner/2,DIA3outer/2,
	DIA2inner/2,DIA2outer/2,shield_thickness/2,0.*deg,360.*deg);
    logicOCInnerCone = new G4LogicalVolume(solidOCInnerCone,ShieldingMaterial,
	"OCInnerCone");
    physiOCInnerCone = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
	logicOCInnerCone,"OCInnerCone",logicCollimator,false,0);
  }

//
// Magnetic Field Chamber
//
solidBChamber=0;logicBChamber=0;physiBChamber=0;
if(ptAssembly_length > 0.)
  { G4double zdepth = (CollimatorHeight/2 - shield_thickness - FeReturnYoke_length/2);
    solidBChamber = new G4Tubs("BChamber",0.*mm,CollimatorDia/2,
	(FeReturnYoke_length)/2, 0.*deg,360.*deg);
    logicBChamber = new G4LogicalVolume(solidBChamber,defaultMaterial,"BChamber");
    physiBChamber = new G4PVPlacement(0,G4ThreeVector(0.,0.,zdepth),
	logicBChamber,"Magnetic Field Chamber",logicCollimatorChamber,false,0);
  }

// Sm-Co Magnets (revised per manufacturing specifications)
solidSmCoMagnet1=0;	logicSmCoMagnet1=0;physiSmCoMagnet1=0;
solidSmCoMagnet2=0;	logicSmCoMagnet2=0;physiSmCoMagnet2=0;
if (ptAssembly_length > 0.)
  { solidSmCoMagnet1 = new G4Tubs("Sm-Co Magnet", 0.*mm, 
	SmCoMagnet_height/2, SmCoMagnet_thickness/2, 0.*deg, 360.*deg);
    logicSmCoMagnet1 = new G4LogicalVolume(solidSmCoMagnet1,
	MagnetMaterial,	"Sm-Co Magnet 1");
    physiSmCoMagnet1 = new G4PVPlacement(MagTrapRotation,
	G4ThreeVector( (0.450)*25.4*mm - SmCoMagnet_thickness/2, 0., 0.),
	logicSmCoMagnet1,"Sm-Co Magnet 1",logicBChamber,false,0);


    solidSmCoMagnet2 = new G4Tubs("Sm-Co Magnet", 0.*mm,
	SmCoMagnet_height/2, SmCoMagnet_thickness/2, 0.*deg,360.*deg);
    logicSmCoMagnet2 = new G4LogicalVolume(solidSmCoMagnet2,
	MagnetMaterial,"Sm-Co Magnet 2");
    physiSmCoMagnet2 = new G4PVPlacement(MagTrapRotation,
	G4ThreeVector( -(0.450)*25.4*mm + SmCoMagnet_thickness/2, 0., 0.),
	logicSmCoMagnet2,"Sm-Co Magnet 2",logicBChamber,false,0);
  }

// Fe Return Yoke
solidReturnYokeA=0;logicReturnYokeA=0;physiReturnYokeA=0;
solidReturnYokeB=0;logicReturnYokeB=0;physiReturnYokeB=0;
solidReturnYokeCD=0;
	logicReturnYokeC=0;physiReturnYokeC=0;
	logicReturnYokeD=0;physiReturnYokeD=0;
if (ptAssembly_length > 0.)
  { solidReturnYokeA = new G4Tubs("Fe Return Yoke (A)",FeReturnYoke_idia/2,
	FeReturnYoke_odia/2,FeReturnYoke_length/2,
	30.*deg,120.*deg);
    logicReturnYokeA = new G4LogicalVolume(solidReturnYokeA,
	ReturnYokeMaterial,"Fe Return Yoke (A)");
    physiReturnYokeA = new G4PVPlacement(0,G4ThreeVector(),
	logicReturnYokeA,"Fe Return Yoke (A)",logicBChamber,false,0);

    solidReturnYokeB = new G4Tubs("Fe Return Yoke (B)",FeReturnYoke_idia/2,
	FeReturnYoke_odia/2,FeReturnYoke_length/2,
	210.*deg,120.*deg);
    logicReturnYokeB = new G4LogicalVolume(solidReturnYokeB,
	ReturnYokeMaterial,"Fe Return Yoke (B)");
    physiReturnYokeB = new G4PVPlacement(0,G4ThreeVector(),
	logicReturnYokeB,"Fe Return Yoke (B)",logicBChamber,false,0);

    solidReturnYokeCD = new G4Box("Fe Return Yoke (CD)",
	(0.100/2)*25.4*mm, (0.260)*25.4*mm,
	FeReturnYoke_length/2 );

    logicReturnYokeC = new G4LogicalVolume(solidReturnYokeCD,
	ReturnYokeMaterial,"Fe Return Yoke (C)");
    physiReturnYokeC = new G4PVPlacement(0,G4ThreeVector(0.5*2.54*cm,0.,0.),
	logicReturnYokeC,"Fe Return Yoke C",logicBChamber,false,0);

    logicReturnYokeD = new G4LogicalVolume(solidReturnYokeCD,
	ReturnYokeMaterial,"Fe Return Yoke (D)");
    physiReturnYokeD = new G4PVPlacement(0,G4ThreeVector(-0.5*2.54*cm,0.,0.),
	logicReturnYokeD,"Fe Return Yoke D",logicBChamber,false,0);
  }

// Aluminium Baffles
solidAlBaffle1=0;logicAlBaffle1=0;physiAlBaffle1=0;
solidAlBaffle2=0;logicAlBaffle2=0;physiAlBaffle2=0;
solidAlBaffle3=0;logicAlBaffle3=0;physiAlBaffle3=0;
if (ptAssembly_length > 0.)
  { solidAlBaffle1 = new G4Tubs("Baffle 1",AlColl1_dia/2,
	AlColl_odia/2,
	AlColl_thickness/2, 0.*deg, 360.*deg);
    logicAlBaffle1 = new G4LogicalVolume(solidAlBaffle1,HousingMaterial,
	"Baffle 1");
    physiAlBaffle1 = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(AlColl_thickness/2 + AlColl_spacing + AlColl_thickness/2)),
	logicAlBaffle1,"Baffle 1",logicBChamber,false,0);

     solidAlBaffle2 = new G4Tubs("Baffle 2",AlColl2_dia/2,
	AlColl_odia/2,
	AlColl_thickness/2, 0.*deg, 360.*deg);
    logicAlBaffle2 = new G4LogicalVolume(solidAlBaffle2,HousingMaterial,
	"Baffle 2");
    physiAlBaffle2 = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
	logicAlBaffle2,"Baffle 2",logicBChamber,false,0);

    solidAlBaffle3 = new G4Tubs("Baffle 3",AlColl3_dia/2,
	AlColl_odia/2,
	AlColl_thickness/2, 0.*deg, 360.*deg);
    logicAlBaffle3 = new G4LogicalVolume(solidAlBaffle3,HousingMaterial,
	"Baffle 3");
    physiAlBaffle3 = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	-(AlColl_thickness/2 + AlColl_spacing + AlColl_thickness/2)),
	logicAlBaffle3,"Baffle 3",logicBChamber,false,0);
  }

// Misc. Shielding
//  Front Cover
solidFrontCover=0;logicFrontCover=0;physiFrontCover=0;
if (ptAssembly_length > 0.)
  { G4double zdepth = (ptAssembly_length/2 - frontCover_thickness/2);
    solidFrontCover = new G4Tubs("Front Cover",DIA1/2,ptAssembly_rad,
	frontCover_thickness/2, 0.*deg,360.*deg);
    logicFrontCover = new G4LogicalVolume(solidFrontCover,
	HousingMaterial,"Front Cover");
    physiFrontCover = new G4PVPlacement(0,G4ThreeVector(0.,0.,zdepth),
	logicFrontCover,"Front Cover",logic_ptAssembly,false,0);
  }

//  Partition
solidPartition=0;logicPartition=0;physiPartition=0;
if (ptAssembly_length > 0.)
  { G4double zdepth = (ptAssembly_length/2 - frontCover_thickness
	- CollimatorSep + partition_thickness/2);
    solidPartition = new G4Cons("Partition",DIA5/2,DeteHousing_dia/2,
	DIA4/2,DeteHousing_dia/2,
	partition_thickness/2,0.*deg,360.*deg);
    logicPartition = new G4LogicalVolume(solidPartition,HousingMaterial,
	"Partition");
    physiPartition = new G4PVPlacement(0,G4ThreeVector(0.,0.,zdepth),
	logicPartition,"Partition",logic_ptAssembly,false,0);
  }

// Rear Cover
solidRearCover=0;logicRearCover=0;physiRearCover=0;
if (ptAssembly_length > 0.)
  { G4double zdepth = (-ptAssembly_length/2 + rearCover_thickness/2);
    solidRearCover = new G4Tubs("Rear Cover",0.*mm,DeteHousing_dia/2,
	rearCover_thickness/2, 0.*deg, 360.*deg);
    logicRearCover = new G4LogicalVolume(solidRearCover,
	HousingMaterial,"Rear Cover");
    physiRearCover = new G4PVPlacement(0,G4ThreeVector(0.,0.,zdepth),
	logicRearCover,"Rear Cover",logic_ptAssembly,false,0);
  }

//
// Detector Housing
//
 solidDeteHousing=0;logicDeteHousing=0;physiDeteHousing=0;

if (ptAssembly_length > 0.)
  { G4double zdepth = (-ptAssembly_length/2 + rearCover_thickness
	+ DeteHousing_length/2 ); 
    solidDeteHousing = new G4Tubs("Detector Housing",0.*mm, DeteHousing_dia/2,
	DeteHousing_length/2,0.*deg,360.*deg);
    logicDeteHousing = new G4LogicalVolume(solidDeteHousing,
	defaultMaterial,"Detector Housing");
    physiDeteHousing = new G4PVPlacement(0,G4ThreeVector(0.,0.,zdepth),
	logicDeteHousing,"Detector Housing",logic_ptAssembly,false,0);
  }

//
// Front Shield / Innermost Collimator
//
 solidInnerCollimator=0;logicInnerCollimator=0;physiInnerCollimator=0;

if (ptAssembly_length > 0.)
  { solidInnerCollimator = new G4Cons("Inner Collimator",DIA6/2,
	DeteHousing_dia/2,DIA5/2,DeteHousing_dia/2,
	shield_thickness/2,0.*deg,360.*deg);
    logicInnerCollimator = new G4LogicalVolume(solidInnerCollimator,
	ShieldingMaterial,"Inner Collimator");
    physiInnerCollimator = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(DeteHousing_length/2 - shield_thickness/2)), // placement
	logicInnerCollimator,"Innermost Collimator",logicDeteHousing,false,0);
  }

//
// Inner Housing
//
solidInnerHousing=0;logicInnerHousing=0;physiInnerHousing=0;

if (ptAssembly_length > 0.)
  { solidInnerHousing = new G4Tubs("InnerHousing",0.*mm,DeteHousing_dia/2,
	(0.477/2)*2.54*cm,0.*deg,360.*deg);
    logicInnerHousing = new G4LogicalVolume(solidInnerHousing,
	defaultMaterial,"Inner Housing");
    physiInnerHousing = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
	logicInnerHousing,"Inner Housing",logicDeteHousing,false,0);
  }

//
// Detector Stack 
//
solidDStack=0;logicDStack=0;physiDStack=0;

if (ptAssembly_length > 0.)
  { solidDStack = new G4Tubs("Detector Stack",0.*mm,DetectorRadius,
	(0.476 - 2*0.050)/2*2.54*cm,0.*deg,360.*deg);
    logicDStack = new G4LogicalVolume(solidDStack,defaultMaterial,"Detector Stack");
    physiDStack = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
	logicDStack,"Detector Stack",logicInnerHousing,false,0);
  }

//
// Belville Washer 1 (really flat cone, so ~ a tube!)
//
solidBelville1=0;logicBelville1=0;physiBelville1=0;

if (ptAssembly_length > 0.)
  { solidBelville1 = new G4Tubs("Belville 1",BelvilleInnerRadius,BelvilleOuterRadius,
	BelvilleThickness/2,0.*deg,360.*deg);
    logicBelville1 = new G4LogicalVolume(solidBelville1,BelvilleMaterial,"Belville 1");
    physiBelville1 = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.182*2.54*cm),
	logicBelville1,"Belville 1",logicDStack,false,0);
  }

//
// Contact 1
//
solidContact1=0;logicContact1=0;physiContact1=0;

if (ptAssembly_length > 0.)
  { solidContact1 = new G4Tubs("Contact 1",ContactInnerRadius,ContactOuterRadius,
	ContactThickness/2,0.*deg,360.*deg);
    logicContact1 = new G4LogicalVolume(solidContact1,ContactMaterial,"Contact 1");
    physiContact1 = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(ContactThickness+DetectorThickness+(ContactThickness)/2)), // rev3
	logicContact1,"Contact 1",logicDStack,false,0);
  }

//
// Detector 1 
//
solid_Detector1=0;logic_Detector1=0;physi_Detector1=0;

if (ptAssembly_length  > 0.)
  { solid_Detector1 = new G4Tubs("Detector 1",0.*mm,DetectorRadius,
	DetectorThickness/2, 0.*deg,360.*deg);
    logic_Detector1 = new G4LogicalVolume(solid_Detector1,defaultMaterial,"D1");
    physi_Detector1 = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(ContactThickness + DetectorThickness/2)),
	logic_Detector1,"D1",logicDStack,false,0);
  }

// Ring Mount (front)
solid_D1RingMntA=0;logic_D1RingMntA=0;physi_D1RingMntA=0;
if (ptAssembly_length > 0.)
  { solid_D1RingMntA = new  G4Tubs("D1 Ring Mount A",(0.422/2)*25.4*mm,
	DetectorRadius,((0.059+0.031)/2)*25.4*mm,
	0.*deg,360.*deg);
    logic_D1RingMntA = new G4LogicalVolume(solid_D1RingMntA,EpoxyMaterial,
	"D1 Ring Mount A");
    physi_D1RingMntA = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(DetectorThickness/2 - ( (0.059 + 0.031)/2 )*25.4*mm)),
	logic_D1RingMntA,"D1 Ring Mount A",logic_Detector1,false,0);
  }

// Ring Mount (back)
solid_D1RingMntB=0;logic_D1RingMntB=0;physi_D1RingMntB=0;
if (ptAssembly_length > 0.)
  { solid_D1RingMntB = new  G4Tubs("D1 Ring Mount B",(0.236/2)*25.4*mm,
	DetectorRadius,((0.055)/2)*25.4*mm,
	0.*deg,360.*deg);
    logic_D1RingMntB = new G4LogicalVolume(solid_D1RingMntB,EpoxyMaterial,
	"D1 Ring Mount B");
    physi_D1RingMntB = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(-DetectorThickness/2 + ( (0.055)/2 )*25.4*mm)),
	logic_D1RingMntB,"D1 Ring Mount B",logic_Detector1,false,0);
  }

// D1 Surface Deposit (Al)
solid_D1SurfDeposit=0;logic_D1SurfDeposit=0;physi_D1SurfDeposit=0;
if (ptAssembly_length > 0.)
  { solid_D1SurfDeposit = new G4Tubs("D1 Surface Film",0.*mm,D1Silicon_rad,
	(ContactDeposit/AlDensity)/2,0.*deg,360.*deg);
    logic_D1SurfDeposit = new G4LogicalVolume(solid_D1SurfDeposit,D1CoatingMaterial,
	"D1 Surface Deposit (Al)");
    physi_D1SurfDeposit = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(-DetectorThickness/2 + 0.055*25.4*mm + DActiveLayer
	+ (ContactDeposit/AlDensity/2))),
	logic_D1SurfDeposit,"Al Surface Deposit",logic_Detector1,false,0);
  }

//
// D1 Active Layer
//
solid_D1ActiveLayer=0;logic_D1ActiveLayer=0;physi_D1ActiveLayer=0;

if (ptAssembly_length > 0.)
  { solid_D1ActiveLayer = new G4Tubs("D1 Active Layer",0.*mm,D1Silicon_rad,
	DActiveLayer/2,0.*deg,360.*deg);
    logic_D1ActiveLayer = new G4LogicalVolume(solid_D1ActiveLayer,DetectorMaterial,
	"D1 Active Layer");
    physi_D1ActiveLayer = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(-DetectorThickness/2 + 0.055*25.4*mm + DActiveLayer/2)),logic_D1ActiveLayer,
	"D1 Active Layer",logic_Detector1,false,0);
  }

//
// Contact 2
//
solidContact2=0;logicContact2=0;physiContact2=0;

if (ptAssembly_length > 0.)
  { solidContact2 = new G4Tubs("Contact 2",ContactInnerRadius,ContactOuterRadius,
        ContactThickness/2 - 0.001*mm,0.*deg,360.*deg);
    logicContact2 = new G4LogicalVolume(solidContact2,ContactMaterial,"Contact 2");
    physiContact2 = new G4PVPlacement(0,G4ThreeVector(0.,0.,ContactThickness/2),
        logicContact2,"Contact 2",logicDStack,false,0);
  }

// 
// Contact 3
//
solidContact3=0;logicContact3=0;physiContact3=0;

if (ptAssembly_length > 0.)
  { solidContact3 = new G4Tubs("Contact 3",ContactInnerRadius,ContactOuterRadius,
        ContactThickness/2 - 0.001*mm,0.*deg,360.*deg);
    logicContact3 = new G4LogicalVolume(solidContact3,ContactMaterial,"Contact 3");
    physiContact3 = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	-ContactThickness/2),
        logicContact3,"Contact 3",logicDStack,false,0);
  }

//
// Detector 2 
//
solid_Detector2=0;logic_Detector2=0;physi_Detector2=0;

if (ptAssembly_length  > 0.)
  { solid_Detector2 = new G4Tubs("Detector 2",0.*mm,DetectorRadius,
        DetectorThickness/2, 0.*deg,360.*deg);
    logic_Detector2 = new G4LogicalVolume(solid_Detector2,defaultMaterial,"D2");
    physi_Detector2 = new G4PVPlacement(0,G4ThreeVector(0.,0.,-0.088*2.54*cm),
        logic_Detector2,"D2",logicDStack,false,0);
  }

// Ring Mount (front)
solid_D2RingMntA=0;logic_D2RingMntA=0;physi_D2RingMntA=0;
if (ptAssembly_length > 0.)
  { solid_D2RingMntA = new  G4Tubs("D2 Ring Mount A",(0.362/2)*25.4*mm,
	DetectorRadius,((0.059+0.031)/2)*25.4*mm,0.*deg,360.*deg);
    logic_D2RingMntA = new G4LogicalVolume(solid_D2RingMntA,EpoxyMaterial,
	"D2 Ring Mount A");
    physi_D2RingMntA = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(DetectorThickness/2 - ( (0.059 + 0.031)/2 )*25.4*mm)),
	logic_D2RingMntA,"D2 Ring Mount A",logic_Detector2,false,0);
  }

// Ring Mount (back)
solid_D2RingMntB=0;logic_D2RingMntB=0;physi_D2RingMntB=0;
if (ptAssembly_length > 0.)
  { solid_D2RingMntB = new  G4Tubs("D2 Ring Mount B",(0.236/2)*25.4*mm,
	DetectorRadius,((0.055)/2)*25.4*mm,0.*deg,360.*deg);
    logic_D2RingMntB = new G4LogicalVolume(solid_D2RingMntB,EpoxyMaterial,
	"D2 Ring Mount B");
    physi_D2RingMntB = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(-DetectorThickness/2 + ( (0.055)/2 )*25.4*mm)),
	logic_D2RingMntB,"D2 Ring Mount B",logic_Detector2,false,0);
  }

// D2 Surface Deposit (Au)
solid_D2SurfDeposit=0;logic_D2SurfDeposit=0;physi_D2SurfDeposit=0;
if (ptAssembly_length > 0.)
  { solid_D2SurfDeposit = new G4Tubs("D2 Surface Film",0.*mm,D2Silicon_rad,
	(ContactDeposit/AlDensity)/2,0.*deg,360.*deg);
    logic_D2SurfDeposit = new G4LogicalVolume(solid_D2SurfDeposit,D2CoatingMaterial,
	"D2 Surface Deposit (Al)");
    physi_D2SurfDeposit = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(-DetectorThickness/2 + 0.055*25.4*mm + DActiveLayer
	+ (ContactDeposit/AlDensity/2))),
	logic_D2SurfDeposit,"Al (Au) Surface Deposit",logic_Detector2,false,0);
  }

//
// D2 Active Layer
//
solid_D2ActiveLayer=0;logic_D2ActiveLayer=0;physi_D2ActiveLayer=0;

if (ptAssembly_length > 0.)
  { solid_D2ActiveLayer = new G4Tubs("D2 Active Layer",0.*mm,D2Silicon_rad,
        DActiveLayer/2,0.*deg,360.*deg);
    logic_D2ActiveLayer = new G4LogicalVolume(solid_D2ActiveLayer,DetectorMaterial,
        "D2 Active Layer");
    physi_D2ActiveLayer = new G4PVPlacement(0,G4ThreeVector(0.,0.,
        (-DetectorThickness/2 + 0.055*25.4*mm + DActiveLayer/2)),logic_D2ActiveLayer,
        "D2 Active Layer",logic_Detector2,false,0);
  }

//
// Contact 4
//
solidContact4=0;logicContact4=0;physiContact4=0;

if (ptAssembly_length > 0.)
  { solidContact4 = new G4Tubs("Contact 4",ContactInnerRadius,ContactOuterRadius,
        ContactThickness/2,0.*deg,360.*deg);
    logicContact4 = new G4LogicalVolume(solidContact4,ContactMaterial,"Contact 4");
    physiContact4 = new G4PVPlacement(0,G4ThreeVector(0.,0.,
        -(ContactThickness+DetectorThickness+(ContactThickness)/2)), // rev3
        logicContact4,"Contact 4",logicDStack,false,0);
  }

//
// Belville Washer 2 (very flat cones, so ~tubes!))
//
solidBelville2=0;logicBelville2=0;physiBelville2=0;

if (ptAssembly_length > 0.)
  { solidBelville2 = new G4Tubs("Belville 2",BelvilleInnerRadius,BelvilleOuterRadius,
        BelvilleThickness/2,0.*deg,360.*deg);
    logicBelville2 = new G4LogicalVolume(solidBelville2,BelvilleMaterial,"Belville 2");
    physiBelville2 = new G4PVPlacement(0,G4ThreeVector(0.,0.,-0.182*2.54*cm),
        logicBelville2,"Belville 2",logicDStack,false,0);
  }

// Detector Insulator & Holder
//
// i. Front Insulator
solidFrontInsul=0;logicFrontInsul=0;physiFrontInsul=0;
if (ptAssembly_length > 0.)
  { solidFrontInsul = new G4Tubs("Front Insulator",BelvilleInnerRadius,
	DetectorRadius,InsulatorThickness/2,
	0.*deg,360.*deg);
    logicFrontInsul = new G4LogicalVolume(solidFrontInsul,InsulatorMaterial,
	"i. Front");
    physiFrontInsul = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(2*0.015 + 0.146 + 0.012 + 0.050/2)*2.54*cm),logicFrontInsul,
	"i. Front",logicInnerHousing,false,0);
  }

// Longitudinal Insulator
solidLongiInsul=0;logicLongiInsul=0;physiLongiInsul=0;
if (ptAssembly_length > 0.)
  { solidLongiInsul = new G4Tubs("Longitudinal Insulator",DetectorRadius,
	inner_DH_dia/2,(0.477/2)*2.54*cm,0.*deg,360.*deg);
    logicLongiInsul = new G4LogicalVolume(solidLongiInsul,InsulatorMaterial,
	"ii. Longitudinal");
    physiLongiInsul = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
	logicLongiInsul,"ii. Longitudinal",logicInnerHousing,false,0);
  }

// Back Insulator
solidBackInsul=0;logicLongiInsul=0;physiLongiInsul=0;
if (ptAssembly_length > 0.)
  { solidBackInsul = new G4Tubs("Back Insulator",BelvilleInnerRadius,
	DetectorRadius,InsulatorThickness/2,
	0.*deg,360.*deg);
    logicBackInsul = new G4LogicalVolume(solidBackInsul,InsulatorMaterial,
	"iii. Back");
    physiBackInsul = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	-(2*0.015 + 0.146 + 0.012 + 0.001/2 + 0.050/2)*2.54*cm),logicBackInsul,
	"iii. Back",logicInnerHousing,false,0);
  }

//
// Housing Shield
//
solidHouseShielding=0;logicHouseShielding=0;physiHouseShielding=0;

if (ptAssembly_length > 0.)
  { solidHouseShielding = new G4Tubs("House Shielding",inner_DH_dia/2,DeteHousing_dia/2,	
	(0.477/2)*2.54*cm,0.*deg,360.*deg);
    logicHouseShielding = new G4LogicalVolume(solidHouseShielding,ShieldingMaterial,
	"House Shielding");
    physiHouseShielding = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
	logicHouseShielding,"House Shielding",logicInnerHousing,false,0);
  }

// Back Shield
solidBackShield=0;logicBackShield=0;physiBackShield=0;
if (ptAssembly_length > 0.)
  { solidBackShield = new G4Tubs("Back Shield",0.*mm,DeteHousing_dia/2,
	shield_thickness/2,0.*deg,360.*deg);
    logicBackShield = new G4LogicalVolume(solidBackShield,
	ShieldingMaterial,"Shield, Back");
    physiBackShield = new G4PVPlacement(0,G4ThreeVector(0.,0.,
	(-DeteHousing_length/2 + shield_thickness/2)),
	logicBackShield,"Shield, Back",logicDeteHousing,false,0);
  }




  PrintModuleParameters();     

  //always return the physical World
  //
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MEPED_DetectorConstruction::PrintModuleParameters()
{
  G4cout 
//	<< "\n%-----------------------------------------------------------"
        << "\n%--> MEPED module Proton Telescope: 2 Totally Depleted "
        << DActiveLayer/mm << " mm Silicon Surface Barrier Detectors"  << ". "; 
//      << "\n%-----------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MEPED_DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) DetectorMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetShieldingMaterial(G4String materialChoice)
{
  // search by name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) ShieldingMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetHousingMaterial(G4String materialChoice)
{
  // search by name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) HousingMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetContactMaterial(G4String materialChoice)
{
  // search by name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) ContactMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetD1CoatingMaterial(G4String materialChoice)
{ G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) D1CoatingMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetD2CoatingMaterial(G4String materialChoice)
{ G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) D2CoatingMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetEpoxyMaterial(G4String materialChoice)
{ G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) EpoxyMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetBelvilleMaterial(G4String materialChoice)
{
  // search by name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) BelvilleMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetInsulatorMaterial(G4String materialChoice)
{
  // search by name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) InsulatorMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetMagnetMaterial(G4String materialChoice)
{ G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) MagnetMaterial = pttoMaterial;
}
void MEPED_DetectorConstruction::SetReturnYokeMaterial(G4String materialChoice)
{ G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) ReturnYokeMaterial = pttoMaterial;
}

void MEPED_DetectorConstruction::SetDetectorThickness(G4double val)
{
  // change Detector thickness
  DetectorThickness = val;
}

//  G4RunManager::GetRunManager()->GeometryHasBeenModified();
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void MEPED_DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructModule());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

