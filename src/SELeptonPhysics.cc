//
////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Title:  Lepton physics for Space Electronics Applications             //
//  Date:   22 March 2005                                                 //
//  Author: D.H. Wright (SLAC)                                            //
//  									  //
//  Last Updated: 20 April 2010, K. Yando (Dartmouth)			  //
//	Reason: Migrate from G4MultipleScattering to 			  //
//		G4eMultipleScattering for electrons and to		  //
//		G4MuMultipleScattering for muons / tauons		  //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
//

#include "SELeptonPhysics.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// processes
#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"
#include "G4hIonisation.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"

// models
#include "G4ElectroNuclearReaction.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"


SELeptonPhysics::SELeptonPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{;}


SELeptonPhysics::~SELeptonPhysics()
{;}


void SELeptonPhysics::ConstructParticle()
{ 
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();
}


void SELeptonPhysics::ConstructProcess()
{
  G4ProcessManager* pManager = 0;

  // Model for e+/e- nuclear reactions
   
  G4ElectroNuclearReaction* theElectronReaction =
                                   new G4ElectroNuclearReaction();
  
  // Electron physics

  pManager = G4Electron::Electron()->GetProcessManager();
  pManager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4eIonisation(),        -1, 2, 2);
  pManager->AddProcess(new G4eBremsstrahlung(),    -1,-1, 3);  

  G4ElectronNuclearProcess* theElectronNuclearProcess =
                                   new G4ElectronNuclearProcess();
  theElectronNuclearProcess->RegisterMe(theElectronReaction);
  pManager->AddDiscreteProcess(theElectronNuclearProcess);

  //Positron physics

  pManager = G4Positron::Positron()->GetProcessManager(); 
  pManager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4eIonisation(),        -1, 2, 2);
  pManager->AddProcess(new G4eBremsstrahlung(),    -1,-1, 3);  
  pManager->AddProcess(new G4eplusAnnihilation(),   0,-1, 4);

  G4PositronNuclearProcess* thePositronNuclearProcess =
                                   new G4PositronNuclearProcess();
  thePositronNuclearProcess->RegisterMe(theElectronReaction);
  pManager->AddDiscreteProcess(thePositronNuclearProcess);

  // Muon-

  pManager = G4MuonMinus::MuonMinus()->GetProcessManager(); 
  pManager->AddProcess(new G4MuMultipleScattering(),    -1, 1, 1);
  pManager->AddProcess(new G4MuIonisation(),          -1, 2, 2);
  pManager->AddProcess(new G4MuBremsstrahlung(),      -1,-1, 3);  
  pManager->AddProcess(new G4MuPairProduction(),      -1,-1, 4);
  pManager->AddProcess(new G4MuonMinusCaptureAtRest(), 0,-1,-1);
  
  // Muon+

  pManager = G4MuonPlus::MuonPlus()->GetProcessManager(); 
  pManager->AddProcess(new G4MuMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4MuIonisation(),       -1, 2, 2);
  pManager->AddProcess(new G4MuBremsstrahlung(),   -1,-1, 3);  
  pManager->AddProcess(new G4MuPairProduction(),   -1,-1, 4);

  // Tau-

  pManager = G4TauMinus::TauMinus()->GetProcessManager();
  pManager->AddProcess(new G4MuMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);
 
  // Tau+

  pManager = G4TauPlus::TauPlus()->GetProcessManager();
  pManager->AddProcess(new G4MuMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

}
