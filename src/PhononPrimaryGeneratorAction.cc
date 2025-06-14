/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/PhononPrimaryGeneratorAction.cc
/// \brief Implementation of the PhononPrimaryGeneratorAction class
//
// $Id: e75f788b103aef810361fad30f75077829192c13 $
//
// 20140519  Allow the user to specify phonon type by name in macro; if
//	     "geantino" is set, use random generator to select.

#include "PhononPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononLong.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

PhononPrimaryGeneratorAction::PhononPrimaryGeneratorAction() { 
    // just to initialize the G4Particle Gun
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);    // need the number of particles

  // default particle kinematics ("geantino" triggers random phonon choice)
  fParticleGun->SetParticleDefinition(G4Geantino::Definition()); // particle type
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection()); // unit vector direction
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0)); // pos
  fParticleGun->SetParticleEnergy(0.0075*eV);   // energy
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


PhononPrimaryGeneratorAction::~PhononPrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 
void PhononPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // picks the actual phonon type and then generates the event
  if (fParticleGun->GetParticleDefinition() == G4Geantino::Definition()) {
    G4double selector = G4UniformRand();
    // probability of each branch of phonon
    // need to change these for different initial energies? and materials?

    // now changed to Si values, but this is the generated phonon
    if (selector < 0.531) {
        fParticleGun->SetParticleDefinition(G4PhononTransSlow::Definition());
    }
    else if (selector < 0.907) {
        fParticleGun->SetParticleDefinition(G4PhononTransFast::Definition());
    }
    else {
        fParticleGun->SetParticleDefinition(G4PhononLong::Definition());
    }
  }

  // probably unnecessary? 
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());

  // a G4PrimaryVertex is where and when the primary particles are created
  // (x,y,z,t) and linkedlist of primary particles that will be fired
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


