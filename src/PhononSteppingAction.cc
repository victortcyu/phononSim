#include "PhononSteppingAction.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononLong.hh"

// Constructor: open file and write header line
PhononSteppingAction::PhononSteppingAction() {
    fout_.open("phonon_tracking.csv");
    fout_ << "trackID,x/nm,y/nm,z/nm,time_ns,energy_eV\n";
}

// Destructor: close file
PhononSteppingAction::~PhononSteppingAction() {
    if (fout_.is_open()) fout_.close();
}

void PhononSteppingAction::UserSteppingAction(const G4Step* step) {
    // 1) Only keep phonon tracks
    auto track = step->GetTrack();
    if (track->GetDefinition() != G4PhononLong::Definition() && 
        track->GetDefinition() != G4PhononTransSlow::Definition() &&
        track->GetDefinition() != G4PhononTransFast::Definition())
        return;

    // 2) Get pre-step information
    auto prePoint = step->GetPreStepPoint();
    auto pos = prePoint->GetPosition();       // G4StepPoint::GetPosition() :contentReference[oaicite:3]{index=3}
    auto time = prePoint->GetGlobalTime();     // G4StepPoint::GetGlobalTime() :contentReference[oaicite:4]{index=4}
    auto energy = prePoint->GetKineticEnergy();  // G4StepPoint::GetKineticEnergy() :contentReference[oaicite:5]{index=5}

    // 3) Write out: trackID, x/nm, y/nm, z/nm, time/ns, energy/eV
    fout_ << track->GetTrackID() << ","
        << pos.x() / nm << "," << pos.y() / nm << "," << pos.z() / nm << ","
        << time / ns << "," << energy / eV << "\n";
}
