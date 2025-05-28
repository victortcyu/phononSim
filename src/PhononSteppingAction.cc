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
    fout_ << "trackID, stepNumber, x/nm, y/nm, z/nm, time_ns, energy_meV\n";
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
    auto postPoint = step->GetPostStepPoint();
    auto prevPhysVol = prePoint->GetPhysicalVolume();
    if (!prevPhysVol // undefined pointer
        || postPoint->GetPhysicalVolume()->GetName() != "TrackingRegion" // not in the correct region
        || prePoint->GetStepStatus() != fGeomBoundary // prevent double-counting by checking it came from outside 
        )
        return;
    auto pos = postPoint->GetPosition();       // G4StepPoint::GetPosition() :contentReference[oaicite:3]{index=3}
    auto time = postPoint->GetGlobalTime();     // G4StepPoint::GetGlobalTime() :contentReference[oaicite:4]{index=4}
    auto energy = postPoint->GetKineticEnergy();  // G4StepPoint::GetKineticEnergy() :contentReference[oaicite:5]{index=5}

    // 3) Write out: trackID, x/nm, y/nm, z/nm, time/ns, energy/eV
    // note that for such geometric crossings, our post-step point will always be on the boundary
    // thus the z value will not be interesting. However, the step will now be in the new volume
    fout_ << track->GetTrackID() << "," << track->GetCurrentStepNumber() << ","
        << pos.x() / nm << "," << pos.y() / nm << "," << pos.z() / nm << ","
        << time / ns << "," << energy / eV * 1e3 << "\n";
}
