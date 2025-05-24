#ifndef PHONONSTEPPINGACTION_H
#define PHONONSTEPPINGACTION_H

#include "G4UserSteppingAction.hh"
#include <fstream>

/// SteppingAction to record every phonon step into a CSV file.
class PhononSteppingAction : public G4UserSteppingAction {
public:
    PhononSteppingAction();
    ~PhononSteppingAction() override;

    /// Called by Geant4 after each step.
    /// @param step pointer to the current G4Step (see G4UserSteppingAction docs :contentReference[oaicite:2]{index=2}).
    void UserSteppingAction(const G4Step* step) override;

private:
    std::ofstream fout_;
};

#endif // PHONONSTEPPINGACTION_H
