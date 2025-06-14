#include "PhononDetectorConstruction.hh"
#include "PhononSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhononDetectorConstruction::PhononDetectorConstruction()
    : fGalactic(0), fSi(0), 
    fWorldPhys(0), fSiSlab(0), fabsorbLayer(0), siVacuumSides(0), siVacuumBottom(0),
    fConstructed(false) {
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhononDetectorConstruction::~PhononDetectorConstruction() {
    delete siVacuumSides;
    delete siVacuumBottom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// deals with the physical geometry, clearing it if necessary, then setting up
G4VPhysicalVolume* PhononDetectorConstruction::Construct()
{
    if (fConstructed) {
        if (!G4RunManager::IfGeometryHasBeenDestroyed()) {
            // Run manager hasn't cleaned volume stores. This code shouldn't execute
            G4GeometryManager::GetInstance()->OpenGeometry();
            G4PhysicalVolumeStore::GetInstance()->Clean();
            G4LogicalVolumeStore::GetInstance()->Clean();
            G4SolidStore::GetInstance()->Clean();
        }
        // Have to completely remove all lattices to avoid warning on reconstruction
        G4LatticeManager::GetLatticeManager()->Reset();
        // Clear all LogicalSurfaces
        // NOTE: No need to redefine the G4CMPSurfaceProperties
        G4CMPLogicalBorderSurface::CleanSurfaceTable();
    }

    DefineMaterials();
    SetupGeometry();
    fConstructed = true;

    return fWorldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhononDetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();

    fGalactic = nist->FindOrBuildMaterial("G4_Galactic");
    fSi = nist->FindOrBuildMaterial("G4_Si");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhononDetectorConstruction::SetupGeometry()
{
    auto SetColour = [&](G4LogicalVolume* lv, G4double p_r, G4double p_g, G4double p_b) {
        auto vis = new G4VisAttributes(G4Colour(p_r, p_g, p_b));
        vis->SetVisibility(true);
        lv->SetVisAttributes(vis);
        };

    G4bool checkOverlaps = true;
    const G4double slab_X = 1200 * um;
    const G4double slab_Y = 4800 * um;
    const G4double slab_half_X = slab_X / 2;
    const G4double slab_half_Y = slab_Y / 2;
    const G4double slab_thk_from_hetero = 2*nm + 50 * nm + 3 * nm
        + 225 * nm;
    const G4double slab_thk_from_subs = 50 * um; 
    /* ^^ Change me if simulating different thicknesses of substrate*/
    const G4double slab_half_thk = (slab_thk_from_hetero + slab_thk_from_subs) / 2.0;

    const G4double absorb_thickness = 10 * nm;
    const G4double absorb_half_thk = absorb_thickness / 2.0;

    //const G4double trackHX = 25 * nm, trackHY = 25 * nm, trackHZ = 1.5 * nm;
    const G4double trackHX = slab_half_X, trackHY = slab_half_Y, trackHZ = 1.5 * nm;

    const G4double total_half_thk = slab_half_thk + absorb_half_thk;
    auto solidWorld = new G4Box("World",                           // its name
        1.2 * slab_half_X, 1.2 * slab_half_Y, 1.2 * total_half_thk);  // its size

    auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
        fGalactic,                                       // its material
        "World");                                        // its name

    fWorldPhys = new G4PVPlacement(nullptr,  // no rotation
        G4ThreeVector(),                           // at (0,0,0)
        logicWorld,                                // its logical volume
        "World",                                   // its name
        nullptr,                                   // its mother  volume
        false,                                     // no boolean operation
        0,                                         // copy number
        checkOverlaps);                            // overlaps checking

    G4Box* siSlab = new G4Box("SiSlab", slab_half_X,
        slab_half_Y, slab_half_thk);
    G4LogicalVolume* logicSlab = new G4LogicalVolume(siSlab, fSi, "SiSlab");
    fSiSlab = new G4PVPlacement(nullptr, G4ThreeVector(), logicSlab, "SiSlab", logicWorld, false, 0);
    // blue
    SetColour(logicSlab, 0.2, 0.2, 0.8);

    G4Box* solidTrack = new G4Box("TrackingRegion", trackHX, trackHY, trackHZ);
    G4LogicalVolume* logicTrack = new G4LogicalVolume(solidTrack, fSi, "TrackingRegion");
    fSolidTrack = new G4PVPlacement(nullptr, 
        G4ThreeVector(0, 0, 0), // centered within logicQW
        logicTrack, "TrackingRegion", logicSlab, false, 1, checkOverlaps);
    SetColour(logicTrack, 0.9, 0.1, 0.1);

    G4Box* absorbLayer = new G4Box("bottomBox", slab_half_X,
        slab_half_Y, absorb_half_thk);
    G4LogicalVolume* absorbLogic = new G4LogicalVolume(absorbLayer, fGalactic, "bottomBox");
    fabsorbLayer = new G4PVPlacement(nullptr,
        G4ThreeVector(0, 0, -absorb_half_thk - slab_half_thk),
        absorbLogic, "bottomBox", logicWorld, false, 0);
    // green
    SetColour(absorbLogic, 0.2, 0.8, 0.2);

    ////
    ////Silicon lattice information
    ////

    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* SiLogical = LM->LoadLattice(fSi, "Si");

    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
    // aligns the lattice normal [1,0,0] with the +Z direction
    SiPhysical->SetMillerOrientation(1, 0, 0);
    LM->RegisterLattice(fSolidTrack, SiPhysical);
    LM->RegisterLattice(fSiSlab, SiPhysical);


    //
    // surfaces determine phonon reflection/absorption
    //
    if (!fConstructed) {
        const G4double GHz = 1e9 * hertz;
        
        // from file:///C:/Users/victo/Downloads/BF00693457.pdf, for low temp
        // only for silicon, need to find sige and al separately.
        // (5/22) APPROX'd the anharmonic reflection with the bulk value
        // (5/29) decided to change to no anharmonic scattering at vacuum boundary
        const std::vector<G4double> SiVacAnhCoeffs = { 0, 0, 0, 0, 0, 0/*1.2e-10*/ };
        const std::vector<G4double> SiGeVacAnhCoeffs = { 0, 0, 0, 0, 0, 0/* NO VALUE */ };

        // ASSUME rms roughness is 2 nm and silicon speed of sound 6350 
        // using effective sound velocity https://en.wikipedia.org/wiki/Debye_model 
        // the Ziman becomes e^(-1.56650674*10^(-5) * f^2)

        // Silicon-vacuum boundary
        const std::vector<G4double> SiVacSpecCoeffs =
        { 1, 0.0002302027389270810, -2.1735704309859052e-5, 
            5.2926139131064460e-8, 3.6540831265189605e-11 };
        const std::vector<G4double> SiVacDiffCoeffs =
        { 0, -0.0002302027389270810, 2.1735704309859052e-5,
            -5.2926139131064460e-8, -3.6540831265189605e-11/*, -1.2e-10*/ };
        // the reflCutoff is chosen because then the Ziman fit is e^(-6), vanishing
        // the anhCutoff is the Debye frequency
        const G4double SiVacAnhCutoff = 15000., SiVacReflCutoff = 618.;   // GHz units

        // 
        // Silicon boundaries
        // 
        // no absorption
        siVacuumSides = new G4CMPSurfaceProperty("siVacuumSides", 
            1.0, 0.0, 0.0, 0.0,   // q absorption, q refl, e min k (to absorb), hole min k
            0.0, 1.0, 0.0 /*set to be vanishing because it's ONLY USED
                           at high frequency where we have no specular refl
                           see G4CMPSurfaceProperty.cc documentation*/,
                          /* 6/2: not sure about the above comment 
                          VARY PARAMETER FOR TESTING/SENSITIVITY!
                          */
            0.0);  // phonon abs, phonon refl (implying transmission), ph specular, p min k
        siVacuumSides->AddScatteringProperties(SiVacAnhCutoff, SiVacReflCutoff, SiVacAnhCoeffs,
            SiVacDiffCoeffs, SiVacSpecCoeffs, GHz, GHz, GHz);
        // AttachPhononSensor(siVacuumSides);

        // high absorption
        siVacuumBottom = new G4CMPSurfaceProperty("siVacuumBottom", 
            1.0, 0.0, 0.0, 0.0,
            0.9 /*can vary this parameter to vary the back-reflection*/, 0.0, 0.0, 0.0);
        siVacuumBottom->AddScatteringProperties(SiVacAnhCutoff, SiVacReflCutoff, SiVacAnhCoeffs,
            SiVacDiffCoeffs, SiVacSpecCoeffs, GHz, GHz, GHz);
        // AttachPhononSensor(siVacuumBottom);

        // tracking-slab interface
        trivialProp = new G4CMPSurfaceProperty("trivial",
            0.0, 0.0, 0.0, 0.0,   // qAbs, qRefl, eMinK, hMinK
            0.0, 0.0, 0.0, 0.0); // pAbs, pRefl, pSpec, pMinK
        std::vector<G4double> zeroA(1, 0.0), zeroD(1, 0.0), zeroS(1, 0.0);
        trivialProp->AddScatteringProperties(1e6, 1e6, zeroA, zeroD, zeroS, GHz, GHz, GHz);
    }

    //
    // Surface attachment
    //
    new G4CMPLogicalBorderSurface("siWall", fSiSlab, fWorldPhys,
        siVacuumSides);
    new G4CMPLogicalBorderSurface("siBottom", fSiSlab, fabsorbLayer,
        siVacuumBottom);
    new G4CMPLogicalBorderSurface("tracking", fSiSlab, fSolidTrack,
        trivialProp);
    new G4CMPLogicalBorderSurface("tracking", fSolidTrack, fSiSlab,
        trivialProp);



    //                                        
    // Visualization attributes
    //
    logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    simpleBoxVisAtt->SetVisibility(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Attach material properties and electrode/sensor handler to surface

void PhononDetectorConstruction::
AttachPhononSensor(G4CMPSurfaceProperty* surfProp) {
    if (!surfProp) return;		// No surface, nothing to do

    // Specify properties of aluminum sensor, same on both detector faces
    // See G4CMPPhononElectrode.hh or README.md for property keys

    // Properties must be added to existing surface-property table
    auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
    sensorProp->AddConstProperty("filmAbsorption", 0.20);    // True sensor area
    sensorProp->AddConstProperty("filmThickness", 600. * nm);
    sensorProp->AddConstProperty("gapEnergy", 173.715e-6 * eV);
    sensorProp->AddConstProperty("lowQPLimit", 3.);
    sensorProp->AddConstProperty("phononLifetime", 242. * ps);
    sensorProp->AddConstProperty("phononLifetimeSlope", 0.29);
    sensorProp->AddConstProperty("vSound", 3.26 * km / s);
    sensorProp->AddConstProperty("subgapAbsorption", 0.1);

    // Attach electrode object to handle KaplanQP interface
    surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
}

