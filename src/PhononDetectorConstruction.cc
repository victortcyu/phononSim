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
    : fGalactic(0), fSi(0), fSiGe(0), fAl(0),
    fWorldPhys(0), fSiTest(0), siVacuumSides(0), siVacuumBottom(0),
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
    G4double densSiGe = (0.7 * 2.33 + 0.3 * 5.32) * g / cm3;
    fSiGe = new G4Material("SiGe_70_30", densSiGe, 2);
    fSiGe->AddElement(nist->FindOrBuildElement("Si"), 0.7);
    fSiGe->AddElement(nist->FindOrBuildElement("Ge"), 0.3);
    fAl = nist->FindOrBuildMaterial("G4_Al");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhononDetectorConstruction::SetupGeometry()
{
    /* color setting code */
    auto SetColour = [&](G4LogicalVolume* lv, G4double p_r, G4double p_g, G4double p_b) {
        auto vis = new G4VisAttributes(G4Colour(p_r, p_g, p_b));
        vis->SetVisibility(true);
        lv->SetVisAttributes(vis);
        };
    //// 


    G4bool checkOverlaps = true;
    const G4double heterostr_X = 1200 * um;
    const G4double heterostr_Y = 4800 * um;
    const G4double heterostr_half_X = heterostr_X / 2;
    const G4double heterostr_half_Y = heterostr_Y / 2;
    const G4double cap_thickness = 2 * nm;
    const G4double spacer_thickness = 50 * nm;
    const G4double qw_thickness = 3 * nm;
    const G4double buffer_thickness = 225 * nm;
    const G4double heterostr_total_thickness = cap_thickness + spacer_thickness + qw_thickness
        + buffer_thickness;
    const G4double heterostr_total_half_thickness = heterostr_total_thickness / 2.0;

    const G4double absorb_layer_thickness = 10 * nm;
    const G4double absorb_layer_half_thickness = absorb_layer_thickness / 2.0;

    const G4double substrate_thickness = 50 * um;
    const G4double total_thickness = heterostr_total_half_thickness + absorb_layer_half_thickness
        + substrate_thickness;
    const G4double total_half_thickness = total_thickness / 2.0;



    auto solidWorld = new G4Box("World",                           // its name
        1.2 * heterostr_half_X, 1.2 * heterostr_half_Y, 1.2 * total_half_thickness);  // its size

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

    // G4double z0 = -total_half_thickness;

    G4Box* siliconTest = new G4Box("SiBox", heterostr_half_X, 
        heterostr_half_Y, heterostr_total_half_thickness);
    G4LogicalVolume* logicTest = new G4LogicalVolume(siliconTest, fSi, "SiBox");
    fSiTest = new G4PVPlacement(nullptr, G4ThreeVector(), logicTest, "SiBox", logicWorld, false, 0);
    SetColour(logicTest, 0.2, 0.2, 0.8);

    G4Box* absorbLayer = new G4Box("bottomBox", heterostr_half_X,
        heterostr_half_Y, absorb_layer_half_thickness);
    G4LogicalVolume* absorbLogic = new G4LogicalVolume(absorbLayer, fGalactic, "bottomBox");
    fabsorbLayer = new G4PVPlacement(nullptr, 
        G4ThreeVector(0,0,-absorb_layer_half_thickness - heterostr_total_half_thickness), 
        absorbLogic, "bottomBox", logicWorld, false, 0);
    SetColour(absorbLogic, 0.2, 0.8, 0.2);


    /*G4Box* solidBuffer = new G4Box("SiGeBuffer", halfX, halfY, tBuffer / 2);
    G4LogicalVolume* logicBuffer = new G4LogicalVolume(solidBuffer, matSiGe, "SiGeBuffer");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z0 + tBuffer / 2), logicBuffer, "SiGeBuffer", logicWorld, false, 0);
    SetColour(logicBuffer, 0.2, 0.2, 0.8);
    z0 += tBuffer;*/

    //G4Box* solidQW = new G4Box("SiQW", halfX, halfY, tQW / 2);
    //G4LogicalVolume* logicQW = new G4LogicalVolume(solidQW, matSi, "SiQW");
    //new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z0 + tQW / 2), logicQW, "SiQW", logicWorld, false, 0);
    //SetColour(logicQW, 0.8, 0.8, 0.2);
    //z0 += tQW;

    //G4Box* solidSpacer = new G4Box("SiGeSpacer", halfX, halfY, tSpacer / 2);
    //G4LogicalVolume* logicSpacer = new G4LogicalVolume(solidSpacer, matSiGe, "SiGeSpacer");
    //new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z0 + tSpacer / 2), logicSpacer, "SiGeSpacer", logicWorld, false, 0);
    //SetColour(logicSpacer, 0.2, 0.8, 0.2);
    //z0 += tSpacer;

    //G4Box* solidCap = new G4Box("SiCap", halfX, halfY, tCap / 2);
    //G4LogicalVolume* logicCap = new G4LogicalVolume(solidCap, matSi, "SiCap");
    //new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z0 + tCap / 2), logicCap, "SiCap", logicWorld, false, 0);
    //SetColour(logicCap, 0.8, 0.2, 0.2);
    
    //
    // World
    //
    //G4VSolid* worldSolid = new G4Box("World",16.*cm,16.*cm,16.*cm);
    //G4LogicalVolume* worldLogical =
    //  new G4LogicalVolume(worldSolid,fLiquidHelium,"World");
    //worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
    //fWorldPhys = new G4PVPlacement(
    //    0, // rotation
    //    G4ThreeVector(),     // translation from mother's origin, but no mother. 
    //                         // center of world = coord system
    //    worldLogical, // logical volume to place
    //    "World", //name 
    //    0, // no mother 
    //    false,0);

    //                               
    //// Germanium cylinder - this is the volume in which we will propagate phonons
    ////  
    //G4VSolid* fGermaniumSolid = new G4Tubs("fGermaniumSolid",0.*cm,3.81*cm,
    //                                       1.27*cm, 0.*deg, 360.*deg);
    //G4LogicalVolume* fGermaniumLogical =
    //  new G4LogicalVolume(fGermaniumSolid,fGermanium,"fGermaniumLogical");
    //G4VPhysicalVolume* GePhys =
    //  new G4PVPlacement(0,G4ThreeVector(),fGermaniumLogical,"fGermaniumPhysical",
    //                    worldLogical,false,0);

    ////
    ////Silicon lattice information
    ////

    G4LatticeManager* LM = G4LatticeManager::GetLatticeManager();
    G4LatticeLogical* SiLogical = LM->LoadLattice(fSi, "Si");

    // G4LatticePhysical assigns G4LatticeLogical a physical orientation
    G4LatticePhysical* SiPhysical = new G4LatticePhysical(SiLogical);
    // aligns the lattice normal [1,0,0] with the +Z direction
    SiPhysical->SetMillerOrientation(1, 0, 0);
    LM->RegisterLattice(fSiTest, SiPhysical);

    ////
    //// Aluminum - crystal end caps. This is where phonon hits are registered
    ////
    //G4VSolid* fAluminumSolid = new G4Tubs("aluminiumSolid",0.*cm,3.81*cm,0.01*cm,
    //                                      0.*deg, 360.*deg);
    //G4LogicalVolume* fAluminumLogical =
    //  new G4LogicalVolume(fAluminumSolid,fAluminum,"fAluminumLogical");
    //G4VPhysicalVolume* aluminumTopPhysical = new G4PVPlacement(0,
    //  G4ThreeVector(0.,0.,1.28*cm), fAluminumLogical, "fAluminumPhysical",
    //  worldLogical,false,0);
    //G4VPhysicalVolume* aluminumBotPhysical = new G4PVPlacement(0,
    //  G4ThreeVector(0.,0.,-1.28*cm), fAluminumLogical, "fAluminumPhysical",
    //  worldLogical,false,1);

    //
    // detector -- Note : "sensitive detector" is attached to Germanium crystal
    //
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if (!electrodeSensitivity)
        electrodeSensitivity = new PhononSensitivity("PhononElectrode");
    SDman->AddNewDetector(electrodeSensitivity);
    logicTest->SetSensitiveDetector(electrodeSensitivity);

    //
    // surface between Al and Ge determines phonon reflection/absorption
    //
    if (!fConstructed) {
        const G4double GHz = 1e9 * hertz;

        //the following coefficients and cutoff values are not well-motivated
        //the code below is used only to demonstrate how to set these values.
        
        // from file:///C:/Users/victo/Downloads/BF00693457.pdf, for low temp
        // only for silicon, need to find sige and al separately.
        // APPROX'd the anharmonic reflection with the bulk value
        const std::vector<G4double> anhCoeffs = { 0, 0, 0, 0, 0, 1.2e-10 };

        // ASSUME rms roughness is 2 nm and silicon speed of sound 6350 
        // using effective sound velocity https://en.wikipedia.org/wiki/Debye_model 
        // the Ziman becomes e^(-1.56650674*10^(-5) * f^2)

        // these are all for a vacuum-silicon boundary
        // SiGe will need to use a different fit because of the different speed of sound
        const std::vector<G4double> specCoeffs =
        { 1, 0.0002302027389270810, -2.1735704309859052e-5, 
            5.2926139131064460e-8, 3.6540831265189605e-11 };
        const std::vector<G4double> diffCoeffs =
        { 0, -0.0002302027389270810, 2.1735704309859052e-5,
            -5.2926139131064460e-8, -3.6540831265189605e-11, -1.2e-10 };
        // the reflCutoff is chosen because then the Ziman fit is e^(-6), vanishing
        const G4double anhCutoff = 15000., reflCutoff = 618.;   // Units external

        // no absorption
        siVacuumSides = new G4CMPSurfaceProperty("siVacuumSides", 
            1.0, 0.0, 0.0, 0.0,   // q absorption, q refl, e min k (to absorb), hole min k
            0.0, 1.0, 0.0 /*vanishing because it's only used at high frequency
                           where we have no specular refl
                           see G4CMPSurfaceProperty.cc documentation*/,
            0.0);  // phonon abs, phonon refl (implying transmission), ph specular, p min k
        siVacuumSides->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
            diffCoeffs, specCoeffs, GHz, GHz, GHz);
        // AttachPhononSensor(siVacuumSides);

        // high absorption
        siVacuumBottom = new G4CMPSurfaceProperty("siVacuumBottom", 
            1.0, 0.0, 0.0, 0.0,
            1.0 /*can vary this parameter to vary the back-reflection*/, 0.0, 0.0, 0.0);
        siVacuumBottom->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
            diffCoeffs, specCoeffs, GHz, GHz, GHz);
        // AttachPhononSensor(siVacuumBottom);

    }
    //
    // Separate surfaces for the bottom vs other sides
    //
    new G4CMPLogicalBorderSurface("siWall", fSiTest, fWorldPhys,
        siVacuumSides);
    new G4CMPLogicalBorderSurface("siBottom", fSiTest, fabsorbLayer,
        siVacuumBottom);


    //                                        
    // Visualization attributes
    //
    logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    simpleBoxVisAtt->SetVisibility(true);
    logicTest->SetVisAttributes(simpleBoxVisAtt);
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

