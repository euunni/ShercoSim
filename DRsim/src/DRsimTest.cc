#include "DRsimDetectorConstruction.hh"
#include "DRsimCellParameterisation.hh"
#include "DRsimFilterParameterisation.hh"
#include "DRsimMirrorParameterisation.hh"
#include "DRsimSiPMSD.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4IntersectionSolid.hh"
#include "G4SDManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"

#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

using namespace std;

G4ThreadLocal DRsimMagneticField* DRsimDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* DRsimDetectorConstruction::fFieldMgr = 0;

int DRsimDetectorConstruction::fNofRow = 1;
int DRsimDetectorConstruction::fNofModules = fNofRow * fNofRow;

