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

int DRsimDetectorConstruction::fNofRow = 5;
int DRsimDetectorConstruction::fNofCol = 8;
int DRsimDetectorConstruction::fNofModules = fNofRow * fNofCol;

DRsimDetectorConstruction::DRsimDetectorConstruction()
: G4VUserDetectorConstruction(), fMessenger(0), fMaterials(NULL) {
  DefineCommands();
  DefineMaterials();

  clad_C_rMin = 0.49*mm;
  clad_C_rMax = 0.50*mm;
  clad_C_Dz   = 2.5*m;
  clad_C_Sphi = 0.;
  clad_C_Dphi = 2.*M_PI;

  core_C_rMin = 0.*mm;
  core_C_rMax = 0.49*mm;
  core_C_Dz   = 2.5*m;
  core_C_Sphi = 0.;
  core_C_Dphi = 2.*M_PI;

  clad_S_rMin = 0.485*mm;
  clad_S_rMax = 0.50*mm;
  clad_S_Dz   = 2.5*m;
  clad_S_Sphi = 0.;
  clad_S_Dphi = 2.*M_PI;

  core_S_rMin = 0.*mm;
  core_S_rMax = 0.485*mm;
  core_S_Dz   = 2.5*m;
  core_S_Sphi = 0.;
  core_S_Dphi = 2.*M_PI;

  PMTT = 0.3*mm;
  filterT = 0.01*mm;
  reflectorT = 0.03*mm;

  fVisAttrOrange = new G4VisAttributes(G4Colour(1.0,0.5,0.,1.0));
  fVisAttrOrange->SetVisibility(true);
  fVisAttrBlue = new G4VisAttributes(G4Colour(0.,0.,1.0,1.0));
  fVisAttrBlue->SetVisibility(true);
  fVisAttrGray = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.3));
  fVisAttrGray->SetVisibility(true);
  fVisAttrGreen = new G4VisAttributes(G4Colour(0.3,0.7,0.3));
  fVisAttrGreen->SetVisibility(true);
  fVisAttrCyan = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  fVisAttrCyan->SetVisibility(true);
  fVisAttrYellow = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
  fVisAttrYellow->SetVisibility(true);
  fVisAttrMagenta = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
  fVisAttrMagenta->SetVisibility(true);
}



DRsimDetectorConstruction::~DRsimDetectorConstruction() {
  delete fMessenger;
  delete fMaterials;

  delete fVisAttrOrange;
  delete fVisAttrBlue;
  delete fVisAttrGray;
  delete fVisAttrGreen;
  delete fVisAttrCyan;
  delete fVisAttrYellow;
  delete fVisAttrMagenta;
}

void DRsimDetectorConstruction::DefineMaterials() {
  fMaterials = DRsimMaterials::GetInstance();
}

G4VPhysicalVolume* DRsimDetectorConstruction::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  checkOverlaps = false;

  G4VSolid* worldSolid             = new G4Box("worldBox",10.*m,10.*m,10.*m);
  worldLogical                     = new G4LogicalVolume(worldSolid,FindMaterial("G4_Galactic"),"worldLogical");
  G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0,checkOverlaps);

  float moduleUnitDimension = 37.5;

  fFrontL     = 0.;
  fTowerDepth = 100.; 
  fModuleH    = moduleUnitDimension;
  fModuleW    = moduleUnitDimension;
  fFiberUnitH = 1.;


  // fRandomSeed = 1;

  doFiber     = true;
  doReflector = false;
  doPMT       = true;

  fiberUnit   = new G4Box("fiber_SQ", (fFiberUnitH/2) *mm, (1./2) *mm, (fTowerDepth/2) *mm);
  fiberClad   = new G4Tubs("fiber",  0, clad_C_rMax, fTowerDepth/2., 0 *deg, 360. *deg);   // S is the same
  fiberCoreC  = new G4Tubs("fiberC", 0, core_C_rMax, fTowerDepth/2., 0 *deg, 360. *deg);
  fiberCoreS  = new G4Tubs("fiberS", 0, core_S_rMax, fTowerDepth/2., 0 *deg, 360. *deg);

  dimCalc = new dimensionCalc();
  dimCalc->SetFrontL(fFrontL);
  dimCalc->SetTower_height(fTowerDepth);
  dimCalc->SetPMTT(PMTT+filterT);
  dimCalc->SetReflectorT(reflectorT);
  dimCalc->SetNofModules(fNofModules);
  dimCalc->SetNofRow(fNofRow);
  dimCalc->SetNofCol(fNofCol);
  dimCalc->SetModuleHeight(fModuleH);
  dimCalc->SetModuleWidth(fModuleW);

  ModuleBuild(ModuleLogical,PMTGLogical,PMTfilterLogical,PMTcellLogical,PMTcathLogical,ReflectorMirrorLogical,fiberUnitIntersection,fiberCladIntersection,fiberCoreIntersection,fModuleProp);

  delete dimCalc;
  return worldPhysical;
}

void DRsimDetectorConstruction::ConstructSDandField() {
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SiPMName = "SiPMSD";

  // ! Not a memory leak - SDs are deleted by G4SDManager. Deleting them manually will cause double delete!
  for (int i = 0; i < fNofModules; i++) {
    DRsimSiPMSD* SiPMSDmodule = new DRsimSiPMSD("Module_"+std::to_string(i), "ModuleHC_"+std::to_string(i), fModuleProp.at(i));
    SDman->AddNewDetector(SiPMSDmodule);
    PMTcathLogical[i]->SetSensitiveDetector(SiPMSDmodule);
  }
}

void DRsimDetectorConstruction::ModuleBuild(G4LogicalVolume* ModuleLogical_[], 
                                            G4LogicalVolume* PMTGLogical_[], 
                                            G4LogicalVolume* PMTfilterLogical_[], 
                                            G4LogicalVolume* PMTcellLogical_[], 
                                            G4LogicalVolume* PMTcathLogical_[], 
                                            G4LogicalVolume* ReflectorMirrorLogical_[],
                                            std::vector<G4LogicalVolume*> fiberUnitIntersection_[], 
                                            std::vector<G4LogicalVolume*> fiberCladIntersection_[], 
                                            std::vector<G4LogicalVolume*> fiberCoreIntersection_[], 
                                            std::vector<DRsimInterface::DRsimModuleProperty>& ModuleProp_) {
  

  std::vector<G4PVPlacement*> tModulePhyVec;

  for ( int nModule = 0; nModule < fNofModules; nModule++ ) {
    bool isCeren = DRsimInterface::IsCerenkov(nModule);
    
    DRsimInterface::DRsimModuleProperty ModulePropSingle;
    ModulePropSingle.towerXY   = std::make_pair(1, 1);
    ModulePropSingle.ModuleNum = nModule;
    ModuleProp_.push_back(ModulePropSingle);

    if ( isCeren == true ) {
      
      auto tModule = new G4Box("Module", (1008. / 2.) *mm, (37.5 / 2.) *mm, (37.5 / 2.) *mm);
      auto tModuleLog = new G4LogicalVolume(tModule, FindMaterial("Iron"), "tLogModule_" + std::to_string(nModule));
      
      auto tAlHousing = new G4Box("tAlHousing", (997. / 2.) *mm, (25.5 / 2.) *mm, (25.5 / 2.) *mm);
      auto tAlHousingLog = new G4LogicalVolume(tAlHousing, FindMaterial("Aluminum"), "tAlHousingLog");
      auto tAlHousingPhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(0.5,0.,0.), tAlHousingLog, "tAlHousingPhy", tModuleLog, false, 0, false);

      auto tActiveMat = new G4Box("tActiveMat", (996.75 / 2.) *mm, (25. / 2.) *mm, (25. / 2.) *mm);
      auto tActiveMatLog_C = new G4LogicalVolume(tActiveMat, FindMaterial("Water"), "tActiveMatLog_C");
      auto tActiveMatPhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(0.125,0.,0.), tActiveMatLog_C, "tActiveMatPhy", tAlHousingLog, false, 0, false);
      
      new G4LogicalSkinSurface("AlSurf", tAlHousingLog, FindSurface("AlSurf"));
      
      auto tPMTGlass = new G4Box("tPMTGlass", (1. / 2.) *mm, (25.5 / 2.) *mm, (25.5 / 2.) *mm);
      auto tPMTGlassLog = new G4LogicalVolume(tPMTGlass, FindMaterial("Glass"), "tPMTGlassLog");
      auto tPMTGlassPhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(499.5,0.,0.), tPMTGlassLog, "tPMTGlassPhy", tModuleLog, false, 0, false);
    
      auto tOpCookie = new G4Box("tOpCookie", (3. / 2.) *mm, (25.5 / 2.) *mm, (25.5 / 2.) *mm);
      auto tOpCookieLog = new G4LogicalVolume(tOpCookie, FindMaterial("Gelatin"), "tOpCookieLog");
      auto tOpCookiePhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(501.5,0.,0.), tOpCookieLog, "tOpCookiePhy", tModuleLog, false, 0, false);
    
      auto tPMT = new G4Box("tPMT", (1. / 2.) *mm, (25.5 / 2.) *mm, (25.5 / 2.) *mm);
      PMTcathLogical_[nModule] = new G4LogicalVolume(tPMT, FindMaterial("Silicon"), "tPMTLog");
      auto tPMTPhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(503.5,0.,0.), PMTcathLogical_[nModule], "tPMTPhy", tModuleLog, false, 0, false);
      
      new G4LogicalSkinSurface("PMTsurf", PMTcathLogical_[nModule], FindSurface("PMTsurf"));
      
      tModulePhyVec.push_back(new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.,1.,0.), 90. *deg), dimCalc->GetOrigin(nModule), tModuleLog, "Module_" + std::to_string(nModule), worldLogical, true, nModule, false));

      tActiveMatLog_C->SetVisAttributes(fVisAttrBlue);
      tPMTGlassLog->SetVisAttributes(fVisAttrCyan);
      tAlHousingLog->SetVisAttributes(fVisAttrGray);
      tOpCookieLog->SetVisAttributes(fVisAttrYellow);
      PMTcathLogical_[nModule]->SetVisAttributes(fVisAttrGreen); 
    } else {
      
      auto tModule = new G4Box("Module", (1008. / 2.) *mm, (37.5 / 2.) *mm, (37.5 / 2.) *mm);
      auto tModuleLog = new G4LogicalVolume(tModule, FindMaterial("Iron"), "tLogModule_" + std::to_string(nModule));
      
      auto tAlHousing = new G4Box("tAlHousing", (997. / 2.) *mm, (25.5 / 2.) *mm, (25.5 / 2.) *mm);
      auto tAlHousingLog = new G4LogicalVolume(tAlHousing, FindMaterial("Aluminum"), "tAlHousingLog");
      auto tAlHousingPhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(0.5,0.,0.), tAlHousingLog, "tAlHousingPhy", tModuleLog, false, 0, false);

      auto tActiveMat = new G4Box("tActiveMat", (996.75 / 2.) *mm, (25. / 2.) *mm, (25. / 2.) *mm);
      auto tActiveMatLog_S = new G4LogicalVolume(tActiveMat, FindMaterial("LS"), "tActiveMatLog_S");
      auto tActiveMatPhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(0.125,0.,0.), tActiveMatLog_S, "tActiveMatPhy", tAlHousingLog, false, 0, false);
      
      new G4LogicalSkinSurface("AlSurf", tAlHousingLog, FindSurface("AlSurf"));
      
      auto tPMTGlass = new G4Box("tPMTGlass", (1. / 2.) *mm, (25.5 / 2.) *mm, (25.5 / 2.) *mm);
      auto tPMTGlassLog = new G4LogicalVolume(tPMTGlass, FindMaterial("Glass"), "tPMTGlassLog");
      auto tPMTGlassPhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(499.5,0.,0.), tPMTGlassLog, "tPMTGlassPhy", tModuleLog, false, 0, false);
    
      auto tOpCookie = new G4Box("tOpCookie", (3. / 2.) *mm, (25.5 / 2.) *mm, (25.5 / 2.) *mm);
      auto tOpCookieLog = new G4LogicalVolume(tOpCookie, FindMaterial("Gelatin"), "tOpCookieLog");
      auto tOpCookiePhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(501.5,0.,0.), tOpCookieLog, "tOpCookiePhy", tModuleLog, false, 0, false);
    
      auto tPMT = new G4Box("tPMT", (1. / 2.) *mm, (25.5 / 2.) *mm, (25.5 / 2.) *mm);
      PMTcathLogical_[nModule] = new G4LogicalVolume(tPMT, FindMaterial("Silicon"), "tPMTLog");
      auto tPMTPhy = new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(503.5,0.,0.), PMTcathLogical_[nModule], "tPMTPhy", tModuleLog, false, 0, false);

      new G4LogicalSkinSurface("PMTsurf", PMTcathLogical_[nModule], FindSurface("PMTsurf"));
      
      tModulePhyVec.push_back(new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.,1.,0.), 90. *deg), dimCalc->GetOrigin(nModule), tModuleLog, "Module_" + std::to_string(nModule), worldLogical, true, nModule, false));    
      
      tActiveMatLog_S->SetVisAttributes(fVisAttrOrange);
      tPMTGlassLog->SetVisAttributes(fVisAttrCyan);
      tAlHousingLog->SetVisAttributes(fVisAttrGray);
      tOpCookieLog->SetVisAttributes(fVisAttrYellow);
      PMTcathLogical_[nModule]->SetVisAttributes(fVisAttrGreen); 
    } 
  }


  // for (int i = 0; i < fNofModules; i++) {    
  //   moduleName = setModuleName(i);
    
  //   dimCalc->SetisModule(true);
  //   module = new G4Box("Mudule", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (fTowerDepth/2.) *mm );
  //   ModuleLogical_[i] = new G4LogicalVolume(module,FindMaterial("Copper"),moduleName);
  //   // G4VPhysicalVolume* modulePhysical = new G4PVPlacement(0,dimCalc->GetOrigin(i),ModuleLogical_[i],moduleName,worldLogical,false,0,checkOverlaps);
  //   new G4PVPlacement(0,dimCalc->GetOrigin(i),ModuleLogical_[i],moduleName,worldLogical,false,0,checkOverlaps);

  //   if ( doPMT ) {
  //     dimCalc->SetisModule(false);  
  //     pmtg = new G4Box("PMTG", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (PMTT+filterT)/2. *mm );
  //     PMTGLogical_[i]  = new G4LogicalVolume(pmtg,FindMaterial("G4_AIR"),moduleName);
  //     new G4PVPlacement(0,dimCalc->GetOrigin_PMTG(i),PMTGLogical_[i],moduleName,worldLogical,false,0,checkOverlaps);
  //   }

  //   FiberImplement(i,,fiberUnitIntersection_,fiberCladIntersection_,fiberCoreIntersection_);

  //   DRsimInterface::DRsimModuleProperty ModulePropSingle;
  //   ModulePropSingle.towerXY   = fTowerXY;
  //   ModulePropSingle.ModuleNum = i;
  //   ModuleProp_.push_back(ModulePropSingle);ModuleLogical_

  //   if ( doPMT ) {
  //     G4VSolid* SiPMlayerSolid = new G4Box("SiPMlayerSolid", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (PMTT/2.) *mm );
  //     G4LogicalVolume* SiPMlayerLogical = new G4LogicalVolume(SiPMlayerSolid,FindMaterial("G4_AIR"),"SiPMlayerLogical");
  //     new G4PVPlacement(0,G4ThreeVector(0.,0.,filterT/2.),SiPMlayerLogical,"SiPMlayerPhysical",PMTGLogical_[i],false,0,checkOverlaps);

  //     G4VSolid* filterlayerSolid = new G4Box("filterlayerSolid", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (filterT/2.) *mm );
  //     G4LogicalVolume* filterlayerLogical = new G4LogicalVolume(filterlayerSolid,FindMaterial("Glass"),"filterlayerLogical");
  //     new G4PVPlacement(0,G4ThreeVector(0.,0.,-PMTT/2.),filterlayerLogical,"filterlayerPhysical",PMTGLogical_[i],false,0,checkOverlaps);

  //     G4VSolid* PMTcellSolid = new G4Box("PMTcellSolid", 1.2/2. *mm, 1.2/2. *mm, PMTT/2. *mm );
  //     PMTcellLogical_[i] = new G4LogicalVolume(PMTcellSolid,FindMaterial("Glass"),"PMTcellLogical_");

  //     // DRsimCellParameterisation* PMTcellParam = new DRsimCellParameterisation(fTowerXY.first,fTowerXY.second);
  //     DRsimCellParameterisation* PMTcellParam = new DRsimCellParameterisation(fFiberX, fFiberY, fFiberWhich);
  //     G4PVParameterised* PMTcellPhysical = new G4PVParameterised("PMTcellPhysical",PMTcellLogical_[i],SiPMlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second,PMTcellParam);

  //     G4VSolid* PMTcathSolid = new G4Box("PMTcathSolid", 1.2/2. *mm, 1.2/2. *mm, filterT/2. *mm );
  //     PMTcathLogical_[i] = new G4LogicalVolume(PMTcathSolid,FindMaterial("Silicon"),"PMTcathLogical_");
  //     new G4PVPlacement(0,G4ThreeVector(0.,0.,(PMTT-filterT)/2.*mm),PMTcathLogical_[i],"PMTcathPhysical",PMTcellLogical_[i],false,0,checkOverlaps);
  //     new G4LogicalSkinSurface("Photocath_surf",PMTcathLogical_[i],FindSurface("SiPMSurf"));

  //     G4VSolid* filterSolid = new G4Box("filterSolid", 1.2/2. *mm, 1.2/2. *mm, filterT/2. *mm );
  //     PMTfilterLogical_[i] = new G4LogicalVolume(filterSolid,FindMaterial("Gelatin"),"PMTfilterLogical_");

  //     int filterNo = (int)(fTowerXY.first * fTowerXY.second) / 2;
  //     if ( fTowerXY.first % 2 == 1)
  //       filterNo++;

  //     // DRsimFilterParameterisation* filterParam = new DRsimFilterParameterisation(fTowerXY.first,fTowerXY.second);
  //     DRsimFilterParameterisation* filterParam = new DRsimFilterParameterisation(fFiberX, fFiberY, fFiberWhich);
  //     G4PVParameterised* filterPhysical = new G4PVParameterised("filterPhysical",PMTfilterLogical_[i],filterlayerLogical,kXAxis,filterNo,filterParam);
  //     new G4LogicalBorderSurface("filterSurf",filterPhysical,PMTcellPhysical,FindSurface("FilterSurf"));
          
  //     PMTcathLogical_[i]->SetVisAttributes(fVisAttrGreen);
  //     PMTfilterLogical_[i]->SetVisAttributes(fVisAttrOrange);
  //   }

    // if ( doReflector ) {
    //   G4VSolid* ReflectorlayerSolid = new G4Box("ReflectorlayerSolid", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (reflectorT/2.) *mm );
    //   G4LogicalVolume* ReflectorlayerLogical = new G4LogicalVolume(ReflectorlayerSolid,FindMaterial("G4_Galactic"),"ReflectorlayerLogical");
    //   new G4PVPlacement(0,dimCalc->GetOrigin_Reflector(i),ReflectorlayerLogical,"ReflectorlayerPhysical",worldLogical,false,0,checkOverlaps);

    //   G4VSolid* mirrorSolid = new G4Box("mirrorSolid", 1.2/2. *mm, 1.2/2. *mm, reflectorT/2. *mm );
    //   ReflectorMirrorLogical_[i] = new G4LogicalVolume(mirrorSolid,FindMaterial("Aluminum"),"ReflectorMirrorLogical_");

    //   // DRsimMirrorParameterisation* mirrorParam = new DRsimMirrorParameterisation(fTowerXY.first,fTowerXY.second);
    //   DRsimMirrorParameterisation* mirrorParam = new DRsimMirrorParameterisation(fFiberX, fFiberY, fFiberWhich);
    //   G4PVParameterised* mirrorPhysical = new G4PVParameterised("mirrorPhysical",ReflectorMirrorLogical_[i],ReflectorlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second/2,mirrorParam);
    //   // new G4LogicalBorderSurface("MirrorSurf",mirrorPhysical,modulePhysical,FindSurface("MirrorSurf"));
    //   new G4LogicalSkinSurface("MirrorSurf",ReflectorMirrorLogical_[i],FindSurface("MirrorSurf"));

    //   ReflectorMirrorLogical_[i]->SetVisAttributes(fVisAttrGray);
  //   }
  // }
}

void DRsimDetectorConstruction::DefineCommands() {}

void DRsimDetectorConstruction::FiberImplement(G4int i, G4LogicalVolume* ModuleLogical__[], 
                                              std::vector<G4LogicalVolume*> fiberUnitIntersection__[], std::vector<G4LogicalVolume*> fiberCladIntersection__[], 
                                              std::vector<G4LogicalVolume*> fiberCoreIntersection__[]) {

  fFiberX.clear();
  fFiberY.clear();
  fFiberWhich.clear();

  int NofFiber = (int)(fModuleW / 1.5);   
  int NofPlate = (int)(fModuleH / 1.5);
  fBottomEdge = fmod(fModuleW, 1.5) / 2.;
  fLeftEdge = fmod(fModuleH, 1.5) / 2.;

  std::cout << "NofFiber : " << NofFiber << " | NofPlate : " << NofPlate << std::endl;

  double randDeviation = 0.; //  double randDeviation = fFiberUnitH - 1.;
  fTowerXY = std::make_pair(NofPlate,NofFiber);
  
  G4bool fWhich = false;  
  for (int k = 0; k < NofPlate; k++) {
    for (int j = 0; j < NofFiber; j++) { 
      /*
        ? fX : # of plate , fY : # of fiber in the plate
      */
      G4float fX = -fModuleH*mm/2 + k*1.5*mm + 0.75*mm + fBottomEdge*mm;
      G4float fY = -fModuleW*mm/2 + j*1.5*mm + 0.75*mm + fLeftEdge*mm;
      fWhich = !fWhich;
      fFiberX.push_back(fX);
      fFiberY.push_back(fY);
      fFiberWhich.push_back(fWhich);
    }
    if ( NofFiber%2==0 ) { fWhich = !fWhich; }   
  }
  
  if ( doFiber ) {
    for (unsigned int j = 0; j<fFiberX.size(); j++) {

      if ( !fFiberWhich.at(j) ) { //c fibre

        tfiberCladIntersection = new G4IntersectionSolid("fiberClad",fiberClad,module,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
        fiberCladIntersection__[i].push_back(new G4LogicalVolume(tfiberCladIntersection,FindMaterial("FluorinatedPolymer"),name));
        new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),fiberCladIntersection__[i].at(j),name,ModuleLogical__[i],false,j,checkOverlaps);

        tfiberCoreIntersection = new G4IntersectionSolid("fiberCore",fiberCoreC,module,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
        fiberCoreIntersection__[i].push_back(new G4LogicalVolume(tfiberCoreIntersection,FindMaterial("PMMA"),name));
        new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberCoreIntersection__[i].at(j),name,fiberCladIntersection__[i].at(j),false,j,checkOverlaps);

        fiberCladIntersection__[i].at(j)->SetVisAttributes(fVisAttrGray);
        fiberCoreIntersection__[i].at(j)->SetVisAttributes(fVisAttrBlue);
      } else { // s fibre

        tfiberCladIntersection = new G4IntersectionSolid("fiberClad",fiberClad,module,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
        fiberCladIntersection__[i].push_back(new G4LogicalVolume(tfiberCladIntersection,FindMaterial("PMMA"),name));
        new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),fiberCladIntersection__[i].at(j),name,ModuleLogical__[i],false,j,checkOverlaps);

        tfiberCoreIntersection = new G4IntersectionSolid("fiberCore",fiberCoreS,module,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
        fiberCoreIntersection__[i].push_back(new G4LogicalVolume(tfiberCoreIntersection,FindMaterial("Polystyrene"),name));
        new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberCoreIntersection__[i].at(j),name,fiberCladIntersection__[i].at(j),false,j,checkOverlaps);

        fiberCladIntersection__[i].at(j)->SetVisAttributes(fVisAttrGray);
        fiberCoreIntersection__[i].at(j)->SetVisAttributes(fVisAttrOrange);
      }
    }
  }
}


