#include "DRsimMaterials.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>


DRsimMaterials* DRsimMaterials::fInstance = 0;

DRsimMaterials::DRsimMaterials() {
  fNistMan = G4NistManager::Instance();
  CreateMaterials();
}

DRsimMaterials::~DRsimMaterials() {}

DRsimMaterials* DRsimMaterials::GetInstance() {
  if (fInstance==0) fInstance = new DRsimMaterials();

  return fInstance;
}

G4Material* DRsimMaterials::GetMaterial(const G4String matName) {
  G4Material* mat = fNistMan->FindOrBuildMaterial(matName);

  if (!mat) mat = G4Material::GetMaterial(matName);
  if (!mat) {
    std::ostringstream o;
    o << "Material " << matName << " not found!";
    G4Exception("DRsimMaterials::GetMaterial","",FatalException,o.str().c_str());
  }

  return mat;
}

G4OpticalSurface* DRsimMaterials::GetOpticalSurface(const G4String surfName) {
  if (surfName=="SiPMSurf") return fSiPMSurf;
  else if (surfName=="FilterSurf") return fFilterSurf;
  else if (surfName=="MirrorSurf") return fMirrorSurf;
  else {
    std::ostringstream o;
    o << "OpticalSurface " << surfName << " not found!";
    G4Exception("DRsimMaterials::GetOpticalSurface","",FatalException,o.str().c_str());
  }

  return nullptr;
}

void DRsimMaterials::CreateMaterials() {
  fNistMan->FindOrBuildMaterial("G4_Galactic");
  fNistMan->FindOrBuildMaterial("G4_AIR");

  G4String symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z=1., a=1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z=6., a=12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z=7., a=14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z=8., a=16.00*g/mole);
  G4Element* F  = new G4Element("Fluorine",symbol="F" , z=9., a=18.9984*g/mole);

  // Lab
  G4int fNofC;
  G4int fNofH;
  char fName[20];

  for ( int i = 0; i < 6; i++ ) {
    fNofC = i + 15;
    fNofH = 2 * (i + 9) + 6;
    fLAB[i] = new G4Material(fName, density=0.863*g/cm3, 2);
    fLAB[i]->AddElement(C, fNofC);
    fLAB[i]->AddElement(H, fNofH);
    fLAB[i]->SetChemicalFormula("AROMATIC");
    
    G4double fLABMol = C->GetA()*fNofC + H->GetA()*fNofH;
    G4MaterialPropertiesTable* mpLAB = new G4MaterialPropertiesTable();
    mpLAB->AddConstProperty("LABMol", fLABMol/g);
    fLAB[i]->SetMaterialPropertiesTable(mpLAB);
  }

  // PPO
  fPPO = new G4Material("PPO", density=1.094*g/cm3, 4);
  fPPO->SetChemicalFormula("FLUOR");
  fPPO->AddElement(C, 15);
  fPPO->AddElement(H, 11);
  fPPO->AddElement(N, 1);
  fPPO->AddElement(O, 1);

  G4double fPPOMol = C->GetA()*15 + H->GetA()*11 + N->GetA()*1 + O->GetA()*1;
  G4MaterialPropertiesTable* mpPPO = new G4MaterialPropertiesTable();
  mpPPO->AddConstProperty("PPOMol", fPPOMol/g);
  fPPO->SetMaterialPropertiesTable(mpPPO);

  // Bis-MSB
  fBisMSB = new G4Material("Bis-MSB", density=1.3*g/cm3, 2);
  fBisMSB->SetChemicalFormula("WLS"); // Wavelength Shifter
  fBisMSB->AddElement(C, 24);
  fBisMSB->AddElement(H, 22);

  G4double fBisMol = C->GetA()*24 + H->GetA()*22;
  G4MaterialPropertiesTable* mpBisMSB = new G4MaterialPropertiesTable();
  mpBisMSB->AddConstProperty("BisMol", fBisMol/g);
  fBisMSB->SetMaterialPropertiesTable(mpBisMSB);

  // LS
  G4double fLSdensity =  0.865*g/cm3;
  fLS = new G4Material("LS", fLSdensity, 8);
  G4double fPPOFrac = 3*g / (1e3*cm3 * fLSdensity);
  G4double fBisFrac = 0.03*g / (1e3*cm3 * fLSdensity);
  
  fLS->AddMaterial(fLAB[0], 0.0047 / (1.0 + fPPOFrac + fBisFrac));
  fLS->AddMaterial(fLAB[1], 0.097 / (1.0 + fPPOFrac + fBisFrac));
  fLS->AddMaterial(fLAB[2], 0.3385 / (1.0 + fPPOFrac + fBisFrac));
  fLS->AddMaterial(fLAB[3], 0.3472 / (1.0 + fPPOFrac + fBisFrac));
  fLS->AddMaterial(fLAB[4], 0.2083 / (1.0 + fPPOFrac + fBisFrac));
  fLS->AddMaterial(fLAB[5], 0.0043 / (1.0 + fPPOFrac + fBisFrac));
  fLS->AddMaterial(fBisMSB, fBisFrac / (1.0 + fPPOFrac + fBisFrac));

  // Water 
  fWater = new G4Material("Water", density=1*g/cm3, 2, kStateLiquid);
  fWater->AddElement(H, 2);
  fWater->AddElement(O, 1);  

  fCu = new G4Material("Copper"  , z = 29., a = 63.546 * g/mole, density = 8.96  * g/cm3);
  fW  = new G4Material("Tungsten", z = 74., a = 183.84 * g/mole, density = 19.30 * g/cm3);
  fFe = new G4Material("Iron"    , z = 26., a = 55.845 * g/mole, density = 7.874 * g/cm3);
  fPb = new G4Material("Lead"    , z = 82., a = 207.2  * g/mole, density = 11.35 * g/cm3);

  fSi = new G4Material("Silicon", z=14., a=28.09*g/mole, density=2.33*g/cm3);
  fAl = new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.699*g/cm3);

  fVacuum = G4Material::GetMaterial("G4_Galactic");
  fAir = G4Material::GetMaterial("G4_AIR");

  fFluoPoly = new G4Material("FluorinatedPolymer", density=1.43*g/cm3, ncomponents=2);
  fFluoPoly->AddElement(C, 2);
  fFluoPoly->AddElement(F, 2);

  fGlass = new G4Material("Glass", density=1.032*g/cm3, 2);
  fGlass->AddElement(C, 91.533*perCent);
  fGlass->AddElement(H, 8.467*perCent);

  fPS = new G4Material("Polystyrene", density=1.05*g/cm3, ncomponents=2);
  fPS->AddElement(C, natoms=8);
  fPS->AddElement(H, natoms=8);

  fPMMA = new G4Material("PMMA", density= 1.19*g/cm3, ncomponents=3);
  fPMMA->AddElement(C, natoms=5);
  fPMMA->AddElement(H, natoms=8);
  fPMMA->AddElement(O, natoms=2);

  fGelatin = new G4Material("Gelatin", density=1.27*g/cm3, ncomponents=4);
  fGelatin->AddElement(C, natoms=102);
  fGelatin->AddElement(H, natoms=151);
  fGelatin->AddElement(N, natoms=31);
  fGelatin->AddElement(O, natoms=39);

  G4MaterialPropertiesTable* mpAir;
  G4MaterialPropertiesTable* mpPS;
  G4MaterialPropertiesTable* mpPMMA;
  G4MaterialPropertiesTable* mpFluoPoly;
  G4MaterialPropertiesTable* mpGlass;
  G4MaterialPropertiesTable* mpSiPM;
  G4MaterialPropertiesTable* mpFilter;
  G4MaterialPropertiesTable* mpFilterSurf;
  G4MaterialPropertiesTable* mpMirror;
  G4MaterialPropertiesTable* mpMirrorSurf;
  G4MaterialPropertiesTable* mpLS;
  G4MaterialPropertiesTable* mpWater;
  G4MaterialPropertiesTable* mpAl;


  // LS Refractive Index
  G4double opEn_RI_LS[] = { // from 800nm to 200nm with 100nm step
    1.54980*eV, 1.77120*eV, 2.06640*eV, 2.47968*eV, 3.09960*eV, 4.13281*eV, 6.19921*eV
  };

  const G4int RIEnt_LS = sizeof(opEn_RI_LS) / sizeof(G4double);

  G4double RI_LS[RIEnt_LS] = { 
    1.47571, 1.47871, 1.48334, 1.49107, 1.50541, 1.53694, 1.63112
  };

  // LS Absorption Length
  G4double waveLen_LS, AbsLen_LS_tmp;
  std::vector<G4double> opEn_Abs_LS;
  std::vector<G4double> AbsLen_LS;

  std::ifstream in;
  in.open("AbsLength_LS.txt", std::ios::in);
  
  while (true) { // wavelength[nm] * opEn[eV] = 1,239.84
    in >> waveLen_LS >> AbsLen_LS_tmp;
    
    if ( !in.good() )
      break;

    opEn_Abs_LS.push_back((1239.84/waveLen_LS)*eV); // Unit : eV
    AbsLen_LS.push_back((AbsLen_LS_tmp/1000.)*m); // Unit : m 
  }
  in.close();

  std::reverse(opEn_Abs_LS.begin(), opEn_Abs_LS.end());
  std::reverse(AbsLen_LS.begin(), AbsLen_LS.end());

  const G4int AbsEnt_LS = AbsLen_LS.size();

  mpLS = new G4MaterialPropertiesTable();
  mpLS->AddProperty("RINDEX",opEn_RI_LS,RI_LS,RIEnt_LS);
  mpLS->AddProperty("ABSLENGTH",&(opEn_Abs_LS[0]),&(AbsLen_LS[0]),AbsEnt_LS);
  mpLS->AddConstProperty("SCINTILLATIONYIELD",9.656/keV);
  mpLS->AddConstProperty("RESOLUTIONSCALE",1.0);
  fLS->SetMaterialPropertiesTable(mpLS);
  fLS->GetIonisation()->SetBirksConstant(0.117*mm/MeV);

  // LS : dy_dwavelength, Otical scattering fraction, REEMISSION_PROB

  // Water Refractive Index
  G4double opEn_RI_Water[] = { // from 800nm to 200nm with 50nm step
    1.54980*eV, 1.65312*eV, 1.77120*eV, 1.90745*eV, 2.06640*eV, 2.25426*eV,
    2.47968*eV, 2.75520*eV, 3.09960*eV, 3.54241*eV, 4.13281*eV, 4.95937*eV, 6.19921*eV
  };
  
  const G4int RIEnt_Water = sizeof(opEn_RI_Water) / sizeof(G4double);

  G4double RI_Water[RIEnt_Water] = { 
    1.3292, 1.32986, 1.33065, 1.33165, 1.33293, 1.33458,
    1.33676, 1.3397, 1.34378, 1.34978, 1.35942, 1.37761, 1.42516
  };

  // Water Absorption length
  G4double opEn_Abs_Water[] = {1.54980*eV, 6.19921*eV};
  const G4int AbsEnt_Water = sizeof(opEn_Abs_Water) / sizeof(G4double);
  G4double AbsLen_Water[AbsEnt_Water] = {10.*m, 10.*m};

  mpWater = new G4MaterialPropertiesTable();
  mpWater->AddProperty("RINDEX",opEn_RI_Water,RI_Water,RIEnt_Water);
  mpWater->AddProperty("ABSLENGTH",&(opEn_Abs_Water[0]),&(AbsLen_Water[0]),AbsEnt_Water);
  fWater->SetMaterialPropertiesTable(mpWater);

  // Aluminium Reflection
  G4double AlRef[AbsEnt_LS]; std::fill_n(AlRef, AbsEnt_LS, 0.95);

  mpAl = new G4MaterialPropertiesTable();
  mpAl->AddProperty("REFLECTIVITY",&(opEn_Abs_LS[0]),AlRef,AbsEnt_LS);
  fAl->SetMaterialPropertiesTable(mpAl);




  G4double opEn[] = { // from 900nm to 300nm with 25nm step
    1.37760*eV, 1.41696*eV, 1.45864*eV, 1.50284*eV, 1.54980*eV, 1.59980*eV, 1.65312*eV, 1.71013*eV,
    1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV, 2.15625*eV, 2.25426*eV, 2.36160*eV,
    2.47968*eV, 2.61019*eV, 2.75520*eV, 2.91728*eV, 3.09960*eV, 3.30625*eV, 3.54241*eV, 3.81490*eV, 4.13281*eV
  };

  const G4int nEnt = sizeof(opEn) / sizeof(G4double);

  G4double RI_Air[nEnt]; std::fill_n(RI_Air,nEnt,1.0);
  mpAir = new G4MaterialPropertiesTable();
  mpAir->AddProperty("RINDEX",opEn,RI_Air,nEnt);
  fAir->SetMaterialPropertiesTable(mpAir);
  
  G4double RI_PMMA[nEnt] = {
    1.48329, 1.48355, 1.48392, 1.48434, 1.48467, 1.48515, 1.48569, 1.48628,
    1.48677, 1.48749, 1.48831, 1.48899, 1.49000, 1.49119, 1.49219, 1.49372,
    1.49552, 1.49766, 1.49953, 1.50252, 1.50519, 1.51000, 1.51518, 1.52182, 1.53055
  };
  G4double AbsLen_PMMA[nEnt] = {
    0.414*m, 0.543*m, 0.965*m, 2.171*m, 2.171*m, 3.341*m, 4.343*m, 1.448*m,
    4.343*m, 14.48*m, 21.71*m, 8.686*m, 28.95*m, 54.29*m, 43.43*m, 48.25*m,
    54.29*m, 48.25*m, 43.43*m, 28.95*m, 21.71*m, 4.343*m, 2.171*m, 0.869*m, 0.434*m
  };

  mpPMMA = new G4MaterialPropertiesTable();
  mpPMMA->AddProperty("RINDEX",opEn,RI_PMMA,nEnt);
  mpPMMA->AddProperty("ABSLENGTH",opEn,AbsLen_PMMA,nEnt);
  fPMMA->SetMaterialPropertiesTable(mpPMMA);

  G4double RI_FluoPoly[nEnt]; std::fill_n(RI_FluoPoly, nEnt, 1.42);
  mpFluoPoly = new G4MaterialPropertiesTable();
  mpFluoPoly->AddProperty("RINDEX",opEn,RI_FluoPoly,nEnt);
  fFluoPoly->SetMaterialPropertiesTable(mpFluoPoly);

  G4double RI_PS[nEnt] = {
    1.57483, 1.57568, 1.57644, 1.57726, 1.57817, 1.57916, 1.58026, 1.58148,
    1.58284, 1.58435, 1.58605, 1.58796, 1.59013, 1.59328, 1.59621, 1.59960,
    1.60251, 1.60824, 1.61229, 1.62032, 1.62858, 1.63886, 1.65191, 1.66888, 1.69165
  };
  G4double AbsLen_PS[nEnt] = {
    2.714*m, 3.102*m, 3.619*m, 4.343*m, 5.791*m, 7.896*m, 4.343*m, 7.896*m,
    5.429*m, 36.19*m, 17.37*m, 36.19*m, 5.429*m, 28.95*m, 21.71*m, 14.48*m,
    12.41*m, 8.686*m, 7.238*m, 1.200*m, 0.200*m, 0.500*m, 0.200*m, 0.100*m, 0.100*m
  };
  G4double scintFast_PS[nEnt] = {
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.03, 0.07, 0.13,
    0.33, 0.63, 1.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00
  };
  mpPS = new G4MaterialPropertiesTable();
  mpPS->AddProperty("RINDEX",opEn,RI_PS,nEnt);
  mpPS->AddProperty("ABSLENGTH",opEn,AbsLen_PS,nEnt);
  mpPS->AddProperty("FASTCOMPONENT",opEn,scintFast_PS,nEnt);
  mpPS->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  mpPS->AddConstProperty("RESOLUTIONSCALE",1.0);
  mpPS->AddConstProperty("FASTTIMECONSTANT",2.8*ns);
  fPS->SetMaterialPropertiesTable(mpPS);
  fPS->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4double RI_Glass[nEnt]; std::fill_n(RI_Glass, nEnt, 1.52);
  G4double Abslength_Glass[nEnt]; std::fill_n(Abslength_Glass, nEnt, 420.*cm);
  mpGlass = new G4MaterialPropertiesTable();
  mpGlass->AddProperty("RINDEX",opEn,RI_Glass,nEnt);
  mpGlass->AddProperty("ABSLENGTH",opEn,Abslength_Glass,nEnt);
  fGlass->SetMaterialPropertiesTable(mpGlass);

  G4double refl_SiPM[nEnt]; std::fill_n(refl_SiPM, nEnt, 0.);
  G4double eff_SiPM[nEnt] = {
    0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
    0.11, 0.13, 0.15, 0.17, 0.19, 0.20, 0.22, 0.23,
    0.24, 0.25, 0.24, 0.23, 0.21, 0.20, 0.17, 0.14, 0.10
  };
  mpSiPM = new G4MaterialPropertiesTable();
  mpSiPM->AddProperty("REFLECTIVITY",opEn,refl_SiPM,nEnt);
  mpSiPM->AddProperty("EFFICIENCY",opEn,eff_SiPM,nEnt);
  fSiPMSurf = new G4OpticalSurface("SiPMSurf",glisur,polished,dielectric_metal);
  fSiPMSurf->SetMaterialPropertiesTable(mpSiPM);

  G4double filterEff[nEnt] = {
    0.913, 0.913, 0.913, 0.913, 0.913, 0.913, 0.913, 0.913, 
    0.913, 0.912, 0.910, 0.907, 0.904, 0.899, 0.884, 0.692, 
    0.015, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000
  };
  
  G4double filterRef[nEnt]; std::fill_n(filterRef,nEnt,0.);
  G4double RI_gel[nEnt]; std::fill_n(RI_gel,nEnt,1.52);
  mpFilter = new G4MaterialPropertiesTable();
  mpFilter->AddProperty("RINDEX",opEn,RI_gel,nEnt);
  fGelatin->SetMaterialPropertiesTable(mpFilter);
  mpFilterSurf = new G4MaterialPropertiesTable();
  mpFilterSurf->AddProperty("TRANSMITTANCE",opEn,filterEff,nEnt);
  mpFilterSurf->AddProperty("REFLECTIVITY",opEn,filterRef,nEnt);
  fFilterSurf = new G4OpticalSurface("FilterSurf",glisur,polished,dielectric_dielectric);
  fFilterSurf->SetMaterialPropertiesTable(mpFilterSurf);

  G4double MirrorRef[nEnt]; std::fill_n(MirrorRef, nEnt, 0.9);
  G4double MirrorEff[nEnt]; std::fill_n(MirrorEff, nEnt, 0.);

  mpMirrorSurf = new G4MaterialPropertiesTable();
  mpMirrorSurf->AddProperty("TRANSMITTANCE",opEn,MirrorEff,nEnt);
  mpMirrorSurf->AddProperty("REFLECTIVITY",opEn,MirrorRef,nEnt);
  fMirrorSurf = new G4OpticalSurface("MirrorSurf",glisur,polished,dielectric_metal);
  fMirrorSurf->SetMaterialPropertiesTable(mpMirrorSurf);
}
