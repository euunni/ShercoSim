#include "DRsimCellParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DRsimCellParameterisation::DRsimCellParameterisation(std::vector<G4float>& x, std::vector<G4float>& y, std::vector<G4bool>& which)
: G4VPVParameterisation()
{
  for (G4int copyNo = 0; copyNo < which.size(); copyNo++ ) {
    
    fXCell.push_back( x.at(copyNo) );
    fYCell.push_back( y.at(copyNo) );
  }
  // fNumx = numx;
  // fNumy = numy;
}

DRsimCellParameterisation::~DRsimCellParameterisation()
{}

void DRsimCellParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const {
  physVol->SetTranslation(G4ThreeVector(fXCell[copyNo],fYCell[copyNo],0.));
}
