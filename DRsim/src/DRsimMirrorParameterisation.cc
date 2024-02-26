#include "RecoInterface.h"
#include "DRsimMirrorParameterisation.hh"
#include "DRsimCellParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DRsimMirrorParameterisation::DRsimMirrorParameterisation(std::vector<G4float>& x, std::vector<G4float>& y, std::vector<G4bool>& which)
: G4VPVParameterisation()
{
  for (G4int copyNo = 0; copyNo < which.size(); copyNo++ ) {

    if ( !which.at(copyNo) ) {
      fXMirror.push_back( x.at(copyNo) );
      fYMirror.push_back( y.at(copyNo) );
    }
  }
}

DRsimMirrorParameterisation::~DRsimMirrorParameterisation() {}

void DRsimMirrorParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const {
  physVol->SetTranslation(G4ThreeVector(fXMirror[copyNo],fYMirror[copyNo],0.));
}
