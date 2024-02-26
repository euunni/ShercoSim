#ifndef DRsimMirrorParameterisation_h
#define DRsimMirrorParameterisation_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VisAttributes.hh"
#include <vector>

class G4VPhysicalVolume;

class DRsimMirrorParameterisation : public G4VPVParameterisation {
public:
  DRsimMirrorParameterisation(std::vector<G4float>& x, std::vector<G4float>& y, std::vector<G4bool>& which);
  virtual ~DRsimMirrorParameterisation();

  virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;

private:
  std::vector<G4double> fXMirror;
  std::vector<G4double> fYMirror;
  G4int fNumx;
  G4int fNumy;
};

#endif
