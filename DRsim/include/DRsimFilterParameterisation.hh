#ifndef DRsimFilterParameterisation_h
#define DRsimFilterParameterisation_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VisAttributes.hh"
#include <vector>

class G4VPhysicalVolume;

class DRsimFilterParameterisation : public G4VPVParameterisation {
public:
  DRsimFilterParameterisation(std::vector<G4float>& x, std::vector<G4float>& y, std::vector<G4bool>& which);
  virtual ~DRsimFilterParameterisation();

  virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;

private:
  std::vector<G4double> fXFilter;
  std::vector<G4double> fYFilter;
  G4int fNumx;
  G4int fNumy;
};

#endif
