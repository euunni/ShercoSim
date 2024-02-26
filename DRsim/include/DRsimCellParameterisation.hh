#ifndef DRsimCellParameterisation_h
#define DRsimCellParameterisation_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VisAttributes.hh"
#include <vector>
class G4VPhysicalVolume;

class DRsimCellParameterisation : public G4VPVParameterisation {
public:
  DRsimCellParameterisation(std::vector<G4float>& x, std::vector<G4float>& y, std::vector<G4bool>& which);
  virtual ~DRsimCellParameterisation();

  virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;

private:
  std::vector<G4double> fXCell;
  std::vector<G4double> fYCell;
  G4int fNumx;
  G4int fNumy;
};

#endif
