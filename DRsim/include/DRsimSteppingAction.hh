#ifndef DRsimSteppingAction_h
#define DRsimSteppingAction_h 1

#include "DRsimInterface.h"
#include "DRsimEventAction.hh"

#include "G4UserSteppingAction.hh"
#include "G4LogicalVolume.hh"
#include "G4Step.hh"

using namespace std;

class DRsimSteppingAction : public G4UserSteppingAction {
public:
  DRsimSteppingAction(DRsimEventAction* eventAction);
  virtual ~DRsimSteppingAction();
  virtual void UserSteppingAction(const G4Step*);

private:
  DRsimEventAction* fEventAction;
  DRsimInterface::DRsimEdepData fEdep;
  DRsimInterface::DRsimLeakageData fLeak;
  
  G4VPhysicalVolume* GetMotherTower(G4TouchableHandle touchable) { return touchable->GetVolume(touchable->GetHistoryDepth()-1); }

  G4int GetModuleNum(G4String towerName) {
    
    if (towerName.find("0") != std::string::npos) return 0;
    if (towerName.find("1") != std::string::npos) return 1;
    if (towerName.find("2") != std::string::npos) return 2;
    if (towerName.find("3") != std::string::npos) return 3;
    if (towerName.find("4") != std::string::npos) return 4;
    if (towerName.find("5") != std::string::npos) return 5;
    if (towerName.find("6") != std::string::npos) return 6;
    if (towerName.find("7") != std::string::npos) return 7;
    if (towerName.find("8") != std::string::npos) return 8;
    if (towerName.find("9") != std::string::npos) return 9;

    return 9999;
  }
};

#endif
