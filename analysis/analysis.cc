#include "RootInterface.h"
#include "RecoInterface.h"
#include "DRsimInterface.h"
#include "functions.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <tuple>


int main(int argc, char* argv[]) {

  float high = 1.;
  float low = 0.;

  gStyle->SetOptFit(1);

  TH1F* tEdep = new TH1F("Total_Edep","Total Energy deposit;MeV;Evt",100,low,high*3000);
  tEdep->Sumw2(); tEdep->SetLineColor(kBlack); tEdep->SetLineWidth(2);
  TH1F* tCtime = new TH1F("Total_C_Time","Total timing of Cerenkov ch.;ns;Evt",150,0,30);
  tCtime->Sumw2(); tCtime->SetLineColor(kBlue); tCtime->SetLineWidth(2);
  TH1F* tStime = new TH1F("Total_S_Time","Total timing of Scintillation ch.;ns;Evt",150,0,30);
  tStime->Sumw2(); tStime->SetLineColor(kRed); tStime->SetLineWidth(2);
  TH1I* tChit = new TH1I("Total_C_Hit","Total hits of Cerenkov ch",100,0.,1500.);
  tChit->Sumw2(); tChit->SetLineColor(kBlue); tChit->SetLineWidth(2);
  TH1I* tShit = new TH1I("Total_S_Hit","Total hits of Scintillation ch",100,0.,30000.);
  tShit->Sumw2(); tShit->SetLineColor(kRed); tShit->SetLineWidth(2);

  // TH1F* Edep = new TH1F("Edep","Energy deposit;MeV;Evt",100,low*1000.,high*1000.);
  // Edep->Sumw2(); Edep->SetLineColor(kBlack); Edep->SetLineWidth(2);
  // TH1F* Ctime = new TH1F("C_Time","timing of Cerenkov ch.;ns;Evt",600,10,70);
  // Ctime->Sumw2(); Ctime->SetLineColor(kBlue); Ctime->SetLineWidth(2);
  // TH1F* Stime = new TH1F("S_Time","timing of Scintillation ch.;ns;Evt",600,10,70);
  // Stime->Sumw2(); Stime->SetLineColor(kRed); Stime->SetLineWidth(2);
  // TH1I* Chit = new TH1I("C_Hit","hits of Cerenkov ch",100,0.,3000.);
  // Chit->Sumw2(); Chit->SetLineColor(kBlue); Chit->SetLineWidth(2);
  // TH1I* Shit = new TH1I("S_Hit","hits of Scintillation ch",100,0.,40000.);
  // Shit->Sumw2(); Shit->SetLineColor(kRed); Shit->SetLineWidth(2);

  // RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>("/d0/scratch/haeun/Sherco/output/v240501/mu_20/test.root", 1);
  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>("/u/user/haeun/Sherco/v240501/ShercoSim/install/input/240626/pi_4/pi_4.root", 1);
  drInterface->set("DRsim","DRsimEventData");

  unsigned int entries = drInterface->entries();
  while (drInterface->numEvt() < entries) {
    if (drInterface->numEvt() % 100 == 0) printf("Analyzing %dth event ...\n", drInterface->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface->read(drEvt);

    float ftEdep = 0.;

    for (auto edepItr = drEvt.Edeps.begin(); edepItr != drEvt.Edeps.end(); ++edepItr) {
      auto edep = *edepItr;
      ftEdep += edep.Edep;
    }

    int ftC_hits = 0; int ftS_hits = 0;

    for (auto towerItr = drEvt.towers.begin(); towerItr != drEvt.towers.end(); ++towerItr) {
      auto sipmItr = *towerItr;
      std::vector<DRsimInterface::DRsimSiPMData> sipmData = sipmItr.SiPMs;
      int nModule = sipmItr.ModuleNum;

      for (int i = 0; i < sipmData.size(); i++) {

        DRsimInterface::DRsimTimeStruct timeItr = sipmData[i].timeStruct;

        for(auto TmpItr = timeItr.begin(); TmpItr != timeItr.end(); ++TmpItr) {
          auto timeData = *TmpItr;
          if(DRsimInterface::IsCerenkov(nModule)) {
            tCtime->Fill((timeData.first.first + timeData.first.second)/2., timeData.second);
            ftC_hits += timeData.second;
          } else {
            tStime->Fill((timeData.first.first + timeData.first.second)/2., timeData.second);
            ftS_hits += timeData.second;
          }
        }
      }
    }

    tEdep->Fill(ftEdep);
    tChit->Fill(ftC_hits);
    tShit->Fill(ftS_hits);
  }

  TCanvas* c = new TCanvas("c","");

  tEdep->Draw("Hist"); c->SaveAs("_TotalEdep.png");
  tChit->Draw("Hist"); c->SaveAs("_TotalChit.png");
  tShit->Draw("Hist"); c->SaveAs("_TotalShit.png");
  tCtime->Draw("Hist"); c->SaveAs("_TotalCtime.png");
  tStime->Draw("Hist"); c->SaveAs("_TotalStime.png");
}
