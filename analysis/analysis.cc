#include "RootInterface.h"
#include "RecoInterface.h"
#include "DRsimInterface.h"
#include "functions.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <tuple>
#include "Riostream.h"


int main(int argc, char* argv[]) {

  TString particle = argv[1]; // ele, mu, pi, pro
  TString energy = argv[2];
  TString filename = argv[3];

  float high = 40.;
  float low = 0.;

  const int row = 27; 
  const int col = 27;
  int numModule = row * col;

  gStyle->SetOptFit(1);

  TH1F* tEdep = new TH1F("Total_Edep","Total Energy deposit;MeV;Evt",100,0.,50000.);
  tEdep->Sumw2(); tEdep->SetLineColor(kBlack); tEdep->SetLineWidth(2);
  TH1F* tCtime = new TH1F("Total_C_Time","Total timing of Cerenkov ch.;ns;Evt",150,0,30);
  tCtime->Sumw2(); tCtime->SetLineColor(kBlue); tCtime->SetLineWidth(2);
  TH1F* tStime = new TH1F("Total_S_Time","Total timing of Scintillation ch.;ns;Evt",150,0,30);
  tStime->Sumw2(); tStime->SetLineColor(kRed); tStime->SetLineWidth(2);
  TH1I* tChit = new TH1I("Total_C_Hit","Total hits of Cerenkov ch",100,0.,10000.) ;
  tChit->Sumw2(); tChit->SetLineColor(kBlue); tChit->SetLineWidth(2);
  TH1I* tShit = new TH1I("Total_S_Hit","Total hits of Scintillation ch",100,0.,1000000.);
  tShit->Sumw2(); tShit->SetLineColor(kRed); tShit->SetLineWidth(2);
  TH1F* tP_leak = new TH1F("Pleak","Momentum leak;MeV;Evt",100,1000.*low,1000.*high);
  tP_leak->Sumw2(); tP_leak->SetLineWidth(2);
  TH1F* tP_leak_nu = new TH1F("Pleak_nu","Neutrino energy leak;MeV;Evt",100,0.,1000.*high);
  tP_leak_nu->Sumw2(); tP_leak_nu->SetLineWidth(2);
  
  TH2F* tEdep_2D = new TH2F("Edep_2D",";;", col, 0, col, row, 0, row);
  TH1F* tEdep_Towers[numModule];
  TString name;
  for (int i = 0; i < numModule; i++) {
    name = std::to_string(i) + "_Tower_Edep";
    tEdep_Towers[i] = new TH1F(name, ";MeV;Evt", 1000, 0., 1000000.);
  }

  // RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>("/u/user/haeun/Sherco/v240628/ShercoSim/install/input/241128/27by27/" + std::string(particle) + "/" + std::string(energy) + "/" + std::string(filename) + ".root", 1);
  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>("/u/user/haeun/Sherco/v240628/ShercoSim/install/input/241128/27by27/ele/" + std::string(filename) + ".root", 1);

  drInterface->set("DRsim","DRsimEventData");

  unsigned int entries = drInterface->entries();
  while (drInterface->numEvt() < entries) {

    if (drInterface->numEvt() % 100 == 0) printf("Analyzing %dth event ...\n", drInterface->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface->read(drEvt);

    float ftEdep = 0.;
    float Edep_Towers[numModule] = {0};
    // int hits_Towers[numModule] = {0};

    for (auto edepItr = drEvt.Edeps.begin(); edepItr != drEvt.Edeps.end(); ++edepItr) {
      auto edep = *edepItr;
      ftEdep += edep.Edep;

      int moduleNum = edep.ModuleNum;
      Edep_Towers[moduleNum] += edep.Edep;
    }

    float Pleak = 0.;
    float muak_nu = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      if ( std::abs(leak.pdgId)==12 || std::abs(leak.pdgId)==14 || std::abs(leak.pdgId)==16 ) {
        muak_nu += leak4vec.P();
      } else {
        Pleak += leak4vec.P();
      }
    }
    tP_leak->Fill(Pleak);
    tP_leak_nu->Fill(muak_nu);

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

    for (int i = 0; i < numModule; i++) {
      tEdep_Towers[i]->Fill(Edep_Towers[i]);
    }

    for (int i = 0; i < numModule; i++) {
      int xIdx = (i / row);
      int yIdx = (i % row);
      tEdep_2D->Fill(xIdx, yIdx, Edep_Towers[i]/entries);
      // std::cout << "Module" << i << " (" << xIdx << ", " << yIdx << ") " << "-> Edep : " << Edep_Towers[i]/entries << std::endl;
    }
  } // End of event loop

  std::vector<std::pair<int, double>> data;
  for (int i = 0; i < numModule; i++ ) {
      data.push_back(std::make_pair(i, tEdep_Towers[i]->GetMean()));
  }

  std::sort(data.begin(), data.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.first < b.first;
  });

  std::ofstream out;
  out.open("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_Edep.csv", std::ios::out | std::ios::app);
  out << "Total Edep : " << tEdep->GetMean() << " MeV" << std::endl;
  for (const auto& itr : data) {
    out << "Module_" << (itr.first)+1 << " " << itr.second << std::endl;
  }

  TCanvas* c = new TCanvas("c","");

  c->SetLogy(1);
  tP_leak->Draw("Hist"); c->SaveAs("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_Pleak.png");
  tP_leak_nu->Draw("Hist"); c->SaveAs("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_Pleak_nu.png");
  c->SetLogy(0);

  tEdep->Draw("Hist"); c->SaveAs("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_TotalEdep.png");
  tChit->Draw("Hist"); c->SaveAs("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_TotalChit.png");
  tShit->Draw("Hist"); c->SaveAs("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_TotalShit.png");
  tCtime->Draw("Hist"); c->SaveAs("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_TotalCtime.png");
  tStime->Draw("Hist"); c->SaveAs("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_TotalStime.png");

  gStyle->SetPaintTextFormat("4.1f");
  c->cd();
  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.15);
  tEdep_2D->GetXaxis()->SetLabelFont(42);
  tEdep_2D->GetYaxis()->SetLabelFont(42);

  // For 4 by 5
  // c->SetCanvasSize(1200,1200);
  // tEdep_2D->SetMarkerSize(1.2);

  // For 27 by 27
  c->SetCanvasSize(1800,1400);
  tEdep_2D->SetMarkerSize(0.4);
  tEdep_2D->GetXaxis()->SetLabelSize(0.025);
  tEdep_2D->GetYaxis()->SetLabelSize(0.025);

  for (int i = 1; i <= col; i++) {
    tEdep_2D->GetXaxis()->SetBinLabel(i, std::to_string(i).c_str());
  }
  for (int i = 1; i <= row; i++) {
    tEdep_2D->GetYaxis()->SetBinLabel(i, std::to_string(i).c_str());
  }

  tEdep_2D->Draw("COL0Z text");
  tEdep_2D->SetStats(0);
  c->SaveAs("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_Edep2D.pdf");

  c->SetLogz(1);
  tEdep_2D->Draw("COL0Z TEXT"); 
  tEdep_2D->SetStats(0);
  c->SaveAs("./plot/241128/27by27/" + particle + "/" + + energy + "/" + filename + "_Edep2D_Log.pdf"); 
  c->SetLogz(0);
}
