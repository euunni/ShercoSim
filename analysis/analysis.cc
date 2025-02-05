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
#include "TFile.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <tuple>
#include "Riostream.h"


int main(int argc, char* argv[]) {

  TString filename = argv[1];
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);

  const int row = 27; 
  const int col = 27;
  int numModule = row * col;

  gStyle->SetOptFit(1);

  TH1F* tEdep = new TH1F("Total_Edep","Total Energy deposit;MeV;Evt",100,low,1000.*high);
  tEdep->Sumw2(); tEdep->SetLineColor(kBlack); tEdep->SetLineWidth(2);
  TH1F* tCtime = new TH1F("Total_C_Time","Total timing of Cerenkov.;ns;Evt",150,0,30);
  tCtime->Sumw2(); tCtime->SetLineColor(kBlue); tCtime->SetLineWidth(2);
  TH1F* tStime = new TH1F("Total_S_Time","Total timing of Scintillation.;ns;Evt",150,0,30);
  tStime->Sumw2(); tStime->SetLineColor(kRed); tStime->SetLineWidth(2);
  TH1I* tChit = new TH1I("Total_C_Hit","Total hits of Cerenkov",100,1000.,10000.) ;
  tChit->Sumw2(); tChit->SetLineColor(kBlue); tChit->SetLineWidth(2);
  TH1I* tShit = new TH1I("Total_S_Hit","Total hits of Scintillation",100,100000.,400000.);
  tShit->Sumw2(); tShit->SetLineColor(kRed); tShit->SetLineWidth(2);
  TH1F* tP_leak = new TH1F("Pleak","Momentum leak;MeV;Evt",100,0.,1000.*high);
  tP_leak->Sumw2(); tP_leak->SetLineWidth(2);
  TH1F* tP_leak_nu = new TH1F("Pleak_nu","Neutrino energy leak;MeV;Evt",100,0.,1000.*high);
  tP_leak_nu->Sumw2(); tP_leak_nu->SetLineWidth(2);
  
  TH2F* tEdep_2D = new TH2F("Edep_2D",";;", col, 0, col, row, 0, row);
  TH2F* tHits_2D = new TH2F("Hits_2D",";;", col, 0, col, row, 0, row);
  TH2F* tE_2D = new TH2F("Energy_2D",";;", col, 0, col, row, 0, row);

  TH1F* tEdep_Towers[numModule];
  TH1F* tHits_Towers[numModule];
  TString nameEdep;
  TString nameHits;
  for (int i = 0; i < numModule; i++) {
    nameEdep = std::to_string(i) + "_Tower_Edep";
    nameHits = std::to_string(i) + "_Tower_Hits";
    tEdep_Towers[i] = new TH1F(nameEdep, ";MeV;Evt", 1000, 0., 1000000.);
    tHits_Towers[i] = new TH1F(nameHits, ";Npe;Evt", 1000, 0., 1000000.);
  }

  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>("/ ** Path for yours ** /" + std::string(filename) + ".root", 1);
  drInterface->set("DRsim","DRsimEventData");

  unsigned int entries = drInterface->entries();
  while (drInterface->numEvt() < entries) {

    if (drInterface->numEvt() % 100 == 0) printf("Analyzing %dth event ...\n", drInterface->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface->read(drEvt);

    float ftEdep = 0.;
    float Edep_Towers[numModule] = {0};
    float hits_Towers[numModule] = {0};

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

    int ftC_hits = 0;
    int ftS_hits = 0;

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
            hits_Towers[nModule] += timeData.second;
          } else {
            tStime->Fill((timeData.first.first + timeData.first.second)/2., timeData.second);
            ftS_hits += timeData.second;
            hits_Towers[nModule] += timeData.second;
          }
        }
      }
    }

    tEdep->Fill(ftEdep);
    tChit->Fill(ftC_hits);
    tShit->Fill(ftS_hits);

    for (int i = 0; i < numModule; i++) {
      int xIdx = (i / row);
      int yIdx = (i % row);
    
      tEdep_Towers[i]->Fill(Edep_Towers[i]);
      tHits_Towers[i]->Fill(hits_Towers[i]);
      tEdep_2D->Fill(xIdx, yIdx, Edep_Towers[i]/entries);
      tHits_2D->Fill(xIdx, yIdx, hits_Towers[i]/entries);
    }
  } // End of event loop


  std::vector<std::pair<int, double>> dataEdep;
  for (int i = 0; i < numModule; i++ ) {
    dataEdep.push_back(std::make_pair(i, tEdep_Towers[i]->GetMean()));
  }

  std::vector<std::pair<int, double>> dataHits;
  for (int i = 0; i < numModule; i++ ) {
    dataHits.push_back(std::make_pair(i, tHits_Towers[i]->GetMean()));
  }

  std::sort(dataEdep.begin(), dataEdep.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.first < b.first;
  });

  std::sort(dataHits.begin(), dataHits.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.first < b.first;
  });

  std::ofstream outEdep;
  outEdep.open("/ ** Path for yours ** /" + filename + "_Edep.csv", std::ios::out | std::ios::app);
  outEdep << "Total Edep : " << tEdep->GetMean() << " MeV" << std::endl;
  for (const auto& itr : dataEdep) {
    outEdep << "Module_" << (itr.first) << " " << itr.second << std::endl;
  }

  std::ofstream outHits;
  outHits.open("/ ** Path for yours ** /" + filename + "_Hits.csv", std::ios::out | std::ios::app);
  outHits << "Total Chits : " << tChit->GetMean() << std::endl;
  outHits << "Total Shits : " << tShit->GetMean() << std::endl;
  for (const auto& itr : dataHits) {
    if (itr.first % 2 == 0) {
      outHits << "C " << "Module_" << (itr.first) << " " << itr.second << std::endl;
    } else {
      outHits << "S " << "Module_" << (itr.first) << " " << itr.second << std::endl;
    }
  }

  TCanvas* c = new TCanvas("c","");

  c->SetLogy(1);
  tP_leak->Draw("Hist"); c->SaveAs("/ ** Path for yours ** /" + filename + "_Pleak.png");
  tP_leak_nu->Draw("Hist"); c->SaveAs("/ ** Path for yours ** /" + filename + "_Pleak_nu.png");
  
  c->SetLogy(0);
  tEdep->Draw("Hist"); c->SaveAs("/ ** Path for yours ** /" + filename + "_TotalEdep.png");
  tChit->Draw("Hist"); c->SaveAs("/ ** Path for yours ** /" + filename + "_TotalChit.png");
  tShit->Draw("Hist"); c->SaveAs("/ ** Path for yours ** /" + filename + "_TotalShit.png");
  tCtime->Draw("Hist"); c->SaveAs("/ ** Path for yours ** /" + filename + "_TotalCtime.png");
  tStime->Draw("Hist"); c->SaveAs("/ ** Path for yours ** /" + filename + "_TotalStime.png");

  gStyle->SetPaintTextFormat("4.1f");
  c->cd();
  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.15);
  c->SetCanvasSize(1800,1400);

  tEdep_2D->SetMarkerSize(0.4);
  tEdep_2D->GetXaxis()->SetLabelFont(42);
  tEdep_2D->GetYaxis()->SetLabelFont(42);
  tEdep_2D->GetXaxis()->SetLabelSize(0.025);
  tEdep_2D->GetYaxis()->SetLabelSize(0.025);

  tHits_2D->SetMarkerSize(0.4);
  tHits_2D->GetXaxis()->SetLabelFont(42);
  tHits_2D->GetYaxis()->SetLabelFont(42);
  tHits_2D->GetXaxis()->SetLabelSize(0.025);
  tHits_2D->GetYaxis()->SetLabelSize(0.025);

  for (int i = 1; i <= col; i++) {
    tEdep_2D->GetXaxis()->SetBinLabel(i, std::to_string(i).c_str());
    tHits_2D->GetXaxis()->SetBinLabel(i, std::to_string(i).c_str());
    tE_2D->GetXaxis()->SetBinLabel(i, std::to_string(i).c_str());
  }
  for (int i = 1; i <= row; i++) {
    tEdep_2D->GetYaxis()->SetBinLabel(i, std::to_string(i).c_str());
    tHits_2D->GetYaxis()->SetBinLabel(i, std::to_string(i).c_str());
    tE_2D->GetYaxis()->SetBinLabel(i, std::to_string(i).c_str());
  }

  tEdep_2D->Draw("COL0Z TEXT");
  tEdep_2D->SetStats(0);
  c->SaveAs("/ ** Path for yours ** /" + filename + "_Edep2D.pdf");

  c->SetLogz(1);

  tEdep_2D->Draw("COL0Z TEXT"); 
  tEdep_2D->SetStats(0);
  c->SaveAs("/ ** Path for yours ** /" + filename + "_Edep2D_Log.pdf"); 

  tHits_2D->Draw("COL0Z TEXT"); 
  tHits_2D->SetStats(0);
  c->SaveAs("/ ** Path for yours ** /" + filename + "_Hits2D_Log.pdf"); 
}

