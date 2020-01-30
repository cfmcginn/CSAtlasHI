//Author: Chris McGinn (2019.12.06)

//cpp
#include <fstream>
#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/ncollFunctions_5TeV.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int deriveCentWeights(std::string inFileName, std::string jzWeightsName = "")
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, ".root")) return 1;
  const bool doJZWeights = check.checkFileExt(jzWeightsName, "txt");
  std::map<std::string, double> weightMap;

  std::cout << "JZWEIGHTSNAME: " << jzWeightsName << std::endl;
  if(doJZWeights){
    std::ifstream inFile(jzWeightsName.c_str());
    std::string tempStr;
    while(std::getline(inFile, tempStr)){
      if(tempStr.size() == 0) continue;
      std::vector<std::string> tempVect = commaSepStringToVect(tempStr);
      if(tempVect.size() == 0) continue;
      if(tempVect[0].size() == 0) continue;
      if(tempVect[0].substr(0, 1).find("#") != std::string::npos) continue;

      weightMap[tempVect[0]] = std::stod(tempVect[1]);
    }
    inFile.close();

    std::cout << "WEIGHTS: " << std::endl;
    for(auto const& iter : weightMap){
      std::cout << " " << iter.first << ", " << iter.second << std::endl;
    }
  }

  const std::string dateStr = getDateStr();

  kirchnerPalette kPal;

  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);

  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  std::string outFileName = inFileName.substr(0, inFileName.rfind(".root"));
  while(outFileName.find("/") != std::string::npos){
    outFileName.replace(0, outFileName.find("/")+1, "");
  }
  outFileName = "output/" + dateStr + "/" + outFileName + "_CentWeights_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* cent_h = new TH1D("cent_h", ";Centrality (%);Effective Counts", 100, -0.5, 99.5);
  TH1D* cent_JZWeights_h = new TH1D("cent_JZWeights_h", ";Centrality (%);Effective Counts (JZ Weighted)", 100, -0.5, 99.5);
  TH1D* cent_FullWeights_h = new TH1D("cent_FullWeights_h", ";Centrality (%);Effective Counts (Full Weighted)", 100, -0.5, 99.5);

  centerTitles({cent_h, cent_JZWeights_h, cent_FullWeights_h});

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("clusterJetsCS");

  Int_t jzVal_;
  Float_t cent_;

  inTree_p->SetBranchStatus("*", 0);
  if(doJZWeights) inTree_p->SetBranchStatus("jzVal", 1);
  inTree_p->SetBranchStatus("cent", 1);

  if(doJZWeights) inTree_p->SetBranchAddress("jzVal", &jzVal_);
  inTree_p->SetBranchAddress("cent", &cent_);

  const Int_t nEntries = inTree_p->GetEntries();

  for(Int_t entry = 0; entry < nEntries; ++entry){
    inTree_p->GetEntry(entry);

    Double_t weight = 1.0;
    if(doJZWeights){
      std::string jzValStr = "JZ" + std::to_string(jzVal_);
      if(weightMap.count(jzValStr) == 0) std::cout << "WARNING CANNOT FIND WEIGHT FOR JTVAL \'" << jzValStr << "\'" << std::endl;
      weight = weightMap[jzValStr];
    }

    cent_h->Fill(cent_);
    cent_JZWeights_h->Fill(cent_, weight);    
  }

  std::vector<double> centWeights;
  for(Int_t cI = 0; cI < 100; ++cI){
    double ave = (findNcoll_Renorm(cI*2) + findNcoll_Renorm(cI*2 + 1))/2.;
    centWeights.push_back(ave/cent_JZWeights_h->GetBinContent(cI+1));
  }

  double maxWeight = centWeights[0];
  for(Int_t cI = 0; cI < 100; ++cI){
    centWeights[cI] /= maxWeight;
  }


  for(Int_t entry = 0; entry < nEntries; ++entry){
    inTree_p->GetEntry(entry);

    Int_t centInt = (Int_t)cent_;
    Double_t weight = centWeights[centInt];
    if(doJZWeights){
      std::string jzValStr = "JZ" + std::to_string(jzVal_);
      if(weightMap.count(jzValStr) == 0) std::cout << "WARNING CANNOT FIND WEIGHT FOR JTVAL \'" << jzValStr << "\'" << std::endl;
      weight *= weightMap[jzValStr];
    }


    cent_FullWeights_h->Fill(cent_, weight);
  }

  outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = outFileName + ".txt";
  std::ofstream outFile(outFileName.c_str());
  outFile << "#LowVal,HighVal,Weight" << std::endl;
  for(Int_t cI = 0; cI < 100; ++cI){
    outFile << cI << "," << cI+1 << "," << centWeights[cI] << std::endl;
  }
  outFile.close();

  inFile_p->Close();
  delete inFile_p;
  
  outFile_p->cd();

  TLegend* leg_p = new TLegend(0.2, 0.2, 0.4, 0.4);
  leg_p->SetFillStyle(0);
  leg_p->SetFillColor(0);
  leg_p->SetBorderSize(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(16);

  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.14);
  canv_p->SetBottomMargin(0.14);

  cent_JZWeights_h->SetMarkerSize(1);
  cent_JZWeights_h->SetMarkerStyle(24);
  cent_JZWeights_h->SetMarkerColor(1);
  cent_JZWeights_h->SetLineColor(1);

  cent_FullWeights_h->SetMarkerSize(1);
  cent_FullWeights_h->SetMarkerStyle(24);
  cent_FullWeights_h->SetMarkerColor(kPal.getColor(1));
  cent_FullWeights_h->SetLineColor(kPal.getColor(1));

  cent_JZWeights_h->SetMinimum(TMath::Min(getMinGTZero(cent_FullWeights_h), getMinGTZero(cent_JZWeights_h))/2.);

  cent_JZWeights_h->DrawCopy("HIST E1 P");
  cent_FullWeights_h->DrawCopy("HIST E1 P SAME");

  gStyle->SetOptStat(0);
  gPad->SetLogy();

  leg_p->AddEntry(cent_JZWeights_h, "Sample Weighted", "P L");
  leg_p->AddEntry(cent_FullWeights_h, "Sample+Centrality Weighted", "P L");

  leg_p->Draw("SAME");

  std::string saveName = "pdfDir/" + dateStr + "/centWeights_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;
  delete leg_p;

  cent_h->Write("", TObject::kOverwrite);
  cent_JZWeights_h->Write("", TObject::kOverwrite);
  cent_FullWeights_h->Write("", TObject::kOverwrite);

  delete cent_h;
  delete cent_JZWeights_h;
  delete cent_FullWeights_h;

  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 3){
    std::cout << "Usage: ./bin/deriveCentWeights.exe <inFileName> <jzWeightsName-default=\'\'>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += deriveCentWeights(argv[1]);
  else if(argc == 3) retVal += deriveCentWeights(argv[1], argv[2]);
  return retVal;
}
