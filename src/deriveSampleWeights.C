//AUTHOR: Chris McGinn (2019.12.06)

//cpp
#include <fstream>
#include <iostream>
#include <map>
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
#include "include/getLinBins.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int deriveSampleWeights(std::string inXSectionFileName, std::string inEntriesFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inXSectionFileName, ".txt")) return 1;

  bool entriesFileTXT = true;
  if(!check.checkFileExt(inEntriesFileName, ".txt")){
    entriesFileTXT = false;
    if(!check.checkFileExt(inEntriesFileName, ".root")) return 1;
  }

  kirchnerPalette kPal;

  const std::string dateStr = getDateStr();

  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);

  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::map<std::string, std::vector<double> > jzMapToXSec;
  std::map<std::string, double> jzMapToEntries;
  std::map<std::string, double> jzMapToWeights;
  std::ifstream inXSectionFile(inXSectionFileName.c_str());
  std::string tempStr;
  while(std::getline(inXSectionFile, tempStr)){
    if(tempStr.size() == 0) continue;
    std::vector<std::string> tempVect = commaSepStringToVect(tempStr);
    if(tempVect.size() == 0) continue;
    if(tempVect[0].size() == 0) continue;
    if(tempVect[0].substr(0,1).find("#") != std::string::npos) continue;

    jzMapToXSec[tempVect[0]] = {std::stod(tempVect[1]), std::stod(tempVect[2])};
    jzMapToEntries[tempVect[0]] = 0.0;
  }
  inXSectionFile.close();

  std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(entriesFileTXT){
    std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    std::ifstream inEntriesFile(inEntriesFileName.c_str());
    while(std::getline(inEntriesFile, tempStr)){
      if(tempStr.size() == 0) continue;
      std::vector<std::string> tempVect = commaSepStringToVect(tempStr);
      if(tempVect.size() == 0) continue;
      if(tempVect[0].size() == 0) continue;
      if(tempVect[0].substr(0,1).find("#") != std::string::npos) continue;
      
      jzMapToEntries[tempVect[0]] = std::stod(tempVect[1]);
    }
    inEntriesFile.close();
    std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  }
  else{
    std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    TFile* inFile_p = new TFile(inEntriesFileName.c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("clusterJetsCS");
    Int_t jzVal_;

    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("jzVal", 1);

    inTree_p->SetBranchAddress("jzVal", &jzVal_);

    const Int_t nEntries = inTree_p->GetEntries();

    for(Int_t entry = 0; entry < nEntries; ++entry){
      inTree_p->GetEntry(entry);

      std::string jzValStr = "JZ" + std::to_string(jzVal_);
      if(jzMapToEntries.count(jzValStr) == 0) std::cout << "WARNING MAP IS MISSING VALUE " << jzValStr << std::endl;
      ++(jzMapToEntries[jzValStr]);      
    }

    inFile_p->Close();
    delete inFile_p;

    std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  }
  
  std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  double maxWeight = -1.0;
  for(auto const& iter : jzMapToXSec){
    if(jzMapToEntries[iter.first] <= 0) continue;

    std::cout << iter.first << ", " << iter.second[0] << ", " << iter.second[1] << ", " << jzMapToEntries[iter.first] << std::endl;

    jzMapToWeights[iter.first] = iter.second[0]*iter.second[1]/jzMapToEntries[iter.first];
    if(maxWeight < jzMapToWeights[iter.first]) maxWeight = jzMapToWeights[iter.first];

    std::cout << " " << jzMapToWeights[iter.first] << std::endl;
  }
  
  std::cout << "Maxweight: " << maxWeight << std::endl;

  std::cout << "RENORMALIZED WEIGHTS: " << std::endl;
  for(auto const& iter : jzMapToXSec){
    std::cout << iter.first << ": " << jzMapToWeights[iter.first]/maxWeight << " (" << jzMapToEntries[iter.first] << ")" << std::endl;
    jzMapToWeights[iter.first] /= maxWeight;
  }

  std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::string outFileName = "";
  if(!entriesFileTXT){
    outFileName = inEntriesFileName.substr(0, inEntriesFileName.rfind(".root"));
    while(outFileName.find("/") != std::string::npos){
      outFileName.replace(0, outFileName.find("/")+1, "");
    }
    outFileName = "output/" + dateStr + "/" + outFileName + "_DerivedWeights_" + dateStr + ".root";
    
    const Int_t nJtPtBins = 110;
    Float_t jtPtLow = 20.;
    Float_t jtPtHigh = 460.;
    Double_t jtPtBins[nJtPtBins+1];
    getLinBins(jtPtLow, jtPtHigh, nJtPtBins, jtPtBins);
    
    TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
    TH1D* jtpt_h = new TH1D("jtpt_h", ";Jet p_{T} [GeV];Effective Counts", nJtPtBins, jtPtBins);
    TH1D* jtpt_Weighted_h = new TH1D("jtpt_Weighted_h", ";Jet p_{T} [GeV];Effective Counts (Weighted)", nJtPtBins, jtPtBins);
    centerTitles({jtpt_h, jtpt_Weighted_h});

    TFile* inFile_p = new TFile(inEntriesFileName.c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("clusterJetsCS");
    Int_t jzVal_;

    const Int_t nMaxJets = 500;
    Int_t njtTruth_;
    Float_t jtptTruth_[nMaxJets];

    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("jzVal", 1);
    inTree_p->SetBranchStatus("njtTruth", 1);
    inTree_p->SetBranchStatus("jtptTruth", 1);

    inTree_p->SetBranchAddress("jzVal", &jzVal_);
    inTree_p->SetBranchAddress("njtTruth", &njtTruth_);
    inTree_p->SetBranchAddress("jtptTruth", jtptTruth_);

    const Int_t nEntries = inTree_p->GetEntries();

    for(Int_t entry = 0; entry < nEntries; ++entry){
      inTree_p->GetEntry(entry);

      std::string jzValStr = "JZ" + std::to_string(jzVal_);
      if(jzMapToEntries.count(jzValStr) == 0) std::cout << "WARNING MAP IS MISSING VALUE " << jzValStr << std::endl;

      Double_t weight = jzMapToWeights[jzValStr];      

      for(Int_t jI = 0; jI < njtTruth_; ++jI){
	jtpt_h->Fill(jtptTruth_[jI]);
	jtpt_Weighted_h->Fill(jtptTruth_[jI], weight);
      }
    }

    inFile_p->Close();
    delete inFile_p;

    outFile_p->cd();

    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.14);
    canv_p->SetBottomMargin(0.14);

    jtpt_h->SetMarkerColor(1);
    jtpt_h->SetLineColor(1);
    jtpt_h->SetMarkerStyle(24);
    jtpt_h->SetMarkerSize(1);

    jtpt_Weighted_h->SetMarkerColor(kPal.getColor(1));
    jtpt_Weighted_h->SetLineColor(kPal.getColor(1));
    jtpt_Weighted_h->SetMarkerStyle(25);
    jtpt_Weighted_h->SetMarkerSize(1);

    jtpt_h->SetMinimum(getMinGTZero(jtpt_Weighted_h)/2.);

    jtpt_h->DrawCopy("HIST E1 P");
    jtpt_Weighted_h->DrawCopy("HIST E1 P SAME");

    gPad->SetLogy();
    gStyle->SetOptStat(0);

    TLegend* leg_p = new TLegend(0.2, 0.2, 0.4, 0.4);
    leg_p->SetFillStyle(0);
    leg_p->SetFillColor(0);
    leg_p->SetBorderSize(0);
    leg_p->SetTextFont(43);
    leg_p->SetTextSize(16);

    leg_p->AddEntry(jtpt_h, "Unweighted", "P L");
    leg_p->AddEntry(jtpt_Weighted_h, "Weighted", "P L");

    leg_p->Draw("SAME");

    std::string saveName = "pdfDir/" + dateStr + "/sampleWeights_" + dateStr + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;
    delete leg_p;

    jtpt_h->Write("", TObject::kOverwrite);
    delete jtpt_h;

    jtpt_Weighted_h->Write("", TObject::kOverwrite);
    delete jtpt_Weighted_h;

    outFile_p->Close();
    delete outFile_p;
  }
  else{
    outFileName = inEntriesFileName.substr(0, inEntriesFileName.rfind(".txt"));
    while(outFileName.find("/") != std::string::npos){
      outFileName.replace(0, outFileName.find("/")+1, "");
    }
    outFileName = "output/" + dateStr + "/" + outFileName + "_DerivedWeights_" + dateStr + ".root";
  }

  std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  outFileName.replace(outFileName.rfind(".root"), 5, "");
  outFileName = outFileName + ".txt";
  
  std::ofstream outFile(outFileName.c_str());

  for(auto const& iter : jzMapToWeights){
    outFile << iter.first << "," << iter.second << " (" << jzMapToEntries[iter.first] << ")" << std::endl;
  }

  outFile.close();

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/deriveSampleWeights.exe <inXSectionFileName> <inEntriesFileName>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += deriveSampleWeights(argv[1], argv[2]);
  return retVal;
}
