//Author: Chris McGinn (2019.12.06)

//cpp
#include <fstream>
#include <iostream>
#include <string>

//ROOT
#include "TDatime.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/ncollFunctions_5TeV.h"
#include "include/stringUtil.h"

int deriveCentWeights(std::string inFileName, std::string jzWeightsName = "")
{
  if(!checkFileExt(inFileName, ".root")) return 1;
  const bool doJZWeights = checkFileExt(jzWeightsName, "txt");
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

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);

  std::string outFileName = inFileName.substr(0, inFileName.rfind(".root"));
  while(outFileName.find("/") != std::string::npos){
    outFileName.replace(0, outFileName.find("/")+1, "");
  }
  outFileName = "output/" + dateStr + "/" + outFileName + "_CentWeights_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* cent_h = new TH1D("cent_h", ";Centrality (%);Counts", 100, -0.5, 99.5);
  TH1D* cent_JZWeights_h = new TH1D("cent_JZWeights_h", ";Centrality (%);Counts (JZ Weighted)", 100, -0.5, 99.5);
  TH1D* cent_FullWeights_h = new TH1D("cent_FullWeights_h", ";Centrality (%);Counts (Full Weighted)", 100, -0.5, 99.5);
  
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
