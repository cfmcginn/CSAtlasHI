//Author: Chris McGinn (2020.02.05)

//cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TFile.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/stringUtil.h"

int validateRho(std::string rhoFileName, std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(rhoFileName, ".root")) return 1;
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);
  
  Int_t run_, evt_;
  UInt_t lumi_;

  std::vector<float>* etaMin_p=nullptr;
  std::vector<float>* etaMax_p=nullptr;
  std::vector<float>* et_p=nullptr;

  std::vector<float>* towers_pt_p=nullptr;
  std::vector<float>* towers_eta_p=nullptr;
  std::vector<float>* towers_phi_p=nullptr;

  std::map<std::string, ULong64_t> runLumiEvtToRhoEntry;
  
  TFile* rhoFile_p = new TFile(rhoFileName.c_str(), "READ");
  TTree* rhoTree_p = (TTree*)rhoFile_p->Get("rhoTree");

  rhoTree_p->SetBranchStatus("*", 0);
  rhoTree_p->SetBranchStatus("runNumber", 1);
  rhoTree_p->SetBranchStatus("lumiBlock", 1);
  rhoTree_p->SetBranchStatus("evtNumber", 1);

  rhoTree_p->SetBranchAddress("runNumber", &run_);
  rhoTree_p->SetBranchAddress("lumiBlock", &lumi_);
  rhoTree_p->SetBranchAddress("evtNumber", &evt_);

  const ULong64_t nRhoEntries = rhoTree_p->GetEntries();
  for(ULong64_t entry = 0; entry < nRhoEntries; ++entry){
    rhoTree_p->GetEntry(entry);

    std::string runLumiEvtStr = std::to_string(run_) + "_" + std::to_string(lumi_) + "_" + std::to_string(evt_);

    runLumiEvtToRhoEntry[runLumiEvtStr] = entry;
  }  

  rhoTree_p->SetBranchStatus("*", 0);
  rhoTree_p->SetBranchStatus("etaMin", 1);
  rhoTree_p->SetBranchStatus("etaMax", 1);
  rhoTree_p->SetBranchStatus("et", 1);

  rhoTree_p->SetBranchAddress("etaMin", &etaMin_p);
  rhoTree_p->SetBranchAddress("etaMax", &etaMax_p);
  rhoTree_p->SetBranchAddress("et", &et_p);

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("jetTree");

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("runNumber", 1);
  inTree_p->SetBranchStatus("lumiBlock", 1);
  inTree_p->SetBranchStatus("eventNumber", 1);
  inTree_p->SetBranchStatus("towers_pt", 1);
  inTree_p->SetBranchStatus("towers_phi", 1);
  inTree_p->SetBranchStatus("towers_eta", 1);
  
  inTree_p->SetBranchAddress("runNumber", &run_);
  inTree_p->SetBranchAddress("lumiBlock", &lumi_);
  inTree_p->SetBranchAddress("eventNumber", &evt_);
  inTree_p->SetBranchAddress("towers_pt", &towers_pt_p);
  inTree_p->SetBranchAddress("towers_phi", &towers_phi_p);
  inTree_p->SetBranchAddress("towers_eta", &towers_eta_p);

  std::vector<float> perEtaEnergySum;  
  std::vector<float> fullEtaBins;
  inTree_p->GetEntry(0);
  for(unsigned int eI = 0; eI < etaMin_p->size(); ++eI){
    perEtaEnergySum.push_back(0.0);
    fullEtaBins.push_back(etaMin_p->at(eI));
  }
  fullEtaBins.push_back(etaMax_p->at(etaMax_p->size()-1));

  std::cout << "ETABINS: ";
  for(unsigned int eI = 0; eI < fullEtaBins.size(); ++eI){
    std::cout << fullEtaBins[eI] << ", ";
  }
  std::cout << std::endl;
  
  const ULong64_t nEntries = inTree_p->GetEntries();
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    inTree_p->GetEntry(entry);

    std::string runLumiEvtStr = std::to_string(run_) + "_" + std::to_string(lumi_) + "_" + std::to_string(evt_);

    if(runLumiEvtToRhoEntry.count(runLumiEvtStr) == 0){
      std::cout << "CANNOT FIND ENTRY (Run_Lumi_Evt) " << runLumiEvtStr << ", continue" << std::endl;
      continue;
    }
    ULong64_t rhoEntry = runLumiEvtToRhoEntry[runLumiEvtStr];

    rhoTree_p->GetEntry(rhoEntry);

    for(unsigned int eI = 0; eI < etaMin_p->size(); ++eI){perEtaEnergySum[eI] = 0.0;}

    for(unsigned int tI = 0; tI < towers_pt_p->size(); ++tI){
      
    }
  }
  
  inFile_p->Close();
  delete inFile_p;
  
  rhoFile_p->Close();
  delete rhoFile_p;

  std::cout << "VALIDATE RHO COMPLETE. return 0" << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/validateRho.exe <rhoFileName> <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += validateRho(argv[1], argv[2]);
  return retVal;
}
