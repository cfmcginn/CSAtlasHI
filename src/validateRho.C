//Author: Chris McGinn (2020.02.05)

//c+cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int validateRho(std::string rhoFileName, std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(rhoFileName, ".root")) return 1;
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  //DEBUG BOOL FROM ENV VAR
  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);
  
  const std::string outFileName = "output/" + dateStr + "/validateRho_" + dateStr + ".root";
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p = new TTree("validateRhoTree", "");  

  std::vector<float>* etRecalc_p = new std::vector<float>;
  std::vector<float>* etATLAS_p = new std::vector<float>;

  outTree_p->Branch("etRecalc", &etRecalc_p);
  outTree_p->Branch("etATLAS", &etATLAS_p);
  
  Int_t run_, evt_;
  UInt_t lumi_;

  std::vector<float>* etaMin_p=nullptr;
  std::vector<float>* etaMax_p=nullptr;
  std::vector<float>* et_p=nullptr;

  std::vector<float>* towers_pt_p=nullptr;
  std::vector<float>* towers_eta_p=nullptr;
  std::vector<float>* towers_phi_p=nullptr;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
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

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  std::vector<float> fullEtaBins;
  rhoTree_p->GetEntry(0);
  for(unsigned int eI = 0; eI < etaMin_p->size(); ++eI){
    etRecalc_p->push_back(0.0);
    etATLAS_p->push_back(0.0);
    fullEtaBins.push_back(etaMin_p->at(eI));
  }
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  fullEtaBins.push_back(etaMax_p->at(etaMax_p->size()-1));

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "ETABINS: ";
  for(unsigned int eI = 0; eI < fullEtaBins.size(); ++eI){
    std::cout << fullEtaBins[eI] << ", ";
  }
  std::cout << std::endl;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const ULong64_t nEntries = inTree_p->GetEntries();
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    inTree_p->GetEntry(entry);

    std::string runLumiEvtStr = std::to_string(run_) + "_" + std::to_string(lumi_) + "_" + std::to_string(evt_);
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(runLumiEvtToRhoEntry.count(runLumiEvtStr) == 0){
      std::cout << "CANNOT FIND ENTRY (Run_Lumi_Evt) " << runLumiEvtStr << ", continue" << std::endl;
      continue;
    }
    ULong64_t rhoEntry = runLumiEvtToRhoEntry[runLumiEvtStr];

    rhoTree_p->GetEntry(rhoEntry);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(unsigned int eI = 0; eI < etaMin_p->size(); ++eI){
      etATLAS_p->at(eI) = et_p->at(eI);
      etRecalc_p->at(eI) = 0.0;
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(unsigned int tI = 0; tI < towers_pt_p->size(); ++tI){
      int etaPos = ghostPos(fullEtaBins, towers_eta_p->at(tI));
      etRecalc_p->at(etaPos) += towers_pt_p->at(tI);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    outTree_p->Fill();
  }
  
  inFile_p->Close();
  delete inFile_p;
  
  rhoFile_p->Close();
  delete rhoFile_p;

  outFile_p->cd();

  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;

  //store some config info in file
  std::string etaMinStr = "";
  std::string etaMaxStr = "";
  for(unsigned int eI = 0; eI < etaMin_p->size(); ++eI){
    etaMinStr = etaMinStr + prettyString(etaMin_p->at(eI), 1, false) + ",";
    etaMaxStr = etaMaxStr + prettyString(etaMax_p->at(eI), 1, false) + ",";
  }

  TEnv configEnv;
  configEnv.SetValue("ETAMINBINS", etaMinStr.c_str());
  configEnv.SetValue("ETAMAXBINS", etaMaxStr.c_str());
  configEnv.Write("config", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
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
