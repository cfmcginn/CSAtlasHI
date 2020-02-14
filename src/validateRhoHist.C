//Author: Chris McGinn (2020.02.06)

//cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/getLinBins.h"
#include "include/histDefUtility.h"
#include "include/stringUtil.h"

int validateRhoHist(std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("validateRhoTree");  

  std::vector<float>* etATLAS_p=nullptr;
  std::vector<float>* etRecalc_p=nullptr;

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("etATLAS", 1);
  inTree_p->SetBranchStatus("etRecalc", 1);

  inTree_p->SetBranchAddress("etATLAS", &etATLAS_p);
  inTree_p->SetBranchAddress("etRecalc", &etRecalc_p);

  TEnv* config_p = (TEnv*)inFile_p->Get("config");
  configParser configPar(config_p);
  std::map<std::string, std::string> configMap = configPar.GetConfigMap();
  std::vector<std::string> essentialParams = {"ETAMINBINS", "ETAMAXBINS"};
  for(auto const & param : essentialParams){
    if(configMap.count(param) == 0){
      std::cout << "VALIDATERHOHIST ERROR: Essential parameter \'" << param << "\' not found in config taken from TEnv in file \'" << inFileName << "\'. Please fix. return 1" << std::endl;
      return 1;
    }
  }

  std::vector<float> etaMinBins = strToVectF(configMap["ETAMINBINS"]);
  std::vector<float> etaMaxBins = strToVectF(configMap["ETAMAXBINS"]);

  const Int_t nMaxEtaBins = 500;
  const Int_t nEtaBins = etaMinBins.size();
  if(nEtaBins > nMaxEtaBins){
    std::cout << "VALIDATERHOHIST ERROR: Necessary number of etabins from config, \'" << nEtaBins << "\' is greater than max, \'" << nMaxEtaBins << "\'. Please fix. return 1" << std::endl;
    return 1;
  }
  Double_t etaBins[nMaxEtaBins+1];
  for(unsigned int eI = 0; eI < etaMinBins.size(); ++eI){
    etaBins[eI] = etaMinBins[eI];
  }
  etaBins[nEtaBins] = etaMaxBins[etaMaxBins.size()-1];

  const Int_t nDeltaRhoBins = 100;
  const Float_t rhoLow = -20;
  const Float_t rhoHigh = 20;
  Double_t deltaRhoBins[nDeltaRhoBins+1];
  getLinBins(rhoLow, rhoHigh, nDeltaRhoBins, deltaRhoBins);  

  std::string outFileName = "output/" + dateStr + "/validateRhoHist_" + dateStr + ".root";
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* deltaEt_p = new TH1D("deltaEt_h", ";#Sigma(E_{T})_{Recalc.} - #Sigma(E_{T})_{ATLAS} [GeV];Counts", nDeltaRhoBins, deltaRhoBins);
  TH2D* deltaEtVEta_p = new TH2D("deltaEtVEta_h", ";#eta;#Sigma(E_{T})_{Recalc.} - #Sigma(E_{T})_{ATLAS} [GeV]", nEtaBins, etaBins, nDeltaRhoBins, deltaRhoBins);
  centerTitles(deltaEt_p);
  centerTitles(deltaEtVEta_p);
  
  const ULong64_t nEntries = inTree_p->GetEntries();
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    inTree_p->GetEntry(entry);

    if(etRecalc_p->size() != etaMinBins.size()){
      std::cout << "VALIDATERHOHIST ERROR: etaMinBins size \'" << etaMinBins.size() << "\' differs from etRecalc_p size, \'" << etRecalc_p->size() << "\'. Please fix. return 1" << std::endl;
      return 1;
    }

    for(unsigned int eI = 0; eI < etRecalc_p->size(); ++eI){      
      double cent = (etaMinBins.at(eI) + etaMaxBins.at(eI))/2.;
      double val = etRecalc_p->at(eI) - etATLAS_p->at(eI);      
      if(val <= rhoLow) val = (deltaRhoBins[0] + deltaRhoBins[1])/2.;
      if(val >= rhoHigh) val = (deltaRhoBins[nDeltaRhoBins-1] + deltaRhoBins[nDeltaRhoBins])/2.;

      deltaEt_p->Fill(val);
      deltaEtVEta_p->Fill(cent, val);
    }
  }  

  outFile_p->cd();

  deltaEt_p->Write("", TObject::kOverwrite);
  delete deltaEt_p;

  deltaEtVEta_p->Write("", TObject::kOverwrite);
  delete deltaEtVEta_p;

  outFile_p->Close();
  delete outFile_p;

  
  inFile_p->Close();
  delete inFile_p;
  
  std::cout << "VALIDATE RHO HIST COMPLETE. return 0" << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/validateRhoHist.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += validateRhoHist(argv[1]);
  return retVal;
}
