//Author: Chris McGinn (2020.07.05)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/ghostUtil.h"

int analyzeTowers(std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const double delta = 0.0001;
  std::vector<double> etaVals, phiVals;
  std::map<unsigned int, ULong64_t> etaValCounts, phiValCounts;

  std::vector<float> etaBounds = {-5.0};
  std::vector<std::vector<int> > countsPerEta;
  while(etaBounds[etaBounds.size()-1] < 5.0){
    etaBounds.push_back(etaBounds[etaBounds.size()-1] + 0.1);
    countsPerEta.push_back({});
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* gammaJetTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

  std::vector<float>* tower_pt_p=nullptr;
  std::vector<float>* tower_eta_p=nullptr;
  std::vector<float>* tower_phi_p=nullptr;

  gammaJetTree_p->SetBranchStatus("*", 0);
  gammaJetTree_p->SetBranchStatus("tower_pt", 1);
  gammaJetTree_p->SetBranchStatus("tower_eta", 1);
  gammaJetTree_p->SetBranchStatus("tower_phi", 1);

  gammaJetTree_p->SetBranchAddress("tower_pt", &tower_pt_p);
  gammaJetTree_p->SetBranchAddress("tower_eta", &tower_eta_p);
  gammaJetTree_p->SetBranchAddress("tower_phi", &tower_phi_p);

  const ULong64_t nEntries = TMath::Min((ULong64_t)100, (ULong64_t)gammaJetTree_p->GetEntries());
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    gammaJetTree_p->GetEntry(entry);

    for(unsigned int eI = 0; eI < countsPerEta.size(); ++eI){
      countsPerEta[eI].push_back(0);
    }
    
    for(unsigned int tI = 0; tI < tower_eta_p->size(); ++tI){
      bool isFound = false;
      for(unsigned int eI = 0; eI < etaVals.size(); ++eI){
	if(TMath::Abs(etaVals[eI] - tower_eta_p->at(tI)) < delta){
	  isFound = true;
	  ++(etaValCounts[eI]);
	  break;
	}
      }
     
      if(!isFound){
	etaVals.push_back(tower_eta_p->at(tI));
	etaValCounts[etaVals.size()-1] = 1;
      }

      isFound = false;
      for(unsigned int eI = 0; eI < phiVals.size(); ++eI){
	if(TMath::Abs(phiVals[eI] - tower_phi_p->at(tI)) < delta){
	  isFound = true;
	  ++(phiValCounts[eI]);
	  break;
	}
      }
     
      if(!isFound){
	phiVals.push_back(tower_phi_p->at(tI));
	phiValCounts[phiVals.size()-1] = 1;
      }

      Int_t etaPos = ghostPos(etaBounds, tower_eta_p->at(tI));
      ++(countsPerEta[etaPos].at(countsPerEta[etaPos].size()-1));
    }
  }

  inFile_p->Close();
  delete inFile_p;
  
  std::vector<double> reducedEtaVals;
  for(unsigned int eI = 0; eI < etaVals.size(); ++eI){
    if(etaValCounts[eI] > 200) reducedEtaVals.push_back(etaVals[eI]);
  }

  std::vector<double> reducedPhiVals;
  for(unsigned int eI = 0; eI < phiVals.size(); ++eI){
    if(phiValCounts[eI] > 200) reducedPhiVals.push_back(phiVals[eI]);
  }

  std::vector<double> phiTestVals;
  const Int_t nPhiBins = 64;
  for(unsigned int pI = 0; pI < nPhiBins+1; ++pI){
    phiTestVals.push_back(-TMath::Pi() + pI*2.*TMath::Pi()/(Double_t)nPhiBins);
  }


  std::cout << reducedEtaVals.size() << " UNIQUE ETA VALUES: ";
  std::sort(reducedEtaVals.begin(), reducedEtaVals.end());
  for(unsigned int eI = 0; eI < reducedEtaVals.size(); ++eI){
    std::cout << reducedEtaVals[eI] << ", ";
  }
  std::cout << std::endl;
  std::cout << reducedPhiVals.size() << " UNIQUE PHI VALUES: ";
  std::sort(reducedPhiVals.begin(), reducedPhiVals.end());
  for(unsigned int eI = 0; eI < reducedPhiVals.size(); ++eI){
    std::cout << reducedPhiVals[eI] << ", ";
  }
  std::cout << std::endl;
  std::cout << phiTestVals.size() << " UNIQUE PHI VALUES2: ";
  for(unsigned int eI = 0; eI < phiTestVals.size()-1; ++eI){
    double midVal = (phiTestVals[eI] + phiTestVals[eI+1])/2.;
    std::cout << midVal << ", ";
  }
  std::cout << std::endl;

  std::cout << "N TOWERS PER ETA BOUNDARIES: " << std::endl;
  for(unsigned int eI = 0; eI < etaBounds.size()-1; ++eI){
    std::vector<int> tempCounts = countsPerEta[eI];
    std::sort(std::begin(tempCounts), std::end(tempCounts));

    std::cout << " " << etaBounds[eI] << "-" << etaBounds[eI+1] << ": " << tempCounts[0] << "-" << tempCounts[tempCounts.size()-1] << std::endl;
  }

  std::cout << "ANALYZETOWERS COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/analyzeTowers.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += analyzeTowers(argv[1]);
  return retVal;
}
