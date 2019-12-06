//cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TDatime.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TNamed.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/etaPhiFunc.h"
#include "include/getLogBins.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/returnRootFileContentsList.h"

int makeClusterHist(std::string inFileName, std::string jzWeightsName = "")
{
  if(!checkFileExt(inFileName, "root")) return 1;
  const bool doJZWeights = checkFileExt(jzWeightsName, "txt");
  std::map<std::string, double> weightMap;
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

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> tnamedList = returnRootFileContentsList(inFile_p, "TNamed");
  std::map<std::string, std::string> tnamedMap;
  for(auto const& iter : tnamedList){
    TNamed* tempName_p = (TNamed*)inFile_p->Get(iter.c_str());
    tnamedMap[tempName_p->GetName()] = tempName_p->GetTitle();
  }
  inFile_p->Close();
  delete inFile_p;

  for(auto const & iter : tnamedMap){
    std::cout << iter.first << ", " << iter.second << std::endl;
  }

  if(tnamedMap.count("nEvents") == 0){
    std::cout << "Parameter \'nEvents\' not found in input paramMap. return 1" << std::endl;
    return 1;
  }
  else if(tnamedMap.count("nJtAlgo") == 0){
    std::cout << "Parameter \'nJtAlgo\' not found in input paramMap. return 1" << std::endl;
    return 1;
  }
  else if(tnamedMap.count("jtAlgos") == 0){
    std::cout << "Parameter \'jtAlgos\' not found in input paramMap. return 1" << std::endl;
    return 1;
  }

  const std::string caloTrkStr = tnamedMap["caloTrackStr"];
  
  const bool doATLAS = std::stoi(tnamedMap["doATLAS"]);
  const bool doTruth = std::stoi(tnamedMap["doTruth"]);
  
  const Int_t nMaxJtAlgo = 8;
  const Int_t nJtAlgo = std::stoi(tnamedMap["nJtAlgo"]);
  std::vector<std::string> jtAlgos = commaSepStringToVect(tnamedMap["jtAlgos"]);

  if(nJtAlgo != (Int_t)jtAlgos.size()){
    std::cout << "Mismatch between nJtAlgo \'" << nJtAlgo << "\' and jtAlgos.size() \'" << jtAlgos.size() << "\'. return 1" << std::endl;
    return 1;
  }
  
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  
  std::string outFileName = inFileName.substr(0, inFileName.find(".root"));
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName = "output/" + dateStr + "/" + outFileName + "_HIST_" + dateStr + ".root";
  
  const Int_t nCentBins = 5;
  const Int_t centBinsLow[nCentBins] = {0, 10, 20, 30, 40};
  const Int_t centBinsHigh[nCentBins] = {10, 20, 30, 40, 60};
  Int_t nEventPerCent[nCentBins];
  std::vector<std::string> centBinsStr;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    centBinsStr.push_back("Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHigh[cI]));
    nEventPerCent[cI] = 0;
  }

  const Float_t maxJtAbsEta = 2.4;//std::stod(tnamedMap["maxJtAbsEta"]);
  
  const Int_t nJtPtBins = 5;
  Float_t jtPtLow = 60;
  Float_t jtPtHigh = 200;
  if(isStrSame(caloTrkStr, "trk")){
    jtPtLow *= 0.7;
    jtPtHigh *= 0.7;
  }
  Double_t jtPtBins[nJtPtBins+1];
  getLogBins(jtPtLow, jtPtHigh, nJtPtBins, jtPtBins);
  std::vector<std::string> jtPtBinsStrVect;
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* spectra_p[nMaxJtAlgo+2][nCentBins];
  TH1D* matchedATLASSpectra_p[nMaxJtAlgo][nCentBins];
  TH1D* matchedTruthSpectra_p[nMaxJtAlgo][nCentBins];
  TH1D* recoOverGen_VPt_p[nMaxJtAlgo][nCentBins][nJtPtBins];
  TH1D* recoGen_DeltaEta_p[nMaxJtAlgo][nCentBins][nJtPtBins];
  TH1D* recoGen_DeltaPhi_p[nMaxJtAlgo][nCentBins][nJtPtBins];
  
  for(Int_t jI = 0; jI < nJtAlgo+2; ++jI){
    std::string algo = "Truth";
    if(jI < nJtAlgo) algo = jtAlgos[jI];
    else if(jI == nJtAlgo) algo = "ATLAS";

    if(!doATLAS && jI > nJtAlgo) break;
    if(!doTruth && jI > nJtAlgo+1) break;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string nameStr = algo + "_" + centBinsStr[cI];
      
      spectra_p[jI][cI] = new TH1D(("spectra_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];Counts", nJtPtBins, jtPtBins);

      if(doTruth){
	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	  std::string ptStr = "JtPt" + prettyString(jtPtBins[jI2], 1, true) + "to" + prettyString(jtPtBins[jI2+1], 1, true);
	  jtPtBinsStrVect.push_back(ptStr);
	  recoOverGen_VPt_p[jI][cI][jI2] = new TH1D(("recoOverGen_VPt_" + nameStr + "_" + ptStr + "_h").c_str(), ";Reco./Gen.;Counts", 21, 0.0, 2.0);
	  recoGen_DeltaEta_p[jI][cI][jI2] = new TH1D(("recoGen_DeltaEta_" + nameStr + "_" + ptStr + "_h").c_str(), ";#eta_{Reco.} - #eta_{Gen.};Counts", 21, -0.3, 0.3);
	  recoGen_DeltaPhi_p[jI][cI][jI2] = new TH1D(("recoGen_DeltaPhi_" + nameStr + "_" + ptStr + "_h").c_str(), ";#phi_{Reco.} - #phi_{Gen.};Counts", 21, -0.3, 0.3);

	  centerTitles({recoOverGen_VPt_p[jI][cI][jI2], recoGen_DeltaEta_p[jI][cI][jI2], recoGen_DeltaPhi_p[jI][cI][jI2]});
	}
      }
      
      if(jI < nJtAlgo && doATLAS){
	matchedATLASSpectra_p[jI][cI] = new TH1D(("matchedATLASSpectra_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];Counts w/ ATLAS jet match", nJtPtBins, jtPtBins);
	if(doTruth){
	  matchedTruthSpectra_p[jI][cI] = new TH1D(("matchedTruthSpectra_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];Counts w/ Truth jet match", nJtPtBins, jtPtBins);
	  centerTitles(matchedTruthSpectra_p[jI][cI]);
	}
	centerTitles(matchedATLASSpectra_p[jI][cI]);
      }
     
      centerTitles(spectra_p[jI][cI]);
    }
  }

  Int_t jtVal_;
  Float_t cent_;

  const Int_t nMaxJets = 500;
  Int_t njt_[nJtAlgo];
  Float_t jtpt_[nJtAlgo][nMaxJets];
  Float_t jteta_[nJtAlgo][nMaxJets];
  Float_t jtphi_[nJtAlgo][nMaxJets];
  Int_t atlasmatchpos_[nJtAlgo][nMaxJets];
  Int_t truthmatchpos_[nJtAlgo][nMaxJets];

  Int_t njtATLAS_;
  Float_t jtptATLAS_[nMaxJets];
  Float_t jtetaATLAS_[nMaxJets];
  Float_t jtphiATLAS_[nMaxJets];

  Int_t njtTruth_;
  Float_t jtptTruth_[nMaxJets];
  Float_t jtchgptTruth_[nMaxJets];
  Float_t jtetaTruth_[nMaxJets];
  Float_t jtphiTruth_[nMaxJets];

  inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* csTree_p = (TTree*)inFile_p->Get("clusterJetsCS");

  csTree_p->SetBranchStatus("*", 0);
  csTree_p->SetBranchStatus("cent", 1);
  if(doJZWeights) csTree_p->SetBranchStatus("jtVal", 1);

  csTree_p->SetBranchAddress("cent", &cent_);
  if(doJZWeights) csTree_p->SetBranchAddress("jtVal", &jtVal_);

  for(Int_t jI = 0; jI < nJtAlgo; ++jI){
    csTree_p->SetBranchStatus(("njt" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("jtpt" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("jteta" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("jtphi" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("atlasmatchpos" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("truthmatchpos" + jtAlgos[jI]).c_str(), 1);

    csTree_p->SetBranchAddress(("njt" + jtAlgos[jI]).c_str(), &njt_[jI]);
    csTree_p->SetBranchAddress(("jtpt" + jtAlgos[jI]).c_str(), jtpt_[jI]);
    csTree_p->SetBranchAddress(("jteta" + jtAlgos[jI]).c_str(), jteta_[jI]);
    csTree_p->SetBranchAddress(("jtphi" + jtAlgos[jI]).c_str(), jtphi_[jI]);
    csTree_p->SetBranchAddress(("atlasmatchpos" + jtAlgos[jI]).c_str(), atlasmatchpos_[jI]);
    csTree_p->SetBranchAddress(("truthmatchpos" + jtAlgos[jI]).c_str(), truthmatchpos_[jI]);
  }

  csTree_p->SetBranchStatus("njtATLAS", 1);
  csTree_p->SetBranchStatus("jtptATLAS", 1);
  csTree_p->SetBranchStatus("jtetaATLAS", 1);
  csTree_p->SetBranchStatus("jtphiATLAS", 1);
  csTree_p->SetBranchStatus("njtTruth", 1);
  csTree_p->SetBranchStatus("jtptTruth", 1);
  csTree_p->SetBranchStatus("jtchgptTruth", 1);
  csTree_p->SetBranchStatus("jtetaTruth", 1);
  csTree_p->SetBranchStatus("jtphiTruth", 1);
  
  csTree_p->SetBranchAddress("njtATLAS", &njtATLAS_);
  csTree_p->SetBranchAddress("jtptATLAS", jtptATLAS_);
  csTree_p->SetBranchAddress("jtetaATLAS", jtetaATLAS_);
  csTree_p->SetBranchAddress("jtphiATLAS", jtphiATLAS_);
  csTree_p->SetBranchAddress("njtTruth", &njtTruth_);
  csTree_p->SetBranchAddress("jtptTruth", jtptTruth_);
  csTree_p->SetBranchAddress("jtchgptTruth", jtchgptTruth_);
  csTree_p->SetBranchAddress("jtetaTruth", jtetaTruth_);
  csTree_p->SetBranchAddress("jtphiTruth", jtphiTruth_);

  const Int_t nEntries = csTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/50);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    csTree_p->GetEntry(entry);

    Int_t centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(cent_ >= centBinsLow[cI] && cent_ < centBinsHigh[cI]){
	centPos = cI;
	break;
      }
    }
    if(centPos < 0) continue;

    Double_t weight = 1.0;
    if(doJZWeights){
      std::string jtValStr = "JZ" + std::to_string(jtVal_);
      if(weightMap.count(jtValStr) == 0) std::cout << "WARNING CANNOT FIND WEIGHT FOR JTVAL \'" << jtValStr << "\'" << std::endl;
      weight = weightMap[jtValStr];
    }

    ++(nEventPerCent[centPos]);
    
    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      for(Int_t jI = 0; jI < njt_[aI]; ++jI){
	if(TMath::Abs(jteta_[aI][jI]) > maxJtAbsEta) continue;
	if(jtpt_[aI][jI] < jtPtLow) continue;
	if(jtpt_[aI][jI] >= jtPtHigh) continue;
	
	spectra_p[aI][centPos]->Fill(jtpt_[aI][jI], weight);
      }    
    }

    if(doATLAS){
      for(Int_t jI = 0; jI < njtATLAS_; ++jI){
	if(TMath::Abs(jtetaATLAS_[jI]) > maxJtAbsEta) continue;
	if(jtptATLAS_[jI] < jtPtLow) continue;
	if(jtptATLAS_[jI] >= jtPtHigh) continue;
	
	spectra_p[nJtAlgo][centPos]->Fill(jtptATLAS_[jI], weight);
	
	for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	  //	  bool isFilled = false;
	  
	  for(Int_t jI2 = 0; jI2 < njt_[aI]; ++jI2){
	    if(atlasmatchpos_[aI][jI2] == jI){
	      matchedATLASSpectra_p[aI][centPos]->Fill(jtpt_[aI][jI2], weight);
	      //      isFilled = true;
	      break;
	    }
	  }

	  /*	  
	  if(!isFilled && jtptATLAS_[jI] > 80.){
	    std::cout << "NO MATCH IN ALGO " << jtAlgos[aI] << " (entry, ptATLAS): " << entry << ", " << jtptATLAS_[jI] << std::endl;
	  }
	  */
	}
      }

      if(doTruth){
	if(isStrSame(caloTrkStr, "trk")){
	  for(Int_t jI = 0; jI < njtTruth_; ++jI){
	    jtptTruth_[jI] = jtchgptTruth_[jI];
	  }
	}

	for(Int_t jI = 0; jI < njtTruth_; ++jI){
	  if(TMath::Abs(jtetaTruth_[jI]) > maxJtAbsEta) continue;
	  if(jtptTruth_[jI] < jtPtLow) continue;
	  if(jtptTruth_[jI] >= jtPtHigh) continue;

	  Int_t jtPos = -1;
	  for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	    if(jtptTruth_[jI] >= jtPtBins[jI2] && jtptTruth_[jI] < jtPtBins[jI2+1]){
	      jtPos = jI2;
	      break;
	    }
	  }
	  if(jtPos == -1 && jtptTruth_[jI] == jtPtBins[nJtPtBins]) jtPos = nJtPtBins-1;
	  
	  spectra_p[nJtAlgo+1][centPos]->Fill(jtptTruth_[jI], weight);
	  	  
	  for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	    bool isFilled = false;
	    
	    for(Int_t jI2 = 0; jI2 < njt_[aI]; ++jI2){
	      if(truthmatchpos_[aI][jI2] == jI){
		matchedTruthSpectra_p[aI][centPos]->Fill(jtptTruth_[jI], weight);
		recoOverGen_VPt_p[aI][centPos][jtPos]->Fill(jtpt_[aI][jI2]/jtptTruth_[jI], weight);
		recoGen_DeltaEta_p[aI][centPos][jtPos]->Fill(jteta_[aI][jI2] - jtetaTruth_[jI], weight);
		recoGen_DeltaPhi_p[aI][centPos][jtPos]->Fill(getDPHI(jtphi_[aI][jI2], jtphiTruth_[jI]), weight);
		isFilled = true;
		break;
	      }
	    }
	    
	    if(!isFilled && jtptTruth_[jI] > 80.){
	      std::cout << "NO MATCH IN ALGO " << jtAlgos[aI] << " (entry, ptTruth): " << entry << ", " << jtptTruth_[jI] << std::endl;
	    }
	  }
	}
      }
    }
  }
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t aI = 0; aI < nJtAlgo+2; ++aI){
    if(!doATLAS && aI > nJtAlgo) break;
    if(!doTruth && aI > nJtAlgo+1) break;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      spectra_p[aI][cI]->Write("", TObject::kOverwrite);
      delete spectra_p[aI][cI];

      if(aI < nJtAlgo){
	if(doATLAS){
	  matchedATLASSpectra_p[aI][cI]->Write("", TObject::kOverwrite);
	  delete matchedATLASSpectra_p[aI][cI];

	  if(doTruth){
	    matchedTruthSpectra_p[aI][cI]->Write("", TObject::kOverwrite);
	    delete matchedTruthSpectra_p[aI][cI];

	    for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	      recoOverGen_VPt_p[aI][cI][jI]->Write("", TObject::kOverwrite);
	      delete recoOverGen_VPt_p[aI][cI][jI];

	      recoGen_DeltaEta_p[aI][cI][jI]->Write("", TObject::kOverwrite);
	      delete recoGen_DeltaEta_p[aI][cI][jI];

	      recoGen_DeltaPhi_p[aI][cI][jI]->Write("", TObject::kOverwrite);
	      delete recoGen_DeltaPhi_p[aI][cI][jI];
	    }
	  }
	}
      }
    }
  }

  TDirectoryFile* params_p = (TDirectoryFile*)outFile_p->mkdir("params");
  params_p->cd();

  std::string jtPtBinsStr = "";
  std::string jtPtBinsStr2 = "";
  for(Int_t jI = 0; jI < nJtPtBins+1; ++jI){
    jtPtBinsStr = jtPtBinsStr + prettyString(jtPtBins[jI], 1, false) + ",";
    if(jI != nJtPtBins) jtPtBinsStr2 = jtPtBinsStr2 + jtPtBinsStrVect[jI] + ",";
  }
  
  std::string centBinsLowStr = "";
  std::string centBinsHighStr = "";
  std::string centBinsStr2 = "";
  std::string nEventPerCentStr = "";
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    centBinsLowStr = centBinsLowStr + std::to_string(centBinsLow[cI]) + ",";
    centBinsHighStr = centBinsHighStr + std::to_string(centBinsHigh[cI]) + ",";
    centBinsStr2 = centBinsStr2 + centBinsStr[cI] + ",";
    nEventPerCentStr = nEventPerCentStr + std::to_string(nEventPerCent[cI]) + ",";
  }

  std::string jtAlgos2 = "";
  for(Int_t aI = 0; aI < nJtAlgo+1; ++aI){
    std::string algo = "ATLAS";
    if(aI < nJtAlgo) algo = jtAlgos[aI];

    jtAlgos2 = jtAlgos2 + algo + ",";
  }
  
  std::map<std::string, std::string> paramMap;
  paramMap["nJtPtBins"] = std::to_string(nJtPtBins);
  paramMap["jtPtBins"] = jtPtBinsStr;
  paramMap["jtPtBinsStr"] = jtPtBinsStr2;

  paramMap["nCentBins"] = std::to_string(nCentBins);
  paramMap["centBinsLow"] = centBinsLowStr;
  paramMap["centBinsHigh"] = centBinsHighStr;
  paramMap["centBinsStr"] = centBinsStr2;
  paramMap["nEventPerCent"] = nEventPerCentStr;
  paramMap["maxJtAbsEta"] = prettyString(maxJtAbsEta, 2, false);

  for(auto const & iter : tnamedMap){paramMap[iter.first] = iter.second;}
  paramMap["nJtAlgo"] = std::to_string(nJtAlgo);
  paramMap["jtAlgos"] = jtAlgos2;

  for(auto const & iter : paramMap){
    TNamed tempName(iter.first.c_str(), iter.second.c_str());
    tempName.Write("", TObject::kOverwrite);
  }
  
  outFile_p->Close();
  delete outFile_p;    
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 3){
    std::cout << "Usage: ./bin/makeClusterHist.exe <inFileName> <jzWeightsName-default=\'\'>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 2) retVal += makeClusterHist(argv[1]);
  else if(argc == 3) retVal += makeClusterHist(argv[1], argv[2]);
  return retVal;
}
