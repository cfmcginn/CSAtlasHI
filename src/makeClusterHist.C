//cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TDirectoryFile.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/etaPhiFunc.h"
#include "include/getLogBins.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/ncollFunctions_5TeV.h"
#include "include/plotUtilities.h"
#include "include/returnRootFileContentsList.h"
#include "include/sharedFunctions.h"
#include "include/stringUtil.h"

int makeClusterHist(std::string inConfigFileName)
{
  globalDebugHandler gDebug;
  const bool doDebug = gDebug.GetDoGlobalDebug();

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"INFILENAME",
					"DOJZWEIGHTS",
					"DOCENTWEIGHTS"};

  if(!checkConfigContainsParams(inConfig_p, reqParams)) return 1;

  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  const bool doJZWeights = inConfig_p->GetValue("DOJZWEIGHTS", 0);
  const bool doCentWeights = inConfig_p->GetValue("DOCENTWEIGHTS", 0);

  TEnv* outEnv_p = new TEnv();
  outEnv_p->SetValue("DOJZWEIGHTS", doJZWeights);
  outEnv_p->SetValue("DOCENTWEIGHTS", doCentWeights);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::map<int, double> jzWeightMap;
  std::map<int, double> centWeightMap;
  if(doCentWeights){
    double centNorm = findAvgNColl_Cent(0,1);
    for(int i = 0; i < 100; ++i){
      double weight = findAvgNColl_Cent(i, i+1);
      centWeightMap[i] = weight/centNorm;
    }
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* fileConfig_p = (TEnv*)inFile_p->Get("config");

  std::vector<std::string> reqParamsFile = {"NJTALGO",
					    "JTALGOS",
					    "ISMC",
					    "RECOJTMINPT"};

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nMaxJtAlgo = 20;
  const Int_t nJtAlgo = fileConfig_p->GetValue("NJTALGO", 0);
  std::string jtAlgosStr = fileConfig_p->GetValue("JTALGOS", "");
  std::vector<std::string> jtAlgos = commaSepStringToVect(jtAlgosStr);
  const bool isMC = fileConfig_p->GetValue("ISMC", 0);

  const double minJtPt = fileConfig_p->GetValue("RECOJTMINPT", 0.0);

  outEnv_p->SetValue("NJTALGO", nJtAlgo);
  outEnv_p->SetValue("JTALGOS", jtAlgosStr.c_str());
  outEnv_p->SetValue("ISMC", isMC);
 
  std::vector<double> uniqueXSec;
  std::vector<double> uniqueFilterEff;
  std::map<int, unsigned long long> xSecCounter;
  std::map<int, unsigned long long> centCounter;
  for(int i = 0; i < 100; ++i){
    centCounter[i] = 0;
  }

  Float_t cent_;
  Float_t xSectionNB_;
  Float_t filterEff_;
  
  TTree* csTree_p = (TTree*)inFile_p->Get("clusterJetsCS");
  const Int_t nEntries = csTree_p->GetEntries();
  if(isMC && (doJZWeights || doCentWeights)){
    csTree_p->SetBranchStatus("*", 0);

    if(doCentWeights) csTree_p->SetBranchStatus("cent", 1);
    if(doJZWeights){
      csTree_p->SetBranchStatus("xSectionNB", 1);
      csTree_p->SetBranchStatus("filterEff", 1);
    }

    if(doCentWeights) csTree_p->SetBranchAddress("cent", &cent_);
    if(doJZWeights){
      csTree_p->SetBranchAddress("xSectionNB", &xSectionNB_);
      csTree_p->SetBranchAddress("filterEff", &filterEff_);
    }

    for(Int_t entry = 0; entry < nEntries; ++entry){
      csTree_p->GetEntry(entry);
      
      if(doJZWeights){
	int xsecPos = -1;
	for(unsigned int xI = 0; xI < uniqueXSec.size(); ++xI){
	  if(TMath::Abs(uniqueXSec[xI] - xSectionNB_) < 0.1){
	    if(TMath::Abs(uniqueFilterEff[xI] - filterEff_) < 0.00001){
	      xsecPos = xI;
	      break;
	    }
	  }
	}
	
	if(xsecPos >= 0) ++(xSecCounter[xsecPos]);
	else{
	  xSecCounter[uniqueXSec.size()] = 1;
	  uniqueXSec.push_back(xSectionNB_);
	  uniqueFilterEff.push_back(filterEff_);
	}
      }

      if(doCentWeights){
	int centInt = (int)cent_;
	++(centCounter[centInt]);
      }
    }

    if(doJZWeights){
      std::cout << "JZ WEIGHTS: " << std::endl;
      for(unsigned int xI = 0; xI < uniqueXSec.size(); ++xI){
	jzWeightMap[xI] = uniqueXSec[xI]*uniqueFilterEff[xI]/xSecCounter[xI];
	
	std::cout << " " << xI << ": " << uniqueXSec[xI] << "*" << uniqueFilterEff[xI] << "/" << xSecCounter[xI] << "=" << jzWeightMap[xI] << std::endl; 
      }
    }
  }

  inFile_p->Close();
  delete inFile_p;  

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


  if(nJtAlgo != (Int_t)jtAlgos.size()){
    std::cout << "Mismatch between nJtAlgo \'" << nJtAlgo << "\' and jtAlgos.size() \'" << jtAlgos.size() << "\'. return 1" << std::endl;
    return 1;
  }

  const std::string dateStr = getDateStr();

  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::string outFileName = inFileName.substr(0, inFileName.find(".root"));
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName = "output/" + dateStr + "/" + outFileName + "_HIST_" + dateStr + ".root";

  //Re-worked binning - don't want to go above 90 because of lack of stats  
  const Int_t nMaxCentBins = 10;
  const Int_t nCentBins = 7;
  const Int_t centBinsLow[nMaxCentBins] = {0, 10, 20, 30, 40, 50, 70};
  const Int_t centBinsHigh[nMaxCentBins] = {10, 20, 30, 40, 50, 70, 90};
  Int_t nEventPerCent[nMaxCentBins];
  std::vector<std::string> centBinsStr;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    centBinsStr.push_back("Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHigh[cI]));
    nEventPerCent[cI] = 0;
  }

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Float_t maxJtAbsEta = 2.4;
  
  const Int_t nJtPtBins = 12;
  Float_t jtPtLow = 20;
  Float_t jtPtHigh = 400;
  Double_t jtPtBins[nJtPtBins+1];
  getLogBins(jtPtLow, jtPtHigh, nJtPtBins, jtPtBins);
  std::vector<std::string> jtPtBinsStrVect;

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* cent_p = new TH1D("cent_h", ";Centrality (%);Weighted Counts", 100, -0.5, 99.5);
  TH1D* cent_CentWeightOnly_p = new TH1D("cent_CentWeightOnly_h", ";Centrality (%);Weighted Counts", 100, -0.5, 99.5);
  TH1D* cent_Unweighted_p = new TH1D("cent_Unweighted_h", ";Centrality (%);Unweighted Counts", 100, -0.5, 99.5);
  centerTitles({cent_p, cent_CentWeightOnly_p, cent_Unweighted_p});

  TH1D* spectra_p[nMaxJtAlgo+1][nMaxCentBins];
  TH1D* spectraUnmatched_p[nMaxJtAlgo][nMaxCentBins];
  TH1D* spectraChg_p[nMaxCentBins];
  TH1D* matchedTruthSpectra_p[nMaxJtAlgo][nMaxCentBins];
  TH1D* spectra_Unweighted_p[nMaxCentBins];
  TH1D* recoOverGen_VPt_p[nMaxJtAlgo][nMaxCentBins][nJtPtBins];
  TH1D* recoOverGenM_VPt_p[nMaxJtAlgo][nMaxCentBins][nJtPtBins];
  TH1D* recoOverGenMOverPt_VPt_p[nMaxJtAlgo][nMaxCentBins][nJtPtBins];
  TH1D* recoGen_DeltaEta_p[nMaxJtAlgo][nMaxCentBins][nJtPtBins];
  TH1D* recoGen_DeltaPhi_p[nMaxJtAlgo][nMaxCentBins][nJtPtBins];
  
  for(Int_t jI = 0; jI < nJtAlgo+1; ++jI){
    std::string algo = "Truth";
    if(jI < nJtAlgo) algo = jtAlgos[jI];
    if(!isMC && jI == nJtAlgo) break;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string nameStr = algo + "_" + centBinsStr[cI];
      
      spectra_p[jI][cI] = new TH1D(("spectra_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];Counts", nJtPtBins, jtPtBins);
      if(jI == nJtAlgo){
	spectraChg_p[cI] = new TH1D(("spectraChg_" + nameStr + "_h").c_str(), ";Charge Jet p_{T} [GeV];Counts", nJtPtBins, jtPtBins);
	centerTitles(spectraChg_p[cI]);
      }
      else{
	spectraUnmatched_p[jI][cI] = new TH1D(("spectraUnmatched_" + nameStr + "_h").c_str(), ";Charge Jet p_{T} [GeV];Counts", nJtPtBins, jtPtBins);
	centerTitles(spectraUnmatched_p[jI][cI]);
      }

      if(isMC && jI == nJtAlgo){
	spectra_Unweighted_p[cI] = new TH1D(("spectra_Unweighted_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];Counts", nJtPtBins, jtPtBins);
	centerTitles(spectra_Unweighted_p[cI]);
      }


      if(isMC){
	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	  std::string ptStr = "JtPt" + prettyString(jtPtBins[jI2], 1, true) + "to" + prettyString(jtPtBins[jI2+1], 1, true);
	  jtPtBinsStrVect.push_back(ptStr);
	  recoOverGen_VPt_p[jI][cI][jI2] = new TH1D(("recoOverGen_VPt_" + nameStr + "_" + ptStr + "_h").c_str(), ";Reco./Gen.;Counts", 51, 0.0, 2.0);
	  recoOverGenM_VPt_p[jI][cI][jI2] = new TH1D(("recoOverGenM_VPt_" + nameStr + "_" + ptStr + "_h").c_str(), ";Reco. Mass/Gen. Mass;Counts", 21, 0.0, 2.0);
	  recoOverGenMOverPt_VPt_p[jI][cI][jI2] = new TH1D(("recoOverGenMOverPt_VPt_" + nameStr + "_" + ptStr + "_h").c_str(), ";(Reco. M/p_{T})/(Gen. M/p_{T});Counts", 21, 0.0, 2.0);
	  recoGen_DeltaEta_p[jI][cI][jI2] = new TH1D(("recoGen_DeltaEta_" + nameStr + "_" + ptStr + "_h").c_str(), ";#eta_{Reco.} - #eta_{Gen.};Counts", 21, -0.3, 0.3);
	  recoGen_DeltaPhi_p[jI][cI][jI2] = new TH1D(("recoGen_DeltaPhi_" + nameStr + "_" + ptStr + "_h").c_str(), ";#phi_{Reco.} - #phi_{Gen.};Counts", 21, -0.3, 0.3);

	  centerTitles({recoOverGen_VPt_p[jI][cI][jI2], recoOverGenM_VPt_p[jI][cI][jI2], recoOverGenMOverPt_VPt_p[jI][cI][jI2], recoGen_DeltaEta_p[jI][cI][jI2], recoGen_DeltaPhi_p[jI][cI][jI2]});
	}

	matchedTruthSpectra_p[jI][cI] = new TH1D(("matchedTruthSpectra_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];Counts w/ Truth jet match", nJtPtBins, jtPtBins);
	centerTitles({matchedTruthSpectra_p[jI][cI]});      
      }     
      centerTitles(spectra_p[jI][cI]);
    }
  }

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nMaxJets = 500;
  Int_t njt_[nMaxJtAlgo];
  Float_t jtpt_[nMaxJtAlgo][nMaxJets];
  Float_t jteta_[nMaxJtAlgo][nMaxJets];
  Float_t jtphi_[nMaxJtAlgo][nMaxJets];
  Float_t jtm_[nMaxJtAlgo][nMaxJets];
  Int_t atlasmatchpos_[nMaxJtAlgo][nMaxJets];
  Int_t truthmatchpos_[nMaxJtAlgo][nMaxJets];
  Int_t chgtruthmatchpos_[nMaxJtAlgo][nMaxJets];

  Int_t njtTruth_;
  Float_t jtptTruth_[nMaxJets];
  Float_t jtetaTruth_[nMaxJets];
  Float_t jtphiTruth_[nMaxJets];
  Float_t jtmTruth_[nMaxJets];
  Int_t jtmatchposTruth_[nMaxJtAlgo][nMaxJets];

  Int_t nchgjtTruth_;
  Float_t chgjtptTruth_[nMaxJets];
  Float_t chgjtetaTruth_[nMaxJets];
  Float_t chgjtphiTruth_[nMaxJets];
  Float_t chgjtmTruth_[nMaxJets];
  Int_t chgjtmatchposTruth_[nMaxJtAlgo][nMaxJets];

  inFile_p = new TFile(inFileName.c_str(), "READ");
  csTree_p = (TTree*)inFile_p->Get("clusterJetsCS");

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  csTree_p->SetBranchStatus("*", 0);
  csTree_p->SetBranchStatus("cent", 1);
  if(doJZWeights && isMC){
    csTree_p->SetBranchStatus("xSectionNB", 1);
    csTree_p->SetBranchStatus("filterEff", 1);
  }

  csTree_p->SetBranchAddress("cent", &cent_);
  if(doJZWeights && isMC){
    csTree_p->SetBranchAddress("xSectionNB", &xSectionNB_);
    csTree_p->SetBranchAddress("filterEff", &filterEff_);
  }

  for(Int_t jI = 0; jI < nJtAlgo; ++jI){
    csTree_p->SetBranchStatus(("njt" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("jtpt" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("jteta" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("jtphi" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("jtm" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("atlasmatchpos" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("truthmatchpos" + jtAlgos[jI]).c_str(), 1);
    csTree_p->SetBranchStatus(("chgtruthmatchpos" + jtAlgos[jI]).c_str(), 1);

    csTree_p->SetBranchAddress(("njt" + jtAlgos[jI]).c_str(), &njt_[jI]);
    csTree_p->SetBranchAddress(("jtpt" + jtAlgos[jI]).c_str(), jtpt_[jI]);
    csTree_p->SetBranchAddress(("jteta" + jtAlgos[jI]).c_str(), jteta_[jI]);
    csTree_p->SetBranchAddress(("jtphi" + jtAlgos[jI]).c_str(), jtphi_[jI]);
    csTree_p->SetBranchAddress(("jtm" + jtAlgos[jI]).c_str(), jtm_[jI]);
    csTree_p->SetBranchAddress(("atlasmatchpos" + jtAlgos[jI]).c_str(), atlasmatchpos_[jI]);
    csTree_p->SetBranchAddress(("truthmatchpos" + jtAlgos[jI]).c_str(), truthmatchpos_[jI]);
    csTree_p->SetBranchAddress(("chgtruthmatchpos" + jtAlgos[jI]).c_str(), chgtruthmatchpos_[jI]);
  }

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  csTree_p->SetBranchStatus("njtTruth", 1);
  csTree_p->SetBranchStatus("jtptTruth", 1);
  csTree_p->SetBranchStatus("jtetaTruth", 1);
  csTree_p->SetBranchStatus("jtphiTruth", 1);
  csTree_p->SetBranchStatus("jtmTruth", 1);
  for(Int_t aI = 0; aI < nJtAlgo; ++aI){
    csTree_p->SetBranchStatus(("jtmatchpos" + jtAlgos[aI] + "Truth").c_str(), 1);
  }

  csTree_p->SetBranchStatus("nchgjtTruth", 1);
  csTree_p->SetBranchStatus("chgjtptTruth", 1);
  csTree_p->SetBranchStatus("chgjtetaTruth", 1);
  csTree_p->SetBranchStatus("chgjtphiTruth", 1);
  csTree_p->SetBranchStatus("chgjtmTruth", 1);
  for(Int_t aI = 0; aI < nJtAlgo; ++aI){
    csTree_p->SetBranchStatus(("chgjtmatchpos" + jtAlgos[aI] + "Truth").c_str(), 1);
  }


  csTree_p->SetBranchAddress("njtTruth", &njtTruth_);
  csTree_p->SetBranchAddress("jtptTruth", jtptTruth_);
  csTree_p->SetBranchAddress("jtetaTruth", jtetaTruth_);
  csTree_p->SetBranchAddress("jtphiTruth", jtphiTruth_);
  csTree_p->SetBranchAddress("jtmTruth", jtmTruth_);
  for(Int_t aI = 0; aI < nJtAlgo; ++aI){
    csTree_p->SetBranchAddress(("jtmatchpos" + jtAlgos[aI] + "Truth").c_str(), jtmatchposTruth_[aI]);
  }

  csTree_p->SetBranchAddress("nchgjtTruth", &nchgjtTruth_);
  csTree_p->SetBranchAddress("chgjtptTruth", chgjtptTruth_);
  csTree_p->SetBranchAddress("chgjtetaTruth", chgjtetaTruth_);
  csTree_p->SetBranchAddress("chgjtphiTruth", chgjtphiTruth_);
  csTree_p->SetBranchAddress("chgjtmTruth", chgjtmTruth_);
  for(Int_t aI = 0; aI < nJtAlgo; ++aI){
    csTree_p->SetBranchAddress(("chgjtmatchpos" + jtAlgos[aI] + "Truth").c_str(), chgjtmatchposTruth_[aI]);
  }

  const Int_t nDiv = TMath::Max(1, nEntries/50);

  /*
  std::vector<int> entries;
  std::vector<double> weights;
  std::vector<double> jzWeights;
  std::vector<double> centWeights;
  std::vector<double> cents;
  */
  
  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    csTree_p->GetEntry(entry);

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(cent_ > 98.1) continue;

    Int_t centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(cent_ >= centBinsLow[cI] && cent_ < centBinsHigh[cI]){
	centPos = cI;
	break;
      }
    }
    if(centPos < 0) continue;

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    Double_t weight = 1.0;
    Double_t centWeight = 1.0;
    if(doJZWeights && isMC){
      int xsecPos = -1;
      for(unsigned int xI = 0; xI < uniqueXSec.size(); ++xI){
	if(TMath::Abs(uniqueXSec[xI] - xSectionNB_) < 0.1){
	  if(TMath::Abs(uniqueFilterEff[xI] - filterEff_) < 0.00001){
	    xsecPos = xI;
	    break;
	  }
	}
      }

      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      if(xsecPos < 0){
	std::cout << "Couldnt find weight for x-sec/filter eff: " << xSectionNB_ << "/" << filterEff_ << std::endl;
	return 1;
      }
      weight *= jzWeightMap[xsecPos];
    }

    if(doCentWeights){
      int centInt_ = (int)cent_;
      centWeight = centWeightMap[centInt_]/centCounter[centInt_];
      weight *= centWeight;
    }

    if(entry == 804801) std::cout << "JZ WEIGHT, CENT WEIGHT: " << weight/centWeight << ", " << centWeight << std::endl;

    cent_p->Fill(cent_, weight);
    cent_CentWeightOnly_p->Fill(cent_, centWeight);
    cent_Unweighted_p->Fill(cent_);

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ++(nEventPerCent[centPos]);
    
    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      for(Int_t jI = 0; jI < njt_[aI]; ++jI){
	if(TMath::Abs(jteta_[aI][jI]) > maxJtAbsEta) continue;
	if(jtpt_[aI][jI] < jtPtLow) continue;
	if(jtpt_[aI][jI] >= jtPtHigh) continue;
	
	spectra_p[aI][centPos]->Fill(jtpt_[aI][jI], weight);

	if(isMC){
	  if(jtAlgos[aI].find("Trk") != std::string::npos){
	    if(chgtruthmatchpos_[aI][jI] < 0) spectraUnmatched_p[aI][centPos]->Fill(jtpt_[aI][jI], weight);	  
	  }
	  else{
	    if(truthmatchpos_[aI][jI] < 0) spectraUnmatched_p[aI][centPos]->Fill(jtpt_[aI][jI], weight);	  
	  }
	}
      }
    }
    
    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    if(isMC){
      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      for(Int_t jI = 0; jI < njtTruth_; ++jI){
	if(TMath::Abs(jtetaTruth_[jI]) > maxJtAbsEta) continue;
	if(jtptTruth_[jI] < jtPtLow) continue;
	if(jtptTruth_[jI] >= jtPtHigh) continue;
	
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	Int_t jtPos = -1;
	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	  if(jtptTruth_[jI] >= jtPtBins[jI2] && jtptTruth_[jI] < jtPtBins[jI2+1]){
	    jtPos = jI2;
	    break;
	  }
	}
	if(jtPos == -1 && jtptTruth_[jI] == jtPtBins[nJtPtBins]) jtPos = nJtPtBins-1;
	  
	
	spectra_p[nJtAlgo][centPos]->Fill(jtptTruth_[jI], weight);
	spectra_Unweighted_p[centPos]->Fill(jtptTruth_[jI]);


	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      	
	for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	  if(jtAlgos[aI].find("Trk") != std::string::npos) continue;
	  int pos = jtmatchposTruth_[aI][jI];
	  if(pos >= 0){	    
	    matchedTruthSpectra_p[aI][centPos]->Fill(jtptTruth_[jI], weight);
	    
	    recoOverGen_VPt_p[aI][centPos][jtPos]->Fill(jtpt_[aI][pos]/jtptTruth_[jI], weight);
	    recoOverGenM_VPt_p[aI][centPos][jtPos]->Fill(jtm_[aI][pos]/jtmTruth_[jI], weight);
	    recoOverGenMOverPt_VPt_p[aI][centPos][jtPos]->Fill(jtm_[aI][pos]*jtptTruth_[jI]/(jtmTruth_[jI]*jtpt_[aI][pos]), weight);

	    recoGen_DeltaEta_p[aI][centPos][jtPos]->Fill(jteta_[aI][pos] - jtetaTruth_[jI], weight);
	    recoGen_DeltaPhi_p[aI][centPos][jtPos]->Fill(getDPHI(jtphi_[aI][pos], jtphiTruth_[jI]), weight);

	    /*
	    if(cent_ >= 80 && jtAlgos[aI].find("TowerCSGlobalAlpha1IterRho0") != std::string::npos){
	      if(jtptTruth_[jI] >= 33.0 && jtptTruth_[jI] < 42.3){
		if(jtpt_[aI][pos]/jtptTruth_[jI] > 0.35 && jtpt_[aI][pos]/jtptTruth_[jI] < 0.46){	
		  //		if(jtpt_[aI][pos]/jtptTruth_[jI] > 0.46 && jtpt_[aI][pos]/jtptTruth_[jI] < 0.57){
	  
		  weights.push_back(weight);
		  jzWeights.push_back(weight/centWeight);
		  centWeights.push_back(centWeight);
		  cents.push_back(cent_);
		  entries.push_back(entry);
		}
	      }
	    }
	    */
	  }
	}
      }
    
      for(Int_t jI = 0; jI < nchgjtTruth_; ++jI){
	if(TMath::Abs(chgjtetaTruth_[jI]) > maxJtAbsEta) continue;
	if(chgjtptTruth_[jI] < jtPtLow) continue;
	if(chgjtptTruth_[jI] >= jtPtHigh) continue;
	
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	Int_t jtPos = -1;
	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	  if(chgjtptTruth_[jI] >= jtPtBins[jI2] && chgjtptTruth_[jI] < jtPtBins[jI2+1]){
	    jtPos = jI2;
	    break;
	  }
	}
	if(jtPos == -1 && chgjtptTruth_[jI] == jtPtBins[nJtPtBins]) jtPos = nJtPtBins-1;


	spectraChg_p[centPos]->Fill(chgjtptTruth_[jI], weight);
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  if(jtAlgos[aI].find("Trk") == std::string::npos) continue;
	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  int pos = chgjtmatchposTruth_[aI][jI];
	  if(pos >= 0){	    
	    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    matchedTruthSpectra_p[aI][centPos]->Fill(chgjtptTruth_[jI], weight);
	    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    recoOverGen_VPt_p[aI][centPos][jtPos]->Fill(jtpt_[aI][pos]/chgjtptTruth_[jI], weight);
	    recoOverGenM_VPt_p[aI][centPos][jtPos]->Fill(jtm_[aI][pos]/chgjtmTruth_[jI], weight);
	    recoOverGenMOverPt_VPt_p[aI][centPos][jtPos]->Fill(jtm_[aI][pos]*chgjtptTruth_[jI]/(chgjtmTruth_[jI]*jtpt_[aI][pos]), weight);

	    recoGen_DeltaEta_p[aI][centPos][jtPos]->Fill(jteta_[aI][pos] - chgjtetaTruth_[jI], weight);
	    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    recoGen_DeltaPhi_p[aI][centPos][jtPos]->Fill(getDPHI(jtphi_[aI][pos], chgjtphiTruth_[jI]), weight);
	    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  }
	}
      
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      }
    }
  }

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  /*
  std::cout << "FINAL FILLS: " << std::endl;
  for(unsigned int eI = 0; eI < entries.size(); ++eI){
    std::cout << " " << eI << ": " << entries[eI] << ", " << weights[eI] << ", " << jzWeights[eI] << ", " << centWeights[eI] << ", " << cents[eI] << std::endl;
  }
  */

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  cent_p->Write("", TObject::kOverwrite);
  delete cent_p;

  cent_CentWeightOnly_p->Write("", TObject::kOverwrite);
  delete cent_CentWeightOnly_p;

  cent_Unweighted_p->Write("", TObject::kOverwrite);
  delete cent_Unweighted_p;

  for(Int_t aI = 0; aI < nJtAlgo+1; ++aI){
    if(!isMC && aI == nJtAlgo) break;

    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      spectra_p[aI][cI]->Write("", TObject::kOverwrite);
      delete spectra_p[aI][cI];

      
      if(aI == nJtAlgo){
	spectraChg_p[cI]->Write("", TObject::kOverwrite);
	delete spectraChg_p[cI];
      }
      else{
	spectraUnmatched_p[aI][cI]->Write("", TObject::kOverwrite);
	delete spectraUnmatched_p[aI][cI];
      }

      if(aI == nJtAlgo){
	spectra_Unweighted_p[cI]->Write("", TObject::kOverwrite);
	delete spectra_Unweighted_p[cI];
      }

      if(aI < nJtAlgo && isMC){
	matchedTruthSpectra_p[aI][cI]->Write("", TObject::kOverwrite);
	delete matchedTruthSpectra_p[aI][cI];


	for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	  recoOverGen_VPt_p[aI][cI][jI]->Write("", TObject::kOverwrite);
	  delete recoOverGen_VPt_p[aI][cI][jI];

	  recoOverGenM_VPt_p[aI][cI][jI]->Write("", TObject::kOverwrite);
	  delete recoOverGenM_VPt_p[aI][cI][jI];

	  recoOverGenMOverPt_VPt_p[aI][cI][jI]->Write("", TObject::kOverwrite);
	  delete recoOverGenMOverPt_VPt_p[aI][cI][jI];
	  
	  recoGen_DeltaEta_p[aI][cI][jI]->Write("", TObject::kOverwrite);
	  delete recoGen_DeltaEta_p[aI][cI][jI];
	  
	  recoGen_DeltaPhi_p[aI][cI][jI]->Write("", TObject::kOverwrite);
	  delete recoGen_DeltaPhi_p[aI][cI][jI];
	}
      }
    }
  }


  std::string jtPtBinsStr = "";
  std::string jtPtBinsStr2 = "";
  for(Int_t jI = 0; jI < nJtPtBins+1; ++jI){
    jtPtBinsStr = jtPtBinsStr + prettyString(jtPtBins[jI], 1, false) + ",";
    if(jI != nJtPtBins) jtPtBinsStr2 = jtPtBinsStr2 + jtPtBinsStrVect[jI] + ",";
  }

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
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
  
  outEnv_p->SetValue("NJTPTBINS", nJtPtBins);
  outEnv_p->SetValue("JTPTBINS", jtPtBinsStr.c_str());
  outEnv_p->SetValue("JTPTBINSSTR", jtPtBinsStr2.c_str());

  outEnv_p->SetValue("NCENTBINS", nCentBins);
  outEnv_p->SetValue("CENTBINSLOW", centBinsLowStr.c_str());
  outEnv_p->SetValue("CENTBINSHIGH", centBinsHighStr.c_str());
  outEnv_p->SetValue("CENTBINSSTR", centBinsStr2.c_str());
  outEnv_p->SetValue("NEVENTPERCENT", nEventPerCentStr.c_str());
  outEnv_p->SetValue("MAXJTABSETA", maxJtAbsEta);
  outEnv_p->SetValue("MINJTPT", minJtPt);

  outEnv_p->Write("config", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;    
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/makeClusterHist.exe <inConfigFileName>" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 2) retVal += makeClusterHist(argv[1]);
  return retVal;
}
