//AUTHOR: Chris McGinn (2019.12.06)

//cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

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
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/returnRootFileContentsList.h"
#include "include/stringUtil.h"

int deriveSampleWeights(std::string inXSectionFileName, std::string inEntriesFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inXSectionFileName, ".txt")) return 1;

  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  const std::string rootExt = ".root";
  const std::string txtExt = ".txt";

  bool entriesFileTXT = true;
  if(!check.checkFileExt(inEntriesFileName, txtExt)){
    entriesFileTXT = false;
    if(!check.checkFileExt(inEntriesFileName, rootExt)) return 1;
  }

  kirchnerPalette kPal;

  const std::string dateStr = getDateStr();

  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);

  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  std::map<std::string, std::vector<double> > jzMapToXSec;
  std::map<std::string, double> jzMapToEntries;
  std::map<std::string, double> jzMapToWeights;
  std::map<std::string, int> fileNameToJZVal;
  std::map<std::string, std::string> fileNameToTTreeStr;

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

  std::vector<std::vector<std::string> > globalTxtVect;
  bool inIsROOT = false;
  if(entriesFileTXT){
    std::ifstream inEntriesFile(inEntriesFileName.c_str());

    while(std::getline(inEntriesFile, tempStr)){
      if(tempStr.size() == 0) continue;
      std::vector<std::string> tempVect = commaSepStringToVect(tempStr);
      if(tempVect.size() == 0) continue;
      if(tempVect[0].size() == 0) continue;
      if(tempVect[0].substr(0,1).find("#") != std::string::npos) continue;
      
      globalTxtVect.push_back(tempVect);
      if(tempVect[1].rfind(".") != std::string::npos){
	if(isStrSame(rootExt, tempVect[1].substr(tempVect[1].rfind("."), rootExt.size()))) inIsROOT = true;
      }
    }
    inEntriesFile.close();

    for(unsigned int gI = 0; gI < globalTxtVect.size(); ++gI){
      if(!inIsROOT) jzMapToEntries[globalTxtVect[gI][0]] = std::stod(globalTxtVect[gI][1]); 
      else{
	if(!check.checkFileExt(globalTxtVect[gI][1], rootExt)) continue;
	TFile* inFile_p = new TFile(globalTxtVect[gI][1].c_str(), "READ");
	std::vector<std::string> treeList = returnRootFileContentsList(inFile_p, "TTree");
	
	if(treeList.size() == 0){
	  std::cout << "DERIVESAMPLEWEIGHTS WARNING: FILE \'" << globalTxtVect[gI][1] << "\' contains zero TTree. continue" << std::endl;
	  inFile_p->Close();
	  delete inFile_p;
	  continue;
	}
	else if(treeList.size() > 1){
	  if(globalTxtVect[gI].size() < 3){
	    std::cout << "DERIVESAMPLEWEIGHTS WARNING: FILE \'" << globalTxtVect[gI][1] << "\' contains more than 1 TTree and input txt file doesnt not disambiguate (specify as third arg in comma separated list). continue" << std::endl;
	    inFile_p->Close();
	    delete inFile_p;
	    continue;
	  }
	  else{
	    if(!vectContainsStr(globalTxtVect[gI][2], &treeList)){
	      std::cout << "DERIVESAMPLEWEIGHTS ERROR: Specified TTree \'" << globalTxtVect[gI][2] << "\' not found in corresponse TFile \'" << globalTxtVect[gI][1] << "\'. return 1" << std::endl;
	      inFile_p->Close();
	      delete inFile_p;
	      return 1;
	    }
	    else{
	      TTree* tempTree_p = (TTree*)inFile_p->Get(globalTxtVect[gI][2].c_str());
	      fileNameToTTreeStr[globalTxtVect[gI][1]] = globalTxtVect[gI][2];
	      jzMapToEntries[globalTxtVect[gI][0]] += (double)tempTree_p->GetEntries();
	    }	  
	  }
	}
	else{
	  TTree* tempTree_p = (TTree*)inFile_p->Get(treeList[0].c_str());
	  fileNameToTTreeStr[globalTxtVect[gI][1]] = treeList[0];
	  jzMapToEntries[globalTxtVect[gI][0]] += tempTree_p->GetEntries();
	}
	
	if(isStrSame(globalTxtVect[gI][0], "JZ1")) fileNameToJZVal[globalTxtVect[gI][1]] = 1;
	else if(isStrSame(globalTxtVect[gI][0], "JZ2")) fileNameToJZVal[globalTxtVect[gI][1]] = 2;
	else if(isStrSame(globalTxtVect[gI][0], "JZ3")) fileNameToJZVal[globalTxtVect[gI][1]] = 3;
	else if(isStrSame(globalTxtVect[gI][0], "JZ4")) fileNameToJZVal[globalTxtVect[gI][1]] = 4;
	else if(isStrSame(globalTxtVect[gI][0], "JZ5")) fileNameToJZVal[globalTxtVect[gI][1]] = 5;
	else if(isStrSame(globalTxtVect[gI][0], "JZ6")) fileNameToJZVal[globalTxtVect[gI][1]] = 6;
	else{
	  std::cout << "DERIVESAMPLEWEIGHTS ERROR: FILE NAME \'" << globalTxtVect[gI][1] << "\' contains no valid JZX. return 1" << std::endl;
	  inFile_p->Close();
	  delete inFile_p;
	  return 1;
	}
	
	inFile_p->Close();
	delete inFile_p;
      }
    }   
  }
  else{
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
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

      fileNameToTTreeStr[inEntriesFileName] = jzVal_;
      fileNameToTTreeStr[inEntriesFileName] = "clusterJetsCS";
    }

    inFile_p->Close();
    delete inFile_p;

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  }
  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "SAMPLE, X-SECTION, FILTER EFF., NUMBER EVENTS" << std::endl;
  double maxWeight = -1.0;
  for(auto const& iter : jzMapToXSec){
    if(jzMapToEntries[iter.first] <= 0) continue;

    std::cout << iter.first << ", " << iter.second[0] << ", " << iter.second[1] << ", " << prettyString(jzMapToEntries[iter.first], 1, false) << std::endl;

    jzMapToWeights[iter.first] = iter.second[0]*iter.second[1]/jzMapToEntries[iter.first];
    if(maxWeight < jzMapToWeights[iter.first]) maxWeight = jzMapToWeights[iter.first];

    std::cout << " " << jzMapToWeights[iter.first] << std::endl;
  }
  
  std::cout << "Maxweight: " << maxWeight << std::endl;

  std::cout << "RENORMALIZED WEIGHTS: " << std::endl;
  for(auto const& iter : jzMapToXSec){
    std::cout << iter.first << ": " << jzMapToWeights[iter.first]/maxWeight << " (N_Evt=" << prettyString(jzMapToEntries[iter.first], 1, false) << ")" << std::endl;
    jzMapToWeights[iter.first] /= maxWeight;
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::string outFileName = "";
  if(!entriesFileTXT || inIsROOT){
    if(!entriesFileTXT) outFileName = inEntriesFileName.substr(0, inEntriesFileName.rfind(rootExt));
    else outFileName = inEntriesFileName.substr(0, inEntriesFileName.rfind(txtExt));

    while(outFileName.find("/") != std::string::npos){
      outFileName.replace(0, outFileName.find("/")+1, "");
    }
    outFileName = "output/" + dateStr + "/" + outFileName + "_DerivedWeights_" + dateStr + rootExt;

    std::cout << "BUILDING WEIGHTS CANVAS" << std::endl;
    
    const Int_t nJtPtBins = 110;
    Float_t jtPtLow = 20.;
    Float_t jtPtHigh = 460.;
    Double_t jtPtBins[nJtPtBins+1];
    getLinBins(jtPtLow, jtPtHigh, nJtPtBins, jtPtBins);
    
    TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
    TH1D* jtpt_h = new TH1D("jtpt_h", ";Jet p_{T} [GeV];Effective Counts", nJtPtBins, jtPtBins);
    TH1D* jtpt_Weighted_h = new TH1D("jtpt_Weighted_h", ";Jet p_{T} [GeV];Effective Counts (Weighted)", nJtPtBins, jtPtBins);
    centerTitles({jtpt_h, jtpt_Weighted_h});

    //EDIT HERE CFM
    TFile* inFile_p = nullptr;
    std::vector<std::string> fileList;
    if(!entriesFileTXT) fileList.push_back(inEntriesFileName);
    else{
      //fileNameToJZStr mpa

      for(auto const & iter : fileNameToJZVal){
	fileList.push_back(iter.first);
      }
    }

    std::cout << "POST FILE LIST REWORK" << std::endl;

    const Int_t nMaxJets = 500;
    Int_t njtTruth_;
    Float_t jtptTruth_[nMaxJets];

    std::vector<float>* jtptTruth_p=nullptr;

    for(unsigned int fI = 0; fI < fileList.size(); ++fI){
      inFile_p = new TFile(fileList[fI].c_str(), "READ");
      TTree* inTree_p = (TTree*)inFile_p->Get(fileNameToTTreeStr[fileList[fI]].c_str());
      Int_t jzVal_;  
      
      inTree_p->SetBranchStatus("*", 0);
      if(!entriesFileTXT){ 
	inTree_p->SetBranchStatus("jzVal", 1);
	inTree_p->SetBranchStatus("njtTruth", 1);
	inTree_p->SetBranchStatus("jtptTruth", 1);
      
	inTree_p->SetBranchAddress("jzVal", &jzVal_);
	inTree_p->SetBranchAddress("njtTruth", &njtTruth_);
	inTree_p->SetBranchAddress("jtptTruth", jtptTruth_);
      }
      else{
	jzVal_ = fileNameToJZVal[fileList[fI]];

	inTree_p->SetBranchStatus("akt4_truth_jet_pt", 1);

	inTree_p->SetBranchAddress("akt4_truth_jet_pt", &jtptTruth_p);
      }

      const Int_t nEntries = inTree_p->GetEntries();
      
      for(Int_t entry = 0; entry < nEntries; ++entry){
	inTree_p->GetEntry(entry);

	if(entriesFileTXT){
	  njtTruth_ = 0;
	  for(unsigned int jI = 0; jI < jtptTruth_p->size(); ++jI){
	    jtptTruth_[njtTruth_] = (*jtptTruth_p)[jI];
	    ++njtTruth_;
	  }
	}
	
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
    }

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
    outFileName = inEntriesFileName.substr(0, inEntriesFileName.rfind(txtExt));
    while(outFileName.find("/") != std::string::npos){
      outFileName.replace(0, outFileName.find("/")+1, "");
    }
    outFileName = "output/" + dateStr + "/" + outFileName + "_DerivedWeights_" + dateStr + rootExt;
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  outFileName.replace(outFileName.rfind(rootExt), 5, "");
  outFileName = outFileName + txtExt;
  
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
