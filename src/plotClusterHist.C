//cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TBox.h"
#include "TCanvas.h"
#include "TEnv.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/returnRootFileContentsList.h"
#include "include/sharedFunctions.h"
#include "include/stringUtil.h"

const int nColors=15;
const int colors[nColors]={1,kRed-4,kAzure-3,kGreen+2,kMagenta+2,kOrange+2,kRed-4,kAzure-3,kGreen+2,kMagenta+2,kOrange+2,kCyan+3,28,41,kGray};

const int nStyles=15;
const int styles[nStyles]={20,21,33,34,29,24,25,27,28,30,23,20,21,33,34};

const int nSizes = 15;
const float sizes[nSizes]={1,1,1.6,1.2,1.6,1,1,1,1.6,1,1,1,1,1.6,1.2};

void drawWhiteBoxNDC(TCanvas* canv_p, Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  canv_p->cd();
  TBox* box_p = new TBox();
  box_p->SetFillColor(0);
  box_p->DrawBox(x1, y1, x2, y2);
  delete box_p;

  return;
}

void configHist(TH1* inHist_p, Int_t pos, bool doDebug=false)
{
  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  inHist_p->SetMarkerSize(sizes[pos%nSizes]);
  inHist_p->SetMarkerStyle(styles[pos%nStyles]);
  inHist_p->SetMarkerColor(colors[pos%nColors]);
  inHist_p->SetLineColor(colors[pos%nColors]);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  inHist_p->GetXaxis()->SetTitleFont(42);
  inHist_p->GetYaxis()->SetTitleFont(42);
  inHist_p->GetXaxis()->SetLabelFont(42);
  inHist_p->GetYaxis()->SetLabelFont(42);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  inHist_p->GetXaxis()->SetTitleSize(0.035);
  inHist_p->GetYaxis()->SetTitleSize(0.035);
  inHist_p->GetXaxis()->SetLabelSize(0.03);
  inHist_p->GetYaxis()->SetLabelSize(0.03);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  inHist_p->GetYaxis()->SetTitleOffset(1.9);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  return;
}

void configHist(TGraph* inHist_p, Int_t pos, bool doDebug=false)
{
  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pos << ", " << inHist_p << std::endl;

  inHist_p->SetMarkerSize(1);
  inHist_p->SetMarkerStyle(styles[pos%nStyles]);
  inHist_p->SetMarkerColor(colors[pos%nColors]);
  inHist_p->SetLineColor(colors[pos%nColors]);
 
  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pos << ", " << inHist_p << std::endl;
 
  inHist_p->GetXaxis()->SetTitleFont(42);
  inHist_p->GetYaxis()->SetTitleFont(42);
  inHist_p->GetXaxis()->SetLabelFont(42);
  inHist_p->GetYaxis()->SetLabelFont(42);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pos << ", " << inHist_p << std::endl;
  
  inHist_p->GetXaxis()->SetTitleSize(0.035);
  inHist_p->GetYaxis()->SetTitleSize(0.035);
  inHist_p->GetXaxis()->SetLabelSize(0.03);
  inHist_p->GetYaxis()->SetLabelSize(0.03);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pos << ", " << inHist_p << std::endl;
  
  inHist_p->GetYaxis()->SetTitleOffset(1.9);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << pos << ", " << inHist_p << std::endl;
  
  return;
}

void plotResponseSet(TEnv* fileConfig_p, TEnv* inConfig_p, std::vector<TH1*> histSet, std::vector<std::string> setLabels, std::string globalNameStr, std::string globalLabelStr, std::string dateStr)
{
  const int globalFont = 42;
  const double globalSize = 0.035;

  std::string frontStr = "";
  if(globalNameStr.find("recoGenSigma_DeltaEta") != std::string::npos) frontStr = "DELTAETASIGMA";
  else if(globalNameStr.find("recoGenSigma_DeltaPhi") != std::string::npos) frontStr = "DELTAPHISIGMA";
  else if(globalNameStr.find("recoOverGenSigmaOverMean") != std::string::npos) frontStr = "RECOOVERGENSIGMA";
  else if(globalNameStr.find("recoOverGenSigma") != std::string::npos) frontStr = "RECOOVERGENSIGMA";
  else if(globalNameStr.find("recoOverGenMean") != std::string::npos) frontStr = "RECOOVERGENMEAN";
  else if(globalNameStr.find("recoOverGenMSigmaOverMean") != std::string::npos) frontStr = "RECOOVERGENMSIGMA";
  else if(globalNameStr.find("recoOverGenMSigma") != std::string::npos) frontStr = "RECOOVERGENMSIGMA";
  else if(globalNameStr.find("recoOverGenMMean") != std::string::npos) frontStr = "RECOOVERGENMMEAN";

  const double max = inConfig_p->GetValue((frontStr + "MAX").c_str(), 1.1);
  const double min = inConfig_p->GetValue((frontStr + "MIN").c_str(), 0.0);
  const double legX = inConfig_p->GetValue((frontStr + "LEGX").c_str(), 0.6);
  const double legY = inConfig_p->GetValue((frontStr + "LEGY").c_str(), 0.6);
  const double labelX = inConfig_p->GetValue((frontStr + "LABELX").c_str(), 0.2);
  const double labelY = inConfig_p->GetValue((frontStr + "LABELY").c_str(), 0.6);
  const bool labelAlignRight = inConfig_p->GetValue((frontStr + "LABELALIGNRIGHT").c_str(), 0);
  
  const bool doLogX = inConfig_p->GetValue((frontStr + "DOLOGX").c_str(), 0);
  const bool doAllLeg = inConfig_p->GetValue((frontStr + "DOALLLEG").c_str(), 1);

  std::map<std::string, std::string> algoNameSwaps;
  for(Int_t i = 0; i < 20; ++i){
    std::string tempStr = inConfig_p->GetValue(("ALGONAMESUB." + std::to_string(i)).c_str(), "");
    if(tempStr.size() == 0) continue;
    std::vector<std::string> tempVect = strToVect(tempStr);
    algoNameSwaps[tempVect[0]] = tempVect[1];
  }

  std::map<std::string, bool> permaLeg;
  if(!doAllLeg){
    for(Int_t eI = 0; eI < 20; ++eI){
      std::string tempStr = inConfig_p->GetValue((frontStr + "PERMALEG." + std::to_string(eI)).c_str(), "");
      if(tempStr.size() == 0) continue;
      permaLeg[tempStr] = true;
    }
  }

  std::string caloTrackStr = fileConfig_p->GetValue("CALOTRACKSTR", "");

  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.16);
  canv_p->SetBottomMargin(0.14);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(globalFont);
  label_p->SetTextSize(globalSize);
  
  double maxLegStr = 0.0;
  int nLeg = 0;
  for(unsigned int hI = 0; hI < histSet.size(); ++hI){
    ++nLeg;
    
    if(setLabels[hI].find("Cent") != std::string::npos){
      setLabels[hI].replace(setLabels[hI].find("Cent"), 4, "");
      setLabels[hI].replace(setLabels[hI].find("to"), 2, "-");
      setLabels[hI] = setLabels[hI] + "%";
    }

    label_p->SetText(0.1, 0.1, setLabels[hI].c_str());
    if(label_p->GetXsize() > maxLegStr) maxLegStr = label_p->GetXsize();
  }

  TLegend* leg_p = new TLegend(legX, legY - 0.04*nLeg, legX + 0.8*maxLegStr, legY);
  leg_p->SetTextFont(globalFont);
  leg_p->SetTextSize(globalSize);
  
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  for(unsigned int hI = 0; hI < histSet.size(); ++hI){
    histSet[hI]->SetMaximum(max);
    histSet[hI]->SetMinimum(min);

    if(hI == 0) histSet[hI]->DrawCopy("HIST E1 P");
    else histSet[hI]->DrawCopy("HIST E1 P SAME");

    std::string newLabel = setLabels[hI];
    if(newLabel.find("Cent") != std::string::npos){
      newLabel.replace(newLabel.find("Cent"), 4, "");
      newLabel.replace(newLabel.find("to"), 2, "-");
      newLabel = newLabel + "%";
    }
    leg_p->AddEntry(histSet[hI], newLabel.c_str(), "P L");
  }

  if(doAllLeg || permaLeg.count(globalLabelStr) != 0) leg_p->Draw("SAME");

  std::string globalLabelStr2 = globalLabelStr;
  if(globalLabelStr2.find("Cent") != std::string::npos){
    globalLabelStr2.replace(globalLabelStr2.find("Cent"), 4, "");
    globalLabelStr2.replace(globalLabelStr2.find("to"), 2, "-");
    globalLabelStr2 = globalLabelStr2 + "%";
    globalLabelStr2 = "#bf{#color[" + std::to_string(kRed) + "]{" + globalLabelStr2 + "}}";
  }
  else{
    if(algoNameSwaps.count(globalLabelStr2) != 0) globalLabelStr2 = algoNameSwaps[globalLabelStr2];
  }

  std::string minJtPtStr = fileConfig_p->GetValue("MINJTPT", "");
  std::string maxJtAbsEtaStr = fileConfig_p->GetValue("MAXJTABSETA", "");

  while(minJtPtStr.find(".") != std::string::npos && minJtPtStr.find(".") != minJtPtStr.size()-2){
    if(minJtPtStr.substr(minJtPtStr.size()-1, 1).find("0") != std::string::npos) minJtPtStr.replace(minJtPtStr.size()-1, 1, "");
    else break;
  }

  while(maxJtAbsEtaStr.find(".") != std::string::npos && maxJtAbsEtaStr.find(".") != maxJtAbsEtaStr.size()-2){
    if(maxJtAbsEtaStr.substr(maxJtAbsEtaStr.size()-1, 1).find("0") != std::string::npos) maxJtAbsEtaStr.replace(maxJtAbsEtaStr.size()-1, 1, "");
    else break;
  }

  std::vector<std::string> globalLabels;
  for(Int_t i = 0; i < 20; ++i){
    std::string tempStr = inConfig_p->GetValue(("GLOBALLABEL." + std::to_string(i)).c_str(), "");
    if(tempStr.size() == 0) continue;
    
    while(tempStr.find("\"") != std::string::npos){
      tempStr.replace(tempStr.find("\""), 1, "");
    }

    globalLabels.push_back(tempStr);
  }

  std::vector<std::string> labels = globalLabels;
  labels.push_back(globalLabelStr2);
  labels.push_back(caloTrackStr + " jet");
  labels.push_back("|#eta_{Jet}| < " + maxJtAbsEtaStr + "; p_{T,Reco.} > " + minJtPtStr);
  //  labels.push_back("p_{T,Reco.} > " + minJtPtStr);

  double maxLabelStr = 0.0;
  if(labelAlignRight){
    for(unsigned int lI = 0; lI < labels.size(); ++lI){
      label_p->SetText(labelX, labelY, labels[lI].c_str());
      if(label_p->GetXsize() > maxLabelStr) maxLabelStr = label_p->GetXsize();
    }
  }

  if(!doAllLeg && permaLeg.count(globalLabelStr) == 0){
    labels.clear();
    for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
      labels.push_back("");
    }
    labels.push_back(globalLabelStr2);
  }
  
  for(unsigned int lI = 0; lI < labels.size(); ++lI){
    if(labelAlignRight){
      label_p->SetText(labelX, labelY, labels[lI].c_str());
      while(label_p->GetXsize() < maxLabelStr){
	labels[lI] = " " + labels[lI];
	label_p->SetText(labelX, labelY, labels[lI].c_str());
      }
    }
    
    label_p->DrawLatex(labelX, labelY - 0.045*lI, labels[lI].c_str());
  }

  
  if(doLogX) gPad->SetLogy();
    
  for(Int_t eI = 0; eI < 20; ++eI){
    std::string tempStr = inConfig_p->GetValue((frontStr + "WHITEBOX." + std::to_string(eI)).c_str(), "");
    if(tempStr.size() == 0) continue;
    
    std::vector<float> pos = strToVectF(tempStr);
    if(pos.size() == 4) drawWhiteBoxNDC(canv_p, pos[0], pos[1], pos[2], pos[3]);
  }
  
  std::string saveName = "pdfDir/" + dateStr + "/" + globalNameStr + "_" + globalLabelStr + "_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;
  delete leg_p;
  delete label_p;
  
  return;
}

int plotClusterHist(std::string inConfigFileName)
{
  globalDebugHandler gDebug;
  const bool doDebug = gDebug.GetDoGlobalDebug();

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const int globalFont = 42;
  const double globalSize = 0.035;

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"INFILENAME",
					"GLOBALTAG",
					"CALOTRACKSTR"};

  if(!checkConfigContainsParams(inConfig_p, reqParams)) return 1;                         
  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;
  const std::string globalStr = inConfig_p->GetValue("GLOBALTAG", "");
  const std::string caloTrackStr = inConfig_p->GetValue("CALOTRACKSTR", "");
  

  std::string globalStr2 = globalStr;
  if(globalStr2.size() == 0) globalStr2 = "NoGlobalStr";
  while(globalStr2.find(" ") != std::string::npos){
    globalStr2.replace(globalStr2.find(" "), 1, "");
  }

  const std::string dateStr = getDateStr();

  std::vector<std::string> globalLabels;
  for(Int_t i = 0; i < 20; ++i){
    std::string tempStr = inConfig_p->GetValue(("GLOBALLABEL." + std::to_string(i)).c_str(), "");
    if(tempStr.size() == 0) continue;
    
    while(tempStr.find("\"") != std::string::npos){
      tempStr.replace(tempStr.find("\""), 1, "");
    }

    globalLabels.push_back(tempStr);
  }

  std::vector<std::string> algosToSkip;
  for(Int_t i = 0; i < 20; ++i){
    std::string tempStr = inConfig_p->GetValue(("ALGOTOSKIP." + std::to_string(i)).c_str(), "");
    if(tempStr.size() == 0) continue;
    algosToSkip.push_back(tempStr);
  }

  std::map<std::string, std::string> algoNameSwaps;
  for(Int_t i = 0; i < 20; ++i){
    std::string tempStr = inConfig_p->GetValue(("ALGONAMESUB." + std::to_string(i)).c_str(), "");
    if(tempStr.size() == 0) continue;
    std::vector<std::string> tempVect = strToVect(tempStr);
    algoNameSwaps[tempVect[0]] = tempVect[1];
  }

  kirchnerPalette kPal;

  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const Int_t nMaxJtPtBins = 50;
  
  const Int_t nMaxCentBins = 20;
  const Int_t nMaxJtAlgo = 20;
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* fileConfig_p = (TEnv*)inFile_p->Get("config");

  std::vector<std::string> fileParams = {"ISMC",
					 "NJTPTBINS",
					 "JTPTBINS",
					 "JTPTBINSSTR",
					 "NCENTBINS",
					 "CENTBINSHIGH",
					 "CENTBINSSTR", 
					 "NEVENTPERCENT",
					 "MAXJTABSETA",
					 "MINJTPT", 
					 "JTALGOS"};  

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(!checkConfigContainsParams(fileConfig_p, fileParams)) return 1;

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const Bool_t isMC = fileConfig_p->GetValue("ISMC", 0);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  fileConfig_p->SetValue("CALOTRACKSTR", caloTrackStr.c_str());
  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const Int_t nJtPtBins = fileConfig_p->GetValue("NJTPTBINS", 0);
  std::vector<std::string> jtPtBinsStr = commaSepStringToVect(fileConfig_p->GetValue("JTPTBINS", ""));

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  Double_t jtPtBins[nMaxJtPtBins+1];
  for(unsigned int jI = 0; jI < jtPtBinsStr.size(); ++jI){
    jtPtBins[jI] = std::stod(jtPtBinsStr[jI]);
  }
  jtPtBinsStr = commaSepStringToVect(fileConfig_p->GetValue("JTPTBINSSTR", ""));

  
  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const Int_t nCentBins = fileConfig_p->GetValue("NCENTBINS", 0);
  std::vector<std::string> centBinsStr = commaSepStringToVect(fileConfig_p->GetValue("CENTBINSSTR", ""));
  std::vector<std::string> centBinsHighStr = commaSepStringToVect(fileConfig_p->GetValue("CENTBINSHIGH", ""));
  std::vector<std::string> nEventPerCentStr = commaSepStringToVect(fileConfig_p->GetValue("NEVENTPERCENT", ""));
  const std::string jtAbsEtaMaxStr = fileConfig_p->GetValue("MAXJTABSETA", "");
  const double minJtPt = fileConfig_p->GetValue("MINJTPT", 0.0);
  std::string minJtPtStr = prettyString(minJtPt, 1, false);
  std::vector<Double_t> nEventPerCent;
  for(unsigned int cI = 0; cI < nEventPerCentStr.size(); ++cI){
    nEventPerCent.push_back(std::stod(nEventPerCentStr[cI]));
  }

  std::string maxJtAbsEtaStr = jtAbsEtaMaxStr;

  while(minJtPtStr.find(".") != std::string::npos && minJtPtStr.find(".") != minJtPtStr.size()-2){
    if(minJtPtStr.substr(minJtPtStr.size()-1, 1).find("0") != std::string::npos) minJtPtStr.replace(minJtPtStr.size()-1, 1, "");
    else break;
  }

  while(maxJtAbsEtaStr.find(".") != std::string::npos && maxJtAbsEtaStr.find(".") != maxJtAbsEtaStr.size()-2){
    if(maxJtAbsEtaStr.substr(maxJtAbsEtaStr.size()-1, 1).find("0") != std::string::npos) maxJtAbsEtaStr.replace(maxJtAbsEtaStr.size()-1, 1, "");
    else break;
  }

  Double_t periphMostVal = -1;
  Int_t periphMostBin = -1;
  for(unsigned int cI = 0; cI < centBinsHighStr.size(); ++cI){
    double centVal = std::stoi(centBinsHighStr[cI]);
    if(centVal >= 70 && centVal >= periphMostVal){
      periphMostBin = cI;
      periphMostVal = centVal;
    }
  }

  std::vector<std::string> jtAlgos = commaSepStringToVect(fileConfig_p->GetValue("JTALGOS", ""));
  unsigned int pos = 0;
  while(pos < jtAlgos.size()){
    if(jtAlgos[pos].find(caloTrackStr) == std::string::npos) jtAlgos.erase(jtAlgos.begin() + pos);
    else if(caloTrackStr.find("Tower") != std::string::npos && jtAlgos[pos].find("4GeVCut") != std::string::npos) jtAlgos.erase(jtAlgos.begin() + pos);    
    else if(vectContainsStr(jtAlgos[pos], &algosToSkip)) jtAlgos.erase(jtAlgos.begin() + pos);
    else ++pos;
  }

  const Int_t nJtAlgo = jtAlgos.size();

  if(nCentBins > nMaxCentBins){
    std::cout << "nCentBins \'" << nCentBins << "\' exceeds maximum \'" << nMaxCentBins << "\'. return 1" << std::endl;
    return 1;
  }

  if(nJtAlgo > nMaxJtAlgo){
    std::cout << "nJtAlgo \'" << nJtAlgo << "\' exceeds maximum \'" << nMaxJtAlgo << "\'. return 1" << std::endl;
    return 1;
  }

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  TH1D* spectra_p[nMaxJtAlgo][nMaxCentBins];
  TH1D* spectraUnmatched_p[nMaxJtAlgo][nMaxCentBins];
  TH1D* spectraTruth_p[nMaxCentBins];
  TH1D* spectraChgTruth_p[nMaxCentBins];

  TH1D* recoOverGen_VPt_p[nMaxJtAlgo][nMaxCentBins][nMaxJtPtBins];
  TH1D* recoOverGenM_VPt_p[nMaxJtAlgo][nMaxCentBins][nMaxJtPtBins];
  TH1D* recoOverGenMOverPt_VPt_p[nMaxJtAlgo][nMaxCentBins][nMaxJtPtBins];
  TF1* recoOverGenFit_p[nMaxJtAlgo][nMaxCentBins][nMaxJtPtBins];
  TH1D* recoGen_DeltaEta_p[nMaxJtAlgo][nMaxCentBins][nMaxJtPtBins];
  TH1D* recoGen_DeltaPhi_p[nMaxJtAlgo][nMaxCentBins][nMaxJtPtBins];
  TH1D* matchedTruthSpectra_p[nMaxJtAlgo][nMaxCentBins];
  TGraphAsymmErrors* truthEff_p[nMaxJtAlgo][nMaxCentBins];
  TGraphAsymmErrors* fakeRate_p[nMaxJtAlgo][nMaxCentBins];

  for(Int_t aI = 0; aI < nJtAlgo; ++ aI){
    for(Int_t cI = 0; cI < nCentBins; ++ cI){
      std::string nameStr = jtAlgos[aI] + "_" + centBinsStr[cI];
      
      spectra_p[aI][cI] = (TH1D*)inFile_p->Get(("spectra_" + nameStr + "_h").c_str());   
      if(isMC){
	spectraUnmatched_p[aI][cI] = (TH1D*)inFile_p->Get(("spectraUnmatched_" + nameStr + "_h").c_str());   
	fakeRate_p[aI][cI] = new TGraphAsymmErrors();
	fakeRate_p[aI][cI]->BayesDivide(spectraUnmatched_p[aI][cI], spectra_p[aI][cI]);
      }
      
      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
      if(isMC){
	if(aI == 0){
	  spectraTruth_p[cI] = (TH1D*)inFile_p->Get(("spectra_Truth_" + centBinsStr[cI] + "_h").c_str());
	  spectraChgTruth_p[cI] = (TH1D*)inFile_p->Get(("spectraChg_Truth_" + centBinsStr[cI] + "_h").c_str());
	}

	if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	  matchedTruthSpectra_p[aI][cI] = (TH1D*)inFile_p->Get(("matchedTruthSpectra_" + nameStr + "_h").c_str());       
	  truthEff_p[aI][cI] = new TGraphAsymmErrors();
	  
	  std::cout << nameStr << std::endl;

	  if(nameStr.find("Trk") == std::string::npos) truthEff_p[aI][cI]->BayesDivide(matchedTruthSpectra_p[aI][cI], spectraTruth_p[cI]);
	  else truthEff_p[aI][cI]->BayesDivide(matchedTruthSpectra_p[aI][cI], spectraChgTruth_p[cI]);
	}
      }

      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      setSumW2(spectra_p[aI][cI]);

      spectra_p[aI][cI]->Scale(1./nEventPerCent[cI]);
      for(Int_t bIX = 0; bIX < spectra_p[aI][cI]->GetXaxis()->GetNbins(); ++bIX){
	double val = spectra_p[aI][cI]->GetBinContent(bIX+1)/spectra_p[aI][cI]->GetBinWidth(bIX+1);
	double err = spectra_p[aI][cI]->GetBinError(bIX+1)/spectra_p[aI][cI]->GetBinWidth(bIX+1);

	spectra_p[aI][cI]->SetBinContent(bIX+1, val);
	spectra_p[aI][cI]->SetBinError(bIX+1, err);
      }

      configHist(spectra_p[aI][cI], aI, doDebug);
      spectra_p[aI][cI]->GetYaxis()->SetTitle("#frac{1}{N_{Event}} #frac{dN_{Jet}}{dp_{T}} [Gev^{-1}]");
      centerTitles(spectra_p[aI][cI]);
    
      if(isMC){
	if(jtAlgos[aI].find("ATLAS") != std::string::npos) continue;
	if(jtAlgos[aI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;
	if(jtAlgos[aI].find("Truth") != std::string::npos) continue;
	
	for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	  recoOverGen_VPt_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoOverGen_VPt_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());
	  recoOverGenM_VPt_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoOverGenM_VPt_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());
	  recoOverGenMOverPt_VPt_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoOverGenMOverPt_VPt_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());

	  recoGen_DeltaEta_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoGen_DeltaEta_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());
	  recoGen_DeltaPhi_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoGen_DeltaPhi_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());


	  configHist(recoOverGen_VPt_p[aI][cI][jI], aI, doDebug);
	  configHist(recoOverGenM_VPt_p[aI][cI][jI], aI, doDebug);
	  configHist(recoOverGenMOverPt_VPt_p[aI][cI][jI], aI, doDebug);
	  configHist(recoGen_DeltaEta_p[aI][cI][jI], aI, doDebug);
	  configHist(recoGen_DeltaPhi_p[aI][cI][jI], aI, doDebug);

	  recoOverGen_VPt_p[aI][cI][jI]->SetMarkerColor(1);
	  recoOverGen_VPt_p[aI][cI][jI]->SetLineColor(1);

	  recoOverGenM_VPt_p[aI][cI][jI]->SetMarkerColor(1);
	  recoOverGenM_VPt_p[aI][cI][jI]->SetLineColor(1);

	  recoOverGenMOverPt_VPt_p[aI][cI][jI]->SetMarkerColor(1);
	  recoOverGenMOverPt_VPt_p[aI][cI][jI]->SetLineColor(1);

	  recoGen_DeltaEta_p[aI][cI][jI]->SetMarkerColor(1);
	  recoGen_DeltaEta_p[aI][cI][jI]->SetLineColor(1);

	  recoGen_DeltaPhi_p[aI][cI][jI]->SetMarkerColor(1);
	  recoGen_DeltaPhi_p[aI][cI][jI]->SetLineColor(1);
	}
      }
    }
  }

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(globalFont);
  label_p->SetTextSize(globalSize);

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.16);
    canv_p->SetBottomMargin(0.14);

    TLegend* leg_p = new TLegend(0.6, 0.75, 0.95, 0.98);
    leg_p->SetTextFont(43);
    leg_p->SetTextSize(16);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);
    
    Double_t globalMin = 1000000;    
    Double_t globalMax = -1;
    
    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      globalMin = TMath::Min(globalMin, getMinGTZero(spectra_p[aI][cI]));
      globalMax = TMath::Max(globalMax, getMax(spectra_p[aI][cI]));
    }
    
    bool isDrawn = false;
    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      if(caloTrackStr.find("Trk") != std::string::npos && jtAlgos[aI].find("ATLAS") != std::string::npos) continue;
      if(jtAlgos[aI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;


      std::string tempStr = jtAlgos[aI];
      if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[aI]];
      leg_p->AddEntry(spectra_p[aI][cI], tempStr.c_str(), "P L");
		      
      spectra_p[aI][cI]->SetMaximum(3.*globalMax);
      spectra_p[aI][cI]->SetMinimum(globalMin/2.);

      if(!isDrawn) spectra_p[aI][cI]->DrawCopy("HIST E1 P");
      else spectra_p[aI][cI]->DrawCopy("HIST E1 P SAME");

      isDrawn = true;
    }

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
    centLabel.replace(centLabel.find("to"), 2, "-");
    centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
    if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";
    
    if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.264, 0.83, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
    if(isStrSame(caloTrackStr, "Tower")) label_p->DrawLatex(0.264, 0.89, "Calo. jets");
    else if(isStrSame(caloTrackStr, "Trk")) label_p->DrawLatex(0.264, 0.89, "Track jets");
    
    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    gPad->SetLogy();
    gStyle->SetOptStat(0);

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    leg_p->Draw("SAME");
    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    std::string saveName = "pdfDir/" + dateStr + "/spectraComp_" + globalStr2 + "_" + caloTrackStr + "_" + centBinsStr[cI] + "_" + dateStr;
    quietSaveAs(canv_p, saveName + ".pdf");
    quietSaveAs(canv_p, saveName + ".png");

    delete leg_p;
    delete canv_p;
  }

  if(isMC){
    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    const double effMax = inConfig_p->GetValue("EFFMAX", 1.1);
    const double effMin = inConfig_p->GetValue("EFFMIN", 0.0);
    const double effLegX = inConfig_p->GetValue("EFFLEGX", 0.6);
    const double effLegY = inConfig_p->GetValue("EFFLEGY", 0.6);
    const double effLabelX = inConfig_p->GetValue("EFFLABELX", 0.2);
    const double effLabelY = inConfig_p->GetValue("EFFLABELY", 0.6);
    const bool effLabelAlignRight = inConfig_p->GetValue("EFFLABELALIGNRIGHT", 0);
    const bool effDoLogX = inConfig_p->GetValue("EFFDOLOGX", 0);

    const bool effDoAllLeg = inConfig_p->GetValue("EFFDOALLLEG", 1);

  
    double maxLabelStr = 0.0;
    if(effLabelAlignRight){
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.16);
      canv_p->SetBottomMargin(0.14);

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	std::vector<std::string> labels = globalLabels;
	std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	centLabel.replace(centLabel.find("to"), 2, "-");
	centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";

	labels.push_back(centLabel);
	labels.push_back("|#eta_{Jet}| < " + maxJtAbsEtaStr + "; p_{T,Reco.} > " + minJtPtStr);
	labels.push_back(caloTrackStr + " jets");

	for(unsigned int lI = 0; lI < labels.size(); ++lI){       
	  label_p->SetText(effLabelX, effLabelY, labels[lI].c_str());
	  if(label_p->GetXsize() > maxLabelStr) maxLabelStr = label_p->GetXsize();
	}
      }
      delete canv_p;
    }

    std::map<std::string, bool> effPermaLeg;
    if(!effDoAllLeg){
      for(Int_t eI = 0; eI < 20; ++eI){
	std::string tempStr = inConfig_p->GetValue(("EFFPERMALEG." + std::to_string(eI)).c_str(), "");
	if(tempStr.size() == 0) continue;
	
	effPermaLeg[tempStr] = true;
      }
    }

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::vector<std::string> labels = globalLabels;

      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.16);
      canv_p->SetBottomMargin(0.14);

      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      double maxLegStr = 0.0;
      int nLeg = 0;
      for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	  if(jtAlgos[aI].find("NoSub") == std::string::npos || cI == periphMostBin){
	    ++nLeg;

	    std::string tempStr = jtAlgos[aI];
	    if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[aI]];
	    label_p->SetText(0.1, 0.1, tempStr.c_str());
	    if(label_p->GetXsize() > maxLegStr) maxLegStr = label_p->GetXsize();
	  }
	}	
      }

      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      TLegend* leg_p = new TLegend(effLegX, effLegY - 0.04*nLeg, effLegX + maxLegStr*0.8, effLegY);
      leg_p->SetTextFont(globalFont);
      leg_p->SetTextSize(globalSize);
      leg_p->SetBorderSize(0);
      leg_p->SetFillColor(0);
      leg_p->SetFillStyle(0);
      
    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      TH1D* dummyHist_p = new TH1D("dummyHist_h", ";Truth Jet p_{T} [GeV];Efficiency", nJtPtBins, jtPtBins);
      centerTitles(dummyHist_p);
      dummyHist_p->SetMaximum(effMax);
      dummyHist_p->SetMinimum(effMin);

      dummyHist_p->DrawCopy("");
      
      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	  if(jtAlgos[aI].find("NoSub") == std::string::npos || cI == periphMostBin){
	    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    configHist(truthEff_p[aI][cI], aI, doDebug);	
	    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	    std::string tempStr = jtAlgos[aI];
	    if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[aI]];
	    
	    leg_p->AddEntry(truthEff_p[aI][cI], tempStr.c_str(), "P L");
	    truthEff_p[aI][cI]->Draw("P");
	  }
	}
      }

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(effMax > 1.0 && effMin < 1.0) line_p->DrawLine(jtPtBins[0], 1.0, jtPtBins[nJtPtBins], 1.0);
      
      std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
      centLabel.replace(centLabel.find("to"), 2, "-");
      centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
      //      if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";

      
      if(effDoAllLeg || effPermaLeg.count(centBinsStr[cI]) != 0){
	labels.push_back(centLabel);
	labels.push_back("|#eta_{jet}| < " + maxJtAbsEtaStr + "; p_{T,Reco.} > " + minJtPtStr);
	labels.push_back(caloTrackStr + " jets");
      }
      else{
	labels.clear();
	for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
	  labels.push_back("");
	}
	labels.push_back(centLabel);
      }
    
      for(unsigned int lI = 0; lI < labels.size(); ++lI){
	if(effLabelAlignRight){
	  label_p->SetText(effLabelX, effLabelY, labels[lI].c_str());
	  while(label_p->GetXsize() < maxLabelStr){
	    labels[lI] = " " + labels[lI];
	    label_p->SetText(effLabelX, effLabelY, labels[lI].c_str());
	  }
	}

	label_p->DrawLatex(effLabelX, effLabelY - 0.045*lI, labels[lI].c_str());
      }

      gStyle->SetOptStat(0);      

      if(effDoAllLeg || effPermaLeg.count(centBinsStr[cI]) != 0) leg_p->Draw("SAME");

      if(effDoLogX) gPad->SetLogx();

      for(Int_t eI = 0; eI < 20; ++eI){
	std::string tempStr = inConfig_p->GetValue(("EFFWHITEBOX." + std::to_string(eI)).c_str(), "");
	if(tempStr.size() == 0) continue;
      
	std::vector<float> pos = strToVectF(tempStr);
	if(pos.size() == 4) drawWhiteBoxNDC(canv_p, pos[0], pos[1], pos[2], pos[3]);
      }

      
      std::string saveName = "pdfDir/" + dateStr + "/eff_" + globalStr2 + "_" + caloTrackStr + "_" + centBinsStr[cI] + "_" + dateStr;
      quietSaveAs(canv_p, saveName + ".pdf");
      //      quietSaveAs(canv_p, saveName + ".png");
      
      delete dummyHist_p;
      delete leg_p;
      delete canv_p;
    }

    const double fakeMax = inConfig_p->GetValue("FAKEMAX", 1.1);
    const double fakeMin = inConfig_p->GetValue("FAKEMIN", 0.0);
    const double fakeLegX = inConfig_p->GetValue("FAKELEGX", 0.6);
    const double fakeLegY = inConfig_p->GetValue("FAKELEGY", 0.6);
    const double fakeLabelX = inConfig_p->GetValue("FAKELABELX", 0.2);
    const double fakeLabelY = inConfig_p->GetValue("FAKELABELY", 0.6);
    const bool fakeLabelAlignRight = inConfig_p->GetValue("FAKELABELALIGNRIGHT", 0);
    const bool fakeDoLogX = inConfig_p->GetValue("FAKEDOLOGX", 0);
    const bool fakeDoLogY = inConfig_p->GetValue("FAKEDOLOGY", 0);

    const bool fakeDoAllLeg = inConfig_p->GetValue("FAKEDOALLLEG", 1);

  
    maxLabelStr = 0.0;
    if(fakeLabelAlignRight){
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.16);
      canv_p->SetBottomMargin(0.14);

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	std::vector<std::string> labels = globalLabels;
	std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	centLabel.replace(centLabel.find("to"), 2, "-");
	centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";

	labels.push_back(centLabel);
	labels.push_back("|#eta_{Jet}| < " + maxJtAbsEtaStr + "; p_{T,Reco.} > " + minJtPtStr);
	labels.push_back(caloTrackStr + " jets");

	for(unsigned int lI = 0; lI < labels.size(); ++lI){       
	  label_p->SetText(fakeLabelX, fakeLabelY, labels[lI].c_str());
	  if(label_p->GetXsize() > maxLabelStr) maxLabelStr = label_p->GetXsize();
	}
      }
      delete canv_p;
    }

    std::map<std::string, bool> fakePermaLeg;
    if(!fakeDoAllLeg){
      for(Int_t eI = 0; eI < 20; ++eI){
	std::string tempStr = inConfig_p->GetValue(("FAKEPERMALEG." + std::to_string(eI)).c_str(), "");
	if(tempStr.size() == 0) continue;
	
	fakePermaLeg[tempStr] = true;
      }
    }

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::vector<std::string> labels = globalLabels;

      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.16);
      canv_p->SetBottomMargin(0.14);

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      double maxLegStr = 0.0;
      int nLeg = 0;
      for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	  if(jtAlgos[aI].find("NoSub") == std::string::npos || cI == periphMostBin){
	    ++nLeg;
	    std::string tempStr = jtAlgos[aI];
	    if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[aI]];
	    label_p->SetText(0.1, 0.1, tempStr.c_str());
	    if(label_p->GetXsize() > maxLegStr) maxLegStr = label_p->GetXsize();
	  }
	}	
      }

      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      TLegend* leg_p = new TLegend(fakeLegX, fakeLegY - 0.04*nLeg, fakeLegX + maxLegStr*0.8, fakeLegY);
      leg_p->SetTextFont(globalFont);
      leg_p->SetTextSize(globalSize);
      leg_p->SetBorderSize(0);
      leg_p->SetFillColor(0);
      leg_p->SetFillStyle(0);
      
    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      TH1D* dummyHist_p = new TH1D("dummyHist_h", ";Reco. Jet p_{T} [GeV];Fake fraction", nJtPtBins, jtPtBins);
      centerTitles(dummyHist_p);
      dummyHist_p->SetMaximum(fakeMax);
      dummyHist_p->SetMinimum(fakeMin);

      dummyHist_p->DrawCopy("");
      
      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	  if(jtAlgos[aI].find("NoSub") == std::string::npos || cI == periphMostBin){
	    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    configHist(fakeRate_p[aI][cI], aI, doDebug);	
	    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    
	    std::string tempStr = jtAlgos[aI];
	    if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[aI]];
	    leg_p->AddEntry(fakeRate_p[aI][cI], tempStr.c_str(), "P L");
	    fakeRate_p[aI][cI]->Draw("P");
	  }
	}
      }

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
      centLabel.replace(centLabel.find("to"), 2, "-");
      centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
      //      if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";

      
      if(fakeDoAllLeg || fakePermaLeg.count(centBinsStr[cI]) != 0){
	labels.push_back(centLabel);
	labels.push_back("|#eta_{jet}| < " + maxJtAbsEtaStr + "; p_{T,Reco.} > " + minJtPtStr);
	labels.push_back(caloTrackStr + " jets");
      }
      else{
	labels.clear();
	for(unsigned int gI = 0; gI < globalLabels.size(); ++gI){
	  labels.push_back("");
	}
	labels.push_back(centLabel);
      }
    
      for(unsigned int lI = 0; lI < labels.size(); ++lI){
	if(fakeLabelAlignRight){
	  label_p->SetText(fakeLabelX, fakeLabelY, labels[lI].c_str());
	  while(label_p->GetXsize() < maxLabelStr){
	    labels[lI] = " " + labels[lI];
	    label_p->SetText(fakeLabelX, fakeLabelY, labels[lI].c_str());
	  }
	}

	label_p->DrawLatex(fakeLabelX, fakeLabelY - 0.045*lI, labels[lI].c_str());
      }

      gStyle->SetOptStat(0);      

      if(fakeDoAllLeg || fakePermaLeg.count(centBinsStr[cI]) != 0) leg_p->Draw("SAME");

      if(fakeDoLogX) gPad->SetLogx();
      if(fakeDoLogY) gPad->SetLogy();

      for(Int_t eI = 0; eI < 20; ++eI){
	std::string tempStr = inConfig_p->GetValue(("FAKEWHITEBOY." + std::to_string(eI)).c_str(), "");
	if(tempStr.size() == 0) continue;
      
	std::vector<float> pos = strToVectF(tempStr);
	if(pos.size() == 4) drawWhiteBoxNDC(canv_p, pos[0], pos[1], pos[2], pos[3]);
      }

      
      std::string saveName = "pdfDir/" + dateStr + "/fake_" + globalStr2 + "_" + caloTrackStr + "_" + centBinsStr[cI] + "_" + dateStr;
      quietSaveAs(canv_p, saveName + ".pdf");
      //      quietSaveAs(canv_p, saveName + ".png");
      
      delete dummyHist_p;
      delete leg_p;
      delete canv_p;
    }


    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    TH1D* recoOverGenMean_p[nMaxJtAlgo][nMaxCentBins];
    TH1D* recoOverGenSigma_p[nMaxJtAlgo][nMaxCentBins];
    TH1D* recoOverGenSigmaOverMean_p[nMaxJtAlgo][nMaxCentBins];

    TH1D* recoOverGenMMean_p[nMaxJtAlgo][nMaxCentBins];
    TH1D* recoOverGenMSigma_p[nMaxJtAlgo][nMaxCentBins];

    TH1D* recoOverGenMOverPtMean_p[nMaxJtAlgo][nMaxCentBins];
    TH1D* recoOverGenMOverPtSigma_p[nMaxJtAlgo][nMaxCentBins];

    TH1D* recoGenSigma_DeltaEta_p[nMaxJtAlgo][nMaxCentBins];

    TH1D* recoGenSigma_DeltaPhi_p[nMaxJtAlgo][nMaxCentBins];

    Int_t nX = 1;
    Int_t nY = 1;
    if(nJtPtBins == 2) nX = 2;
    else if(nJtPtBins == 3) nX = 3;
    else if(nJtPtBins == 4) nX = 4;
    else if(nJtPtBins == 5 || nJtPtBins == 6){nX = 3; nY = 2;}
    else if(nJtPtBins == 7 || nJtPtBins == 8){nX = 4; nY = 2;} 
    else if(nJtPtBins == 9){nX = 3; nY = 3;} 
    else if(nJtPtBins == 10 || nJtPtBins == 11 || nJtPtBins == 12 ){nX = 4; nY = 3;} 
    else if(nJtPtBins == 13 || nJtPtBins == 14 || nJtPtBins == 15 ){nX = 5; nY = 3;} 
 
  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      std::vector<TH1*> centHistMean, centHistSigma, centHistSigmaOverMean, centHistMMean, centHistMSigma, centHistMOverPtMean, centHistMOverPtSigma, centHistSigmaEta, centHistSigmaPhi;
    for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
      if(jtAlgos[jI].find("ATLAS") != std::string::npos) continue;
      if(jtAlgos[jI].find("Truth") != std::string::npos) continue;


      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(jtAlgos[jI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

     	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	std::string nameStr = jtAlgos[jI] + "_" + centBinsStr[cI];
	recoOverGenMean_p[jI][cI] = new TH1D(("recoOverGenMean_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#LTReco./Gen.#GT", nJtPtBins, jtPtBins);       

	recoOverGenSigma_p[jI][cI] = new TH1D(("recoOverGenSigma_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma(Reco./Gen.)", nJtPtBins, jtPtBins);
	recoOverGenSigmaOverMean_p[jI][cI] = new TH1D(("recoOverGenSigmaOverMean_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma(Reco./Gen.)/#LTReco./Gen.#GT", nJtPtBins, jtPtBins);

	recoOverGenMMean_p[jI][cI] = new TH1D(("recoOverGenMMean_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#LTReco. Mass/Gen. Mass#GT", nJtPtBins, jtPtBins);
	recoOverGenMSigma_p[jI][cI] = new TH1D(("recoOverGenMSigma_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma(Reco. Mass/Gen. Mass)", nJtPtBins, jtPtBins);

	recoOverGenMOverPtMean_p[jI][cI] = new TH1D(("recoOverGenMOverPtMean_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#LT(Reco. M/p_{T})/(Gen. M/p_{T})#GT", nJtPtBins, jtPtBins);
	recoOverGenMOverPtSigma_p[jI][cI] = new TH1D(("recoOverGenMOverPtSigma_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma((Reco. M/p_{T})/(Gen. M/p_{T}))", nJtPtBins, jtPtBins);

	recoGenSigma_DeltaEta_p[jI][cI] = new TH1D(("recoGenSigma_DeltaEta_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma(#eta_{Reco.} - #eta_{Gen.})", nJtPtBins, jtPtBins);
	recoGenSigma_DeltaPhi_p[jI][cI] = new TH1D(("recoGenSigma_DeltaPhi_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma(#phi_{Reco.} - #phi_{Gen.})", nJtPtBins, jtPtBins);
	
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    centerTitles({recoOverGenMean_p[jI][cI], recoOverGenSigma_p[jI][cI], recoOverGenSigmaOverMean_p[jI][cI], recoOverGenMMean_p[jI][cI], recoOverGenMSigma_p[jI][cI], recoOverGenMOverPtMean_p[jI][cI], recoOverGenMOverPtSigma_p[jI][cI], recoGenSigma_DeltaEta_p[jI][cI], recoGenSigma_DeltaPhi_p[jI][cI]});

 	TCanvas* canvFit_p = new TCanvas("canvFit_p", "", 450*nX, 450*nY);
	canvFit_p->SetTopMargin(0.01);
	canvFit_p->SetLeftMargin(0.01);
	canvFit_p->SetBottomMargin(0.01);
	canvFit_p->SetRightMargin(0.01);
 	
	canvFit_p->Divide(nX, nY);
	
	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  centerTitles(recoOverGen_VPt_p[jI][cI][jI2]);

	  //Pre-fit to define our fit range

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  if(isStrSame(caloTrackStr, "Tower")){
	    TF1* tempFit_p = new TF1("tempFit_p", "gaus", 0.0, 2.0);
	    recoOverGen_VPt_p[jI][cI][jI2]->Fit(tempFit_p, "Q N M", "", 0.0, 2.0);
	    Double_t mean = tempFit_p->GetParameter(1);
	    Int_t binPos = recoOverGen_VPt_p[jI][cI][jI2]->FindBin(mean);
	    Int_t iter = 1;
	    double integral = 0.0;
	    Double_t minVal = 0.0;
	    Double_t maxVal = 2.0;
	    
	    Int_t minBin = TMath::Max(1, recoOverGen_VPt_p[jI][cI][jI2]->FindBin(minJtPt/jtPtBins[jI2]));
	    Int_t maxBin = recoOverGen_VPt_p[jI][cI][jI2]->GetNbinsX()+1;
	    
	    while(integral < 0.95){       
	      Int_t lowBin = TMath::Max(minBin, binPos - iter);
	      Int_t highBin = TMath::Min(maxBin, binPos + iter);
	      
	      integral = recoOverGen_VPt_p[jI][cI][jI2]->Integral(lowBin, highBin)/recoOverGen_VPt_p[jI][cI][jI2]->Integral(minBin, maxBin);
	      minVal = (recoOverGen_VPt_p[jI][cI][jI2]->GetBinCenter(lowBin)+recoOverGen_VPt_p[jI][cI][jI2]->GetBinLowEdge(lowBin))/2.;
	      maxVal = (recoOverGen_VPt_p[jI][cI][jI2]->GetBinCenter(highBin)+recoOverGen_VPt_p[jI][cI][jI2]->GetBinLowEdge(highBin+1))/2.;	    
	      
	      ++iter;
	    }
	    
	    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	    delete tempFit_p;
	    
	    recoOverGenFit_p[jI][cI][jI2] = new TF1(("recoOverGenFit_" + nameStr + "_" + jtPtBinsStr[jI]).c_str(), "gaus", minVal, maxVal);
	    recoOverGen_VPt_p[jI][cI][jI2]->Fit(recoOverGenFit_p[jI][cI][jI2], "Q N M", "", minVal, maxVal);
	  }

	  Double_t mean = recoOverGen_VPt_p[jI][cI][jI2]->GetMean();
	  Double_t meanErr = recoOverGen_VPt_p[jI][cI][jI2]->GetMeanError();

	  Double_t sigma = recoOverGen_VPt_p[jI][cI][jI2]->GetStdDev();
	  Double_t sigmaErr = recoOverGen_VPt_p[jI][cI][jI2]->GetStdDevError();

	  Double_t relErr = TMath::Sqrt(meanErr/mean + sigmaErr/sigma);
	  Double_t sigmaOverMean = sigma/mean;
	  Double_t sigmaOverMeanErr = sigmaOverMean*relErr;

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  if(isStrSame(caloTrackStr, "Tower")){
	    mean = recoOverGenFit_p[jI][cI][jI2]->GetParameter(1);
	    meanErr = recoOverGenFit_p[jI][cI][jI2]->GetParError(1);
	    
	    sigma = recoOverGenFit_p[jI][cI][jI2]->GetParameter(2);
	    sigmaErr = recoOverGenFit_p[jI][cI][jI2]->GetParError(2);
	    
	    relErr = TMath::Sqrt(meanErr/mean + sigmaErr/sigma);
	    sigmaOverMean = sigma/mean;
	    sigmaOverMeanErr = sigmaOverMean*relErr;
	  }

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  recoOverGenMean_p[jI][cI]->SetBinContent(jI2+1, mean);
	  recoOverGenMean_p[jI][cI]->SetBinError(jI2+1, meanErr);
	  recoOverGenSigma_p[jI][cI]->SetBinContent(jI2+1, sigma);
	  recoOverGenSigma_p[jI][cI]->SetBinError(jI2+1, sigmaErr);

	  recoOverGenMMean_p[jI][cI]->SetBinContent(jI2+1, recoOverGenM_VPt_p[jI][cI][jI2]->GetMean());
	  recoOverGenMMean_p[jI][cI]->SetBinError(jI2+1, recoOverGenM_VPt_p[jI][cI][jI2]->GetMeanError());
	  recoOverGenMSigma_p[jI][cI]->SetBinContent(jI2+1, recoOverGenM_VPt_p[jI][cI][jI2]->GetStdDev());
	  recoOverGenMSigma_p[jI][cI]->SetBinError(jI2+1, recoOverGenM_VPt_p[jI][cI][jI2]->GetStdDevError());

	  recoOverGenMOverPtMean_p[jI][cI]->SetBinContent(jI2+1, recoOverGenMOverPt_VPt_p[jI][cI][jI2]->GetMean());
	  recoOverGenMOverPtMean_p[jI][cI]->SetBinError(jI2+1, recoOverGenMOverPt_VPt_p[jI][cI][jI2]->GetMeanError());
	  recoOverGenMOverPtSigma_p[jI][cI]->SetBinContent(jI2+1, recoOverGenMOverPt_VPt_p[jI][cI][jI2]->GetStdDev());
	  recoOverGenMOverPtSigma_p[jI][cI]->SetBinError(jI2+1, recoOverGenMOverPt_VPt_p[jI][cI][jI2]->GetStdDevError());


	  recoGenSigma_DeltaEta_p[jI][cI]->SetBinContent(jI2+1, recoGen_DeltaEta_p[jI][cI][jI2]->GetStdDev());
	  recoGenSigma_DeltaEta_p[jI][cI]->SetBinError(jI2+1, recoGen_DeltaEta_p[jI][cI][jI2]->GetStdDevError());

	  recoGenSigma_DeltaPhi_p[jI][cI]->SetBinContent(jI2+1, recoGen_DeltaPhi_p[jI][cI][jI2]->GetStdDev());
	  recoGenSigma_DeltaPhi_p[jI][cI]->SetBinError(jI2+1, recoGen_DeltaPhi_p[jI][cI][jI2]->GetStdDevError());

	  recoOverGenSigmaOverMean_p[jI][cI]->SetBinContent(jI2+1, sigmaOverMean);
	  recoOverGenSigmaOverMean_p[jI][cI]->SetBinError(jI2+1, sigmaOverMeanErr);

	  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  canvFit_p->cd();
	  canvFit_p->cd(jI2+1);

	  gPad->SetRightMargin(0.01);
	  gPad->SetTopMargin(0.01);
	  gPad->SetLeftMargin(0.14);
	  gPad->SetBottomMargin(0.12);
	  
	  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  recoOverGen_VPt_p[jI][cI][jI2]->GetXaxis()->SetTitleOffset(2.0);
	  recoOverGen_VPt_p[jI][cI][jI2]->GetYaxis()->SetTitleOffset(3.0);

	  //	  recoOverGen_VPt_p[jI][cI][jI2]->GetXaxis()->SetTitleSize(20);
	  //	  recoOverGen_VPt_p[jI][cI][jI2]->GetYaxis()->SetTitleSize(20);

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  recoOverGen_VPt_p[jI][cI][jI2]->DrawCopy("HIST E1");
	  if(isStrSame(caloTrackStr, "Tower")) recoOverGenFit_p[jI][cI][jI2]->DrawCopy("SAME");	  
	
	  std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	  centLabel.replace(centLabel.find("to"), 2, "-");
	  centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
	  if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";
	
	  std::string tempStr = jtAlgos[jI];
	  if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[jI]];
	  centLabel = tempStr + ", " + centLabel;
	  
	  label_p->DrawLatex(0.20, 0.95, centLabel.c_str());
	  std::string ptStr = prettyString(jtPtBins[jI2], 1, false) + "<p_{T,Jet}<" + prettyString(jtPtBins[jI2+1], 1, false);
	  label_p->DrawLatex(0.20, 0.83, ptStr.c_str());

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.20, 0.77, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
	  if(minJtPtStr.size() != 0) label_p->DrawLatex(0.20, 0.71, ("p_{T,Reco} > " + minJtPtStr).c_str());
	  if(isStrSame(caloTrackStr, "Tower")) label_p->DrawLatex(0.20, 0.89, "Calo. jets");
	  else if(isStrSame(caloTrackStr, "Trk")) label_p->DrawLatex(0.20, 0.89, "Track jets");
	}

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
	configHist(recoOverGenMean_p[jI][cI], cI, doDebug);
	configHist(recoOverGenSigma_p[jI][cI], cI, doDebug);
	configHist(recoOverGenSigmaOverMean_p[jI][cI], cI, doDebug);
	configHist(recoOverGenMMean_p[jI][cI], cI, doDebug);
	configHist(recoOverGenMSigma_p[jI][cI], cI, doDebug);
	configHist(recoOverGenMOverPtMean_p[jI][cI], cI, doDebug);
	configHist(recoOverGenMOverPtSigma_p[jI][cI], cI, doDebug);
	configHist(recoGenSigma_DeltaEta_p[jI][cI], cI, doDebug);
	configHist(recoGenSigma_DeltaPhi_p[jI][cI], cI, doDebug);

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
		
	centHistMean.push_back(recoOverGenMean_p[jI][cI]);
	centHistSigma.push_back(recoOverGenSigma_p[jI][cI]);
	centHistSigmaOverMean.push_back(recoOverGenSigmaOverMean_p[jI][cI]);

	centHistMMean.push_back(recoOverGenMMean_p[jI][cI]);
	centHistMSigma.push_back(recoOverGenMSigma_p[jI][cI]);

	centHistMOverPtMean.push_back(recoOverGenMOverPtMean_p[jI][cI]);
	centHistMOverPtSigma.push_back(recoOverGenMOverPtSigma_p[jI][cI]);

	centHistSigmaEta.push_back(recoGenSigma_DeltaEta_p[jI][cI]);
	centHistSigmaPhi.push_back(recoGenSigma_DeltaPhi_p[jI][cI]);

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	std::string saveName = "pdfDir/" + dateStr + "/recoOverGen_" + nameStr + "_" + dateStr + ".pdf";
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << saveName << ", " << canvFit_p << std::endl;


	quietSaveAs(canvFit_p, saveName);
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	delete canvFit_p;
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      }


	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

   

      for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
	if(jtAlgos[jI].find("NoSub") != std::string::npos && cI != (unsigned int)periphMostBin) continue;
	std::string nameStr = jtAlgos[jI] + "_" + centBinsStr[cI];
	TCanvas* canvFit_p = new TCanvas("canvFit_p", "", 450*nX, 450*nY);
	canvFit_p->SetTopMargin(0.01);
	canvFit_p->SetLeftMargin(0.01);
	canvFit_p->SetBottomMargin(0.01);
	canvFit_p->SetRightMargin(0.01);
	
	canvFit_p->Divide(nX, nY);

	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  centerTitles(recoGen_DeltaEta_p[jI][cI][jI2]);

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  canvFit_p->cd();
	  canvFit_p->cd(jI2+1);

	  gPad->SetRightMargin(0.01);
	  gPad->SetTopMargin(0.01);
	  gPad->SetLeftMargin(0.14);
	  gPad->SetBottomMargin(0.12);

	  recoGen_DeltaEta_p[jI][cI][jI2]->GetXaxis()->SetTitleOffset(2.0);
	  recoGen_DeltaEta_p[jI][cI][jI2]->GetYaxis()->SetTitleOffset(3.0);

	  //	  recoGen_DeltaEta_p[jI][cI][jI2]->GetXaxis()->SetTitleSize(20);
	  //	  recoGen_DeltaEta_p[jI][cI][jI2]->GetYaxis()->SetTitleSize(20);

 	  recoGen_DeltaEta_p[jI][cI][jI2]->SetMinimum(getMinGTZero(recoGen_DeltaEta_p[jI][cI][jI2])/2.);
	  recoGen_DeltaEta_p[jI][cI][jI2]->DrawCopy("HIST E1");
	  gPad->SetLogy();

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	  centLabel.replace(centLabel.find("to"), 2, "-");
	  centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
	  if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";
	
	  std::string tempStr = jtAlgos[jI];
	  if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[jI]];
	  centLabel = tempStr + ", " + centLabel;
	  
	  label_p->DrawLatex(0.20, 0.95, centLabel.c_str());
	  std::string ptStr = prettyString(jtPtBins[jI2], 1, false) + "<p_{T,Jet}<" + prettyString(jtPtBins[jI2+1], 1, false);
	  label_p->DrawLatex(0.20, 0.83, ptStr.c_str());
	  
	  if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.20, 0.77, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
	  if(minJtPtStr.size() != 0) label_p->DrawLatex(0.20, 0.71, ("p_{T,Reco} > " + minJtPtStr).c_str());
	  if(isStrSame(caloTrackStr, "Tower")) label_p->DrawLatex(0.20, 0.89, "Calo. jets");
	  else if(isStrSame(caloTrackStr, "Trk")) label_p->DrawLatex(0.20, 0.89, "Track jets");

	  label_p->DrawLatex(0.20, 0.65, ("#sigma=" + prettyString(recoGen_DeltaEta_p[jI][cI][jI2]->GetStdDev(), 3, false)).c_str());
	}

	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


	std::string saveName = "pdfDir/" + dateStr + "/recoGen_DeltaEta_" + nameStr + "_" + dateStr + ".pdf";
	quietSaveAs(canvFit_p, saveName);
	delete canvFit_p;
      }

      for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
	if(jtAlgos[jI].find("NoSub") != std::string::npos && cI != (unsigned int)periphMostBin) continue;
	std::string nameStr = jtAlgos[jI] + "_" + centBinsStr[cI];

	TCanvas* canvFit_p = new TCanvas("canvFit_p", "", 450*nX, 450*nY);
	canvFit_p->SetTopMargin(0.01);
	canvFit_p->SetLeftMargin(0.01);
	canvFit_p->SetBottomMargin(0.01);
	canvFit_p->SetRightMargin(0.01);
	
	canvFit_p->Divide(nX, nY);

	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	  centerTitles(recoGen_DeltaPhi_p[jI][cI][jI2]);
	  
	  canvFit_p->cd();
	  canvFit_p->cd(jI2+1);

	  gPad->SetRightMargin(0.01);
	  gPad->SetTopMargin(0.01);
	  gPad->SetLeftMargin(0.14);
	  gPad->SetBottomMargin(0.12);
	  
	  recoGen_DeltaPhi_p[jI][cI][jI2]->GetXaxis()->SetTitleOffset(2.0);
	  recoGen_DeltaPhi_p[jI][cI][jI2]->GetYaxis()->SetTitleOffset(3.0);

	  //	  recoGen_DeltaPhi_p[jI][cI][jI2]->GetXaxis()->SetTitleSize(20);
	  //	  recoGen_DeltaPhi_p[jI][cI][jI2]->GetYaxis()->SetTitleSize(20);

 	  recoGen_DeltaPhi_p[jI][cI][jI2]->SetMinimum(getMinGTZero(recoGen_DeltaPhi_p[jI][cI][jI2])/2.);
	  recoGen_DeltaPhi_p[jI][cI][jI2]->DrawCopy("HIST E1");
	
	  gPad->SetLogy();

	  std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	  centLabel.replace(centLabel.find("to"), 2, "-");
	  centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
	  if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";
	
	  std::string tempStr = jtAlgos[jI];
	  if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[jI]];
	  centLabel = tempStr + ", " + centLabel;
	  
	  label_p->DrawLatex(0.20, 0.95, centLabel.c_str());
	  std::string ptStr = prettyString(jtPtBins[jI2], 1, false) + "<p_{T,Jet}<" + prettyString(jtPtBins[jI2+1], 1, false);
	  label_p->DrawLatex(0.20, 0.83, ptStr.c_str());
	  
	  if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.20, 0.77, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
	  if(minJtPtStr.size() != 0) label_p->DrawLatex(0.20, 0.71, ("p_{T,Reco} > " + minJtPtStr).c_str());
	  if(isStrSame(caloTrackStr, "Tower")) label_p->DrawLatex(0.20, 0.89, "Calo. jets");
	  else if(isStrSame(caloTrackStr, "Trk")) label_p->DrawLatex(0.20, 0.89, "Track jets");

	  label_p->DrawLatex(0.20, 0.65, ("#sigma=" + prettyString(recoGen_DeltaPhi_p[jI][cI][jI2]->GetStdDev(), 3, false)).c_str());
	}


	std::string saveName = "pdfDir/" + dateStr + "/recoGen_DeltaPhi_" + nameStr + "_" + dateStr + ".pdf";
	quietSaveAs(canvFit_p, saveName);
	delete canvFit_p;
      }

      plotResponseSet(fileConfig_p, inConfig_p, centHistMean, centBinsStr, "recoOverGenMean_" + globalStr + "_" + caloTrackStr, jtAlgos[jI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, centHistSigma, centBinsStr, "recoOverGenSigma_" + globalStr + "_" + caloTrackStr, jtAlgos[jI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, centHistSigmaOverMean, centBinsStr, "recoOverGenSigmaOverMean_" + globalStr + "_" + caloTrackStr, jtAlgos[jI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, centHistMMean, centBinsStr, "recoOverGenMMean_" + globalStr + "_" + caloTrackStr, jtAlgos[jI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, centHistMSigma, centBinsStr, "recoOverGenMSigma_" + globalStr + "_" + caloTrackStr, jtAlgos[jI], dateStr);

      plotResponseSet(fileConfig_p, inConfig_p, centHistMOverPtMean, centBinsStr, "recoOverGenMOverPtMean_" + globalStr + "_" + caloTrackStr, jtAlgos[jI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, centHistMOverPtSigma, centBinsStr, "recoOverGenMOverPtSigma_" + globalStr + "_" + caloTrackStr, jtAlgos[jI], dateStr);

      plotResponseSet(fileConfig_p, inConfig_p, centHistSigmaEta, centBinsStr, "recoGenSigma_DeltaEta_" + globalStr + "_" + caloTrackStr, jtAlgos[jI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, centHistSigmaPhi, centBinsStr, "recoGenSigma_DeltaPhi_" + globalStr + "_" + caloTrackStr, jtAlgos[jI], dateStr);

      centHistMean.clear();
      centHistSigma.clear();
      centHistSigmaOverMean.clear();
      centHistMMean.clear();
      centHistMSigma.clear();
      centHistMOverPtMean.clear();
      centHistMOverPtSigma.clear();
      centHistSigmaEta.clear();
      centHistSigmaPhi.clear();      
    }
  
    std::vector<TH1*> algoHistMean, algoHistSigma, algoHistMMean, algoHistMSigma, algoHistMOverPtMean, algoHistMOverPtSigma, algoHistSigmaOverMean, algoHistSigmaEta, algoHistSigmaPhi;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::vector<std::string> jtAlgosLabel;
      
      for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
	if(jtAlgos[jI].find("ATLAS") != std::string::npos) continue;
	if(jtAlgos[jI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;
	if(jtAlgos[jI].find("Truth") != std::string::npos) continue;

 	configHist(recoOverGenMean_p[jI][cI], jI, doDebug);
	configHist(recoOverGenSigma_p[jI][cI], jI, doDebug);
	configHist(recoOverGenSigmaOverMean_p[jI][cI], jI, doDebug);
 	configHist(recoOverGenMMean_p[jI][cI], jI, doDebug);
	configHist(recoOverGenMSigma_p[jI][cI], jI, doDebug);
 	configHist(recoOverGenMOverPtMean_p[jI][cI], jI, doDebug);
	configHist(recoOverGenMOverPtSigma_p[jI][cI], jI, doDebug);
	configHist(recoGenSigma_DeltaEta_p[jI][cI], jI, doDebug);
	configHist(recoGenSigma_DeltaPhi_p[jI][cI], jI, doDebug);
	
	algoHistMean.push_back(recoOverGenMean_p[jI][cI]);
	algoHistSigma.push_back(recoOverGenSigma_p[jI][cI]);
	algoHistSigmaOverMean.push_back(recoOverGenSigmaOverMean_p[jI][cI]);
	algoHistMMean.push_back(recoOverGenMMean_p[jI][cI]);
	algoHistMSigma.push_back(recoOverGenMSigma_p[jI][cI]);
	algoHistMOverPtMean.push_back(recoOverGenMOverPtMean_p[jI][cI]);
	algoHistMOverPtSigma.push_back(recoOverGenMOverPtSigma_p[jI][cI]);
	algoHistSigmaEta.push_back(recoGenSigma_DeltaEta_p[jI][cI]);
	algoHistSigmaPhi.push_back(recoGenSigma_DeltaPhi_p[jI][cI]);

	std::string tempStr = jtAlgos[jI];
	if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[jI]];
	jtAlgosLabel.push_back(tempStr);
      }

      plotResponseSet(fileConfig_p, inConfig_p, algoHistMean, jtAlgosLabel, "recoOverGenMean_" + globalStr + "_" + caloTrackStr, centBinsStr[cI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, algoHistSigma, jtAlgosLabel, "recoOverGenSigma_" + globalStr + "_" + caloTrackStr, centBinsStr[cI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, algoHistSigmaOverMean, jtAlgosLabel, "recoOverGenSigmaOverMean_" + globalStr + "_" + caloTrackStr, centBinsStr[cI], dateStr);

      plotResponseSet(fileConfig_p, inConfig_p, algoHistMMean, jtAlgosLabel, "recoOverGenMMean_" + globalStr + "_" + caloTrackStr, centBinsStr[cI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, algoHistMSigma, jtAlgosLabel, "recoOverGenMSigma_" + globalStr + "_" + caloTrackStr, centBinsStr[cI], dateStr);

      plotResponseSet(fileConfig_p, inConfig_p, algoHistMOverPtMean, jtAlgosLabel, "recoOverGenMOverPtMean_" + globalStr + "_" + caloTrackStr, centBinsStr[cI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, algoHistMOverPtSigma, jtAlgosLabel, "recoOverGenMOverPtSigma_" + globalStr + "_" + caloTrackStr, centBinsStr[cI], dateStr);

      plotResponseSet(fileConfig_p, inConfig_p, algoHistSigmaEta, jtAlgosLabel, "recoGenSigma_DeltaEta_" + globalStr + "_" + caloTrackStr, centBinsStr[cI], dateStr);
      plotResponseSet(fileConfig_p, inConfig_p, algoHistSigmaPhi, jtAlgosLabel, "recoGenSigma_DeltaPhi_" + globalStr + "_" + caloTrackStr, centBinsStr[cI], dateStr);

      algoHistMean.clear();
      algoHistSigma.clear();
      algoHistSigmaOverMean.clear();
      algoHistMMean.clear();
      algoHistMSigma.clear();
      algoHistMOverPtMean.clear();
      algoHistMOverPtSigma.clear();
      algoHistSigmaEta.clear();
      algoHistSigmaPhi.clear();      
    }

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    //Special summary plot

    const double summMax = inConfig_p->GetValue("RECOOVERGENSUMMMAX", 60.0);
    const double summMin = inConfig_p->GetValue("RECOOVERGENSUMMMIN", 0.0);
    const double summXMax = inConfig_p->GetValue("RECOOVERGENSUMMXMAX", 60.0);
    const double summXMin = inConfig_p->GetValue("RECOOVERGENSUMMXMIN", -60.0);
    const double summLegX = inConfig_p->GetValue("RECOOVERGENSUMMLEGX", 0.2);
    const double summLegY = inConfig_p->GetValue("RECOOVERGENSUMMLEGY", 0.3);
    const double summLabelX = inConfig_p->GetValue("RECOOVERGENSUMMLABELX", 0.2);
    const double summLabelY = inConfig_p->GetValue("RECOOVERGENSUMMLABELY", 0.9);
    const bool summLabelAlignRight = inConfig_p->GetValue("RECOOVERGENSUMMALIGNRIGHT", 0.9);
    
    const Double_t leftMarg = 0.1;
    const Double_t bottomMarg = 0.13;
    const Double_t topMarg = 0.07;
    const Double_t rightMarg = 0.01;
    TCanvas* summCanv_p = new TCanvas("summCanv_p", "", 450*2, 450);
    summCanv_p->SetTopMargin(0.001);
    summCanv_p->SetBottomMargin(0.001);
    summCanv_p->SetLeftMargin(0.001);
    summCanv_p->SetRightMargin(0.001);

    const Int_t nPt = 4;
    TPad* pads_p[nPt];
    const Double_t padWidths = (1.0 - leftMarg - rightMarg)/(Double_t)nPt;
    for(Int_t pI = 0; pI < nPt; ++pI){
      double x1 = 0.0;
      if(pI != 0) x1 = leftMarg + pI*padWidths;

      std::cout << "PADPOS: " << x1 << ", " << 0.0 << ", " << leftMarg + (pI+1)*padWidths << ", " << 1.0 << std::endl;

      summCanv_p->cd();
      pads_p[pI] = new TPad("pad0", "", x1, 0.0, leftMarg + (pI+1)*padWidths, 1.0);
      pads_p[pI]->SetBottomMargin(bottomMarg);
      pads_p[pI]->SetTopMargin(topMarg);
      if(pI == 0) pads_p[pI]->SetLeftMargin(leftMarg*(1.0 - leftMarg  - padWidths)/(padWidths));
      else pads_p[pI]->SetLeftMargin(0.0001);

      if(pI == nPt-1) pads_p[pI]->SetRightMargin(rightMarg);
      else pads_p[pI]->SetRightMargin(0.0001);

      summCanv_p->cd();
      pads_p[pI]->Draw("SAME");
    }
  
    TGraph* summaryGraph_p[nPt][nMaxJtAlgo];
    for(Int_t pI = 0; pI < nPt; ++pI){
      summCanv_p->cd();
      pads_p[pI]->cd();

      double renorm = padWidths;
      if(pI == 0) renorm += leftMarg;
      TH1F* dummyHist_p = new TH1F("dummyHist_h", ";#LT#Deltap_{T}#GT [GeV];#sigma_{#Deltap_{T}} [GeV]", 10, summXMin, summXMax);
      dummyHist_p->SetMaximum(summMax);
      dummyHist_p->SetMinimum(summMin);
      centerTitles(dummyHist_p);

      dummyHist_p->GetXaxis()->SetTitleSize(padWidths*globalSize*2.6/renorm);
      dummyHist_p->GetYaxis()->SetTitleSize(padWidths*globalSize*2.6/renorm);
      dummyHist_p->GetXaxis()->SetLabelSize(padWidths*globalSize*2.0/renorm);
      dummyHist_p->GetYaxis()->SetLabelSize(padWidths*globalSize*2.0/renorm);

      if(pI == 0){
	dummyHist_p->GetXaxis()->SetLabelOffset(-0.005);
	dummyHist_p->GetXaxis()->SetTitleOffset(0.5*renorm/padWidths);
      }
      else{
	dummyHist_p->GetXaxis()->SetLabelOffset(-0.0225);
	dummyHist_p->GetXaxis()->SetTitleOffset(0.5);
      }

      dummyHist_p->GetXaxis()->SetNdivisions(505);      
      dummyHist_p->DrawCopy("HIST E1 P");

      const Int_t ptPos = TMath::Min(pI*nJtPtBins/nPt, nJtPtBins-1);
      const Double_t ptVal = (jtPtBins[ptPos] + jtPtBins[ptPos+1])/2.;

      for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	if(jtAlgos[aI].find("NoSub") != std::string::npos) summaryGraph_p[pI][aI] = new TGraph(1);
	else summaryGraph_p[pI][aI] = new TGraph(nCentBins);

	if(doDebug) std::cout << pI << ", " << aI << std::endl;
	configHist(summaryGraph_p[pI][aI], aI, doDebug);

	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << ptPos << std::endl;
	
	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << ptPos << std::endl;

	  if(jtAlgos[aI].find("NoSub") != std::string::npos){
	    if(cI == periphMostBin){
	      std::cout << "HIT SUMMARY aI, cI: " << aI << ", " << cI << ", " << summaryGraph_p[pI][aI]->GetN() << std::endl;
	      summaryGraph_p[pI][aI]->SetPoint(0, ptVal*recoOverGenMean_p[aI][cI]->GetBinContent(ptPos+1) - ptVal, ptVal*recoOverGenSigma_p[aI][cI]->GetBinContent(ptPos+1));
	    }
	  }
	  else summaryGraph_p[pI][aI]->SetPoint(cI, ptVal*recoOverGenMean_p[aI][cI]->GetBinContent(ptPos+1) - ptVal, ptVal*recoOverGenSigma_p[aI][cI]->GetBinContent(ptPos+1));

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << ", " << ptPos << std::endl;

	}

	std::cout << "PI, AI: " << pI << ", " << aI << ", " <<  summaryGraph_p[pI][aI]->GetN() << std::endl;
	summaryGraph_p[pI][aI]->Draw("P");
      }

      delete dummyHist_p;
    }

    summCanv_p->cd();
    for(Int_t pI = 0; pI < nPt; ++pI){
      double renorm = padWidths;
      if(pI == 0) renorm += leftMarg;
      const Int_t ptPos = TMath::Min(pI*nJtPtBins/nPt, nJtPtBins-1);
      
      label_p->DrawLatex(leftMarg + padWidths/6. + pI*padWidths, .955, (prettyString(jtPtBins[ptPos], 1, false) + " < p_{T,Gen.} < " + prettyString(jtPtBins[ptPos+1], 1, false)).c_str());
    }

    maxLabelStr = 0.0;
    if(summLabelAlignRight){
      std::vector<std::string> labels = globalLabels;
      labels.push_back("|#eta_{Jet}| < " + maxJtAbsEtaStr + "; p_{T,Reco.} > " + minJtPtStr);
      labels.push_back(caloTrackStr + " jets");

      for(unsigned int lI = 0; lI < labels.size(); ++lI){       
	label_p->SetText(summLabelX, summLabelY, labels[lI].c_str());
	if(label_p->GetXsize() > maxLabelStr) maxLabelStr = label_p->GetXsize();
      }      
    }
    
    std::vector<std::string> labels = globalLabels;
    labels.push_back("|#eta_{Jet}| < " + maxJtAbsEtaStr + "; p_{T,Reco.} > " + minJtPtStr);
    labels.push_back(caloTrackStr + " jets");

    for(unsigned int lI = 0; lI < labels.size(); ++lI){
      if(summLabelAlignRight){
	label_p->SetText(summLabelX, summLabelY, labels[lI].c_str());
	while(label_p->GetXsize() < maxLabelStr){
	  labels[lI] = " " + labels[lI];
	  label_p->SetText(summLabelX, summLabelY, labels[lI].c_str());
	}
      }
      
      label_p->DrawLatex(summLabelX, summLabelY - 0.045*lI, labels[lI].c_str());
    }


    double maxLegStr = 0.0;
    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	std::string tempStr = jtAlgos[aI];
	if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[aI]];

	label_p->SetText(0.1, 0.1, tempStr.c_str());
	if(label_p->GetXsize() > maxLegStr) maxLegStr = label_p->GetXsize();
      }	
    }

    TLegend* leg_p = new TLegend(summLegX, summLegY - 0.04*nJtAlgo, summLegX + maxLegStr*0.8, summLegY);
    leg_p->SetTextFont(globalFont);
    leg_p->SetTextSize(globalSize);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);

    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      std::string tempStr = jtAlgos[aI];
      if(algoNameSwaps.count(tempStr) != 0) tempStr = algoNameSwaps[jtAlgos[aI]];
      leg_p->AddEntry(summaryGraph_p[0][aI], tempStr.c_str(), "P");
    }
    
    leg_p->Draw("SAME");
    
    quietSaveAs(summCanv_p, "pdfDir/" + dateStr + "/recoOverGen_Summary_" + globalStr + "_" + caloTrackStr + "_" + dateStr + ".pdf");

    for(Int_t pI = 0; pI < nPt; ++pI){
      delete pads_p[pI];
    }

    delete summCanv_p;

    for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
      if(jtAlgos[jI].find("ATLAS") != std::string::npos) continue;
      if(jtAlgos[jI].find("Truth") != std::string::npos) continue;
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(jtAlgos[jI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;

	delete recoOverGenMean_p[jI][cI];
	delete recoOverGenSigma_p[jI][cI];
	delete recoOverGenSigmaOverMean_p[jI][cI];

	delete recoOverGenMMean_p[jI][cI];
	delete recoOverGenMSigma_p[jI][cI];

	delete recoOverGenMOverPtMean_p[jI][cI];
	delete recoOverGenMOverPtSigma_p[jI][cI];

	delete recoGenSigma_DeltaEta_p[jI][cI];
	delete recoGenSigma_DeltaPhi_p[jI][cI];

	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	  if(false) delete recoOverGenFit_p[jI][cI][jI2];
	}
      }
    }
  }

  for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
    if(jtAlgos[jI].find("ATLAS") == std::string::npos){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(jtAlgos[jI].find("NoSub") == std::string::npos || cI == periphMostBin){
	  delete truthEff_p[jI][cI];
	}
      }
    }
  }
  
  delete label_p;
  delete line_p;
  
  inFile_p->Close();
  delete inFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/plotClusterHist.exe <inConfigFileName>" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;    
    return 1;
  }
  
  int retVal = 0;
  if(argc == 2) retVal += plotClusterHist(argv[1]);
  return retVal;
}
