//cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TF1.h"
#include "TFile.h"
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

const Int_t nColors = 5;
const Int_t colors[nColors] = {0, 1, 3, 4, 6};

const Int_t nStyles = 4;
const Int_t styles[nStyles] = {24, 25, 28, 46};

void configHist(TH1* inHist_p, Int_t pos)
{
  kirchnerPalette kPal;
  inHist_p->SetMarkerSize(1);
  inHist_p->SetMarkerStyle(styles[pos%nStyles]);
  inHist_p->SetMarkerColor(kPal.getColor(colors[pos%nColors]));
  inHist_p->SetLineColor(kPal.getColor(colors[pos%nColors]));
  
  inHist_p->GetXaxis()->SetTitleFont(43);
  inHist_p->GetYaxis()->SetTitleFont(43);
  inHist_p->GetXaxis()->SetLabelFont(43);
  inHist_p->GetYaxis()->SetLabelFont(43);
  
  inHist_p->GetXaxis()->SetTitleSize(15);
  inHist_p->GetYaxis()->SetTitleSize(15);
  inHist_p->GetXaxis()->SetLabelSize(13);
  inHist_p->GetYaxis()->SetLabelSize(13);
  
  inHist_p->GetYaxis()->SetTitleOffset(1.9);
  
  return;
}

void configHist(TGraph* inHist_p, Int_t pos)
{
  kirchnerPalette kPal;
  inHist_p->SetMarkerSize(1);
  inHist_p->SetMarkerStyle(styles[pos%nStyles]);
  inHist_p->SetMarkerColor(kPal.getColor(colors[pos%nColors]));
  inHist_p->SetLineColor(kPal.getColor(colors[pos%nColors]));
  
  inHist_p->GetXaxis()->SetTitleFont(43);
  inHist_p->GetYaxis()->SetTitleFont(43);
  inHist_p->GetXaxis()->SetLabelFont(43);
  inHist_p->GetYaxis()->SetLabelFont(43);
  
  inHist_p->GetXaxis()->SetTitleSize(15);
  inHist_p->GetYaxis()->SetTitleSize(15);
  inHist_p->GetXaxis()->SetLabelSize(13);
  inHist_p->GetYaxis()->SetLabelSize(13);
  
  inHist_p->GetYaxis()->SetTitleOffset(1.9);
  
  return;
}

void plotResponseSet(TEnv* fileConfig_p, std::vector<TH1*> histSet, std::vector<std::string> setLabels, std::string globalNameStr, std::string globalLabelStr, std::string dateStr, double globalMinOverride = 100000, double globalMaxOverride = -100000)
{
  std::string caloTrackStr = fileConfig_p->GetValue("CALOTRACKSTR", "");

  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.16);
  canv_p->SetBottomMargin(0.14);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);
  
  TLegend* leg_p = new TLegend(0.6, 0.75, 0.95, 0.98);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(16);
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

  double globalMin = globalMinOverride;
  double globalMax = globalMaxOverride;
  if(globalMax < globalMin){
    for(unsigned int hI = 0; hI < histSet.size(); ++hI){
      if(getMin(histSet[hI]) < globalMin) globalMin = getMin(histSet[hI]);
      if(getMax(histSet[hI]) > globalMax) globalMax = getMax(histSet[hI]);
    }

    double interval = globalMax - globalMin;
    globalMin -= interval/10.;
    globalMax += interval/10.;
  }

    
  for(unsigned int hI = 0; hI < histSet.size(); ++hI){
    histSet[hI]->SetMaximum(globalMax);
    histSet[hI]->SetMinimum(globalMin);

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

  leg_p->Draw("SAME");

  std::string globalLabelStr2 = globalLabelStr;
  if(globalLabelStr2.find("Cent") != std::string::npos){
    globalLabelStr2.replace(globalLabelStr2.find("Cent"), 4, "");
    globalLabelStr2.replace(globalLabelStr2.find("to"), 2, "-");
    globalLabelStr2 = globalLabelStr2 + "%";
  }
  label_p->DrawLatex(0.26, 0.93, globalLabelStr2.c_str());
  if(caloTrackStr.find("Tower") != std::string::npos) label_p->DrawLatex(0.26, 0.87, "Calo. jets");
  else label_p->DrawLatex(0.26, 0.87, "Track jets");

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

  label_p->DrawLatex(0.26, 0.81, ("|#eta_{Jet}| < " + maxJtAbsEtaStr).c_str());
  label_p->DrawLatex(0.26, 0.75, ("p_{T,Reco.} > " + minJtPtStr).c_str());

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

  if(!checkConfigContainsParams(fileConfig_p, fileParams)) return 1;

  const Bool_t isMC = fileConfig_p->GetValue("ISMC", 0);

  fileConfig_p->SetValue("CALOTRACKSTR", caloTrackStr.c_str());
  
  const Int_t nJtPtBins = fileConfig_p->GetValue("NJTPTBINS", 0);
  std::vector<std::string> jtPtBinsStr = commaSepStringToVect(fileConfig_p->GetValue("JTPTBINS", ""));

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
  const std::string minJtPtStr = prettyString(minJtPt, 1, false);
  std::vector<Double_t> nEventPerCent;
  for(unsigned int cI = 0; cI < nEventPerCentStr.size(); ++cI){
    nEventPerCent.push_back(std::stod(nEventPerCentStr[cI]));
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
    if(jtAlgos[pos].find(caloTrackStr) != std::string::npos) ++pos;
    else jtAlgos.erase(jtAlgos.begin() + pos);
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
  TH1D* spectraTruth_p[nMaxCentBins];

  TH1D* recoOverGen_VPt_p[nMaxJtAlgo][nCentBins][nMaxJtPtBins];
  TF1* recoOverGenFit_p[nMaxJtAlgo][nCentBins][nMaxJtPtBins];
  TH1D* recoGen_DeltaEta_p[nMaxJtAlgo][nCentBins][nMaxJtPtBins];
  TH1D* recoGen_DeltaPhi_p[nMaxJtAlgo][nCentBins][nMaxJtPtBins];
  TH1D* matchedTruthSpectra_p[nMaxJtAlgo][nMaxCentBins];
  TGraphAsymmErrors* truthEff_p[nMaxJtAlgo][nMaxCentBins];

  for(Int_t aI = 0; aI < nJtAlgo; ++ aI){
    for(Int_t cI = 0; cI < nCentBins; ++ cI){
      std::string nameStr = jtAlgos[aI] + "_" + centBinsStr[cI];
      
      spectra_p[aI][cI] = (TH1D*)inFile_p->Get(("spectra_" + nameStr + "_h").c_str());   
      
      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
      if(isMC){
	if(aI == 0) spectraTruth_p[cI] = (TH1D*)inFile_p->Get(("spectra_Truth_" + centBinsStr[cI] + "_h").c_str());

	if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	  matchedTruthSpectra_p[aI][cI] = (TH1D*)inFile_p->Get(("matchedTruthSpectra_" + nameStr + "_h").c_str());       
	  truthEff_p[aI][cI] = new TGraphAsymmErrors();
	  truthEff_p[aI][cI]->BayesDivide(matchedTruthSpectra_p[aI][cI], spectraTruth_p[cI]);
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

      configHist(spectra_p[aI][cI], aI);
      spectra_p[aI][cI]->GetYaxis()->SetTitle("#frac{1}{N_{Event}} #frac{dN_{Jet}}{dp_{T}} [Gev^{-1}]");
      centerTitles(spectra_p[aI][cI]);
      std::cout << "HIST: " << spectra_p[aI][cI]->GetName() << std::endl;
    
      if(isMC){
	if(jtAlgos[aI].find("ATLAS") != std::string::npos) continue;
	if(jtAlgos[aI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;
	if(jtAlgos[aI].find("Truth") != std::string::npos) continue;
	
	for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	  recoOverGen_VPt_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoOverGen_VPt_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());
	  recoGen_DeltaEta_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoGen_DeltaEta_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());
	  recoGen_DeltaPhi_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoGen_DeltaPhi_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());

	  std::cout << "NAME CHECK: " << recoGen_DeltaEta_p[aI][cI][jI]->GetName() << std::endl;

	  configHist(recoOverGen_VPt_p[aI][cI][jI], aI);
	  configHist(recoGen_DeltaEta_p[aI][cI][jI], aI);
	  configHist(recoGen_DeltaPhi_p[aI][cI][jI], aI);
	}
      }
    }
  }

  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);

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
    
    std::cout << "GLOBALMIN, GLOBALMAX: " << globalMin << ", " << globalMax << std::endl;

    bool isDrawn = false;
    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      if(caloTrackStr.find("Trk") != std::string::npos && jtAlgos[aI].find("ATLAS") != std::string::npos) continue;
      if(jtAlgos[aI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;

      leg_p->AddEntry(spectra_p[aI][cI], jtAlgos[aI].c_str(), "P L");
		      
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
    
    label_p->DrawLatex(0.264, 0.95, centLabel.c_str());
    if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.264, 0.83, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
    if(isStrSame(caloTrackStr, "Tower")) label_p->DrawLatex(0.264, 0.89, "Calo. jets");
    else if(isStrSame(caloTrackStr, "Trk")) label_p->DrawLatex(0.264, 0.89, "Track jets");
    
    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    gPad->SetLogy();
    gStyle->SetOptStat(0);

    if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    leg_p->Draw("SAME");
    
    std::string saveName = "pdfDir/" + dateStr + "/spectraComp_" + globalStr2 + "_" + caloTrackStr + "_" + centBinsStr[cI] + "_" + dateStr;
    quietSaveAs(canv_p, saveName + ".pdf");
    quietSaveAs(canv_p, saveName + ".png");

    delete leg_p;
    delete canv_p;
  }

  if(isMC){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetLeftMargin(0.16);
      canv_p->SetBottomMargin(0.14);
      
      TLegend* leg_p = new TLegend(0.6, 0.25, 0.95, 0.48);
      leg_p->SetTextFont(43);
      leg_p->SetTextSize(16);
      leg_p->SetBorderSize(0);
      leg_p->SetFillColor(0);
      leg_p->SetFillStyle(0);
      
      TH1D* dummyHist_p = new TH1D("dummyHist_h", ";Truth Jet p_{T} [GeV];Efficiency", nJtPtBins, jtPtBins);
      centerTitles(dummyHist_p);
      dummyHist_p->SetMaximum(1.1);
      dummyHist_p->SetMinimum(0.0);

      dummyHist_p->DrawCopy("");
      
      for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	  if(jtAlgos[aI].find("NoSub") == std::string::npos || cI == periphMostBin){
	    configHist(truthEff_p[aI][cI], aI);
	    
	    leg_p->AddEntry(truthEff_p[aI][cI], jtAlgos[aI].c_str(), "P L");
	    truthEff_p[aI][cI]->Draw("P");
	  }
	}
      }

      line_p->DrawLine(jtPtBins[0], 1.0, jtPtBins[nJtPtBins], 1.0);
      
      std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
      centLabel.replace(centLabel.find("to"), 2, "-");
      centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
      if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";
      
      label_p->DrawLatex(0.25, 0.45, centLabel.c_str());
      if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.25, 0.33, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
      if(isStrSame(caloTrackStr, "Tower")) label_p->DrawLatex(0.25, 0.39, "Calo. jets");
      else if(isStrSame(caloTrackStr, "Trk")) label_p->DrawLatex(0.25, 0.39, "Track jets");
      label_p->DrawLatex(0.25, 0.27, ("p_{T,Reco.} > " + prettyString(minJtPt, 1, false)).c_str());
      
      gStyle->SetOptStat(0);      
      leg_p->Draw("SAME");
      
      std::string saveName = "pdfDir/" + dateStr + "/eff_" + globalStr2 + "_" + caloTrackStr + "_" + centBinsStr[cI] + "_" + dateStr;
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
    for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
      if(jtAlgos[jI].find("ATLAS") != std::string::npos) continue;
      if(jtAlgos[jI].find("Truth") != std::string::npos) continue;

      std::vector<TH1*> centHistMean, centHistSigma, centHistSigmaOverMean, centHistSigmaEta, centHistSigmaPhi;

      if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(jtAlgos[jI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;

 	TCanvas* canvFit_p = new TCanvas("canvFit_p", "", 450*nX, 450*nY);
	canvFit_p->SetTopMargin(0.01);
	canvFit_p->SetLeftMargin(0.01);
	canvFit_p->SetBottomMargin(0.01);
	canvFit_p->SetRightMargin(0.01);
	
	canvFit_p->Divide(nX, nY);
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	std::string nameStr = jtAlgos[jI] + "_" + centBinsStr[cI];
	recoOverGenMean_p[jI][cI] = new TH1D(("recoOverGenMean_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#LTReco./Gen.#GT", nJtPtBins, jtPtBins);
	recoOverGenSigma_p[jI][cI] = new TH1D(("recoOverGenSigma_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma(Reco./Gen.)", nJtPtBins, jtPtBins);
	recoOverGenSigmaOverMean_p[jI][cI] = new TH1D(("recoOverGenSigmaOverMean_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma(Reco./Gen.)/#LTReco./Gen.#GT", nJtPtBins, jtPtBins);
	recoGenSigma_DeltaEta_p[jI][cI] = new TH1D(("recoGenSigma_DeltaEta_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma(#eta_{Reco.} - #eta_{Gen.})", nJtPtBins, jtPtBins);
	recoGenSigma_DeltaPhi_p[jI][cI] = new TH1D(("recoGenSigma_DeltaPhi_" + nameStr + "_h").c_str(), ";Truth Jet p_{T} [GeV];#sigma(#phi_{Reco.} - #phi_{Gen.})", nJtPtBins, jtPtBins);
	
	if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  

	centerTitles({recoOverGenMean_p[jI][cI], recoOverGenSigma_p[jI][cI], recoOverGenSigmaOverMean_p[jI][cI], recoGenSigma_DeltaEta_p[jI][cI], recoGenSigma_DeltaPhi_p[jI][cI]});
	
	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  centerTitles(recoOverGen_VPt_p[jI][cI][jI2]);

	  //Pre-fit to define our fit range

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  if(false){
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

	  if(false){
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

	  recoOverGen_VPt_p[jI][cI][jI2]->GetXaxis()->SetTitleSize(20);
	  recoOverGen_VPt_p[jI][cI][jI2]->GetYaxis()->SetTitleSize(20);

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  recoOverGen_VPt_p[jI][cI][jI2]->DrawCopy("HIST E1");
	  if(false) recoOverGenFit_p[jI][cI][jI2]->DrawCopy("SAME");	  
	
	  std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	  centLabel.replace(centLabel.find("to"), 2, "-");
	  centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
	  if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";

	  centLabel = jtAlgos[jI] + ", " + centLabel;
	  
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
	
	configHist(recoOverGenMean_p[jI][cI], cI);
	configHist(recoOverGenSigma_p[jI][cI], cI);
	configHist(recoOverGenSigmaOverMean_p[jI][cI], cI);
	configHist(recoGenSigma_DeltaEta_p[jI][cI], cI);
	configHist(recoGenSigma_DeltaPhi_p[jI][cI], cI);

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
		
	centHistMean.push_back(recoOverGenMean_p[jI][cI]);
	centHistSigma.push_back(recoOverGenSigma_p[jI][cI]);
	centHistSigmaOverMean.push_back(recoOverGenSigmaOverMean_p[jI][cI]);

	centHistSigmaEta.push_back(recoGenSigma_DeltaEta_p[jI][cI]);
	centHistSigmaPhi.push_back(recoGenSigma_DeltaPhi_p[jI][cI]);

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	std::string saveName = "pdfDir/" + dateStr + "/recoOverGen_" + nameStr + "_" + dateStr + ".pdf";
	quietSaveAs(canvFit_p, saveName);
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
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
	  std::cout << "recoGen_DeltaEta_" + nameStr + "_" + jtPtBinsStr[jI] + "_h" << std::endl;
	  std::cout << recoGen_DeltaEta_p[jI][cI][jI2] << std::endl;
	  std::cout << recoGen_DeltaEta_p[jI][cI][jI2]->GetName() << std::endl;
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

	  recoGen_DeltaEta_p[jI][cI][jI2]->GetXaxis()->SetTitleSize(20);
	  recoGen_DeltaEta_p[jI][cI][jI2]->GetYaxis()->SetTitleSize(20);

 	  recoGen_DeltaEta_p[jI][cI][jI2]->SetMinimum(getMinGTZero(recoGen_DeltaEta_p[jI][cI][jI2])/2.);
	  recoGen_DeltaEta_p[jI][cI][jI2]->DrawCopy("HIST E1");
	  gPad->SetLogy();

	  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	  centLabel.replace(centLabel.find("to"), 2, "-");
	  centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
	  if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";

	  centLabel = jtAlgos[jI] + ", " + centLabel;
	  
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

	  recoGen_DeltaPhi_p[jI][cI][jI2]->GetXaxis()->SetTitleSize(20);
	  recoGen_DeltaPhi_p[jI][cI][jI2]->GetYaxis()->SetTitleSize(20);

 	  recoGen_DeltaPhi_p[jI][cI][jI2]->SetMinimum(getMinGTZero(recoGen_DeltaPhi_p[jI][cI][jI2])/2.);
	  recoGen_DeltaPhi_p[jI][cI][jI2]->DrawCopy("HIST E1");
	
	  gPad->SetLogy();

	  std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	  centLabel.replace(centLabel.find("to"), 2, "-");
	  centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
	  if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";

	  centLabel = jtAlgos[jI] + ", " + centLabel;
	  
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


      plotResponseSet(fileConfig_p, centHistMean, centBinsStr, "recoOverGenMean", jtAlgos[jI], dateStr, 0.6, 1.1);
      plotResponseSet(fileConfig_p, centHistSigmaOverMean, centBinsStr, "recoOverGenSigmaOverMean", jtAlgos[jI], dateStr, 0.0, 0.6);
      plotResponseSet(fileConfig_p, centHistSigmaEta, centBinsStr, "recoGenSigma_DeltaEta", jtAlgos[jI], dateStr, 0.0, 0.1);
      plotResponseSet(fileConfig_p, centHistSigmaPhi, centBinsStr, "recoGenSigma_DeltaPhi", jtAlgos[jI], dateStr, 0.0, 0.1);
    }
  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::vector<TH1*> algoHistMean, algoHistSigma, algoHistSigmaOverMean, algoHistSigmaEta, algoHistSigmaPhi;
      std::vector<std::string> jtAlgosLabel;
      
      for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
	if(jtAlgos[jI].find("ATLAS") != std::string::npos) continue;
	if(jtAlgos[jI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;
	if(jtAlgos[jI].find("Truth") != std::string::npos) continue;

	configHist(recoOverGenMean_p[jI][cI], jI);
	configHist(recoOverGenSigma_p[jI][cI], jI);
	configHist(recoOverGenSigmaOverMean_p[jI][cI], jI);
	configHist(recoGenSigma_DeltaEta_p[jI][cI], jI);
	configHist(recoGenSigma_DeltaPhi_p[jI][cI], jI);
	
	algoHistMean.push_back(recoOverGenMean_p[jI][cI]);
	algoHistSigma.push_back(recoOverGenSigma_p[jI][cI]);
	algoHistSigmaOverMean.push_back(recoOverGenSigmaOverMean_p[jI][cI]);
	algoHistSigmaEta.push_back(recoGenSigma_DeltaEta_p[jI][cI]);
	algoHistSigmaPhi.push_back(recoGenSigma_DeltaPhi_p[jI][cI]);

	jtAlgosLabel.push_back(jtAlgos[jI]);
      }

      plotResponseSet(fileConfig_p, algoHistMean, jtAlgosLabel, "recoOverGenMean", centBinsStr[cI], dateStr, 0.6, 1.1);
      plotResponseSet(fileConfig_p, algoHistSigmaOverMean, jtAlgosLabel, "recoOverGenSigmaOverMean", centBinsStr[cI], dateStr, 0.0, 0.6);

      plotResponseSet(fileConfig_p, algoHistSigmaEta, jtAlgosLabel, "recoGenSigma_DeltaEta", centBinsStr[cI], dateStr, 0.0, 0.1);
      plotResponseSet(fileConfig_p, algoHistSigmaPhi, jtAlgosLabel, "recoGenSigma_DeltaPhi", centBinsStr[cI], dateStr, 0.0, 0.1);
    }
  if(doDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
      if(jtAlgos[jI].find("ATLAS") != std::string::npos) continue;
      if(jtAlgos[jI].find("Truth") != std::string::npos) continue;
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(jtAlgos[jI].find("NoSub") != std::string::npos && cI != periphMostBin) continue;

	delete recoOverGenMean_p[jI][cI];
	delete recoOverGenSigma_p[jI][cI];
	delete recoOverGenSigmaOverMean_p[jI][cI];

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
