//cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TNamed.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/returnRootFileContentsList.h"
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

void plotResponseSet(std::map<std::string, std::string> params, std::vector<TH1*> histSet, std::vector<std::string> setLabels, std::string globalNameStr, std::string globalLabelStr, std::string dateStr, double globalMinOverride = 100000, double globalMaxOverride = -100000)
{
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
 label_p->DrawLatex(0.2, 0.93, globalLabelStr2.c_str());
  if(params["caloTrackStr"].find("calo") != std::string::npos) label_p->DrawLatex(0.2, 0.87, "Calo. jets");
  else label_p->DrawLatex(0.2, 0.87, "Track jets");
  
  std::string saveName = "pdfDir/" + dateStr + "/" + globalNameStr + "_" + globalLabelStr + "_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;
  delete leg_p;
  delete label_p;
  
  return;
}

int plotClusterHist(std::string inFileName, std::string globalStr = "")
{
  if(!checkFileExt(inFileName, ".root")) return 1;

  std::string globalStr2 = globalStr;
  if(globalStr2.size() == 0) globalStr2 = "NoGlobalStr";
  while(globalStr2.find(" ") != std::string::npos){
    globalStr2.replace(globalStr2.find(" "), 1, "");
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());  
  delete date;

  kirchnerPalette kPal;

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);

  const Int_t nMaxJtPtBins = 50;
  
  const Int_t nMaxCentBins = 20;
  const Int_t nMaxJtAlgo = 8;
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> nameList = returnRootFileContentsList(inFile_p, "TNamed");
  std::map<std::string, std::string> paramMap;
  for(auto const & iter : nameList){
    TNamed* tempName_p = (TNamed*)inFile_p->Get(iter.c_str());
    std::string name = tempName_p->GetName();
    while(name.find("/") != std::string::npos){name.replace(0, name.find("/")+1, "");}
    paramMap[name] = tempName_p->GetTitle();
  }

  //  const Bool_t doATLAS = std::stoi(paramMap["doATLAS"]);
  const Bool_t doTruth = std::stoi(paramMap["doTruth"]);

  const std::string caloTrackStr = paramMap["caloTrackStr"];
  
  const Int_t nJtPtBins = std::stoi(paramMap["nJtPtBins"]);
  std::vector<std::string> jtPtBinsStr = commaSepStringToVect(paramMap["jtPtBins"]);
  Double_t jtPtBins[nMaxJtPtBins+1];
  for(unsigned int jI = 0; jI < jtPtBinsStr.size(); ++jI){
    jtPtBins[jI] = std::stod(jtPtBinsStr[jI]);
  }
  jtPtBinsStr = commaSepStringToVect(paramMap["jtPtBinsStr"]);
  
  const Int_t nCentBins = std::stoi(paramMap["nCentBins"]);
  std::vector<std::string> centBinsStr = commaSepStringToVect(paramMap["centBinsStr"]);
  std::vector<std::string> nEventPerCentStr = commaSepStringToVect(paramMap["nEventPerCent"]);
  const std::string jtAbsEtaMaxStr = paramMap["maxJtAbsEta"];
  std::vector<Double_t> nEventPerCent;
  for(unsigned int cI = 0; cI < nEventPerCentStr.size(); ++cI){
    nEventPerCent.push_back(std::stod(nEventPerCentStr[cI]));
  }
  
  std::vector<std::string> jtAlgos = commaSepStringToVect(paramMap["jtAlgos"]);
  const Int_t nJtAlgo = jtAlgos.size();

  if(nCentBins > nMaxCentBins){
    std::cout << "nCentBins \'" << nCentBins << "\' exceeds maximum \'" << nMaxCentBins << "\'. return 1" << std::endl;
    return 1;
  }

  if(nJtAlgo > nMaxJtAlgo){
    std::cout << "nJtAlgo \'" << nJtAlgo << "\' exceeds maximum \'" << nMaxJtAlgo << "\'. return 1" << std::endl;
    return 1;
  }

  TH1D* spectra_p[nMaxJtAlgo][nMaxCentBins];
  TH1D* spectraTruth_p[nMaxCentBins];
  TH1D* recoOverGen_VPt_p[nMaxJtAlgo][nCentBins][nMaxJtPtBins];
  TH1D* recoGen_DeltaEta_p[nMaxJtAlgo][nCentBins][nMaxJtPtBins];
  TH1D* recoGen_DeltaPhi_p[nMaxJtAlgo][nCentBins][nMaxJtPtBins];
  TH1D* matchedTruthSpectra_p[nMaxJtAlgo][nMaxCentBins];
  TGraphAsymmErrors* truthEff_p[nMaxJtAlgo][nMaxCentBins];

  for(Int_t aI = 0; aI < nJtAlgo; ++ aI){
    for(Int_t cI = 0; cI < nCentBins; ++ cI){
      std::string nameStr = jtAlgos[aI] + "_" + centBinsStr[cI];
      
      spectra_p[aI][cI] = (TH1D*)inFile_p->Get(("spectra_" + nameStr + "_h").c_str());       
      if(doTruth){
	if(aI == 0) spectraTruth_p[cI] = (TH1D*)inFile_p->Get(("spectra_Truth_" + centBinsStr[cI] + "_h").c_str());

	if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	  matchedTruthSpectra_p[aI][cI] = (TH1D*)inFile_p->Get(("matchedTruthSpectra_" + nameStr + "_h").c_str());       
	  truthEff_p[aI][cI] = new TGraphAsymmErrors();
	  truthEff_p[aI][cI]->Divide(matchedTruthSpectra_p[aI][cI], spectraTruth_p[cI]);
	}
      }

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

      if(doTruth){
	if(jtAlgos[aI].find("ATLAS") != std::string::npos) continue;
	if(jtAlgos[aI].find("Truth") != std::string::npos) continue;
	
	for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	  recoOverGen_VPt_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoOverGen_VPt_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());
	  recoGen_DeltaEta_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoGen_DeltaEta_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());
	  recoGen_DeltaPhi_p[aI][cI][jI] = (TH1D*)inFile_p->Get(("recoGen_DeltaPhi_" + nameStr + "_" + jtPtBinsStr[jI] + "_h").c_str());

	  configHist(recoOverGen_VPt_p[aI][cI][jI], aI);
	  configHist(recoGen_DeltaEta_p[aI][cI][jI], aI);
	  configHist(recoGen_DeltaPhi_p[aI][cI][jI], aI);
	}
      }
    }
  }


  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);

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

    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      leg_p->AddEntry(spectra_p[aI][cI], jtAlgos[aI].c_str(), "P L");
		      
      spectra_p[aI][cI]->SetMaximum(3.*globalMax);
      spectra_p[aI][cI]->SetMinimum(globalMin/2.);

      if(aI == 0) spectra_p[aI][cI]->DrawCopy("HIST E1 P");
      else spectra_p[aI][cI]->DrawCopy("HIST E1 P SAME");
    }

    std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
    centLabel.replace(centLabel.find("to"), 2, "-");
    centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
    if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";
    
    label_p->DrawLatex(0.24, 0.95, centLabel.c_str());
    if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.24, 0.83, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
    if(isStrSame(caloTrackStr, "calo")) label_p->DrawLatex(0.24, 0.89, "Calo. jets");
    else if(isStrSame(caloTrackStr, "trk")) label_p->DrawLatex(0.24, 0.89, "Track jets");
    
    gPad->SetLogy();
    gStyle->SetOptStat(0);

    leg_p->Draw("SAME");
    
    std::string saveName = "pdfDir/" + dateStr + "/spectraComp_" + globalStr2 + "_" + caloTrackStr + "_" + centBinsStr[cI] + "_" + dateStr;
    quietSaveAs(canv_p, saveName + ".pdf");
    quietSaveAs(canv_p, saveName + ".png");

    delete leg_p;
    delete canv_p;
  }

  if(doTruth){
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
      
      TH1D* dummyHist_p = new TH1D("dummyHist_h", ";Jet p_{T} [GeV];Efficiency", nJtPtBins, jtPtBins);
      centerTitles(dummyHist_p);
      dummyHist_p->SetMaximum(1.1);
      dummyHist_p->SetMinimum(0.0);

      dummyHist_p->DrawCopy("");
      
      for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	if(jtAlgos[aI].find("ATLAS") == std::string::npos){
	  configHist(truthEff_p[aI][cI], aI);
	  
	  leg_p->AddEntry(truthEff_p[aI][cI], jtAlgos[aI].c_str(), "P L");
	  truthEff_p[aI][cI]->Draw("P");
	}
      }
      
      std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
      centLabel.replace(centLabel.find("to"), 2, "-");
      centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
      if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";
      
      label_p->DrawLatex(0.44, 0.45, centLabel.c_str());
      if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.44, 0.33, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
      if(isStrSame(caloTrackStr, "calo")) label_p->DrawLatex(0.44, 0.39, "Calo. jets");
      else if(isStrSame(caloTrackStr, "trk")) label_p->DrawLatex(0.44, 0.39, "Track jets");
      
      gStyle->SetOptStat(0);      
      leg_p->Draw("SAME");
      
      std::string saveName = "pdfDir/" + dateStr + "/eff_" + globalStr2 + "_" + caloTrackStr + "_" + centBinsStr[cI] + "_" + dateStr;
      quietSaveAs(canv_p, saveName + ".pdf");
      //      quietSaveAs(canv_p, saveName + ".png");
      
      delete dummyHist_p;
      delete leg_p;
      delete canv_p;
    }


    TH1D* recoOverGenMean_p[nMaxJtAlgo][nMaxCentBins];
    TH1D* recoOverGenSigma_p[nMaxJtAlgo][nMaxCentBins];
    TH1D* recoOverGenSigmaOverMean_p[nMaxJtAlgo][nMaxCentBins];

    Int_t nX = 1;
    Int_t nY = 1;
    if(nJtPtBins == 2) nX = 2;
    else if(nJtPtBins == 3) nX = 3;
    else if(nJtPtBins == 4) nX = 4;
    else if(nJtPtBins == 5 || nJtPtBins == 6){nX = 3; nY = 2;}
    else if(nJtPtBins == 7 || nJtPtBins == 8){nX = 4; nY = 2;} 
 
    for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
      if(jtAlgos[jI].find("ATLAS") != std::string::npos) continue;
      if(jtAlgos[jI].find("Truth") != std::string::npos) continue;

      std::vector<TH1*> centHistMean, centHistSigma, centHistSigmaOverMean;
            
      for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
	TCanvas* canvFit_p = new TCanvas("canvFit_p", "", 450*nX, 450*nY);
	canvFit_p->SetTopMargin(0.01);
	canvFit_p->SetLeftMargin(0.01);
	canvFit_p->SetBottomMargin(0.01);
	canvFit_p->SetRightMargin(0.01);
	
	canvFit_p->Divide(nX, nY);
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	std::string nameStr = jtAlgos[jI] + "_" + centBinsStr[cI];
	recoOverGenMean_p[jI][cI] = new TH1D(("recoOverGenMean_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];#LTReco./Gen.#GT", nJtPtBins, jtPtBins);
	recoOverGenSigma_p[jI][cI] = new TH1D(("recoOverGenSigma_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];#sigma(Reco./Gen.)", nJtPtBins, jtPtBins);
	recoOverGenSigmaOverMean_p[jI][cI] = new TH1D(("recoOverGenSigmaOverMean_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];#sigma(Reco./Gen.)/#LTReco./Gen.#GT", nJtPtBins, jtPtBins);

	centerTitles({recoOverGenMean_p[jI][cI], recoOverGenSigma_p[jI][cI], recoOverGenSigmaOverMean_p[jI][cI]});
	
	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  centerTitles(recoOverGen_VPt_p[jI][cI][jI2]);
	  
	  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  recoOverGenMean_p[jI][cI]->SetBinContent(jI2+1, recoOverGen_VPt_p[jI][cI][jI2]->GetMean());
	  recoOverGenMean_p[jI][cI]->SetBinError(jI2+1, recoOverGen_VPt_p[jI][cI][jI2]->GetMeanError());
	  recoOverGenSigma_p[jI][cI]->SetBinContent(jI2+1, recoOverGen_VPt_p[jI][cI][jI2]->GetStdDev());
	  recoOverGenSigma_p[jI][cI]->SetBinError(jI2+1, recoOverGen_VPt_p[jI][cI][jI2]->GetStdDevError());
	  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  recoOverGenSigmaOverMean_p[jI][cI]->SetBinContent(jI2+1, recoOverGen_VPt_p[jI][cI][jI2]->GetStdDev()/recoOverGen_VPt_p[jI][cI][jI2]->GetMean());
	  recoOverGenSigmaOverMean_p[jI][cI]->SetBinError(jI2+1, recoOverGen_VPt_p[jI][cI][jI2]->GetStdDevError()/recoOverGen_VPt_p[jI][cI][jI2]->GetMean());

	  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  canvFit_p->cd();
	  canvFit_p->cd(jI2+1);

	  gPad->SetRightMargin(0.01);
	  gPad->SetTopMargin(0.01);
	  gPad->SetLeftMargin(0.14);
	  gPad->SetBottomMargin(0.12);
	  
	  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	  recoOverGen_VPt_p[jI][cI][jI2]->DrawCopy("HIST E1");

	
	  std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	  centLabel.replace(centLabel.find("to"), 2, "-");
	  centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
	  if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";

	  centLabel = jtAlgos[jI] + ", " + centLabel;
	  
	  label_p->DrawLatex(0.2, 0.95, centLabel.c_str());
	  std::string ptStr = prettyString(jtPtBins[jI2], 1, false) + "<p_{T,Jet}<" + prettyString(jtPtBins[jI2+1], 1, false);
	  label_p->DrawLatex(0.2, 0.83, ptStr.c_str());
	  
	  if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.2, 0.77, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
	  if(isStrSame(caloTrackStr, "calo")) label_p->DrawLatex(0.2, 0.89, "Calo. jets");
	  else if(isStrSame(caloTrackStr, "trk")) label_p->DrawLatex(0.2, 0.89, "Track jets");
	}

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
	configHist(recoOverGenMean_p[jI][cI], cI);
	configHist(recoOverGenSigma_p[jI][cI], cI);
	configHist(recoOverGenSigmaOverMean_p[jI][cI], cI);

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
		
	centHistMean.push_back(recoOverGenMean_p[jI][cI]);
	centHistSigma.push_back(recoOverGenSigma_p[jI][cI]);
	centHistSigmaOverMean.push_back(recoOverGenSigmaOverMean_p[jI][cI]);

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	std::string saveName = "pdfDir/" + dateStr + "/recoOverGen_" + nameStr + "_" + dateStr + ".pdf";
	quietSaveAs(canvFit_p, saveName);
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	delete canvFit_p;
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      }


      for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
	TCanvas* canvFit_p = new TCanvas("canvFit_p", "", 450*nX, 450*nY);
	canvFit_p->SetTopMargin(0.01);
	canvFit_p->SetLeftMargin(0.01);
	canvFit_p->SetBottomMargin(0.01);
	canvFit_p->SetRightMargin(0.01);
	
	canvFit_p->Divide(nX, nY);

	std::string nameStr = jtAlgos[jI] + "_" + centBinsStr[cI];

	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	  centerTitles(recoGen_DeltaEta_p[jI][cI][jI2]);
	  
	  canvFit_p->cd();
	  canvFit_p->cd(jI2+1);

	  gPad->SetRightMargin(0.01);
	  gPad->SetTopMargin(0.01);
	  gPad->SetLeftMargin(0.14);
	  gPad->SetBottomMargin(0.12);

 	  recoGen_DeltaEta_p[jI][cI][jI2]->SetMinimum(0.0);
	  recoGen_DeltaEta_p[jI][cI][jI2]->DrawCopy("HIST E1");
	
	  std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	  centLabel.replace(centLabel.find("to"), 2, "-");
	  centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
	  if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";

	  centLabel = jtAlgos[jI] + ", " + centLabel;
	  
	  label_p->DrawLatex(0.2, 0.95, centLabel.c_str());
	  std::string ptStr = prettyString(jtPtBins[jI2], 1, false) + "<p_{T,Jet}<" + prettyString(jtPtBins[jI2+1], 1, false);
	  label_p->DrawLatex(0.2, 0.83, ptStr.c_str());
	  
	  if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.2, 0.77, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
	  if(isStrSame(caloTrackStr, "calo")) label_p->DrawLatex(0.2, 0.89, "Calo. jets");
	  else if(isStrSame(caloTrackStr, "trk")) label_p->DrawLatex(0.2, 0.89, "Track jets");

	  label_p->DrawLatex(0.2, 0.71, ("#sigma=" + prettyString(recoGen_DeltaEta_p[jI][cI][jI2]->GetStdDev(), 3, false)).c_str());
	}


	std::string saveName = "pdfDir/" + dateStr + "/recoGen_DeltaEta_" + nameStr + "_" + dateStr + ".pdf";
	quietSaveAs(canvFit_p, saveName);
	delete canvFit_p;
      }

      for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
	TCanvas* canvFit_p = new TCanvas("canvFit_p", "", 450*nX, 450*nY);
	canvFit_p->SetTopMargin(0.01);
	canvFit_p->SetLeftMargin(0.01);
	canvFit_p->SetBottomMargin(0.01);
	canvFit_p->SetRightMargin(0.01);
	
	canvFit_p->Divide(nX, nY);

	std::string nameStr = jtAlgos[jI] + "_" + centBinsStr[cI];

	for(Int_t jI2 = 0; jI2 < nJtPtBins; ++jI2){
	  centerTitles(recoGen_DeltaPhi_p[jI][cI][jI2]);
	  
	  canvFit_p->cd();
	  canvFit_p->cd(jI2+1);

	  gPad->SetRightMargin(0.01);
	  gPad->SetTopMargin(0.01);
	  gPad->SetLeftMargin(0.14);
	  gPad->SetBottomMargin(0.12);
	  
 	  recoGen_DeltaPhi_p[jI][cI][jI2]->SetMinimum(0.0);
	  recoGen_DeltaPhi_p[jI][cI][jI2]->DrawCopy("HIST E1");
	
	  std::string centLabel = centBinsStr[cI].substr(4, centBinsStr[cI].size()) + "%";
	  centLabel.replace(centLabel.find("to"), 2, "-");
	  centLabel = "#bf{#color[" + std::to_string(kRed) + "]{" + centLabel + "}}";
	  if(globalStr.size() != 0) centLabel = centLabel + "; #bf{" + globalStr + "}";

	  centLabel = jtAlgos[jI] + ", " + centLabel;
	  
	  label_p->DrawLatex(0.2, 0.95, centLabel.c_str());
	  std::string ptStr = prettyString(jtPtBins[jI2], 1, false) + "<p_{T,Jet}<" + prettyString(jtPtBins[jI2+1], 1, false);
	  label_p->DrawLatex(0.2, 0.83, ptStr.c_str());
	  
	  if(jtAbsEtaMaxStr.size() != 0) label_p->DrawLatex(0.2, 0.77, ("|#eta_{jet}| < " + jtAbsEtaMaxStr).c_str());
	  if(isStrSame(caloTrackStr, "calo")) label_p->DrawLatex(0.2, 0.89, "Calo. jets");
	  else if(isStrSame(caloTrackStr, "trk")) label_p->DrawLatex(0.2, 0.89, "Track jets");

	  label_p->DrawLatex(0.2, 0.71, ("#sigma=" + prettyString(recoGen_DeltaPhi_p[jI][cI][jI2]->GetStdDev(), 3, false)).c_str());
	}


	std::string saveName = "pdfDir/" + dateStr + "/recoGen_DeltaPhi_" + nameStr + "_" + dateStr + ".pdf";
	quietSaveAs(canvFit_p, saveName);
	delete canvFit_p;
      }


      plotResponseSet(paramMap, centHistMean, centBinsStr, "recoOverGenMean", jtAlgos[jI], dateStr, 0.6, 1.6);
      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      plotResponseSet(paramMap, centHistSigmaOverMean, centBinsStr, "recoOverGenSigmaOverMean", jtAlgos[jI], dateStr, 0.0, 0.6);
      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    }

    //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    
    for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
      std::vector<TH1*> algoHistMean, algoHistSigma, algoHistSigmaOverMean;
      std::vector<std::string> jtAlgosLabel;
      
      for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
	if(jtAlgos[jI].find("ATLAS") != std::string::npos) continue;
	if(jtAlgos[jI].find("Truth") != std::string::npos) continue;

	configHist(recoOverGenMean_p[jI][cI], jI);
	configHist(recoOverGenSigma_p[jI][cI], jI);
	configHist(recoOverGenSigmaOverMean_p[jI][cI], jI);
	
	algoHistMean.push_back(recoOverGenMean_p[jI][cI]);
	algoHistSigma.push_back(recoOverGenSigma_p[jI][cI]);
	algoHistSigmaOverMean.push_back(recoOverGenSigmaOverMean_p[jI][cI]);

	jtAlgosLabel.push_back(jtAlgos[jI]);
      }

      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
      plotResponseSet(paramMap, algoHistMean, jtAlgosLabel, "recoOverGenMean", centBinsStr[cI], dateStr, 0.6, 1.6);
      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      plotResponseSet(paramMap, algoHistSigmaOverMean, jtAlgosLabel, "recoOverGenSigmaOverMean", centBinsStr[cI], dateStr, 0.0, 0.6);
      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    }
    
    //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
      if(jtAlgos[jI].find("ATLAS") != std::string::npos) continue;
      if(jtAlgos[jI].find("Truth") != std::string::npos) continue;
      for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
	delete recoOverGenMean_p[jI][cI];
	delete recoOverGenSigma_p[jI][cI];
	delete recoOverGenSigmaOverMean_p[jI][cI];
      }
    }
  }

  for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
    if(jtAlgos[jI].find("ATLAS") == std::string::npos){
      for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
	delete truthEff_p[jI][cI];
      }
    }
  }
  
  delete label_p;
  
  inFile_p->Close();
  delete inFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 3){
    std::cout << "Usage: ./bin/plotClusterHist.exe <inFileName> <globalStr-default-"">. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 2) retVal += plotClusterHist(argv[1]);
  else if(argc == 3) retVal += plotClusterHist(argv[1], argv[2]);
  return retVal;
}
