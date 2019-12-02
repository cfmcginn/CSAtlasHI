//cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
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

  const Int_t nColors = 5;
  const Int_t colors[nColors] = {0, 1, 3, 4, 6};

  const Int_t nStyles = 4;
  const Int_t styles[nStyles] = {24, 25, 28, 46};

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
  std::vector<std::string> jtPtBinsStr = commaSepStringToVect(paramMap["jtPtBinsStr"]);
  Double_t jtPtBins[nMaxJtPtBins+1];
  for(unsigned int jI = 0; jI < jtPtBinsStr.size(); ++jI){
    jtPtBins[jI] = std::stod(jtPtBinsStr[jI]);
  }
  
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
  //  TH1D* matchedSpectra_p[nMaxJtAlgo][nMaxCentBins];

  for(Int_t aI = 0; aI < nJtAlgo; ++ aI){
    for(Int_t cI = 0; cI < nCentBins; ++ cI){
      spectra_p[aI][cI] = (TH1D*)inFile_p->Get(("spectra_" + jtAlgos[aI] + "_" + centBinsStr[cI] + "_h").c_str());

      setSumW2(spectra_p[aI][cI]);

      spectra_p[aI][cI]->Scale(1./nEventPerCent[cI]);
      for(Int_t bIX = 0; bIX < spectra_p[aI][cI]->GetXaxis()->GetNbins(); ++bIX){
	double val = spectra_p[aI][cI]->GetBinContent(bIX+1)/spectra_p[aI][cI]->GetBinWidth(bIX+1);
	double err = spectra_p[aI][cI]->GetBinError(bIX+1)/spectra_p[aI][cI]->GetBinWidth(bIX+1);

	spectra_p[aI][cI]->SetBinContent(bIX+1, val);
	spectra_p[aI][cI]->SetBinError(bIX+1, err);
      }

      spectra_p[aI][cI]->SetMarkerSize(1);
      spectra_p[aI][cI]->SetMarkerStyle(styles[aI%nStyles]);
      spectra_p[aI][cI]->SetMarkerColor(kPal.getColor(colors[aI%nColors]));
      spectra_p[aI][cI]->SetLineColor(kPal.getColor(colors[aI%nColors]));
      
      spectra_p[aI][cI]->GetXaxis()->SetTitleFont(43);
      spectra_p[aI][cI]->GetYaxis()->SetTitleFont(43);
      spectra_p[aI][cI]->GetXaxis()->SetLabelFont(43);
      spectra_p[aI][cI]->GetYaxis()->SetLabelFont(43);
      
      spectra_p[aI][cI]->GetXaxis()->SetTitleSize(15);
      spectra_p[aI][cI]->GetYaxis()->SetTitleSize(15);
      spectra_p[aI][cI]->GetXaxis()->SetLabelSize(13);
      spectra_p[aI][cI]->GetYaxis()->SetLabelSize(13);

      spectra_p[aI][cI]->GetYaxis()->SetTitleOffset(1.9);
      spectra_p[aI][cI]->GetYaxis()->SetTitle("#frac{1}{N_{Event}} #frac{dN_{Jet}}{dp_{T}} [Gev^{-1}]");
      centerTitles(spectra_p[aI][cI]);
      std::cout << "HIST: " << spectra_p[aI][cI]->GetName() << std::endl;
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
    TH1D* recoOverGenMean_p[nMaxJtAlgo][nMaxCentBins];
    TH1D* recoOverGenSigma_p[nMaxJtAlgo][nMaxCentBins];
    TH1D* recoOverGenSigmaOverMean_p[nMaxJtAlgo][nMaxCentBins];

    for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
      for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
	std::string nameStr = jtAlgos[jI] + "_" + centBinsStr[cI];
	recoOverGenMean_p[jI][cI] = new TH1D(("recoOverGenMean_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];#LTReco./Gen.#GT", nJtPtBins, jtPtBins);
	recoOverGenSigma_p[jI][cI] = new TH1D(("recoOverGenSigma_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];#sigma(Reco./Gen.)", nJtPtBins, jtPtBins);
	recoOverGenSigmaOverMean_p[jI][cI] = new TH1D(("recoOverGenSigmaOverMean_" + nameStr + "_h").c_str(), ";Jet p_{T} [GeV];#sigma(Reco./Gen.)", nJtPtBins, jtPtBins);
      }
    }    



    for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
      for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
	delete recoOverGenMean_p[jI][cI];
	delete recoOverGenSigma_p[jI][cI];
	delete recoOverGenSigmaOverMean_p[jI][cI];
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
