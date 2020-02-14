//Author: Chris McGinn (2020.02.06)

//cpp
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/getLinBins.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int validateRhoPlot(std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TH1D* deltaEt_p = (TH1D*)inFile_p->Get("deltaEt_h");
  TH2D* deltaEtVEta_p = (TH2D*)inFile_p->Get("deltaEtVEta_h");
  TH2D* deltaEtVCent_p = (TH2D*)inFile_p->Get("deltaEtVCent_h");

  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.02);
  canv_p->SetRightMargin(0.02);
  canv_p->SetBottomMargin(0.12);
  canv_p->SetLeftMargin(0.12);

  deltaEt_p->DrawCopy("HIST E1 P");
  gStyle->SetOptStat(0);
  
  std::string saveName = "pdfDir/" + dateStr + "/deltaEt_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;

  canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.03);
  canv_p->SetRightMargin(0.12);
  canv_p->SetBottomMargin(0.12);
  canv_p->SetLeftMargin(0.12);

  deltaEtVEta_p->DrawCopy("COLZ");
  gStyle->SetOptStat(0);
  
  saveName = "pdfDir/" + dateStr + "/deltaEtVEta_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;

  canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.03);
  canv_p->SetRightMargin(0.12);
  canv_p->SetBottomMargin(0.12);
  canv_p->SetLeftMargin(0.12);

  deltaEtVCent_p->DrawCopy("COLZ");
  gStyle->SetOptStat(0);
  
  saveName = "pdfDir/" + dateStr + "/deltaEtVCent_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;

  inFile_p->Close();
  delete inFile_p;
  
  std::cout << "VALIDATE RHO PLOT COMPLETE. return 0" << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/validateRhoPlot.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += validateRhoPlot(argv[1]);
  return retVal;
}
