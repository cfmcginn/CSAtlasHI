#ifndef PLOTUTILITIES_H
#define PLOTUTILITIES_H

//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TBox.h"
#include "TCanvas.h"
#include "TError.h"
#include "TH1.h"
#include "TMath.h"

std::string prettyString(double inVal, const int prec, const bool doDot)
{
  const int maxPrec = 18;
  if(prec > maxPrec){
    std::cout << "PRETTYSTRING ERROR: CANNOT HANDLE REQUESTED PRECISION \'" << prec << "\', max is \'" << maxPrec << "\'. return empty string" << std::endl;
    return "";
  }

  std::string minStr = "";
  if(inVal < 0) minStr = "-";
  inVal = TMath::Abs(inVal);
  std::string retStr = "";
  
  inVal *= TMath::Power(10, prec);

  unsigned long long tempInVal = inVal;
  if(inVal - tempInVal > 0.5) ++tempInVal;
  retStr = std::to_string(tempInVal);
  
  retStr = retStr.substr(0, retStr.size()-prec) + "." + retStr.substr(retStr.size()-prec, prec);

  if(retStr.substr(0,1).find(".") != std::string::npos) retStr = "0" + retStr;
  if(retStr.substr(retStr.size()-1,1).find(".") != std::string::npos) retStr.replace(retStr.size()-1, 1, "");

  if(doDot){
    if(retStr.find(".") != std::string::npos) retStr.replace(retStr.find("."), 1, "p");
    if(minStr.size() != 0) minStr = "Neg";
  }
    
  return minStr + retStr;
}

std::string prettyStringE(const double inVal, const int prec, const bool doDot)
{
  std::string retStr = prettyString(inVal, prec, false);
  int tenScale = retStr.find(".") - 1;
  while(retStr.find(".") != std::string::npos){
    retStr.replace(retStr.find("."), 1, "");
  }
  if(!doDot) retStr = retStr.substr(0,1) + "." + retStr.substr(1,retStr.size());
  else retStr = retStr.substr(0,1) + "p" + retStr.substr(1,retStr.size());

  while(retStr.find(".")+prec+1 < retStr.size()){
    retStr = retStr.substr(0, retStr.size()-1);
  }
  
  retStr = retStr + "E";
  if(tenScale >= 0) retStr = retStr + "+" + std::to_string(std::abs(tenScale));
  else retStr = retStr + "-" + std::to_string(std::abs(tenScale));
  
  return retStr;
}


void prettyCanv(TCanvas* canv_p)
{
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(1.5*canv_p->GetLeftMargin());
  canv_p->SetBottomMargin(canv_p->GetLeftMargin());
  canv_p->SetTopMargin(canv_p->GetLeftMargin()/2.);

  return;
}


void prettyTH1(TH1* hist_p, const double size, const int style, const int col)
{
  hist_p->SetMarkerSize(size);
  hist_p->SetMarkerStyle(style);
  hist_p->SetMarkerColor(col);
  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();

  return;
}


void drawWhiteBox(Double_t x1, Double_t x2, Double_t y1, Double_t y2)
{
  TBox* tempBox_p = new TBox();
  tempBox_p->SetFillColor(0);
  tempBox_p->DrawBox(x1, y1, x2, y2);
  delete tempBox_p;
}


void quietSaveAs(TCanvas* canv_p, const std::string saveName)
{
  Int_t oldLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  canv_p->SaveAs(saveName.c_str());
  gErrorIgnoreLevel = oldLevel;
  
  return;
}


double getNearestFactor10Up(double inVal, UInt_t steps = 0)
{
  double val = 0.00000000000001;
 
  for(UInt_t i = 0; i < 28; ++i){
    //    std::cout << val << std::endl;
    if(inVal >= val && inVal < val*10){
      val *= 10;
      break;
    }
    else val *= 10;
  }

  for(UInt_t i = 0; i < steps; ++i){
    val *= 10;
  }
  
  return val;
}

double getNearestFactor10Down(double inVal, UInt_t steps = 0)
{
  double val = 0.00000000000001;
 
  for(UInt_t i = 0; i < 28; ++i){
    //    std::cout << val << std::endl;
    if(inVal >= val && inVal < val*10){
      break;
    }
    else val *= 10;
  }

  for(UInt_t i = 0; i < steps; ++i){
    val /= 10;
  }
  
  return val;
}

#endif
