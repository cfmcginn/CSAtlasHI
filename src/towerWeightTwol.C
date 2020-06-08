//Author: Chris McGinn (2020.02.14)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//All of this is based on athena tool PhysicsAnalysis/HeavyIonPhys/HIEventUtils/HIEventUtils/HITowerWeightTool.h

//c+cpp
#include <iostream>
#include <vector>

//Class header + Local
#include "include/returnRootFileContentsList.h"
#include "include/stringUtil.h"
#include "include/towerWeightTwol.h"

towerWeightTwol::towerWeightTwol(std::string inWeightFileName)
{
  Init(inWeightFileName);
  return;
}

towerWeightTwol::~towerWeightTwol(){Clean();}

Bool_t towerWeightTwol::Init(std::string inWeightFileName)
{
  m_doGlobalDebug = gDebug.GetDoGlobalDebug();  

  m_weightFileName = inWeightFileName;
  if(!check.checkFileExt(m_weightFileName, "root")){
    std::cout << "towerWeightTwol::Init - Given inWeightFileName \'" << inWeightFileName << "\' is invalid. Init failed, return false" << std::endl;
    Clean();
    return false;
  }

  m_weightFile_p = new TFile(m_weightFileName.c_str(), "READ");
  std::vector<std::string> th3List = returnRootFileContentsList(m_weightFile_p, "TH3F");
  if(!vectContainsStr(m_histName, &th3List)){
    std::cout << "towerWeightTwol::Init - Given inWeightFileName \'" << inWeightFileName << "\' does not contain TH3F \'" << m_histName << "\'. Init failed, return false" << std::endl;
    Clean();
    return false;
  }
  m_h3_w_p = (TH3F*)m_weightFile_p->Get(m_histName.c_str());
  
  return true;
}

//Following function is modelled after:
//https://gitlab.cern.ch/atlas/athena/blob/21.2/PhysicsAnalysis/HeavyIonPhys/HIEventUtils/HIEventUtils/HITowerWeightTool.h
Float_t towerWeightTwol::GetEtaPhiResponse(Float_t eta, Float_t phi, Int_t runNumber)
{
  Float_t weight = 1.0;
  runNumber = 1;
  Int_t binX = m_h3_w_p->GetXaxis()->FindBin(eta);
  Int_t binY = m_h3_w_p->GetYaxis()->FindBin(phi);
  Float_t val = m_h3_w_p->GetBinContent(binX, binY, runNumber);
  if(val > 0) weight = 1./val;
  else if(m_doGlobalDebug) std::cout << "WARNING WEIGHT IS ZERO" << std::endl;
  
  return weight;
}

void towerWeightTwol::Clean()
{
  m_weightFileName = "";
  m_doGlobalDebug = false;
  
  if(m_weightFile_p != nullptr){
    m_weightFile_p->Close();
    delete m_weightFile_p;
    m_weightFile_p = nullptr;
  }
  m_h3_w_p = nullptr;
  
  return;
}
