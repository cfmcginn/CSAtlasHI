//Author: Chris McGinn (2020.02.13)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//All of this is based on athena tool PhysicsAnalysis/HeavyIonPhys/HIEventUtils/HIEventUtils/HITowerWeightTool.h

#ifndef TOWERWEIGHTTWOL_H
#define TOWERWEIGHTTWOL_H

//c+cpp
#include <string>

//ROOT
#include "TFile.h"
#include "TH3F.h"

//Local
#include "include/checkMakeDir.h"
#include "include/globalDebugHandler.h"

class towerWeightTwol{
 public:
  towerWeightTwol(){};
  towerWeightTwol(std::string inWeightFileName);  
  ~towerWeightTwol();
  
  Bool_t Init(std::string inWeightFileName);
  Float_t GetEtaPhiResponse(Float_t eta, Float_t phi, Int_t runNumber);
  void Clean();
  
 private:
  checkMakeDir check;
  globalDebugHandler gDebug;

  bool m_doGlobalDebug = false;
  
  std::string m_weightFileName;
  const std::string m_histName = "h3_eta_phi_response"; 

  TFile* m_weightFile_p = nullptr;
  TH3F* m_h3_w_p = nullptr;
};

#endif
