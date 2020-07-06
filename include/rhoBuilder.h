//class/Author: Chris McGinn (2020.01.28)

#ifndef RHOBUILDER_H
#define RHOBUILDER_H

//cpp
#include <vector>

//Fastjet
#include "fastjet/PseudoJet.hh"

//Local
#include "include/globalDebugHandler.h"

class rhoBuilder{
 public:
  rhoBuilder(){};
  ~rhoBuilder(){};

  rhoBuilder(std::vector<float> inEtaBins);
  rhoBuilder(std::vector<double> inEtaBins);
  bool Init(std::vector<float> inEtaBins);
  bool Init(std::vector<double> inEtaBins);
  bool CalcRhoFromPseudoJet(std::vector<fastjet::PseudoJet>* constituents_p, std::vector<fastjet::PseudoJet>* jets_p=nullptr, bool doTowerExclude=false);
  bool CalcRhoFromPtEtaPhi(std::vector<float>* pt_p, std::vector<float>* eta_p, std::vector<float>* phi_p, std::vector<fastjet::PseudoJet>* jets_p=nullptr, bool doTowerExclude=false);
  bool CalcRhoFromPtEtaPhi(std::vector<double>* pt_p, std::vector<double>* eta_p, std::vector<double>* phi_p, std::vector<fastjet::PseudoJet>* jets_p=nullptr, bool doTowerExclude=false);
  bool SetRho(std::vector<double>* rho_p);
  bool SetRho(std::vector<float>* rho_p);
  bool SetRhoPt(std::vector<double>* rhoPt_p);
  bool SetRhoPt(std::vector<float>* rhoPt_p);
  void Clean();
  void Print();
  
 private:  
  bool m_doGlobalDebug;
  globalDebugHandler gBug;
  
  bool m_isInit;
  std::vector<double> m_etaBins;
  std::vector<double> m_rhoVals;  
  std::vector<double> m_rhoPtVals;  
  std::vector<int> m_nExcluded;

  std::vector<double> m_towerPhiBounds;

  int posInBins(double val, std::vector<double>* bins);
};

#endif 
