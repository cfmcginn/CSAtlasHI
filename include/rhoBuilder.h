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
  bool CalcRhoFromPseudoJet(std::vector<fastjet::PseudoJet>* constituents_p);
  bool CalcRhoFromPtEta(std::vector<float>* pt_p, std::vector<float>* eta_p);
  bool CalcRhoFromPtEta(std::vector<double>* pt_p, std::vector<double>* eta_p);
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

  int posInBins(double val, std::vector<double>* bins);
};

#endif 
