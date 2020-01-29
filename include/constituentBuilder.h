//Author: Chris McGinn (2020.01.23)

//cpp
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//FastJet
#include "fastjet/PseudoJet.hh"

//Local
#include "include/globalDebugHandler.h"

class constituentBuilder
{
 public:
  constituentBuilder(){};
  constituentBuilder(std::vector<float>* pt_p, std::vector<float>* eta_p, std::vector<float>* phi_p, double inUserMinPt = -1);
  ~constituentBuilder(){};
  
  bool InitPtEtaPhi(std::vector<float>* pt_p, std::vector<float>* eta_p, std::vector<float>* phi_p, double inUserMinPt = -1);
  std::vector<fastjet::PseudoJet> GetCleanInputs();
  std::vector<fastjet::PseudoJet> GetGhostedInputs();
  std::vector<fastjet::PseudoJet> GetAllInputs();
  bool IsUserIndexGhosted(unsigned int inUserIndex);
  void Clean();
  
 private:
  globalDebugHandler gBug;
  bool m_doGlobalDebug;

  const double minPt = 0.001;
  
  std::map<unsigned int, double> mapToOrigPt;
  std::map<unsigned int, bool> mapToIsSwapped;
  std::vector<fastjet::PseudoJet> cleanInputs; // This is the collection 'to-be-clustered' at the end of day
  std::vector<fastjet::PseudoJet> ghostedInputs; // This is the collection that is ghosted w/ neg values
};
