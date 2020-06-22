//Author: Chris McG!inn (2020.01.28)

//cpp
#include <cmath>
#include <iostream>

//Local
#include "include/rhoBuilder.h"

//public member functions
rhoBuilder::rhoBuilder(std::vector<float> inEtaBins)
{
  Init(inEtaBins);
  return;
}

rhoBuilder::rhoBuilder(std::vector<double> inEtaBins)
{
  Init(inEtaBins);
  return;
}

bool rhoBuilder::Init(std::vector<float> inEtaBins)
{
  std::vector<double> inEtaBinsD;
  for(unsigned int eI = 0; eI < inEtaBins.size(); ++eI){
    inEtaBinsD.push_back(inEtaBins[eI]);
  }
  return Init(inEtaBinsD);
}

bool rhoBuilder::Init(std::vector<double> inEtaBins)
{
  Clean();//Call clean in case not invoked pre-init
  m_doGlobalDebug = gBug.GetDoGlobalDebug();
  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  m_isInit = true;

  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(unsigned int eI = 0; eI < inEtaBins.size(); ++eI){
    if(eI != 0){
      m_rhoVals.push_back(0.0);
      m_rhoPtVals.push_back(0.0);
      if(inEtaBins[eI] < inEtaBins[eI-1]){
	m_isInit = false;
	break;
      }
    }

    m_etaBins.push_back(inEtaBins[eI]);
  }

  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(!m_isInit){
    std::cout << "ERROR IN RHOBUILDER INIT: Given etaBins input is not ordered. Clean and return false" << std::endl;
    std::cout << " Eta bins given: \'";
    for(unsigned int eI = 0; eI < inEtaBins.size(); ++eI){
      std::cout << inEtaBins[eI] << ",";
    }
    std::cout << "\'" << std::endl;
    
    Clean();
  }

  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE, m_isInit: " << __FILE__ << ", " << __LINE__ << ", " << m_isInit << std::endl;

  return m_isInit;
}


bool rhoBuilder::CalcRhoFromPseudoJet(std::vector<fastjet::PseudoJet>* constituents_p)
{
  std::vector<double> ptD_, etaD_;
  for(unsigned int cI = 0; cI < constituents_p->size(); ++cI){
    ptD_.push_back(constituents_p->at(cI).pt());
    etaD_.push_back(constituents_p->at(cI).eta());
  }

  return CalcRhoFromPtEta(&ptD_, &etaD_);
}

bool rhoBuilder::CalcRhoFromPtEta(std::vector<float>* pt_p, std::vector<float>* eta_p)
{
  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  std::vector<double> ptD_, etaD_;
  for(unsigned int pI = 0; pI < pt_p->size(); ++pI){
    ptD_.push_back((*pt_p)[pI]);
    etaD_.push_back((*eta_p)[pI]);
  }
  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  return CalcRhoFromPtEta(&ptD_, &etaD_);
}

bool rhoBuilder::CalcRhoFromPtEta(std::vector<double>* pt_p, std::vector<double>* eta_p)
{
  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(!m_isInit){
    std::cout << "ERROR IN RHOBUILDER CALCRHOFROMPTETAPHI: rhoBuilder is not initialized! return false" << std::endl;
    return false;
  }

  for(unsigned int rI = 0; rI < m_rhoVals.size(); ++rI){m_rhoVals[rI] = 0.0;}
  for(unsigned int rI = 0; rI < m_rhoPtVals.size(); ++rI){m_rhoPtVals[rI] = 0.0;}
  
  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  for(unsigned int pI = 0; pI < pt_p->size(); ++pI){
    int pos = posInBins((*eta_p)[pI], &m_etaBins);
    if(pos == -1){
      std::cout << "ERROR IN RHOBUILDER CALCRHOFROMPTETAPHI: Eta \'" << (*eta_p)[pI] << "\' not found in given etabins! return false" << std::endl;
      return false;
    }

    //Calculations based on massless assumption
    double E = (*pt_p)[pI]*std::cosh((*eta_p)[pI]);
    m_rhoVals[pos] += E;
    m_rhoPtVals[pos] += (*pt_p)[pI];
  }

  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  for(unsigned int rI = 0; rI < m_rhoVals.size(); ++rI){
    if(m_doGlobalDebug){
      std::cout << "GLOBAL DEBUG FILE, LINE, RI: " << __FILE__ << ", " << __LINE__ << ", " << rI << std::endl;
      std::cout << m_etaBins.size() << std::endl;
      std::cout << m_etaBins[rI] << std::endl;
      std::cout << m_etaBins[rI+1] << std::endl;
      std::cout << M_PI << std::endl;
    }

    
    double area = 2.*M_PI*(m_etaBins[rI+1] - m_etaBins[rI]);
    m_rhoVals[rI] /= area;
    //    m_rhoPtVals[rI] /= area;
    double etaVal = (m_etaBins[rI+1] + m_etaBins[rI])/2.;
    if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE, ETAVAL: " << __FILE__ << ", " << __LINE__ << ", " << etaVal << std::endl;
    m_rhoPtVals[rI] = m_rhoVals[rI]/std::cosh(etaVal);
    if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  }
  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  return true;
}

bool rhoBuilder::SetRho(std::vector<float>* rho_p)
{
  if(!m_isInit){
    std::cout << "ERROR IN RHOBUILDER SETRHO: rhoBuilder is not initialized! return false" << std::endl;
    return false;
  }
  else if(rho_p->size() != m_rhoVals.size()){
    std::cout << "ERROR IN RHOBUILDER SETRHO: Input rho pointer is not same size as internal rho! return false" << std::endl;
    return false;
  } 
  
  for(unsigned int rI = 0; rI < m_rhoVals.size(); ++rI){
    (*rho_p)[rI] = m_rhoVals[rI];
  }  

  return true;
}

bool rhoBuilder::SetRho(std::vector<double>* rho_p)
{
  if(!m_isInit){
    std::cout << "ERROR IN RHOBUILDER SETRHO: rhoBuilder is not initialized! return false" << std::endl;
    return false;
  }
  else if(rho_p->size() != m_rhoVals.size()){
    std::cout << "ERROR IN RHOBUILDER SETRHO: Input rho pointer is not same size as internal rho! return false" << std::endl;
    return false;
  } 
  
  for(unsigned int rI = 0; rI < m_rhoVals.size(); ++rI){
    (*rho_p)[rI] = m_rhoVals[rI];
  }  
  return true;
}

bool rhoBuilder::SetRhoPt(std::vector<float>* rhoPt_p)
{
  if(!m_isInit){
    std::cout << "ERROR IN RHOBUILDER SETRHOPT: rhoBuilder is not initialized! return false" << std::endl;
    return false;
  }
  else if(rhoPt_p->size() != m_rhoPtVals.size()){
    std::cout << "ERROR IN RHOBUILDER SETRHOPT: Input rhoPt pointer is not same size as internal rhoPt! return false" << std::endl;
    return false;
  } 
  
  for(unsigned int rI = 0; rI < m_rhoPtVals.size(); ++rI){
    (*rhoPt_p)[rI] = m_rhoPtVals[rI];
  }  

  return true;
}

bool rhoBuilder::SetRhoPt(std::vector<double>* rhoPt_p)
{
  if(!m_isInit){
    std::cout << "ERROR IN RHOBUILDER SETRHOPT: rhoBuilder is not initialized! return false" << std::endl;
    return false;
  }
  else if(rhoPt_p->size() != m_rhoPtVals.size()){
    std::cout << "ERROR IN RHOBUILDER SETRHOPT: Input rhoPt pointer is not same size as internal rhoPt! return false" << std::endl;
    return false;
  } 
  
  for(unsigned int rI = 0; rI < m_rhoPtVals.size(); ++rI){
    (*rhoPt_p)[rI] = m_rhoPtVals[rI];
  }  
  return true;
}

void rhoBuilder::Clean()
{
  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  m_isInit = false;
  m_etaBins.clear();
  m_rhoVals.clear();

  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  m_doGlobalDebug = false;
  return;
}

void rhoBuilder::Print()
{
  if(!m_isInit){
    std::cout << "ERROR IN RHOBUILDER PRINT: rhoBuilder is not initialized! return" << std::endl;
    return;
  }

  std::cout << "RHOBUILDER PRINT: " << std::endl;
  std::cout << " Eta bins: ";
  for(unsigned int eI = 0; eI < m_etaBins.size()-1; ++eI){
    std::cout << m_etaBins[eI] << ", ";
  }
  std::cout << m_etaBins[m_etaBins.size()-1] << "." << std::endl;
  
  return;
}

//private member functions
int rhoBuilder::posInBins(double val, std::vector<double>* bins)
{
  int pos = -1;

  if(bins->size() == 0 || bins->size() == 1){
    std::cout << "ERROR IN RHOBUILDER POSINBINS: Given bins vector has size zero or 1. return pos -1" << std::endl;
  }
  else if(val < (*bins)[0]){
    pos = 0;    
    std::cout << "WARNING IN RHOBUILDER POSINBINS: Given val \'" << val << "\' is less than lower edge of bins vector, \'" << (*bins)[0] << "\'. return position \'" << pos << "\' corresponding to " << (*bins)[0] << "-" << (*bins)[1] << std::endl;
  }
  else if(val >= (*bins)[bins->size()-1]){
    pos = bins->size()-2;    
    std::cout << "WARNING IN RHOBUILDER POSINBINS: Given val \'" << val << "\' is greater than upper edge of bins vector, \'" << (*bins)[bins->size()] << "\'. return position \'" << pos << "\' corresponding to " << (*bins)[bins->size()-2] << "-" << (*bins)[bins->size()-1] << std::endl;
  }
  else{
    for(unsigned int bI = 0; bI < bins->size()-1; ++bI){
      if(val >= (*bins)[bI] && val < (*bins)[bI+1]){
	pos = bI;
	break;
      }
    }
  }
  
  return pos;
}
