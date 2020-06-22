//Author: Chris McGinn (2020.01.23)

//Local
#include "include/constituentBuilder.h"

constituentBuilder::constituentBuilder(std::vector<float>* pt_p, std::vector<float>* eta_p, std::vector<float>* phi_p, double inUserMinPt)
{
  InitPtEtaPhi(pt_p, eta_p, phi_p, inUserMinPt);
  return;
}

bool constituentBuilder::InitPtEtaPhi(std::vector<float>* pt_p, std::vector<float>* eta_p, std::vector<float>* phi_p, double inUserMinPt)
{
  m_doGlobalDebug = gBug.GetDoGlobalDebug();

  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int pI = 0; pI < pt_p->size(); ++pI){//NOTE WE ARE WORKING FROM A MASSLESS ASSUMPTION FOR NOW - LIKELY DUMB IN PARTICULAR FOR TRACK JETS CHECK ATLAS STANDARD

    if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(inUserMinPt > 0 && (*pt_p)[pI] < inUserMinPt) continue;
    
    mapToOrigPt[pI] = (*pt_p)[pI];
    if((*pt_p)[pI] < minPt){
      (*pt_p)[pI] = minPt;
      mapToIsSwapped[pI] = true;
    }
    else mapToIsSwapped[pI] = false;

    if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    double E = std::cosh(((*eta_p)[pI]))*((*pt_p)[pI]);
    double Px = std::cos(((*phi_p)[pI]))*((*pt_p)[pI]);
    double Py = std::sin(((*phi_p)[pI]))*((*pt_p)[pI]);
    double Pz = std::sinh(((*eta_p)[pI]))*((*pt_p)[pI]);

    if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(mapToIsSwapped[pI]){      
      ghostedInputs.push_back(fastjet::PseudoJet(Px, Py, Pz, E));
      ghostedInputs.at(ghostedInputs.size()-1).set_user_index(pI);
    }
    else{      
      cleanInputs.push_back(fastjet::PseudoJet(Px, Py, Pz, E));
      cleanInputs.at(cleanInputs.size()-1).set_user_index(pI);
    }
  }

  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  return true;
}

bool constituentBuilder::InitPtEtaPhiID(std::vector<float>* pt_p, std::vector<float>* eta_p, std::vector<float>* phi_p, std::vector<bool>* id_p, double inUserMinPt)
{
  m_doGlobalDebug = gBug.GetDoGlobalDebug();

  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<float> pt2, eta2, phi2;
  for(unsigned int pI = 0; pI < pt_p->size(); ++pI){
    if(id_p->at(pI) == false) continue;

    pt2.push_back(pt_p->at(pI));
    eta2.push_back(eta_p->at(pI));
    phi2.push_back(phi_p->at(pI));
  }

  return InitPtEtaPhi(&pt2, &eta2, &phi2, inUserMinPt);
}


std::vector<fastjet::PseudoJet> constituentBuilder::GetCleanInputs(){return cleanInputs;}
std::vector<fastjet::PseudoJet> constituentBuilder::GetGhostedInputs(){return ghostedInputs;}

std::vector<fastjet::PseudoJet> constituentBuilder::GetAllInputs()
{
  std::vector<fastjet::PseudoJet> allInputs;
  allInputs.reserve(cleanInputs.size() + ghostedInputs.size());
  for(auto const & clean : cleanInputs){allInputs.push_back(clean);}
  for(auto const & ghosted : ghostedInputs){allInputs.push_back(ghosted);}
  return allInputs;
}

bool constituentBuilder::IsUserIndexGhosted(unsigned int inUserIndex)
{
  bool retBool = false;
  if(mapToIsSwapped.count(inUserIndex) == 0){
    std::cout << "WARNING CONSTITUENTBUILDER::IsUserIndexGhosted(): Given inUserIndex \'" << inUserIndex << "\' is not found. returning false" << std::endl;
  }
  else retBool = mapToIsSwapped[inUserIndex];
  return retBool;
}

void constituentBuilder::Clean()
{
  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  mapToOrigPt.clear();
  mapToIsSwapped.clear();
  cleanInputs.clear();
  ghostedInputs.clear();

  if(m_doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  m_doGlobalDebug = false;

  return;
}
