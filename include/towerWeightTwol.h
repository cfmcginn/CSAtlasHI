//Author: Chris McGinn (2020.02.13)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//All of this is based on athena tool PhysicsAnalysis/HeavyIonPhys/HIEventUtils/HIEventUtils/HITowerWeightTool.h

#ifndef TOWERWEIGHTTWOL_H
#define TOWERWEIGHTTWOL_H

//c+cpp
#include <iostream>
#include <string>

//ROOT
#include "TFile.h"
#include "TH3D.h"

class towerWeightTwol{
 public:
  towerWeightTwol(){};
  ~towerWeightTwol();
  
 private:
  std::string m_inFileName;
  const std::string m_histName = "h3_w"; 

  TH3D* h3_w_p = nullptr;
};

#endif
