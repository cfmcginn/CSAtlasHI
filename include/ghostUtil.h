#ifndef GHOSTUTIL_H
#define GHOSTUTIL_H

//cpp
#include <iostream>
#include <vector>

//fastjet
#include "fastjet/PseudoJet.hh"

inline int ghostPos(std::vector<float> bins_, double ghostVal)
{
  if(bins_.size() == 0){
    std::cout << "MAKECLUSTERTREE GHOSTPOS ERROR: Given bins have size \'0\'. returning int32 max for absurd result" << std::endl;
    return 2147483647;//lol                                                                  
  }

  int ghostPos = -1;
  if(ghostVal<=bins_.at(0)){
    std::cout << "WARNING MAKECLUSTERTREE GHOSTPOS: Given value \'" << ghostVal << "\' falls below binning (low edge \'" << bins_[0] << "\'). return pos 0" << std::endl;
    ghostPos = 0;
  }
  else if(ghostVal>=bins_.at(bins_.size()-1)){
    ghostPos = bins_.size()-2;//-2 because etabins are 1 greater than rhobins                
    if(ghostVal > bins_.at(bins_.size()-1) + (bins_.at(bins_.size()-1) - bins_.at(bins_.size()-2))) std::cout << "WARNING MAKECLUSTERTREE GHOSTPOS: Given value \'" << ghostVal << "\' is above binning (high edge \'" << bins_[bins_.size()-1] << "\'). return pos " << bins_.size()-2 << ", " << bins_[bins_.size()-2] << "-" << bins_[bins_.size()-1] << std::endl;
  }
  else{
    for(unsigned int ie = 0; ie < bins_.size()-1; ++ie){
      if(ghostVal>=bins_.at(ie) && ghostVal<bins_.at(ie+1)){
        ghostPos = ie;
        break;
      }
    }
  }
  return ghostPos;
}

inline int ghostEtaPos(std::vector<float> etaBins_, fastjet::PseudoJet ghost){return ghostPos(etaBins_, ghost.eta());}

inline int ghostPhiPos(std::vector<float> phiBins_, fastjet::PseudoJet ghost){return ghostPos(phiBins_, ghost.phi_std());}

#endif
