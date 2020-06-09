//Author: Chris McGinn (2020.06.08)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs
#ifndef SHAREDFUNCTIONS_H
#define SHAREDFUNCTIONS_H

//c+cpp
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"

inline bool checkConfigContainsParams(TEnv* inEnv_p, std::vector<std::string> params)
{
  bool retVal = true;
  for(auto const & param : params){
    std::string temp = inEnv_p->GetValue(param.c_str(), "");
    if(temp.size() == 0){
      retVal = false;
      std::cout << "CONFIG missing param \'" << param << "\'. return false" << std::endl;
    }
  }

  return retVal;
}

#endif
