//Author: Chris McGinn (2020.01.24)
#ifndef TTREEUTIL_H
#define TTREEUTIL_H

//cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TObjArray.h"
#include "TTree.h"

//Local
#include "include/stringUtil.h"

bool ttreeContainsBranches(TTree* inTree_p, std::vector<std::string>* listToSearch)
{
  TObjArray* listOfBranchesA = inTree_p->GetListOfBranches();

  std::vector<std::string> listOfBranchesV;
  for(Int_t bI = 0; bI < listOfBranchesA->GetEntries(); ++bI){
    listOfBranchesV.push_back(listOfBranchesA->At(bI)->GetName());
  }

  bool allFound = true;
  for(auto const & val : (*listToSearch)){
    if(!vectContainsStr(val, & listOfBranchesV)){
      std::cout << "\'" << val << "\' not found in list of branches defined by TTree \'" << inTree_p->GetName() << "\'. return false" << std::endl;
      allFound = false;
    }
  }
  
  return allFound;
}
  
#endif
