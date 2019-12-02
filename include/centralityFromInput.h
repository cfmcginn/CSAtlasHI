#ifndef CENTRALITYFROMINPUT_H
#define CENTRALITYFROMINPUT_H

//cpp
#include <fstream>
#include <iostream>
#include <vector>

//Local
#include "include/checkMakeDir.h"

class centralityFromInput
{
 public:
  centralityFromInput(){return;}
  centralityFromInput(std::string inTableFile);
  ~centralityFromInput(){};
  void SetTable(std::string inTableFile);

  double getCent(double inVal);

 private:
  bool isInit;
  bool isDescending;
  std::vector<double> centVals;
};

centralityFromInput::centralityFromInput(std::string inTableFile)
{
  SetTable(inTableFile);
}

void centralityFromInput::SetTable(std::string inTableFile)
{
  if(!checkFile(inTableFile) || inTableFile.find(".txt") == std::string::npos){
    std::cout << "CENTRALITYFROMINPUT: Given table \'" << inTableFile << "\' is invalid. return isInit=false" << std::endl;
    isInit = false;
    return;
  }

  //std::cout << "FILE,LINE: " << __FILE__ << "," << __LINE__ << std::endl;
  
  std::ifstream inFile(inTableFile.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    //std::cout << "FILE,LINE: " << __FILE__ << "," << __LINE__ << std::endl;
    while(tempStr.find(",") != std::string::npos){tempStr.replace(tempStr.find(","), 1, "");}
    if(tempStr.size() == 0) continue;

    centVals.push_back(std::stod(tempStr));

    //std::cout << "FILE,LINE: " << __FILE__ << "," << __LINE__ << std::endl;

    if(centVals.size() > 1){
      //std::cout << "FILE,LINE: " << __FILE__ << "," << __LINE__ << std::endl;
      if(centVals.size() == 2){
	//std::cout << "FILE,LINE: " << __FILE__ << "," << __LINE__ << std::endl;
	if(centVals[centVals.size()-1] <= centVals[centVals.size()-2]) isDescending = true;
	else isDescending = false;

	//std::cout << "FILE,LINE: " << __FILE__ << "," << __LINE__ << std::endl;

      }
      else{
	//std::cout << "FILE,LINE: " << __FILE__ << "," << __LINE__ << std::endl;
	if(isDescending && centVals[centVals.size()-1] >= centVals[centVals.size()-2]){
	  std::cout << "CENTRALITYFROMINPUT: Values in table \'" << inTableFile << "\' are not descending. return isInit=false" << std::endl;	
	  isInit = false;
	  return;
	}
	if(!isDescending && centVals[centVals.size()-1] <= centVals[centVals.size()-2]){
	  std::cout << "CENTRALITYFROMINPUT: Values in table \'" << inTableFile << "\' are not ascending. return isInit=false" << std::endl;	
	  isInit = false;
	  return;
	}
	//std::cout << "FILE,LINE: " << __FILE__ << "," << __LINE__ << std::endl;
      }
    }
  }  

  //std::cout << "FILE,LINE: " << __FILE__ << "," << __LINE__ << std::endl;

  if(centVals.size() != 101){
    std::cout << "CENTRALITYFROMINPUT: Values in table \'" << inTableFile << "\' are not N=101. return isInit=false" << std::endl;	
    isInit = false;
    return;
  }
  
  isInit = true;
  
  return;
}

double centralityFromInput::getCent(double inVal)
{
  double outVal = -1;

  if(!isInit) std::cout << "CENTRALITYFROMINPUT: Initialization failed. getCent call will return -1" << std::endl;
  else{
    for(unsigned int cI = 0; cI < centVals.size()-1; ++cI){
      if(isDescending){
	if(inVal < centVals[cI] && inVal >= centVals[cI+1]){
	  outVal = 99-cI;
	  break;
	}
      }
      else{
	if(inVal > centVals[cI] && inVal <= centVals[cI+1]){
	  outVal = 99-cI;
	  break;
	}
      }
    }
  }

  return outVal;
}

#endif
