//Author: Chris McGinn (2020.02.06)

//cpp
#include <iostream>
#include <string>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"

int validateRhoHist(std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, ".root")) return 1;
  
  std::cout << "VALIDATE RHO HIST COMPLETE. return 0" << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/validateRhoHist.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += validateRhoHist(argv[1]);
  return retVal;
}
