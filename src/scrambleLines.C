//Author: Chris McGinn (2020.07.07)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TRandom3.h"

int scrambleLines(std::string inFileName)
{
  std::ifstream inFile(inFileName.c_str());
  std::string tempStr;
  std::vector<std::string> lines, linesNew;
  while(std::getline(inFile, tempStr)){
    if(tempStr.size() == 0) continue;
    lines.push_back(tempStr);
  } 
  inFile.close();

  TRandom3* randGen_p = new TRandom3(0);

  while(lines.size() > 0){
    int val = randGen_p->Integer(lines.size());
    linesNew.push_back(lines[val]);
    lines.erase(lines.begin() + val);
  }


  std::string outFileName = inFileName.substr(0, inFileName.rfind(".")) + "_SCRAMBLED.txt";
  std::ofstream outFile(outFileName.c_str());
  for(unsigned int lI = 0; lI < linesNew.size(); ++lI){
    outFile << linesNew[lI] << std::endl;
  }
  outFile.close();

  delete randGen_p;

  std::cout << "SCRAMBLELINES COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/scrambleLines.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += scrambleLines(argv[1]);
  return retVal;
}
