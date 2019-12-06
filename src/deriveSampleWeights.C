//AUTHOR: Chris McGinn (2019.12.06)

//cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>

//Local
#include "include/checkMakeDir.h"
#include "include/stringUtil.h"

int deriveSampleWeights(std::string inXSectionFileName, std::string inEntriesFileName)
{
  if(!checkFileExt(inXSectionFileName, ".txt")) return 1;
  if(!checkFileExt(inEntriesFileName, ".txt")) return 1;

  std::map<std::string, std::vector<double> > jzMapToXSec;
  std::map<std::string, double> jzMapToEntries;
  std::map<std::string, double> jzMapToWeights;
  std::ifstream inXSectionFile(inXSectionFileName.c_str());
  std::string tempStr;
  while(std::getline(inXSectionFile, tempStr)){
    if(tempStr.size() == 0) continue;
    std::vector<std::string> tempVect = commaSepStringToVect(tempStr);
    if(tempVect.size() == 0) continue;
    if(tempVect[0].size() == 0) continue;
    if(tempVect[0].substr(0,1).find("#") != std::string::npos) continue;

    jzMapToXSec[tempVect[0]] = {std::stod(tempVect[1]), std::stod(tempVect[2])};
    jzMapToEntries[tempVect[0]] = 0.0;
  }
  inXSectionFile.close();

  std::ifstream inEntriesFile(inEntriesFileName.c_str());
  while(std::getline(inEntriesFile, tempStr)){
    if(tempStr.size() == 0) continue;
    std::vector<std::string> tempVect = commaSepStringToVect(tempStr);
    if(tempVect.size() == 0) continue;
    if(tempVect[0].size() == 0) continue;
    if(tempVect[0].substr(0,1).find("#") != std::string::npos) continue;

    jzMapToEntries[tempVect[0]] = std::stod(tempVect[1]);
  }
  inEntriesFile.close();

  
  double maxWeight = -1.0;
  for(auto const& iter : jzMapToXSec){
    if(jzMapToEntries[iter.first] <= 0) continue;

    std::cout << iter.first << ", " << iter.second[0] << ", " << iter.second[1] << ", " << jzMapToEntries[iter.first] << std::endl;

    jzMapToWeights[iter.first] = iter.second[0]*iter.second[1]/jzMapToEntries[iter.first];
    if(maxWeight < jzMapToWeights[iter.first]) maxWeight = jzMapToWeights[iter.first];

    std::cout << " " << jzMapToWeights[iter.first] << std::endl;
  }
  
  std::cout << "Maxweight: " << maxWeight << std::endl;

  std::cout << "RENORMALIZED WEIGHTS: " << std::endl;
  for(auto const& iter : jzMapToXSec){
    std::cout << iter.first << ": " << jzMapToWeights[iter.first]/maxWeight << std::endl;
  }

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/deriveSampleWeights.exe <inXSectionFileName> <inEntriesFileName>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += deriveSampleWeights(argv[1], argv[2]);
  return retVal;
}
