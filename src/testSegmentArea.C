//Author: Chris McGinn (2020.07.06)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <math.h>
#include <iostream>
#include <string>

//ROOT
#include "TMath.h"
#include "TRandom3.h"

int testSegmentArea()
{
  const double circleEta = 0.1;
  const double circlePhi = 1.0;
  const double circleR = 0.4;

  const double eta1 = 0.5;
  const double eta2 = 0.55;

  TRandom3* randGen_p = new TRandom3(0);

  if(std::fabs(eta1 - circleEta) >= circleR){
    std::cout << "NO AREA OVERLAP" << std::endl;
  }
  else{
    if((circleEta - eta1 < 0 && circleEta - eta2 < 0) || (circleEta - eta1 > 0 && circleEta - eta2 > 0)){
      const double deltaEta1 = std::fabs(eta1 - circleEta);
      const double theta1 = std::acos(deltaEta1/circleR);
      
      double area = theta1*circleR*circleR - deltaEta1*circleR*std::sin(theta1);
      
      const double deltaEta2 = std::fabs(eta2 - circleEta);
      if(circleR > deltaEta2){
	const double theta2 = std::acos(deltaEta2/circleR);
	
	double areaCorrection = theta2*circleR*circleR - deltaEta2*circleR*std::sin(theta2);
	area -= areaCorrection;
      }
      std::cout << "ANALYTIC AREA: " << area << std::endl;
    }
    else{
      const double deltaEta1 = std::fabs(eta1 - circleEta);
      const double deltaEta2 = std::fabs(eta2 - circleEta);
      const double theta1 = std::acos(deltaEta1/circleR);
      const double theta2 = std::acos(deltaEta2/circleR);
      const double theta3 = TMath::Pi() - theta1 - theta2;
      
      double area = deltaEta1*circleR*std::sin(theta1) + deltaEta2*circleR*std::sin(theta2) + theta3*circleR*circleR;
      
      std::cout << "ANALYTIC AREA: " << area << std::endl;      
    }

    int areaCounter = 0;
    int totalCounter = 0;

    while(totalCounter < 10000000){
      double x = randGen_p->Uniform(circleEta - circleR, circleEta + circleR);
      double y = randGen_p->Uniform(circlePhi - circleR, circlePhi + circleR);

      ++totalCounter;
      if(std::sqrt((x - circleEta)*(x - circleEta) + (y - circlePhi)*(y - circlePhi)) < circleR){
	if(x >= eta1 && x < eta2) ++areaCounter;
      }
    }

    std::cout << "NUMERIC AREA: " << circleR*2*circleR*2*((Double_t)areaCounter)/(Double_t)totalCounter << std::endl;
  }

  delete randGen_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += testSegmentArea();
  return retVal;
}
