//Author: Chris McGinn (2020.02.27)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>

//Local
#include "include/sampleHandler.h"
#include "include/stringUtil.h"

sampleHandler::sampleHandler(bool in_isPP, bool in_isMC, bool in_isGamma, int in_year, int in_minPthat)
{
  Init(in_isPP, in_isMC, in_isGamma, in_year, in_minPthat);
  return;
}
 
sampleHandler::~sampleHandler()
{
  Clean();
  return;
}

void sampleHandler::PreInit()
{
  validMinPthatsByYear[2015] = {};
  validMinPthatsByYear[2017] = {0, 20, 35, 50, 60, 70, 140, 160, 280, 400, 800};
  validMinPthatsByYear[2018] = {0, 20, 50, 60, 70, 140, 160, 400, 800};

  //Inc. PYTHIA8 + Overlay
  dataSetNameToIsPP["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.merge.AOD.e4108_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.merge.AOD.e4108_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.merge.AOD.e4108_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.merge.AOD.e4108_d1516_r11439_r11217"] = 0;
  dataSetNameToIsGamma["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.merge.AOD.e4108_d1516_r11439_r11217"] = false;

  dataSetNameToIsPP["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e4108_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e4108_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e4108_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e4108_d1516_r11439_r11217"] = 20;
  dataSetNameToIsGamma["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e4108_d1516_r11439_r11217"] = false;

  dataSetNameToIsPP["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.merge.AOD.e4108_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.merge.AOD.e4108_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.merge.AOD.e4108_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.merge.AOD.e4108_d1516_r11439_r11217"] = 60;
  dataSetNameToIsGamma["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.merge.AOD.e4108_d1516_r11439_r11217"] = false;

  dataSetNameToIsPP["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.merge.AOD.e4108_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.merge.AOD.e4108_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.merge.AOD.e4108_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.merge.AOD.e4108_d1516_r11439_r11217"] = 160;
  dataSetNameToIsGamma["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.merge.AOD.e4108_d1516_r11439_r11217"] = false;

  dataSetNameToIsPP["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.merge.AOD.e4108_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.merge.AOD.e4108_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.merge.AOD.e4108_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.merge.AOD.e4108_d1516_r11439_r11217"] = 400;
  dataSetNameToIsGamma["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.merge.AOD.e4108_d1516_r11439_r11217"] = false;

  dataSetNameToIsPP["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.merge.AOD.e4108_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.merge.AOD.e4108_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.merge.AOD.e4108_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.merge.AOD.e4108_d1516_r11439_r11217"] = 800;
  dataSetNameToIsGamma["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.merge.AOD.e4108_d1516_r11439_r11217"] = false;


  //Gamma+Jet PYTHIA8 + Overlay
  dataSetNameToIsPP["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217"] = 50;
  dataSetNameToIsGamma["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_d1516_r11439_r11217"] = true;
  
  dataSetNameToIsPP["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_d1516_r11439_r11217"] = 70;
  dataSetNameToIsGamma["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_d1516_r11439_r11217"] = true;
  
  dataSetNameToIsPP["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_d1516_r11439_r11217"] = 0;
  dataSetNameToIsMC["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_d1516_r11439_r11217"] = 1;
  dataSetNameToYear["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_d1516_r11439_r11217"] = 2018;
  dataSetNameToMinPthat["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_d1516_r11439_r11217"] = 140;
  dataSetNameToIsGamma["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_d1516_r11439_r11217"] = true;

  //Inc. Jet PYTHIA8
  dataSetNameToIsPP["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.recon.AOD.e4108_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.recon.AOD.e4108_s3238_r11199"] = 0;
  dataSetNameToIsGamma["mc16_5TeV.420010.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.recon.AOD.e4108_s3238_r11199"] = false;

  dataSetNameToIsPP["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e4108_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e4108_s3238_r11199"] = 20;
  dataSetNameToIsGamma["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e4108_s3238_r11199"] = false;

  dataSetNameToIsPP["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e6608_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e6608_s3238_r11199"] = 20;
  dataSetNameToIsGamma["mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e6608_s3238_r11199"] = false;
 
  dataSetNameToIsPP["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e4108_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e4108_s3238_r11199"] = 60;
  dataSetNameToIsGamma["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e4108_s3238_r11199"] = false;

  dataSetNameToIsPP["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e6608_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e6608_s3238_r11199"] = 60;
  dataSetNameToIsGamma["mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e6608_s3238_r11199"] = false;

  dataSetNameToIsPP["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e4108_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e4108_s3238_r11199"] = 160;
  dataSetNameToIsGamma["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e4108_s3238_r11199"] = false;

  dataSetNameToIsPP["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e6608_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e6608_s3238_r11199"] = 160;
  dataSetNameToIsGamma["mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e6608_s3238_r11199"] = false;

  dataSetNameToIsPP["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e4108_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e4108_s3238_r11199"] = 400;
  dataSetNameToIsGamma["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e4108_s3238_r11199"] = false;

  dataSetNameToIsPP["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e6608_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e6608_s3238_r11199"] = 400;
  dataSetNameToIsGamma["mc16_5TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.recon.AOD.e6608_s3238_r11199"] = false;

  dataSetNameToIsPP["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e4108_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e4108_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e4108_s3238_r11199"] = 800;
  dataSetNameToIsGamma["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e4108_s3238_r11199"] = false;

  dataSetNameToIsPP["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToIsMC["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e6608_s3238_r11199"] = 1;
  dataSetNameToYear["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e6608_s3238_r11199"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e6608_s3238_r11199"] = 800;
  dataSetNameToIsGamma["mc16_5TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.recon.AOD.e6608_s3238_r11199"] = false;

  //Gamma+Jet PYTHIA8
  dataSetNameToIsPP["mc16_5TeV.423101.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP35_50.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423101.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP35_50.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423101.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP35_50.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423101.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP35_50.merge.AOD.e5094_s3238_r10441_r10210"] = 35;
  dataSetNameToIsGamma["mc16_5TeV.423101.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP35_50.merge.AOD.e5094_s3238_r10441_r10210"] = true;
 
  dataSetNameToIsPP["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_s3238_r10441_r10210"] = 50;
  dataSetNameToIsGamma["mc16_5TeV.423102.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP50_70.merge.AOD.e5094_s3238_r10441_r10210"] = true;

  dataSetNameToIsPP["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_s3238_r10441_r10210"] = 70;
  dataSetNameToIsGamma["mc16_5TeV.423103.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP70_140.merge.AOD.e5094_s3238_r10441_r10210"] = true;

  dataSetNameToIsPP["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_s3238_r10441_r10210"] = 140;
  dataSetNameToIsGamma["mc16_5TeV.423104.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP140_280.merge.AOD.e5094_s3238_r10441_r10210"] = true;

  dataSetNameToIsPP["mc16_5TeV.423105.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP280_500.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToIsMC["mc16_5TeV.423105.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP280_500.merge.AOD.e5094_s3238_r10441_r10210"] = 1;
  dataSetNameToYear["mc16_5TeV.423105.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP280_500.merge.AOD.e5094_s3238_r10441_r10210"] = 2017;
  dataSetNameToMinPthat["mc16_5TeV.423105.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP280_500.merge.AOD.e5094_s3238_r10441_r10210"] = 280;
  dataSetNameToIsGamma["mc16_5TeV.423105.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP280_500.merge.AOD.e5094_s3238_r10441_r10210"] = true;

  for(int ppI = 0; ppI < 2; ++ppI){
    for(int mcI = 0; mcI < 2; ++mcI){
      for(int gI = 0; gI < 2; ++gI){
	for(auto const & year : validYears){
	  if(ppI == 0 && year == 2017) continue;
	  if(ppI == 1 && (year == 2015 || year == 2018)) continue;

	  for(auto const & min : validMinPthatsByYear[year]){
	    unsigned long long tag = CreateTag(ppI, mcI, gI, year, min);

	    //	    std::cout << "TAG (" << ppI << ", " << mcI << ", " << gI << ", " << year << ", " << min << "): " << tag << std::endl;
	    
	    //Values as extracted in AMI
	    //PYT.8 search: mc16_5TeV.%.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP%.merge.AOD.e5094_s3238_r10441_r10210
	    //PYT.8+Overlay search: mc16_5TeV.%.Pythia8EvtGen_A14NNPDF23LO_gammajet_DP%.merge.AOD.e5094_d1516_r11439_r11217
	    
	    tagToMinPthat[tag] = min;
	  
	    if(tag == 2017011){
	      tagToXSec[tag] = 67890000;//in nanobarn
	      tagToFilterEff[tag] = 0.99713;	    
	    }
	    else if(tag == 202017011){
	      tagToXSec[tag] = 67890000;//in nanobarn
	      tagToFilterEff[tag] = .0028748;	    
	    }
	    else if(tag == 602017011){
	      tagToXSec[tag] = 639960;//in nanobarn
	      tagToFilterEff[tag] = .0042952; 
	    }
	    else if(tag == 1602017011){
	      tagToXSec[tag] = 4719.5;//in nanobarn
	      tagToFilterEff[tag] = .0052994; 
	    }
	    else if(tag == 4002017011){
	      tagToXSec[tag] = 26.602;//in nanobarn
	      tagToFilterEff[tag] = .0045901; 
	    }
	    else if(tag == 8002017011){
	      tagToXSec[tag] = .22476;//in nanobarn
	      tagToFilterEff[tag] = .0021846; 
	    }
	    else if(tag == 352017111){
	      tagToXSec[tag] = 351620;//in nanobarn
	      tagToFilterEff[tag] = 0.000029108;	    
	    }
	    else if(tag == 502017111){
	      tagToXSec[tag] = 85898;//in nanobarn
	      tagToFilterEff[tag] = 0.00003339;
	    }
	    else if(tag == 702017111){
	      tagToXSec[tag] = 21551;//in nanobarn
	      tagToFilterEff[tag] = 0.000045787;
	    }
	    else if(tag == 1402017111){
	      tagToXSec[tag] = 1044;//in nanobarn
	      tagToFilterEff[tag] = 0.000050981;
	    }
	    else if(tag == 2802017111){
	      tagToXSec[tag] = 37.592;//in nanobarn
	      tagToFilterEff[tag] = 0.000043848;
	    }
	    else if(tag == 502018110){
	      tagToXSec[tag] = 85898;//in nanobar
	      tagToFilterEff[tag] = 0.00003339;
	    }
	    else if(tag == 702018110){
	    tagToXSec[tag] = 21551;//in nanobarn
	    tagToFilterEff[tag] = 0.000045787;
	    }
	    else if(tag == 1402018110){
	      tagToXSec[tag] = 1044;//in nanobarn
	      tagToFilterEff[tag] = 0.000050981;
	    }
	    else if(tag == 2018010){
	      tagToXSec[tag] = 0.0;//in nanobar
	      tagToFilterEff[tag] = 0.0;
	    }
	    else if(tag == 202018010){
	      tagToXSec[tag] = 67890000;//in nanobar
	      tagToFilterEff[tag] = 0.0028311;
	    }
	    else if(tag == 602018010){
	      tagToXSec[tag] = 640000;//in nanobar
	      tagToFilterEff[tag] = 0.0042785;
	    }
	    else if(tag == 1602018010){
	      tagToXSec[tag] = 4719.3;//in nanobar
	      tagToFilterEff[tag] = 0.005288;
	    }
	    else if(tag == 4002018010){
	      tagToXSec[tag] = 26.602;//in nanobar
	      tagToFilterEff[tag] =  0.0045851;
	    }
	    else if(tag == 8002018010){
	      tagToXSec[tag] = 0.22476;//in nanobar
	      tagToFilterEff[tag] = 0.0021828;
	    }
	    else{
	      tagToXSec[tag] = 0.0;//in nanobarn
	      tagToFilterEff[tag] = 0.0;
	    }
	  }
	}
      }      
    }
  }

  return;
}

bool sampleHandler::Init(std::string sampleString)
{
  //Pre-Initialize the maps
  PreInit();
     
  m_isInit = false;

  if(dataSetNameToIsPP.count(sampleString) == 0){
    std::cout << "sampleHandler error - Given sample string \'" << sampleString << "\' is not valid. Initialization failed." << std::endl;
    std::cout << "To fix, pick from: " << std::endl;
    for(auto const & sample : dataSetNameToIsPP){
      std::cout << " " << sample.first << std::endl;
    }
    std::cout << "Or modify sampleHandler class. return false" << std::endl;
    Clean();    
    return m_isInit;
  }

  m_isInit = true;

  m_isPP = dataSetNameToIsPP[sampleString];
  m_isMC = dataSetNameToIsMC[sampleString];
  m_isGamma = dataSetNameToIsGamma[sampleString];
  m_year = dataSetNameToYear[sampleString];
  m_minPthat = dataSetNameToMinPthat[sampleString];

  m_tagVal = CreateTag();
  
  return m_isInit;
}

bool sampleHandler::Init(bool in_isPP, bool in_isMC, bool in_isGamma, int in_year, int in_minPthat)
{
  //Pre-Initialize the maps
  PreInit();
     
  m_isInit = false;

  if(!vectContainsInt(in_year, &validYears)){
    std::cout << "sampleHandler error - Given year val \'" << in_year << "\' is not valid. Initialization failed." << std::endl;
    std::cout << " To fix, pick from: ";
    for(auto const & year : validYears){
      std::cout << year << ", ";
    }
    std::cout << " or modify sampleHandler class. return false" << std::endl;
    Clean();
    return m_isInit;
  }
  
  if(in_isMC){
    if(!vectContainsInt(in_minPthat, &(validMinPthatsByYear[in_year]))){
      std::cout << "sampleHandler error - Given year/minPthat vals \'" << in_year << "/" << in_minPthat << "\' is not valid. Initialization failed." << std::endl;
      std::cout << " To fix, pick from: ";
      for(auto const & iter : validMinPthatsByYear){
	std::cout << iter.first << "{";
	std::string goodStuff = "";
	for(auto const & iter2 : iter.second){
	  goodStuff = goodStuff + std::to_string(iter2) + ", ";
	}
	if(goodStuff.find(",") != std::string::npos) goodStuff.replace(goodStuff.rfind(","), goodStuff.size(), "");

	std::cout << goodStuff << "}, ";
      }
      std::cout << " or modify sampleHandler class. return false" << std::endl;
      
      Clean();
      return m_isInit;
    }    
  }

  m_isInit = true;

  m_isPP = in_isPP;
  m_isMC = in_isMC;
  m_isGamma = in_isGamma;
  m_year = in_year;
  m_minPthat = in_minPthat;

  m_tagVal = CreateTag();
  
  return m_isInit;
}

unsigned long long sampleHandler::GetTag()
{
  if(!m_isInit){
    std::cout << "sampleHandler GetTag Error - Called w/o proper initialization. return tag -1" << std::endl;
  }
  
  return m_tagVal;
}

double sampleHandler::GetXSection()
{
  double xsec = 0.0;
  if(m_isInit) xsec = tagToXSec[m_tagVal];
  else std::cout << "SAMPLEHANDLER ERROR - GetXSection() called despite no init. returning 0.0" << std::endl;
  return xsec;
}

double sampleHandler::GetFilterEff()
{
  double eff = 0.0;
  if(m_isInit) eff = tagToFilterEff[m_tagVal];
  else std::cout << "SAMPLEHANDLER ERROR - GetFilterEff() called despite no init. returning 0.0" << std::endl;
  return eff;
}

int sampleHandler::GetMinPthat()
{
  int minPthat = 0;
  if(m_isInit) minPthat = tagToMinPthat[m_tagVal];
  else std::cout << "SAMPLEHANDLER ERROR - GetMinPthat() called despite no init. returning 0" << std::endl;
  return minPthat;
}

void sampleHandler::Clean()
{
  m_isInit = false;

  m_isPP = false;
  m_isMC = false;
  m_isGamma = false;
  m_year = -1;
  m_minPthat = -1;

  m_tagVal = -1;
  
  return;
}

void sampleHandler::PrintTags()
{
  std::cout << "SAMPLEHANDLER - PrintTags: " << std::endl;
  for(auto const & tag : tagToXSec){
    std::cout << " " << tag.first << ": " << tag.second << ", " << tagToFilterEff[tag.first] << std::endl;
  }
  
  return;
}

unsigned long long sampleHandler::CreateTag()
{
  unsigned long long retVal = -1;
  if(m_isInit) retVal = CreateTag(m_isPP, m_isMC, m_isGamma, m_year, m_minPthat);
  else std::cout << "sampleHandler CreateTag Error - Called w/o proper initialization. return tag -1" << std::endl;
  
  return retVal;
}

unsigned long long sampleHandler::CreateTag(bool in_isPP, bool in_isMC, bool in_isGamma, int in_year, int in_minPthat)
{
  unsigned long long isPP = (unsigned long long)in_isPP;
  unsigned long long isMC = (unsigned long long)in_isMC;
  unsigned long long isGamma = (unsigned long long)in_isGamma;
  unsigned long long year = (unsigned long long)in_year;
  unsigned long long minPthat = (unsigned long long)in_minPthat;
  

  unsigned long long retVal = isPP + isMC*10 + isGamma*100 + year*1000 + minPthat*10000000;
  
  return retVal;
}
