//Author: Chris McGinn (2020.01.23)

//cpp
#include <iostream>
#include <map>
#include <math.h>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

//FASTJET
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"

//FASTJET CONTRIB
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/contrib/IterativeConstituentSubtractor.hh"

//Local
#include "include/checkMakeDir.h"
#include "include/centralityFromInput.h"
#include "include/constituentBuilder.h"
#include "include/cppWatch.h"
#include "include/etaPhiFunc.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/pdgToChargeMassClass.h"
#include "include/plotUtilities.h"
#include "include/returnRootFileContentsList.h"
#include "include/rhoBuilder.h"
#include "include/sampleHandler.h"
#include "include/stringUtil.h"
#include "include/ttreeUtil.h"

bool setJet(fastjet::PseudoJet jet, Float_t* jtpt_, Float_t* jteta_, Float_t* jtphi_, Float_t ptMin, Float_t absEtaMax)
{
  if(jet.pt() < ptMin) return false;
  if(TMath::Abs(jet.eta()) >= absEtaMax) return false;

  (*jtpt_) = jet.pt();
  (*jteta_) = jet.eta();
  (*jtphi_) = jet.phi_std();
  
  return true;
}

double calcMass(fastjet::PseudoJet jet)
{
  std::vector<fastjet::PseudoJet> constituents = jet.constituents();
  double E = 0.0;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;

  for(auto const & constituent : constituents){
    E += constituent.e();
    px += constituent.px();
    py += constituent.py();
    pz += constituent.pz();
  }

  double p2 = px*px + py*py + pz*pz;
  double E2 = E*E;
  if(E2 > p2) return TMath::Sqrt(E2 - p2);
  else return TMath::Sqrt(p2 - E2);
}

bool setJet(fastjet::PseudoJet jet, Float_t* jtpt_, Float_t* jteta_, Float_t* jtphi_, Float_t* jtm_, Float_t ptMin, Float_t absEtaMax)
{
  if(jet.pt() < ptMin) return false;
  if(TMath::Abs(jet.eta()) >= absEtaMax) return false;

  (*jtpt_) = jet.pt();
  (*jteta_) = jet.eta();
  (*jtphi_) = jet.phi_std();
  (*jtm_) = calcMass(jet);
  
  return true;
}

void fillArrays(std::vector<fastjet::PseudoJet>* jets, Int_t* njt_, Float_t jtpt_[], Float_t jteta_[], Float_t jtphi_[], Float_t ptMin, Float_t absEtaMax)
{
  (*njt_) = 0;
  for(auto const & jet : (*jets)){
    if(setJet(jet, &(jtpt_[(*njt_)]), &(jteta_[(*njt_)]), &(jtphi_[(*njt_)]), ptMin, absEtaMax)) ++(*njt_);
  }
  
  return;
}

void fillArrays(std::vector<fastjet::PseudoJet>* jets, Int_t* njt_, Float_t jtpt_[], Float_t jteta_[], Float_t jtphi_[], Float_t jtm_[], Float_t ptMin, Float_t absEtaMax)
{
  (*njt_) = 0;
  for(auto const & jet : (*jets)){
    if(setJet(jet, &(jtpt_[(*njt_)]), &(jteta_[(*njt_)]), &(jtphi_[(*njt_)]), &(jtm_[(*njt_)]), ptMin, absEtaMax)) ++(*njt_);
  }
  
  return;
}

void fillArrays(std::vector<float>* jetPts_p, std::vector<float>* jetEtas_p, std::vector<float>* jetPhis_p, Int_t* njt_, Float_t jtpt_[], Float_t jteta_[], Float_t jtphi_[], Float_t ptMin, Float_t absEtaMax)
{
  std::vector<fastjet::PseudoJet> jets;
  for(unsigned int jI = 0; jI < jetPts_p->size(); ++jI){
    double E = std::cosh(((*jetEtas_p)[jI]))*((*jetPts_p)[jI]);
    double Px = std::cos(((*jetPhis_p)[jI]))*((*jetPts_p)[jI]);
    double Py = std::sin(((*jetPhis_p)[jI]))*((*jetPts_p)[jI]);
    double Pz = std::sinh(((*jetEtas_p)[jI]))*((*jetPts_p)[jI]);

    jets.push_back(fastjet::PseudoJet(Px, Py, Pz, E));
  }
  fillArrays(&jets, njt_, jtpt_, jteta_, jtphi_, ptMin, absEtaMax);
  
  return;
}

void rescaleGhosts(std::vector<float> rho_, std::vector<float> etaBins_, std::vector<fastjet::PseudoJet>* ghosts, double rescaleEtaCap = 100.)
{
  for(fastjet::PseudoJet& ighost : (*ghosts)){
    if(TMath::Abs(ighost.eta()) > rescaleEtaCap) continue;
    //    if(isinf(ighost.eta())) continue;
    //    if(TMath::Abs(ighost.eta()) > 10.) continue;

    int ghostPos = ghostEtaPos(etaBins_, ighost);

    double E = (rho_.at(ghostPos))*ighost.area();
    if(E < ighost.e() && ighost.e() < TMath::Power(10, -50)){
      //      std::cout << "Skipping ghost w/ eta \'" << ighost.eta() << "\' in bin \'" << etaBins_[ghostPos] << "-" << etaBins_[ghostPos+1] << "\', rho=\'"  << rho_[ghostPos] << "\', ghostE=" << ighost.e() << ",recalcE=" << E << "." << std::endl;
      continue; // Skip it if its gonna give a nonsense value
    }

    Double_t Et = E/std::cosh(ighost.eta());
    Double_t Px = Et*std::cos(ighost.phi_std());
    Double_t Py = Et*std::sin(ighost.phi_std());
    Double_t Pz = Et*std::sinh(ighost.eta());

    ighost.reset_momentum(Px, Py, Pz, E);
  }
  return;
}

//Main executable - should basically be the only thing here mod a few functions defined above if at all
int makeClusterTree(std::string inConfigFileName)
{
  //DEBUG BOOL FROM ENV VAR
  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  pdgToChargeMass pdgToM;

  //Timing Tools
  cppWatch total, preLoop, mainLoop, postLoop;
  std::vector<cppWatch> subMainLoop;
  unsigned int subMainLoopPos = 0;

  total.start();
  preLoop.start();

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, "config")) return 1; // Check input is valid Config file

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::string inROOTFileName = inConfig_p->GetValue("INFILENAME", "");
  std::string inCentFileName = inConfig_p->GetValue("CENTFILENAME", "");
  
  if(!check.checkFileExt(inROOTFileName, "root")) return 1; // Check input is valid ROOT file
  if(!check.checkFileExt(inCentFileName, "txt")) return 1; // Check centrality table is valid TXT file

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output"); // check output dir exists; if not create
  check.doCheckMakeDir("output/" + dateStr); // check dated output subdir exists; if not create

  centralityFromInput centTable(inCentFileName);

  //Process our config file
  const bool isMC = inConfig_p->GetValue("ISMC", 0);
  const bool doTracks = inConfig_p->GetValue("DOTRACKS", 0);
  const bool doTowers = inConfig_p->GetValue("DOTOWERS", 0);  

  const bool doIterRho = inConfig_p->GetValue("DOITERRHO", 1);  
  Int_t nIterRhoTemp = 1;
  if(doIterRho) ++nIterRhoTemp;
  const Int_t nIterRho = nIterRhoTemp;  

  const std::string trkStr = "Trk";
  const std::string towerStr = "Tower";
  
  if(!doTracks && !doTowers){//No point if we have no inputs
    std::cout << "MAKECLUSTERTREE ERROR: Input config \'" << inConfigFileName << "\' has neither doTracks nor doTowers. Please turn one on. return 1" << std::endl;
    return 1;
  }
  
  const double ghost_area = inConfig_p->GetValue("GHOSTAREA", 0.01);
  
  const double recoJtMinPt = inConfig_p->GetValue("RECOJTMINPT", 10.);
  const double genJtMinPt = inConfig_p->GetValue("GENJTMINPT", 10.);
  const double jtMaxAbsEta = inConfig_p->GetValue("JTMAXABSETA", 2.5);  
  
  TFile* inFile_p = new TFile(inROOTFileName.c_str(), "READ"); 
  TEnv* inFileConfig_p = (TEnv*)inFile_p->Get("config");
  std::string inDataSet = inFileConfig_p->GetValue("INDATASET", "");
  if(inDataSet.size() == 0){
    std::cout << "NO INDATASET FOUND. return 1" << std::endl;
    return 1;
  }

  sampleHandler sHandler;
  sHandler.Init(inDataSet);
  

  std::vector<std::string> treeList = returnRootFileContentsList(inFile_p, "TTree"); //Grab all file ttree names
  std::string treeName = "gammaJetTree_p";
  if(!vectContainsStr(treeName, &treeList)){ //Check contents contain tree we want
    std::cout << "Tree \'" << treeName << "\' is not found in file \'" << inROOTFileName << "\'. return 1" << std::endl;
    return 1;
  }
  TTree* inTree_p = (TTree*)inFile_p->Get(treeName.c_str());

  std::vector<std::string> branchList = {"runNumber",
					 "eventNumber",
					 "lumiBlock",
					 "fcalA_et",
					 "fcalC_et",
					 "akt4hi_em_xcalib_jet_pt",
					 "akt4hi_em_xcalib_jet_eta",
					 "akt4hi_em_xcalib_jet_phi"}; // Define baseline branch names we will want to use

  if(doTracks){//Append track variables according to config
    branchList.push_back("trk_pt");
    branchList.push_back("trk_eta");
    branchList.push_back("trk_phi");
    branchList.push_back("trk_tight_primary");
  }
  if(doTowers){//Append tower variables according to config
    branchList.push_back("tower_pt");
    branchList.push_back("tower_eta");
    branchList.push_back("tower_phi");
  }
  if(isMC){
    branchList.push_back("akt4_truth_jet_pt");
    branchList.push_back("akt4_truth_jet_eta");
    branchList.push_back("akt4_truth_jet_phi");

    branchList.push_back("truth_pt");
    branchList.push_back("truth_eta");
    branchList.push_back("truth_phi");
    branchList.push_back("truth_charge");
    branchList.push_back("truth_pdg");
  }
  
  if(!ttreeContainsBranches(inTree_p, &branchList)) return 1; //Test ttree contains all valid branches

  inTree_p->SetBranchStatus("*", 0); 
  for(auto const & branch : branchList){inTree_p->SetBranchStatus(branch.c_str(), 1);} //Turn back on our branches  

  Int_t runNumber, eventNumber;
  UInt_t lumiBlock;
  Float_t fcalA_et, fcalC_et;

  std::vector<float>* trk_pt_p=nullptr;
  std::vector<float>* trk_eta_p=nullptr;
  std::vector<float>* trk_phi_p=nullptr;
  std::vector<bool>* trk_tight_primary_p=nullptr;

  std::vector<float>* tower_pt_p=nullptr;
  std::vector<float>* tower_eta_p=nullptr;
  std::vector<float>* tower_phi_p=nullptr;

  std::vector<float>* akt4hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_eta_p=nullptr;
  
  std::vector<float>* akt4_truth_jet_pt_p=nullptr;
  std::vector<float>* akt4_truth_jet_eta_p=nullptr;
  std::vector<float>* akt4_truth_jet_phi_p=nullptr;

  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_charge_p=nullptr;
  std::vector<int>* truth_pdg_p=nullptr;

  inTree_p->SetBranchAddress("runNumber", &runNumber);
  inTree_p->SetBranchAddress("eventNumber", &eventNumber);
  inTree_p->SetBranchAddress("lumiBlock", &lumiBlock);
  inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
  inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);

  if(doTracks){
    inTree_p->SetBranchAddress("trk_pt", &trk_pt_p);
    inTree_p->SetBranchAddress("trk_eta", &trk_eta_p);
    inTree_p->SetBranchAddress("trk_phi", &trk_phi_p);
    inTree_p->SetBranchAddress("trk_tight_primary", &trk_tight_primary_p);
  }
  if(doTowers){
    inTree_p->SetBranchAddress("tower_pt", &tower_pt_p);
    inTree_p->SetBranchAddress("tower_eta", &tower_eta_p);
    inTree_p->SetBranchAddress("tower_phi", &tower_phi_p);
  }

  inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_p);
  inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);
  inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);
  
  if(isMC){
    inTree_p->SetBranchAddress("akt4_truth_jet_pt", &akt4_truth_jet_pt_p);
    inTree_p->SetBranchAddress("akt4_truth_jet_eta", &akt4_truth_jet_eta_p);
    inTree_p->SetBranchAddress("akt4_truth_jet_phi", &akt4_truth_jet_phi_p);

    inTree_p->SetBranchAddress("truth_pt", &truth_pt_p);
    inTree_p->SetBranchAddress("truth_eta", &truth_eta_p);
    inTree_p->SetBranchAddress("truth_phi", &truth_phi_p);
    inTree_p->SetBranchAddress("truth_charge", &truth_charge_p);
    inTree_p->SetBranchAddress("truth_pdg", &truth_pdg_p);
  }

  std::string outFileName = "output/" + dateStr + "/" + rootFileNameProc(inConfig_p->GetValue("OUTFILENAME", "outFile"), {"ISMC" + std::to_string(isMC), dateStr}); 
  
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p = new TTree("clusterJetsCS", "");

  unsigned long long sampleTag_ = sHandler.GetTag();
  Float_t xSectionNB_ = sHandler.GetXSection();
  Float_t filterEff_ = sHandler.GetFilterEff();;
  Float_t cent_;

  std::vector<float>* etaBinsOut_p=new std::vector<float>;

  std::vector<float>* trkRhoOut_p=new std::vector<float>;
  std::vector<float>* trkRhoIterOut_p=new std::vector<float>;
  std::vector<float>* trkRhoCorrOut_p=new std::vector<float>;
  std::vector<float>* trkPtRhoOut_p=new std::vector<float>;
  std::vector<float>* trkPtRhoCorrOut_p=new std::vector<float>;

  std::vector<float>* towerRhoOut_p=new std::vector<float>;
  std::vector<float>* towerRhoIterOut_p=new std::vector<float>;
  std::vector<float>* towerRhoCorrOut_p=new std::vector<float>;
  std::vector<float>* towerPtRhoOut_p=new std::vector<float>;
  std::vector<float>* towerPtRhoCorrOut_p=new std::vector<float>;
  
  //Following is set of defined params not supplied in config
  const Int_t nMaxJets = 500;
  const Int_t nMaxJtAlgo = 20; //Number of algos temp hard-coded
  const double rParam = 0.4;
  const double maxGlobalAbsEta = 5.0;
  const fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, rParam, fastjet::E_scheme);

  //Temp. hard-coded etaBins w/ 0.1 spacing
  const Double_t etaWidth = 0.1;
  etaBinsOut_p->push_back(-maxGlobalAbsEta);
  while(etaBinsOut_p->at(etaBinsOut_p->size()-1) < maxGlobalAbsEta - etaWidth/2.){
    etaBinsOut_p->push_back(etaBinsOut_p->at(etaBinsOut_p->size()-1) + etaWidth);
    trkRhoOut_p->push_back(0.0);
    trkRhoIterOut_p->push_back(0.0);
    trkRhoCorrOut_p->push_back(0.0);
    trkPtRhoOut_p->push_back(0.0);
    trkPtRhoCorrOut_p->push_back(0.0);
    towerRhoOut_p->push_back(0.0);
    towerRhoIterOut_p->push_back(0.0);
    towerRhoCorrOut_p->push_back(0.0);
    towerPtRhoOut_p->push_back(0.0);
    towerPtRhoCorrOut_p->push_back(0.0);
  }  
  
  const int active_area_repeats = 1;
  const fastjet::GhostedAreaSpec ghost_spec(maxGlobalAbsEta, active_area_repeats, ghost_area);
  const fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, ghost_spec);
  const std::vector<std::string> baseCS = {"CSJetByJet", "CSGlobal", "CSGlobalIter"};
  const std::vector<int> alphaParams = {1};
  
  std::vector<std::string> jtAlgosNom = {"NoSub", "4GeVCut"};
  std::vector<int> jtAlphasNom = {0, 0};

  std::vector<std::string> jtAlgos, towerTrackStr;
  std::map<std::string, unsigned int> algoToPosMap;
  std::vector<int> jtAlphas;

  if(doTracks) towerTrackStr.push_back(trkStr);
  if(doTowers) towerTrackStr.push_back(towerStr);

  for(unsigned int tI = 0; tI < towerTrackStr.size(); ++tI){//Lets construct some algos
    for(unsigned int jI = 0; jI < jtAlgosNom.size(); ++jI){
      if(towerTrackStr[tI].find("Tower") != std::string::npos){
	if(jtAlgosNom[jI].find("4GeVCut") != std::string::npos){
	  continue;
	}
      }

      jtAlgos.push_back(towerTrackStr[tI] + jtAlgosNom[jI]);
      jtAlphas.push_back(jtAlphasNom[jI]);
    }

    for(int rI = 0; rI < nIterRho; ++rI){
      for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
	for(unsigned int jI = 0; jI < baseCS.size(); ++jI){
	  std::string jtStr = towerTrackStr[tI] + baseCS[jI] + "Alpha" + std::to_string(alphaParams[aI]) + "IterRho" + std::to_string(rI);
	  jtAlgos.push_back(jtStr);
	  jtAlphas.push_back(alphaParams[aI]);
	}
      }
    }
  }

  //build our position map
  for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
    algoToPosMap[jtAlgos[jI]] = jI;
  }
  
  const Int_t nJtAlgo = jtAlgos.size();
  
  if(nJtAlgo > nMaxJtAlgo){
    std::cout << "MAKECLUSTERTREE: nJtAlgo \'" << nJtAlgo << "\' exceends nMaxJtAlgo \'" << nMaxJtAlgo << "\'. return 1" << std::endl;
    return 1;
  }

  if(jtAlphas.size() != jtAlgos.size()){
    std::cout << "MAKECLUSTERTREE: Mismatch between jtAlphas.size() \'" << jtAlphas.size() << "\' and jtAlgos.size() \'" << jtAlgos.size() << "\'. return 1" << std::endl;
    return 1;
  }

  constituentBuilder cBuilder; // We will need this - declare here to avoid constant resizing for memory
  rhoBuilder rBuilder(*etaBinsOut_p);
  rBuilder.Print();
  std::vector<fastjet::PseudoJet> tempInputs, tempJets, globalGhosts, globalGhostsIter, subtracted_particles, subtracted_particles_iter, realJetConst, realJetConstClean, realJetConstDirty, ghostJetConst; // Again, don't want to waste time on resizes so declare all these semi-global
  
  Int_t njt_[nMaxJtAlgo];
  Float_t jtpt_[nMaxJtAlgo][nMaxJets];
  Float_t jteta_[nMaxJtAlgo][nMaxJets];
  Float_t jtphi_[nMaxJtAlgo][nMaxJets];
  Float_t jtm_[nMaxJtAlgo][nMaxJets];
  Int_t atlasmatchpos_[nMaxJtAlgo][nMaxJets];
  Int_t truthmatchpos_[nMaxJtAlgo][nMaxJets];
  Int_t chgtruthmatchpos_[nMaxJtAlgo][nMaxJets];
  
  Int_t njtATLAS_;
  Float_t jtptATLAS_[nMaxJets];
  Float_t jtetaATLAS_[nMaxJets];
  Float_t jtphiATLAS_[nMaxJets];

  Int_t njtTruth_;
  Float_t jtptTruth_[nMaxJets];
  Float_t jtetaTruth_[nMaxJets];
  Float_t jtphiTruth_[nMaxJets];
  Float_t jtmTruth_[nMaxJets];
  Int_t jtmatchChgJtTruth_[nMaxJets];
  Int_t jtmatchposTruth_[nMaxJtAlgo][nMaxJets];

  Int_t nchgjtTruth_;
  Float_t chgjtptTruth_[nMaxJets];
  Float_t chgjtetaTruth_[nMaxJets];
  Float_t chgjtphiTruth_[nMaxJets];
  Float_t chgjtmTruth_[nMaxJets];
  Int_t chgjtmatchJtTruth_[nMaxJets];
  Int_t chgjtmatchposTruth_[nMaxJtAlgo][nMaxJets];

  outTree_p->Branch("sampleTag", &sampleTag_, "sampleTag/l");
  outTree_p->Branch("xSectionNB", &xSectionNB_, "xSectionNB/F");
  outTree_p->Branch("filterEff", &filterEff_, "filterEff/F");
  
  outTree_p->Branch("run", &runNumber, "run/I");
  outTree_p->Branch("lumi", &lumiBlock, "lumi/i");
  outTree_p->Branch("evt", &eventNumber, "evt/I");

  outTree_p->Branch("fcalA_et", &fcalA_et, "fcalA_et/F");
  outTree_p->Branch("fcalC_et", &fcalC_et, "fcalC_et/F");

  outTree_p->Branch("cent", &cent_, "cent/F");
  outTree_p->Branch("etaBins", &etaBinsOut_p);
  outTree_p->Branch("trkRho", &trkRhoOut_p);
  outTree_p->Branch("trkRhoIter", &trkRhoIterOut_p);
  outTree_p->Branch("trkRhoCorr", &trkRhoCorrOut_p);
  outTree_p->Branch("trkPtRho", &trkPtRhoOut_p);
  outTree_p->Branch("trkPtRhoCorr", &trkPtRhoCorrOut_p);
  outTree_p->Branch("towerRho", &towerRhoOut_p);
  outTree_p->Branch("towerRhoIter", &towerRhoIterOut_p);
  outTree_p->Branch("towerRhoCorr", &towerRhoCorrOut_p);
  outTree_p->Branch("towerPtRho", &towerPtRhoOut_p);
  outTree_p->Branch("towerPtRhoCorr", &towerPtRhoCorrOut_p);

  for(Int_t jI = 0; jI < nJtAlgo; ++jI){
    outTree_p->Branch(("njt" + jtAlgos[jI]).c_str(), &(njt_[jI]), ("njt" + jtAlgos[jI] + "/I").c_str());
    outTree_p->Branch(("jtpt" + jtAlgos[jI]).c_str(), jtpt_[jI], ("jtpt" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    outTree_p->Branch(("jteta" + jtAlgos[jI]).c_str(), jteta_[jI], ("jteta" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    outTree_p->Branch(("jtphi" + jtAlgos[jI]).c_str(), jtphi_[jI], ("jtphi" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    outTree_p->Branch(("jtm" + jtAlgos[jI]).c_str(), jtm_[jI], ("jtm" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    outTree_p->Branch(("atlasmatchpos" + jtAlgos[jI]).c_str(), atlasmatchpos_[jI], ("atlasmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());

    if(isMC){
      outTree_p->Branch(("truthmatchpos" + jtAlgos[jI]).c_str(), truthmatchpos_[jI], ("truthmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());
      outTree_p->Branch(("chgtruthmatchpos" + jtAlgos[jI]).c_str(), chgtruthmatchpos_[jI], ("chgtruthmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());
    }
  }
  
  outTree_p->Branch("njtATLAS", &njtATLAS_, "njtATLAS/I");
  outTree_p->Branch("jtptATLAS", jtptATLAS_, "jtptATLAS[njtATLAS]/F");
  outTree_p->Branch("jtetaATLAS", jtetaATLAS_, "jtetaATLAS[njtATLAS]/F");
  outTree_p->Branch("jtphiATLAS", jtphiATLAS_, "jtphiATLAS[njtATLAS]/F");

  if(isMC){
    outTree_p->Branch("njtTruth", &njtTruth_, "njtTruth/I");
    outTree_p->Branch("jtptTruth", jtptTruth_, "jtptTruth[njtTruth]/F");
    outTree_p->Branch("jtetaTruth", jtetaTruth_, "jtetaTruth[njtTruth]/F");
    outTree_p->Branch("jtphiTruth", jtphiTruth_, "jtphiTruth[njtTruth]/F");
    outTree_p->Branch("jtmTruth", jtmTruth_, "jtmTruth[njtTruth]/F");
    outTree_p->Branch("jtmatchChgJtTruth", jtmatchChgJtTruth_, "jtmatchChgJtTruth[njtTruth]/I");

    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      outTree_p->Branch(("jtmatchpos" + jtAlgos[aI] + "Truth").c_str(), jtmatchposTruth_[aI], ("jtmatchpos" + jtAlgos[aI] + "Truth[njtTruth]/I").c_str());      
    }

    outTree_p->Branch("nchgjtTruth", &nchgjtTruth_, "nchgjtTruth/I");
    outTree_p->Branch("chgjtptTruth", chgjtptTruth_, "chgjtptTruth[nchgjtTruth]/F");
    outTree_p->Branch("chgjtetaTruth", chgjtetaTruth_, "chgjtetaTruth[nchgjtTruth]/F");
    outTree_p->Branch("chgjtphiTruth", chgjtphiTruth_, "chgjtphiTruth[nchgjtTruth]/F");
    outTree_p->Branch("chgjtmTruth", chgjtmTruth_, "chgjtmTruth[nchgjtTruth]/F");
    outTree_p->Branch("chgjtmatchJtTruth", chgjtmatchJtTruth_, "chgjtmatchJtTruth[nchgjtTruth]/I");

    for(Int_t aI = 0; aI < nJtAlgo; ++aI){
      outTree_p->Branch(("chgjtmatchpos" + jtAlgos[aI] + "Truth").c_str(), chgjtmatchposTruth_[aI], ("chgjtmatchpos" + jtAlgos[aI] + "Truth[nchgjtTruth]/I").c_str());      
    }
  }
  
  const ULong64_t nEntries = TMath::Min((ULong64_t)200, (ULong64_t)inTree_p->GetEntries());
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::cout << "Processing " << nEntries << " TTree entries..." << std::endl;
  preLoop.stop();
  mainLoop.start();
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    subMainLoopPos = 0;
    bool doSubMain = subMainLoop.size() == 0;
    if(doSubMain) subMainLoop.push_back(cppWatch());
    subMainLoop[subMainLoopPos].start();
    
    if(entry%nDiv == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);

    if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    cent_ = centTable.GetCent(fcalA_et + fcalC_et);
  
    //Pass thru for standard ATLAS reco.
    if(doSubMain) subMainLoop.push_back(cppWatch());
    subMainLoop[subMainLoopPos].stop();
    ++subMainLoopPos;
    subMainLoop[subMainLoopPos].start();
  
    fillArrays(akt4hi_em_xcalib_jet_pt_p, akt4hi_em_xcalib_jet_eta_p, akt4hi_em_xcalib_jet_phi_p, &njtATLAS_, jtptATLAS_, jtetaATLAS_, jtphiATLAS_, recoJtMinPt, jtMaxAbsEta);
    if(isMC){
      fillArrays(akt4_truth_jet_pt_p, akt4_truth_jet_eta_p, akt4_truth_jet_phi_p, &njtTruth_, jtptTruth_, jtetaTruth_, jtphiTruth_, genJtMinPt, jtMaxAbsEta);

      for(Int_t jI = 0; jI < njtTruth_; ++jI){
	for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	  jtmatchposTruth_[aI][jI] = -1;
	}
      }

      std::vector<fastjet::PseudoJet> particles, particlesChg;
      TLorentzVector tL;
      for(unsigned int tI = 0; tI < truth_pt_p->size(); ++tI){
	if(TMath::Abs(truth_pdg_p->at(tI)) == 13) continue; //ATLAS doesn't include muons in jet reco.
	tL.SetPtEtaPhiM(truth_pt_p->at(tI), truth_eta_p->at(tI), truth_phi_p->at(tI), pdgToM.GetMassFromPDG(truth_pdg_p->at(tI)));
	particles.push_back(fastjet::PseudoJet(tL.Px(), tL.Py(), tL.Pz(), tL.E()));

	if(TMath::Abs(truth_charge_p->at(tI)) < 0.1) continue;
	particlesChg.push_back(fastjet::PseudoJet(tL.Px(), tL.Py(), tL.Pz(), tL.E()));
      }
   
      fastjet::ClusterSequence cs(particles, jet_def);
      tempJets = fastjet::sorted_by_pt(cs.inclusive_jets(genJtMinPt));
      std::vector<bool> jetUsed;
      for(unsigned int tI = 0; tI < tempJets.size(); ++tI){
	jetUsed.push_back(false);
      }

      for(Int_t jI = 0; jI < njtTruth_; ++jI){
	for(unsigned int jI2 = 0; jI2 < tempJets.size(); ++jI2){
	  if(jetUsed[jI2]) continue;

	  if(getDR(jtetaTruth_[jI], jtphiTruth_[jI], tempJets[jI2].eta(), tempJets[jI2].phi_std()) < 0.3){
	    jetUsed[jI2] = true;
	    jtmTruth_[jI] = calcMass(tempJets[jI2]);
	    break;
	  }
	}
      }

      fastjet::ClusterSequence csChg(particlesChg, jet_def);
      tempJets = fastjet::sorted_by_pt(csChg.inclusive_jets(genJtMinPt));
      
      fillArrays(&tempJets, &nchgjtTruth_, chgjtptTruth_, chgjtetaTruth_, chgjtphiTruth_, chgjtmTruth_, genJtMinPt, jtMaxAbsEta);   

      for(Int_t jI = 0; jI < nchgjtTruth_; ++jI){
	for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	  chgjtmatchposTruth_[aI][jI] = -1;
	}
      }
    }

    if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    //Now we do our re-clusters; first build our track and calorimeter tower input collections

    if(doSubMain) subMainLoop.push_back(cppWatch());
    subMainLoop[subMainLoopPos].stop();
    ++subMainLoopPos;
    subMainLoop[subMainLoopPos].start();
    
    //Reset all our arrays
    for(Int_t aI = 0; aI < nJtAlgo; ++aI){njt_[aI] = 0;}

    if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    if(doTracks){
      std::vector<std::vector<fastjet::PseudoJet > > jetsToExclude = {{}, {}};

      for(Int_t iI = 0; iI < nIterRho; ++iI){
	cBuilder.Clean();
	cBuilder.InitPtEtaPhiID(trk_pt_p, trk_eta_p, trk_phi_p, trk_tight_primary_p);
	tempInputs = cBuilder.GetAllInputs(); //No ghosted negative inputs needed for tracks, only happens w/ towers
	
	//Do no-sub - this is slow because we run ClusterSequenceArea
	fastjet::ClusterSequenceArea csA(tempInputs, jet_def, area_def);
	tempJets = fastjet::sorted_by_pt(csA.inclusive_jets(recoJtMinPt));
	std::string algo = trkStr + "NoSub";
	if(!vectContainsStr(algo, &jtAlgos)) return 1;
	unsigned int algoPos = algoToPosMap[algo];
	
	if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;	
	
	fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], jtm_[algoPos], recoJtMinPt, jtMaxAbsEta);
	
	if(doSubMain) subMainLoop.push_back(cppWatch());
	subMainLoop[subMainLoopPos].stop();
	++subMainLoopPos;
	subMainLoop[subMainLoopPos].start();
	
	//Build our globalghost collection and run jet-by-jet constituent subtraction
	globalGhosts.clear();
	globalGhostsIter.clear();
	
	if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
	//We need to build our rho
	if(iI == 0){
	  if(!rBuilder.CalcRhoFromPtEtaPhi(trk_pt_p, trk_eta_p, trk_phi_p)) return 1;
	}
	else{
	  if(!rBuilder.CalcRhoFromPtEtaPhi(trk_pt_p, trk_eta_p, trk_phi_p, &(jetsToExclude[0]), 0)) return 1;
	}

	if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	if(!rBuilder.SetRho(trkRhoOut_p)) return 1;
	if(!rBuilder.SetRhoPt(trkPtRhoOut_p)) return 1;
	
	if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
	for(const auto & jet : tempJets){
	  realJetConst.clear();
	  realJetConstClean.clear();
	  realJetConstDirty.clear();
	  ghostJetConst.clear();
	  fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghostJetConst, realJetConst);
	  
	  for(unsigned int rI = 0; rI < realJetConst.size(); ++rI){
	    if(cBuilder.IsUserIndexGhosted(realJetConst[rI].user_index())) realJetConstDirty.push_back(realJetConst[rI]);
	    else realJetConstClean.push_back(realJetConst[rI]);
	  }
	  
	  rescaleGhosts(*trkRhoOut_p, *etaBinsOut_p, &ghostJetConst, 2.5);
	  
	  globalGhosts.insert(std::end(globalGhosts), std::begin(ghostJetConst), std::end(ghostJetConst));
	  const Int_t nRealConst = realJetConstClean.size();
	  if(nRealConst == 0) continue;
	  
	  for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
	    algo = trkStr + "CSJetByJetAlpha" + std::to_string(alphaParams[aI]) + "IterRho" + std::to_string(iI);
	    if(!vectContainsStr(algo, &jtAlgos)) return 1;
	    algoPos = algoToPosMap[algo];
	    
	    fastjet::contrib::ConstituentSubtractor subtractor;
	    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
	    subtractor.set_max_distance(rParam);
	    subtractor.set_alpha(alphaParams[aI]);
	    subtractor.set_remove_all_zero_pt_particles(true);
	    subtractor.set_max_eta(maxGlobalAbsEta);
	    subtracted_particles = subtractor.do_subtraction(realJetConstClean, ghostJetConst);
	    
	    fastjet::PseudoJet subtracted_jet = join(subtracted_particles);
	    if(setJet(subtracted_jet, &(jtpt_[algoPos][njt_[algoPos]]), &(jteta_[algoPos][njt_[algoPos]]), &(jtphi_[algoPos][njt_[algoPos]]), &(jtm_[algoPos][njt_[algoPos]]), recoJtMinPt, jtMaxAbsEta)) ++(njt_[algoPos]);
	  }
	}
	
	if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
	if(doSubMain) subMainLoop.push_back(cppWatch());
	subMainLoop[subMainLoopPos].stop();
	++subMainLoopPos;
	subMainLoop[subMainLoopPos].start();
	
	cBuilder.Clean();
	cBuilder.InitPtEtaPhiID(trk_pt_p, trk_eta_p, trk_phi_p, trk_tight_primary_p, 4.0);
	tempInputs = cBuilder.GetCleanInputs();
	fastjet::ClusterSequence cs4(tempInputs, jet_def);
	tempJets = fastjet::sorted_by_pt(cs4.inclusive_jets(recoJtMinPt));
	algo = trkStr + "4GeVCut";
	if(!vectContainsStr(algo, &jtAlgos)) return 1;
	algoPos = algoToPosMap[algo];
	fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], jtm_[algoPos], recoJtMinPt, jtMaxAbsEta);      
	
	cBuilder.Clean();
	cBuilder.InitPtEtaPhiID(trk_pt_p, trk_eta_p, trk_phi_p, trk_tight_primary_p);
	tempInputs = cBuilder.GetAllInputs(); 
	for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
	  fastjet::contrib::ConstituentSubtractor subtractor;
	  subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
	  subtractor.set_max_distance(rParam);
	  subtractor.set_alpha(alphaParams[aI]);
	  subtractor.set_max_eta(maxGlobalAbsEta);
	  subtractor.set_remove_all_zero_pt_particles(true);
	  //	subtractor.set_keep_original_masses();
	  subtracted_particles = subtractor.do_subtraction(tempInputs, globalGhosts, &globalGhostsIter);
	  
	  fastjet::ClusterSequence cs(subtracted_particles, jet_def);
	  tempJets = fastjet::sorted_by_pt(cs.inclusive_jets(recoJtMinPt));
	  algo = trkStr + "CSGlobalAlpha" + std::to_string(alphaParams[aI]) + "IterRho" + std::to_string(iI);
	  if(!vectContainsStr(algo, &jtAlgos)) return 1;
	  algoPos = algoToPosMap[algo];
	  
	  fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], jtm_[algoPos], recoJtMinPt, jtMaxAbsEta);      	

	  for(unsigned int eI = 0; eI < trkRhoIterOut_p->size(); ++eI){
	    trkRhoIterOut_p->at(eI) = 0.0;
	  }
	  
	  if(!rBuilder.CalcRhoFromPseudoJet(&globalGhostsIter)) return 1;
	  if(!rBuilder.SetRho(trkRhoIterOut_p)) return 1;
	  
	  rescaleGhosts(*trkRhoIterOut_p, *etaBinsOut_p, &globalGhosts, 2.5);
	  
	  subtracted_particles_iter = subtractor.do_subtraction(subtracted_particles, globalGhosts);       	
	  
	  fastjet::ClusterSequence csIter(subtracted_particles_iter, jet_def);
	  tempJets = fastjet::sorted_by_pt(csIter.inclusive_jets(recoJtMinPt));
       	  algo = trkStr + "CSGlobalIterAlpha" + std::to_string(alphaParams[aI]) + "IterRho" + std::to_string(iI);
	  if(!vectContainsStr(algo, &jtAlgos)) return 1;
	  algoPos = algoToPosMap[algo];
	  
	  fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], jtm_[algoPos], recoJtMinPt, jtMaxAbsEta);      		
	}      
      }
    }

    if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    if(doSubMain) subMainLoop.push_back(cppWatch());
    subMainLoop[subMainLoopPos].stop();
    ++subMainLoopPos;
    subMainLoop[subMainLoopPos].start();
  
    if(doTowers){
      std::vector<std::vector<fastjet::PseudoJet > > jetsToExclude = {{}, {}};

      for(Int_t iI = 0; iI < nIterRho; ++iI){
	cBuilder.Clean();
	cBuilder.InitPtEtaPhi(tower_pt_p, tower_eta_p, tower_phi_p);
	tempInputs = cBuilder.GetAllInputs(); 

	//Do no-sub - this is slow because we run ClusterSequenceArea
	fastjet::ClusterSequenceArea csA(tempInputs, jet_def, area_def);
	tempJets = fastjet::sorted_by_pt(csA.inclusive_jets(recoJtMinPt));
	std::string algo = towerStr + "NoSub";
	if(!vectContainsStr(algo, &jtAlgos)) return 1;
	unsigned int algoPos = algoToPosMap[algo];
	
	fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], jtm_[algoPos], recoJtMinPt, jtMaxAbsEta);
	
	if(doSubMain) subMainLoop.push_back(cppWatch());
	subMainLoop[subMainLoopPos].stop();
	++subMainLoopPos;
	subMainLoop[subMainLoopPos].start();
	
	//Build our globalghost collection and run jet-by-jet constituent subtraction
	globalGhosts.clear();
	globalGhostsIter.clear();
	
	if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	
	//We need to build our rho
	if(iI == 0){
	  if(!rBuilder.CalcRhoFromPtEtaPhi(tower_pt_p, tower_eta_p, tower_phi_p)) return 1;
	}
	else{
	  if(!rBuilder.CalcRhoFromPtEtaPhi(tower_pt_p, tower_eta_p, tower_phi_p, &(jetsToExclude[0]), 1)) return 1;
	}

	if(!rBuilder.SetRho(towerRhoOut_p)) return 1;
	if(!rBuilder.SetRhoPt(towerPtRhoOut_p)) return 1;
	
	
	for(const auto & jet : tempJets){
	  realJetConst.clear();
	  realJetConstClean.clear();
	  realJetConstDirty.clear();
	  ghostJetConst.clear();
	  fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghostJetConst, realJetConst);
	  
	  for(unsigned int rI = 0; rI < realJetConst.size(); ++rI){
	    if(cBuilder.IsUserIndexGhosted(realJetConst[rI].user_index())) realJetConstDirty.push_back(realJetConst[rI]);
	    else realJetConstClean.push_back(realJetConst[rI]);
	  }
	  
	  rescaleGhosts(*towerRhoOut_p, *etaBinsOut_p, &ghostJetConst, 5.0);
	  globalGhosts.insert(std::end(globalGhosts), std::begin(ghostJetConst), std::end(ghostJetConst));
	  const Int_t nRealConst = realJetConstClean.size();
	  if(nRealConst == 0) continue;
	  
	  for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
	    algo = towerStr + "CSJetByJetAlpha" + std::to_string(alphaParams[aI]) + "IterRho" + std::to_string(iI);
	    if(!vectContainsStr(algo, &jtAlgos)) return 1;
	    algoPos = algoToPosMap[algo];
	    
	    fastjet::contrib::ConstituentSubtractor subtractor;
	    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
	    subtractor.set_max_distance(rParam);
	    subtractor.set_alpha(alphaParams[aI]);
	    subtractor.set_remove_all_zero_pt_particles(true);
	    subtractor.set_max_eta(maxGlobalAbsEta);
	    subtracted_particles = subtractor.do_subtraction(realJetConstClean, ghostJetConst);
	    
	    fastjet::PseudoJet subtracted_jet = join(subtracted_particles);
	    if(setJet(subtracted_jet, &(jtpt_[algoPos][njt_[algoPos]]), &(jteta_[algoPos][njt_[algoPos]]), &(jtphi_[algoPos][njt_[algoPos]]), &(jtm_[algoPos][njt_[algoPos]]), recoJtMinPt, jtMaxAbsEta)) ++(njt_[algoPos]);
	  }
	}
	
	tempInputs = cBuilder.GetAllInputs(); 
	for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
	  fastjet::contrib::ConstituentSubtractor subtractor;
	  subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
	  subtractor.set_max_distance(rParam);
	  subtractor.set_alpha(alphaParams[aI]);
	  subtractor.set_max_eta(maxGlobalAbsEta);
	  subtractor.set_remove_all_zero_pt_particles(true);
	  //	subtractor.set_keep_original_masses();
	  subtracted_particles = subtractor.do_subtraction(tempInputs, globalGhosts, &globalGhostsIter);
	  
	  fastjet::ClusterSequence cs(subtracted_particles, jet_def);
	  tempJets = fastjet::sorted_by_pt(cs.inclusive_jets(recoJtMinPt));
	  algo = towerStr + "CSGlobalAlpha" + std::to_string(alphaParams[aI]);
	  if(!vectContainsStr(algo, &jtAlgos)) return 1;
	  algoPos = algoToPosMap[algo];
	 
	  if(iI == 0){
	    for(unsigned int tI = 0; tI < tempJets.size(); ++tI){
	      jetsToExclude[0].push_back(tempJets[tI]);
	    }
	  }

	  fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], jtm_[algoPos], recoJtMinPt, jtMaxAbsEta);      	

	  
	  for(unsigned int eI = 0; eI < towerRhoIterOut_p->size(); ++eI){
	    towerRhoIterOut_p->at(eI) = 0.0;
	  }
	  
	  if(iI == 0){
	    if(!rBuilder.CalcRhoFromPseudoJet(&globalGhostsIter)) return 1;
	  }
	  else{
	    if(!rBuilder.CalcRhoFromPseudoJet(&globalGhostsIter, &(jetsToExclude[1]), 1)) return 1;
	  }
	  if(!rBuilder.SetRho(towerRhoIterOut_p)) return 1;
      	  
	  rescaleGhosts(*towerRhoIterOut_p, *etaBinsOut_p, &globalGhosts, 5.0);
	  subtracted_particles = subtractor.do_subtraction(subtracted_particles, globalGhosts);
	  fastjet::ClusterSequence csIter(subtracted_particles, jet_def);
	  tempJets = fastjet::sorted_by_pt(csIter.inclusive_jets(recoJtMinPt));
	  algo = towerStr + "CSGlobalIterAlpha" + std::to_string(alphaParams[aI]) + "IterRho" + std::to_string(iI);
	  if(!vectContainsStr(algo, &jtAlgos)) return 1;
	  algoPos = algoToPosMap[algo];
	  
	  if(iI == 0){
	    for(unsigned int tI = 0; tI < tempJets.size(); ++tI){
	      jetsToExclude[1].push_back(tempJets[tI]);
	    }
	  }

	  fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], jtm_[algoPos], recoJtMinPt, jtMaxAbsEta);      		
	}     
      } 
    }

    if(isMC){
      for(Int_t aI = 0; aI < nJtAlgo; ++aI){
	std::vector<bool> truthJetMatched, chgtruthJetMatched;
	for(Int_t jI = 0; jI < njtTruth_; ++jI){
	  truthJetMatched.push_back(false);
	  jtmatchChgJtTruth_[jI] = -1;	  
	}
	for(Int_t jI = 0; jI < nchgjtTruth_; ++jI){
	  chgtruthJetMatched.push_back(false);
	  chgjtmatchJtTruth_[jI] = -1;	  
	}

	for(Int_t jI = 0; jI < njt_[aI]; ++jI){
	  atlasmatchpos_[aI][jI] = -1;
	  truthmatchpos_[aI][jI] = -1;
	  chgtruthmatchpos_[aI][jI] = -1;
	  
	  for(Int_t jI2 = 0; jI2 < njtTruth_; ++jI2){
	    if(truthJetMatched[jI2]) continue;

	    if(getDR(jteta_[aI][jI], jtphi_[aI][jI], jtetaTruth_[jI2], jtphiTruth_[jI2]) < 0.3){
	      truthmatchpos_[aI][jI] = jI2;
	      jtmatchposTruth_[aI][jI2] = jI;

	      truthJetMatched[jI2] = true;
	      break;
	    }
	  }

	  for(Int_t jI2 = 0; jI2 < nchgjtTruth_; ++jI2){
	    if(chgtruthJetMatched[jI2]) continue;

	    if(getDR(jteta_[aI][jI], jtphi_[aI][jI], chgjtetaTruth_[jI2], chgjtphiTruth_[jI2]) < 0.3){
	      chgtruthmatchpos_[aI][jI] = jI2;
	      chgjtmatchposTruth_[aI][jI2] = jI;
	      chgtruthJetMatched[jI2] = true;
	      break;
	    }
	  }
	}
      }

      for(Int_t jI = 0; jI < njtTruth_; ++jI){
	for(Int_t jI2 = 0; jI2 < nchgjtTruth_; ++jI2){
	  if(chgjtmatchJtTruth_[jI2] >= 0) continue;

	  if(getDR(jtetaTruth_[jI], jtphiTruth_[jI], chgjtetaTruth_[jI2], chgjtphiTruth_[jI2]) < 0.3){
	    chgjtmatchJtTruth_[jI2] = jI;
	    jtmatchChgJtTruth_[jI] = jI2;
	    break;
	  }
	}	
      }
    }
       
    outTree_p->Fill();
    subMainLoop[subMainLoopPos].stop();
  }

  mainLoop.stop();
  postLoop.start();

  if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  etaBinsOut_p->clear();
  trkRhoOut_p->clear();
  trkRhoCorrOut_p->clear();
  trkRhoOut_p->clear();
  trkRhoCorrOut_p->clear();
  towerRhoOut_p->clear();
  towerRhoCorrOut_p->clear();
  towerRhoOut_p->clear();
  towerRhoCorrOut_p->clear();

  delete etaBinsOut_p;
  delete trkRhoOut_p;
  delete trkRhoCorrOut_p;
  delete trkPtRhoOut_p;
  delete trkPtRhoCorrOut_p;
  delete towerRhoOut_p;
  delete towerRhoCorrOut_p;
  delete towerPtRhoOut_p;
  delete towerPtRhoCorrOut_p;

  cBuilder.Clean();
  rBuilder.Clean();  

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;

  std::string jtAlgosStr = "";
  for(auto const & jtAlgo : jtAlgos){
    jtAlgosStr = jtAlgosStr + jtAlgo + ",";
  }
  
  inConfig_p->SetValue("NJTALGO", nJtAlgo);
  inConfig_p->SetValue("JTALGOS", jtAlgosStr.c_str());
  inConfig_p->SetValue("DATE", dateStr.c_str());
  inConfig_p->SetValue("OUTFILENAMEFINAL", outFileName.c_str());
  
  inConfig_p->Write("config", TObject::kOverwrite); 
  outFile_p->Close();
  delete outFile_p;

  delete inConfig_p;

  if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  postLoop.stop();
  total.stop();
  
  double totalWall = total.totalWall();
  double totalCPU = total.totalCPU();
  
  std::cout << "Timing report: " << nEntries << " events in " << totalWall << " seconds, or " << prettyString(totalWall/(double)nEntries, 5, false) << " seconds per event..." << std::endl;
  std::cout << " Pre-eventloop wall (s), CPU: " << preLoop.totalWall() << "/" << totalWall << ", " << preLoop.totalCPU() << "/" << totalCPU << " (" << prettyString(100*preLoop.totalCPU()/totalCPU, 2, false) << "%)" << std::endl;
  std::cout << " Main-eventloop wall (s), CPU: " << mainLoop.totalWall() << "/" << totalWall << ", " << mainLoop.totalCPU() << "/" << totalCPU << " (" << prettyString(100*mainLoop.totalCPU()/totalCPU, 2, false) << "%)" << std::endl;

  for(unsigned int sI = 0; sI < subMainLoop.size(); ++sI){
    std::cout << "POS: " << sI << std::endl;
    std::cout << subMainLoop.size() << ", " << totalWall << ", " << totalCPU << std::endl;

    std::cout << "  Submain-eventloop " << sI << "/" << subMainLoop.size() << " wall (s), CPU: " << subMainLoop[sI].totalWall() << "/" << totalWall << ", " << subMainLoop[sI].totalCPU() << "/" << totalCPU << std::endl; //" (" << prettyString(100*subMainLoop[sI].totalCPU()/totalCPU, 2, false) << "%)" << std::endl;
  }

  if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << " Post-eventloop wall (s), CPU: " << postLoop.totalWall() << "/" << totalWall << ", " << postLoop.totalCPU() << "/" << totalCPU << std::endl;//" (" << prettyString(100*postLoop.totalCPU()/totalCPU, 2, false) << "%)" << std::endl;

  if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/makeClusterTree.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;    
    return 1;
  }
  
  int retVal = 0;
  retVal += makeClusterTree(argv[1]);
  return retVal;
}
