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
#include "include/configParser.h"
#include "include/constituentBuilder.h"
#include "include/cppWatch.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
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

void fillArrays(std::vector<fastjet::PseudoJet>* jets, Int_t* njt_, Float_t jtpt_[], Float_t jteta_[], Float_t jtphi_[], Float_t ptMin, Float_t absEtaMax)
{
  (*njt_) = 0;
  for(auto const & jet : (*jets)){
    if(setJet(jet, &(jtpt_[(*njt_)]), &(jteta_[(*njt_)]), &(jtphi_[(*njt_)]), ptMin, absEtaMax)) ++(*njt_);
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

  //Timing Tools
  cppWatch total, preLoop, mainLoop, postLoop;
  std::vector<cppWatch> subMainLoop;
  unsigned int subMainLoopPos = 0;

  total.start();
  preLoop.start();

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, "txt")) return 1; // Check input is valid Config file
  configParser config(inConfigFileName);//Process input conf(ig
  std::string inROOTFileName = config.GetConfigVal("INFILENAME");
  std::string inCentFileName = config.GetConfigVal("CENTFILENAME");
  
  if(!check.checkFileExt(inROOTFileName, "root")) return 1; // Check input is valid ROOT file
  if(!check.checkFileExt(inCentFileName, "txt")) return 1; // Check centrality table is valid TXT file

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output"); // check output dir exists; if not create
  check.doCheckMakeDir("output/" + dateStr); // check dated output subdir exists; if not create

  centralityFromInput centTable(inCentFileName);

  //Process our config file
  const bool isMC = std::stoi(config.GetConfigVal("ISMC"));
  const bool doTracks = std::stoi(config.GetConfigVal("DOTRACKS"));
  const bool doTowers = std::stoi(config.GetConfigVal("DOTOWERS"));  

  const std::string trkStr = "Trk";
  const std::string towerStr = "Tower";
  
  if(!doTracks && !doTowers){//No point if we have no inputs
    std::cout << "MAKECLUSTERTREE ERROR: Input config \'" << inConfigFileName << "\' has neither doTracks nor doTowers. Please turn one on. return 1" << std::endl;
    return 1;
  }
  
  const double ghost_area = std::stod(config.GetConfigVal("GHOSTAREA"));
  
  const double recoJtMinPt = std::stod(config.GetConfigVal("RECOJTMINPT"));
  const double genJtMinPt = std::stod(config.GetConfigVal("GENJTMINPT"));
  const double jtMaxAbsEta = std::stod(config.GetConfigVal("JTMAXABSETA"));  
  
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

  std::vector<float>* tower_pt_p=nullptr;
  std::vector<float>* tower_eta_p=nullptr;
  std::vector<float>* tower_phi_p=nullptr;

  std::vector<float>* akt4hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_eta_p=nullptr;
  
  std::vector<float>* akt4_truth_jet_pt_p=nullptr;
  std::vector<float>* akt4_truth_jet_eta_p=nullptr;
  std::vector<float>* akt4_truth_jet_phi_p=nullptr;

  inTree_p->SetBranchAddress("runNumber", &runNumber);
  inTree_p->SetBranchAddress("eventNumber", &eventNumber);
  inTree_p->SetBranchAddress("lumiBlock", &lumiBlock);
  inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
  inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);

  if(doTracks){
    inTree_p->SetBranchAddress("trk_pt", &trk_pt_p);
    inTree_p->SetBranchAddress("trk_eta", &trk_eta_p);
    inTree_p->SetBranchAddress("trk_phi", &trk_phi_p);
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
  }

  std::string outFileName = "output/" + dateStr + "/" + rootFileNameProc(config.GetConfigVal("OUTFILENAME"), {"ISMC" + std::to_string(isMC), dateStr}); 
  
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
  const Int_t nMaxJtAlgo = 10; //Number of algos temp hard-coded
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
      jtAlgos.push_back(towerTrackStr[tI] + jtAlgosNom[jI]);
      jtAlphas.push_back(jtAlphasNom[jI]);
    }

    for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
      for(unsigned int jI = 0; jI < baseCS.size(); ++jI){
	std::string jtStr = towerTrackStr[tI] + baseCS[jI] + "Alpha" + std::to_string(alphaParams[aI]);
	jtAlgos.push_back(jtStr);
	jtAlphas.push_back(alphaParams[aI]);
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
  std::vector<fastjet::PseudoJet> tempInputs, tempJets, globalGhosts, globalGhostsIter, subtracted_particles, realJetConst, realJetConstClean, realJetConstDirty, ghostJetConst; // Again, don't want to waste time on resizes so declare all these semi-global
  
  Int_t njt_[nMaxJtAlgo];
  Float_t jtpt_[nMaxJtAlgo][nMaxJets];
  Float_t jteta_[nMaxJtAlgo][nMaxJets];
  Float_t jtphi_[nMaxJtAlgo][nMaxJets];
  Int_t atlasmatchpos_[nMaxJtAlgo][nMaxJets];
  Int_t truthmatchpos_[nMaxJtAlgo][nMaxJets];
  
  Int_t njtATLAS_;
  Float_t jtptATLAS_[nMaxJets];
  Float_t jtetaATLAS_[nMaxJets];
  Float_t jtphiATLAS_[nMaxJets];

  Int_t njtTruth_;
  Float_t jtptTruth_[nMaxJets];
  Float_t jtetaTruth_[nMaxJets];
  Float_t jtphiTruth_[nMaxJets];

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
    outTree_p->Branch(("atlasmatchpos" + jtAlgos[jI]).c_str(), atlasmatchpos_[jI], ("atlasmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());

    if(isMC) outTree_p->Branch(("truthmatchpos" + jtAlgos[jI]).c_str(), truthmatchpos_[jI], ("truthmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());
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
  }
  
  const ULong64_t nEntries = TMath::Min((ULong64_t)1, (ULong64_t)inTree_p->GetEntries());
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
    fillArrays(akt4_truth_jet_pt_p, akt4_truth_jet_eta_p, akt4_truth_jet_phi_p, &njtTruth_, jtptTruth_, jtetaTruth_, jtphiTruth_, genJtMinPt, jtMaxAbsEta);

    if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    //Now we do our re-clusters; first build our track and calorimeter tower input collections

    if(doSubMain) subMainLoop.push_back(cppWatch());
    subMainLoop[subMainLoopPos].stop();
    ++subMainLoopPos;
    subMainLoop[subMainLoopPos].start();

    if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    //Reset all our arrays
    for(Int_t aI = 0; aI < nJtAlgo; ++aI){njt_[aI] = 0;}

    if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    if(doTracks){
      cBuilder.Clean();
      cBuilder.InitPtEtaPhi(trk_pt_p, trk_eta_p, trk_phi_p);
      tempInputs = cBuilder.GetAllInputs(); //No ghosted negative inputs needed for tracks, only happens w/ towers

      if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


      //Do no-sub - this is slow because we run ClusterSequenceArea
      fastjet::ClusterSequenceArea csA(tempInputs, jet_def, area_def);
      tempJets = fastjet::sorted_by_pt(csA.inclusive_jets(recoJtMinPt));
      std::string algo = trkStr + "NoSub";
      if(!vectContainsStr(algo, &jtAlgos)) return 1;
      unsigned int algoPos = algoToPosMap[algo];

      if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


      fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], recoJtMinPt, jtMaxAbsEta);
      
      if(doSubMain) subMainLoop.push_back(cppWatch());
      subMainLoop[subMainLoopPos].stop();
      ++subMainLoopPos;
      subMainLoop[subMainLoopPos].start();

      //Build our globalghost collection and run jet-by-jet constituent subtraction
      globalGhosts.clear();

      if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
      //We need to build our rho
      if(!rBuilder.CalcRhoFromPtEta(trk_pt_p, trk_eta_p)) return 1;
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

	std::cout << "LINE: " << __LINE__ << std::endl;
	rescaleGhosts(*trkRhoOut_p, *etaBinsOut_p, &ghostJetConst, 2.5);

	globalGhosts.insert(std::end(globalGhosts), std::begin(ghostJetConst), std::end(ghostJetConst));
	const Int_t nRealConst = realJetConstClean.size();
	if(nRealConst == 0) continue;

	for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
	  algo = trkStr + "CSJetByJetAlpha" + std::to_string(alphaParams[aI]);
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
	  if(setJet(subtracted_jet, &(jtpt_[algoPos][njt_[algoPos]]), &(jteta_[algoPos][njt_[algoPos]]), &(jtphi_[algoPos][njt_[algoPos]]), recoJtMinPt, jtMaxAbsEta)) ++(njt_[algoPos]);
	}
      }

      if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


      if(doSubMain) subMainLoop.push_back(cppWatch());
      subMainLoop[subMainLoopPos].stop();
      ++subMainLoopPos;
      subMainLoop[subMainLoopPos].start();

      std::cout << "LINE: " << __LINE__ << std::endl;
      cBuilder.Clean();
      cBuilder.InitPtEtaPhi(trk_pt_p, trk_eta_p, trk_phi_p, 4.0);
      tempInputs = cBuilder.GetCleanInputs();
      fastjet::ClusterSequence cs4(tempInputs, jet_def);
      tempJets = fastjet::sorted_by_pt(cs4.inclusive_jets(recoJtMinPt));
      algo = trkStr + "4GeVCut";
      if(!vectContainsStr(algo, &jtAlgos)) return 1;
      algoPos = algoToPosMap[algo];
      
      std::cout << "LINE: " << __LINE__ << std::endl;
      fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], recoJtMinPt, jtMaxAbsEta);      
    
      cBuilder.InitPtEtaPhi(trk_pt_p, trk_eta_p, trk_phi_p);
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
	algo = trkStr + "CSGlobalAlpha" + std::to_string(alphaParams[aI]);
	if(!vectContainsStr(algo, &jtAlgos)) return 1;
	algoPos = algoToPosMap[algo];
      
	fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], recoJtMinPt, jtMaxAbsEta);      	

	for(unsigned int eI = 0; eI < trkRhoIterOut_p->size(); ++eI){
	  trkRhoIterOut_p->at(eI) = 0.0;
	}

	for(fastjet::PseudoJet& ighost : globalGhostsIter){
	  //	  if(isinf(ighost.eta())) continue;
	  //	  if(TMath::Abs(ighost.eta()) > 10.) continue;
	  int ghostPos = ghostEtaPos(*etaBinsOut_p, ighost);
	  trkRhoIterOut_p->at(ghostPos) += ighost.E();
	}

	for(unsigned int rI = 0; rI < trkRhoIterOut_p->size(); ++rI){
	  trkRhoIterOut_p->at(rI) /= 2.*TMath::Pi()*((etaBinsOut_p->at(rI+1) - etaBinsOut_p->at(rI)));
	}
	
	rescaleGhosts(*trkRhoIterOut_p, *etaBinsOut_p, &globalGhosts, 2.5);
	subtracted_particles = subtractor.do_subtraction(subtracted_particles, globalGhosts);
	fastjet::ClusterSequence csIter(subtracted_particles, jet_def);
	tempJets = fastjet::sorted_by_pt(csIter.inclusive_jets(recoJtMinPt));
	algo = trkStr + "CSGlobalIterAlpha" + std::to_string(alphaParams[aI]);
	if(!vectContainsStr(algo, &jtAlgos)) return 1;
	algoPos = algoToPosMap[algo];

	fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], recoJtMinPt, jtMaxAbsEta);      		
      }      

      std::cout << "LINE: " << __LINE__ << std::endl;
    }

    if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    if(doSubMain) subMainLoop.push_back(cppWatch());
    subMainLoop[subMainLoopPos].stop();
    ++subMainLoopPos;
    subMainLoop[subMainLoopPos].start();
  
    if(doTowers){
      cBuilder.Clean();
      cBuilder.InitPtEtaPhi(tower_pt_p, tower_eta_p, tower_phi_p);
      tempInputs = cBuilder.GetAllInputs(); 

      //Do no-sub - this is slow because we run ClusterSequenceArea
      fastjet::ClusterSequenceArea csA(tempInputs, jet_def, area_def);
      tempJets = fastjet::sorted_by_pt(csA.inclusive_jets(recoJtMinPt));
      std::string algo = towerStr + "NoSub";
      if(!vectContainsStr(algo, &jtAlgos)) return 1;
      unsigned int algoPos = algoToPosMap[algo];

      fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], recoJtMinPt, jtMaxAbsEta);
      
      if(doSubMain) subMainLoop.push_back(cppWatch());
      subMainLoop[subMainLoopPos].stop();
      ++subMainLoopPos;
      subMainLoop[subMainLoopPos].start();

      //Build our globalghost collection and run jet-by-jet constituent subtraction
      globalGhosts.clear();

      if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
      //We need to build our rho
      if(!rBuilder.CalcRhoFromPtEta(tower_pt_p, tower_eta_p)) return 1;
      if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
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
	  algo = towerStr + "CSJetByJetAlpha" + std::to_string(alphaParams[aI]);
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
	  if(setJet(subtracted_jet, &(jtpt_[algoPos][njt_[algoPos]]), &(jteta_[algoPos][njt_[algoPos]]), &(jtphi_[algoPos][njt_[algoPos]]), recoJtMinPt, jtMaxAbsEta)) ++(njt_[algoPos]);
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
      
	fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], recoJtMinPt, jtMaxAbsEta);      	

	for(unsigned int eI = 0; eI < towerRhoIterOut_p->size(); ++eI){
	  towerRhoIterOut_p->at(eI) = 0.0;
	}

	for(fastjet::PseudoJet& ighost : globalGhostsIter){
	  //	  if(isinf(ighost.eta())) continue;
	  //	  if(TMath::Abs(ighost.eta()) > 10.) continue;
	  int ghostPos = ghostEtaPos(*etaBinsOut_p, ighost);
	  towerRhoIterOut_p->at(ghostPos) += ighost.E();
	}

	for(unsigned int rI = 0; rI < towerRhoIterOut_p->size(); ++rI){
	  towerRhoIterOut_p->at(rI) /= 2.*TMath::Pi()*((etaBinsOut_p->at(rI+1) - etaBinsOut_p->at(rI)));
	}
	
	rescaleGhosts(*towerRhoIterOut_p, *etaBinsOut_p, &globalGhosts, 5.0);
	subtracted_particles = subtractor.do_subtraction(subtracted_particles, globalGhosts);
	fastjet::ClusterSequence csIter(subtracted_particles, jet_def);
	tempJets = fastjet::sorted_by_pt(csIter.inclusive_jets(recoJtMinPt));
	algo = towerStr + "CSGlobalIterAlpha" + std::to_string(alphaParams[aI]);
	if(!vectContainsStr(algo, &jtAlgos)) return 1;
	algoPos = algoToPosMap[algo];

	fillArrays(&tempJets, &njt_[algoPos], jtpt_[algoPos], jteta_[algoPos], jtphi_[algoPos], recoJtMinPt, jtMaxAbsEta);      		
      }      

    }

    if(doGlobalDebug) std::cout << "DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
       
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
  
  std::map<std::string, std::string> configMap = config.GetConfigMap(); //grab the config in map form
  configMap["NJTALGO"] = std::to_string(nJtAlgo);
  configMap["JTALGOS"] = jtAlgosStr;
  configMap["DATE"] = "dateStr";
  configMap["OUTFILENAMEFINAL"] = outFileName;
  TEnv configEnv; // We will convert to tenv and store in file  
  for(auto const & val : configMap){
    configEnv.SetValue(val.first.c_str(), val.second.c_str()); //Fill out the map
  }
  configEnv.Write("config", TObject::kOverwrite);  
  outFile_p->Close();
  delete outFile_p;

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
