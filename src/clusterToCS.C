//cpp
#include <iostream>
#include <string>
#include <vector>

//OMP
#include <omp.h>

//ROOT
#include "TDatime.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TTree.h"

//FASTJET
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
//#include "fastjet/tools/Subtractor.hh"
//#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

//FASTJET CONTRIB
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/contrib/IterativeConstituentSubtractor.hh"

//Local
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/cppWatch.h"
#include "include/etaPhiFunc.h"
#include "include/stringUtil.h"

//The rewrite of CS is in part a tool to help me understand better the internal workings
//Based on Marta Verweij's work for CMS, here: https://github.com/CmsHI/cmssw/blob/forest_CMSSW_10_3_1/RecoJets/JetProducers/plugins/CSJetProducer.cc

//Bad practice, keep global temporary
const double maxGlobalAbsEta = 5.0;
const double rParam = 0.4;
const fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, rParam, fastjet::E_scheme);
const double ghost_area = 0.01;
const int active_area_repeats = 1;
const fastjet::GhostedAreaSpec ghost_spec(maxGlobalAbsEta, active_area_repeats, ghost_area);
const fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, ghost_spec);
const Float_t minJtPt = 15.;
const Float_t maxJtAbsEta = 3.;
const Int_t nMaxJets = 500;
const Float_t deltaEta = 0.1;

const std::vector<std::string> baseCS = {"CSJetByJet", "CSGlobal", "CSGlobalIter"};
const std::vector<int> alphaParams = {1};

void getJetsFromParticles(std::vector<float> rho_, std::vector<float> etaBins_, std::vector<fastjet::PseudoJet> particles, std::map<std::string, std::vector<fastjet::PseudoJet> >* jets)
{
  if(jets->count("NoSub") == 0){
    std::cout << "NoSub is missing from input jets collection. return" << std::endl;
    return;
  }

  for(auto const & iter : (*jets)){
    (*jets)[iter.first].clear();
  }

  fastjet::ClusterSequenceArea csA(particles, jet_def, area_def);  
  ((*jets)["NoSub"]) = fastjet::sorted_by_pt(csA.inclusive_jets(5.));

  std::vector<fastjet::PseudoJet> globalGhosts, subtracted_particles, subtracted_particles_clean;
  subtracted_particles.reserve(particles.size());
  subtracted_particles_clean.reserve(particles.size());

  for(fastjet::PseudoJet& ijet : ((*jets)["NoSub"]) ) {
    std::vector<fastjet::PseudoJet> realConst, ghosts;
    fastjet::SelectorIsPureGhost().sift(ijet.constituents(), ghosts, realConst);

    for(fastjet::PseudoJet& ighost : ghosts){
      int ghostPos = -1;
      if(ighost.eta()<=etaBins_.at(0)) ghostPos = 0;
      else if(ighost.eta()>=etaBins_.at(etaBins_.size()-1)) ghostPos = rho_.size()-1;
      else{
	for(unsigned int ie = 0; ie < etaBins_.size()-1; ++ie){
	  if(ighost.eta()>=etaBins_.at(ie) && ighost.eta()<etaBins_.at(ie+1)){
	    ghostPos = ie;
	    break;
	  }
	}
      }
      
      double pt = (rho_.at(ghostPos))*ighost.area();
      ighost.reset_momentum_PtYPhiM(pt,ighost.rap(),ighost.phi(),ighost.m());
    }
    globalGhosts.insert(std::end(globalGhosts), std::begin(ghosts), std::end(ghosts));

    const Int_t nRealConst = realConst.size();
    if(nRealConst == 0) continue;

    for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
      fastjet::contrib::ConstituentSubtractor subtractor;
      subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
      subtractor.set_max_distance(rParam);
      subtractor.set_alpha(alphaParams[aI]);
      subtractor.set_remove_all_zero_pt_particles(true);
      subtractor.set_max_eta(maxGlobalAbsEta);
      subtracted_particles = subtractor.do_subtraction(realConst, ghosts);

      for(unsigned int pI = 0; pI < subtracted_particles.size(); ++pI){
	if(subtracted_particles[pI].pt() < 0.1) continue;	
	subtracted_particles_clean.push_back(subtracted_particles[pI]);
      }

      fastjet::PseudoJet subtracted_jet=join(subtracted_particles_clean);
      if(subtracted_jet.pt() < minJtPt) continue;
      if(TMath::Abs(ijet.eta()) >= maxJtAbsEta) continue;

      std::string jtStr = baseCS[0] + "Alpha" + std::to_string(alphaParams[aI]);
      ((*jets)[jtStr]).push_back(subtracted_jet);
      subtracted_particles_clean.clear();
    }
  }

  for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
    fastjet::contrib::ConstituentSubtractor subtractor;
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
    subtractor.set_max_distance(rParam);
    subtractor.set_alpha(alphaParams[aI]);
    subtractor.set_max_eta(maxGlobalAbsEta);
    subtractor.set_remove_all_zero_pt_particles(true);
    
    subtracted_particles = subtractor.do_subtraction(particles, globalGhosts);

    for(unsigned int pI = 0; pI < subtracted_particles.size(); ++pI){
      if(subtracted_particles[pI].pt() < 0.1) continue;

      subtracted_particles_clean.push_back(subtracted_particles[pI]);
    }
    fastjet::ClusterSequence cs(subtracted_particles_clean, jet_def);

    std::string jtStr = baseCS[1] + "Alpha" + std::to_string(alphaParams[aI]);
    ((*jets)[jtStr]) = fastjet::sorted_by_pt(cs.inclusive_jets(5.));
    subtracted_particles_clean.clear();
  }

  for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
    fastjet::contrib::IterativeConstituentSubtractor subtractor;
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
    std::vector<double> max_distances = {0.1, 0.2};
    std::vector<double> alphas = {(double)alphaParams[aI], (double)alphaParams[aI]};
    subtractor.set_parameters(max_distances, alphas);
    subtractor.set_max_eta(maxGlobalAbsEta);
    subtractor.set_remove_all_zero_pt_particles(true);
    
    subtracted_particles = subtractor.do_subtraction(particles, globalGhosts);

    for(unsigned int pI = 0; pI < subtracted_particles.size(); ++pI){
      if(subtracted_particles[pI].pt() < 0.1) continue;    
      subtracted_particles_clean.push_back(subtracted_particles[pI]);
    }
    fastjet::ClusterSequence cs(subtracted_particles_clean, jet_def);

    std::string jtStr = baseCS[2] + "Alpha" + std::to_string(alphaParams[aI]);
    ((*jets)[jtStr]) = fastjet::sorted_by_pt(cs.inclusive_jets(5.));
    subtracted_particles_clean.clear();
  }

  return;
}


int clusterToCS(std::string inFileName, std::string inATLASFileName = "", std::string caloTrackStr = "calo")
{
  if(!checkFileExt(inFileName, ".root")) return 1;
  const std::string centTableStr = "input/centrality_cuts_Gv32_proposed_RCMOD2.txt";
  if(!checkFileExt(centTableStr, "txt")) return 1;
  if(inATLASFileName.size() != 0 && !checkFileExt(inATLASFileName, ".root")) return 1;

  bool doCalo = isStrSame("calo", caloTrackStr);
  std::string treeStr = "caloTree";
  std::string ptEStr = "clusterE";
  std::string phiStr = "clusterPhi";
  std::string etaStr = "clusterEta";
  std::string truthPtStr = "akt4_truth_jet_pt";
  std::string truthEtaStr = "akt4_truth_jet_eta";
  std::string truthPhiStr = "akt4_truth_jet_phi";
  std::string truthChgStr = "";
  if(isStrSame("trk", caloTrackStr)){
    doCalo = false;

    treeStr = "bush";
    ptEStr = "trk_pt";
    phiStr = "trk_phi";
    etaStr = "trk_eta";
    
    truthPtStr = "truth_trk_pt";
    truthEtaStr = "truth_trk_eta";
    truthPhiStr = "truth_trk_phi";
    truthChgStr = "truth_trk_charge";    
  }
  else{
    std::cout << "Input caloTrackStr \'" << caloTrackStr << "\' is not valid. please input \'calo\' or \'trk\'" << std::endl;
    return 1;
  }

  const bool doATLASFile = inATLASFileName.size() != 0;
  const bool sameFileATLAS = isStrSame(inFileName, inATLASFileName);
  bool doTruth = false;
  
  if(doATLASFile){
    TFile* inFile_p = new TFile(inATLASFileName.c_str(), "READ");
    TTree* clusterTree_p = (TTree*)inFile_p->Get("bush");
    TObjArray* branchArr_p = (TObjArray*)clusterTree_p->GetListOfBranches();
    
    for(Int_t bI = 0; bI < branchArr_p->GetEntries(); ++bI){
      std::string branchStr = branchArr_p->At(bI)->GetName();
      if(isStrSame(truthPtStr, branchStr)){
	doTruth = true;
	break;
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }
  
  centralityFromInput centTable(centTableStr);

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);

  cppWatch total, preLoop, preCluster, cluster1, cluster2Sub1, cluster2Sub2, cluster2Sub3, cluster2Sub4, postCluster, postLoop;
  total.start();
  preLoop.start();

  const Int_t nPara = 10;

  std::string outFileName = inFileName.substr(0, inFileName.find(".root"));
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName = "output/" + dateStr + "/" + outFileName + "_NPara" + std::to_string(nPara) + "_CS_" + dateStr + ".root";

  Int_t run_, evt_;
  UInt_t lumi_;

  Int_t runATLAS_, evtATLAS_;
  UInt_t lumiATLAS_;

  std::vector<float>* akt4hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_eta_p=nullptr;

  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<float>* truth_chg_p=nullptr;

  Float_t fcalA_et_, fcalC_et_;
  Float_t cent_;
    
  const Int_t nJtAlgo = 4;
  std::vector<std::string> jtAlgos = {"NoSub"};
  std::vector<int> jtAlphas = {0};
  //Lets construct some algos
  for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
    for(unsigned int jI = 0; jI < baseCS.size(); ++jI){
      std::string jtStr = baseCS[jI] + "Alpha" + std::to_string(alphaParams[aI]);
      jtAlgos.push_back(jtStr);
      jtAlphas.push_back(alphaParams[aI]);
    }
  }  
  
  if(nJtAlgo != (Int_t)jtAlgos.size()){
    std::cout << "Mismatch between nJtAlgo \'" << nJtAlgo << "\' and jtAlgos.size() \'" << jtAlgos.size() << "\'. return 1" << std::endl;
    return 1;
  }
  
  Int_t njt_[nJtAlgo];
  Float_t jtpt_[nJtAlgo][nMaxJets];
  Float_t jteta_[nJtAlgo][nMaxJets];
  Float_t jtphi_[nJtAlgo][nMaxJets];
  Int_t atlasmatchpos_[nJtAlgo][nMaxJets];
  Int_t truthmatchpos_[nJtAlgo][nMaxJets];

  Int_t njtATLAS_;
  Float_t jtptATLAS_[nMaxJets];
  Float_t jtetaATLAS_[nMaxJets];
  Float_t jtphiATLAS_[nMaxJets];

  Int_t njtTruth_;
  Float_t jtptTruth_[nMaxJets];
  Float_t jtetaTruth_[nMaxJets];
  Float_t jtphiTruth_[nMaxJets];

  std::vector<float>* etaBins_p=nullptr;
  std::vector<float>* rho_p=nullptr;

  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  if(!doCalo){
    etaBins_p = new std::vector<float>;
    rho_p = new std::vector<float>;

    Int_t nEtaBins = 2*(maxJtAbsEta + rParam)/deltaEta + 2;
    for(Int_t eI = 0; eI < nEtaBins; ++eI){
      Double_t etaVal = -(maxJtAbsEta + rParam) + eI*deltaEta;
      if(TMath::Abs(etaVal) < deltaEta/2.) etaVal = 0.0;
      etaBins_p->push_back(etaVal);
      if(eI != 0) rho_p->push_back(0);
    }
  }
  
  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* clusterJetsCS_p = new TTree("clusterJetsCS", "");
  clusterJetsCS_p->Branch("run", &run_, "run/I");
  clusterJetsCS_p->Branch("lumi", &lumi_, "lumi/i");
  clusterJetsCS_p->Branch("evt", &evt_, "evt/I");

  clusterJetsCS_p->Branch("fcalA_et", &fcalA_et_, "fcalA_et/F");
  clusterJetsCS_p->Branch("fcalC_et", &fcalC_et_, "fcalC_et/F");

  clusterJetsCS_p->Branch("cent", &cent_, "cent/F");

  for(Int_t jI = 0; jI < nJtAlgo; ++jI){
    clusterJetsCS_p->Branch(("njt" + jtAlgos[jI]).c_str(), &(njt_[jI]), ("njt" + jtAlgos[jI] + "/I").c_str());
    clusterJetsCS_p->Branch(("jtpt" + jtAlgos[jI]).c_str(), jtpt_[jI], ("jtpt" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    clusterJetsCS_p->Branch(("jteta" + jtAlgos[jI]).c_str(), jteta_[jI], ("jteta" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    clusterJetsCS_p->Branch(("jtphi" + jtAlgos[jI]).c_str(), jtphi_[jI], ("jtphi" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    if(doATLASFile){
      clusterJetsCS_p->Branch(("atlasmatchpos" + jtAlgos[jI]).c_str(), atlasmatchpos_[jI], ("atlasmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());
      if(doTruth) clusterJetsCS_p->Branch(("truthmatchpos" + jtAlgos[jI]).c_str(), truthmatchpos_[jI], ("truthmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());

    }      
  }

  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  if(doATLASFile){
    clusterJetsCS_p->Branch("njtATLAS", &njtATLAS_, "njtATLAS/I");
    clusterJetsCS_p->Branch("jtptATLAS", jtptATLAS_, "jtptATLAS[njtATLAS]/F");
    clusterJetsCS_p->Branch("jtetaATLAS", jtetaATLAS_, "jtetaATLAS[njtATLAS]/F");
    clusterJetsCS_p->Branch("jtphiATLAS", jtphiATLAS_, "jtphiATLAS[njtATLAS]/F");      

    if(doTruth){
      clusterJetsCS_p->Branch("njtTruth", &njtTruth_, "njtTruth/I");
      clusterJetsCS_p->Branch("jtptTruth", jtptTruth_, "jtptTruth[njtTruth]/F");
      clusterJetsCS_p->Branch("jtetaTruth", jtetaTruth_, "jtetaTruth[njtTruth]/F");
      clusterJetsCS_p->Branch("jtphiTruth", jtphiTruth_, "jtphiTruth[njtTruth]/F");      
    }
  }

  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* inATLASFile_p = nullptr;
  TTree* atlasTree_p = nullptr;
  std::map<std::string, unsigned int> runLumiEvtStrToEntry;
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* clusterTree_p = (TTree*)inFile_p->Get(treeStr.c_str());

  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<float>* clusterE_p=nullptr;
  std::vector<float>* clusterEta_p=nullptr;
  std::vector<float>* clusterPhi_p=nullptr;

  clusterTree_p->SetBranchStatus("*", 0);
  clusterTree_p->SetBranchStatus("runNumber", 1);
  clusterTree_p->SetBranchStatus("lumiBlock", 1);
  clusterTree_p->SetBranchStatus("eventNumber", 1);
  clusterTree_p->SetBranchStatus("fcalA_et", 1);
  clusterTree_p->SetBranchStatus("fcalC_et", 1);
  clusterTree_p->SetBranchStatus(ptEStr.c_str(), 1);
  clusterTree_p->SetBranchStatus(etaStr.c_str(), 1);
  clusterTree_p->SetBranchStatus(phiStr.c_str(), 1);

  if(doCalo){
    clusterTree_p->SetBranchStatus("etaBins", 1);
    clusterTree_p->SetBranchStatus("rho", 1);
  }
  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  clusterTree_p->SetBranchAddress("runNumber", &run_);
  clusterTree_p->SetBranchAddress("lumiBlock", &lumi_);
  clusterTree_p->SetBranchAddress("eventNumber", &evt_);
  clusterTree_p->SetBranchAddress("fcalA_et", &fcalA_et_);
  clusterTree_p->SetBranchAddress("fcalC_et", &fcalC_et_);
  clusterTree_p->SetBranchAddress(ptEStr.c_str(), &clusterE_p);
  clusterTree_p->SetBranchAddress(etaStr.c_str(), &clusterEta_p);
  clusterTree_p->SetBranchAddress(phiStr.c_str(), &clusterPhi_p);
  if(doCalo){
    clusterTree_p->SetBranchAddress("etaBins", &etaBins_p);
    clusterTree_p->SetBranchAddress("rho", &rho_p);
  }
  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  if(doATLASFile){    
    if(sameFileATLAS){
      inATLASFile_p = inFile_p;
      atlasTree_p = clusterTree_p;
    }
    else{
      inATLASFile_p = new TFile(inATLASFileName.c_str(), "READ");
      atlasTree_p = (TTree*)inATLASFile_p->Get("bush");
      atlasTree_p->SetBranchStatus("*", 0);

      atlasTree_p->SetBranchStatus("runNumber", 1);
      atlasTree_p->SetBranchStatus("lumiBlock", 1);
      atlasTree_p->SetBranchStatus("eventNumber", 1);
      
      atlasTree_p->SetBranchAddress("runNumber", &runATLAS_);
      atlasTree_p->SetBranchAddress("lumiBlock", &lumiATLAS_);
      atlasTree_p->SetBranchAddress("eventNumber", &evtATLAS_);
    
      const Int_t nEntries = atlasTree_p->GetEntries();
      for(Int_t entry = 0; entry < nEntries; ++entry){
	atlasTree_p->GetEntry(entry);
	
	std::string runStr = std::to_string(runATLAS_);
	std::string lumiStr = std::to_string(lumiATLAS_);
	std::string evtStr = std::to_string(evtATLAS_);
	
	while(runStr.size() < 10){runStr = "0" + runStr;}
	while(lumiStr.size() < 10){lumiStr = "0" + lumiStr;}
	while(evtStr.size() < 10){evtStr = "0" + evtStr;}
	
	runLumiEvtStrToEntry[runStr + lumiStr + evtStr] = entry;
      }
    
      atlasTree_p->SetBranchStatus("*", 0);
    }
    
    atlasTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_pt", 1);
    atlasTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_phi", 1);
    atlasTree_p->SetBranchStatus("akt4hi_em_xcalib_jet_eta", 1);

    atlasTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_p);
    atlasTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);
    atlasTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);

    if(doTruth){
      atlasTree_p->SetBranchStatus(truthPtStr.c_str(), 1);
      atlasTree_p->SetBranchStatus(truthEtaStr.c_str(), 1);
      atlasTree_p->SetBranchStatus(truthPhiStr.c_str(), 1);
      if(truthChgStr.size() != 0) atlasTree_p->SetBranchStatus(truthChgStr.c_str(), 1);
      
      atlasTree_p->SetBranchAddress(truthPtStr.c_str(), &truth_pt_p);
      atlasTree_p->SetBranchAddress(truthPhiStr.c_str(), &truth_phi_p);
      atlasTree_p->SetBranchAddress(truthEtaStr.c_str(), &truth_eta_p);
      if(truthChgStr.size() != 0) atlasTree_p->SetBranchAddress(truthChgStr.c_str(), &truth_chg_p);
    }
  }

  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nEntries = 10000;//clusterTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/400);

  clusterJetsCS_p->Branch("run", &run_, "run/I");
  clusterJetsCS_p->Branch("lumi", &lumi_, "lumi/i");
  clusterJetsCS_p->Branch("evt", &evt_, "evt/I");

  clusterJetsCS_p->Branch("fcalA_et", &fcalA_et_, "fcalA_et/F");
  clusterJetsCS_p->Branch("fcalC_et", &fcalC_et_, "fcalC_et/F");

  clusterJetsCS_p->Branch("cent", &cent_, "cent/F");

  std::vector<Int_t> runVect;
  std::vector<UInt_t> lumiVect;
  std::vector<Int_t> evtVect;
  std::vector<float> fcalA_etVect;
  std::vector<float> fcalC_etVect;
  std::vector<float> centVect;
  std::vector<std::vector<float> > rhoVect;
  std::vector<std::vector<fastjet::PseudoJet> > particles;
  std::vector<std::map<std::string, std::vector<fastjet::PseudoJet> > > jets;
  
  std::vector<std::vector<float> > atlasPt, atlasPhi, atlasEta;
  std::vector<std::vector<float> > truthPt, truthPhi, truthEta;
  
  for(Int_t pI = 0; pI < nPara; ++pI){
    particles.push_back({});
    particles[pI].reserve(10000);

    atlasPt.push_back({});
    atlasPhi.push_back({});
    atlasEta.push_back({});

    atlasPt[pI].reserve(500);
    atlasPhi[pI].reserve(500);
    atlasEta[pI].reserve(500);

    truthPt.push_back({});
    truthPhi.push_back({});
    truthEta.push_back({});

    truthPt[pI].reserve(500);
    truthPhi[pI].reserve(500);
    truthEta[pI].reserve(500);

    jets.push_back({});
    for(Int_t jI = 0; jI < nJtAlgo; ++jI){
      jets[pI][jtAlgos[jI]] = {};
      jets[pI][jtAlgos[jI]].reserve(500);
    }
  }
  
  preLoop.stop();  
  TLorentzVector tL;
  
  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    preCluster.start();
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    clusterTree_p->GetEntry(entry);

    for(unsigned int cI = 0; cI < clusterE_p->size(); ++cI){      
      if(TMath::Abs(clusterEta_p->at(cI)) > maxJtAbsEta + rParam) continue;
      Float_t E, ET, Px, Py, Pz;
      E = clusterE_p->at(cI);
      if(doCalo){
	if(E < 0.1) continue;

	ET = E/std::cosh(clusterEta_p->at(cI));      
	if(ET < 0.1) continue;
	Px = ET*std::cos(clusterPhi_p->at(cI));
	Py = ET*std::sin(clusterPhi_p->at(cI));
	Pz = ET*std::sinh(clusterEta_p->at(cI));
      }
      else{
	tL.SetPtEtaPhiM(E, clusterEta_p->at(cI), clusterPhi_p->at(cI), 0.0);
	Px = tL.Px();
	Py = tL.Py();
	Pz = tL.Pz();
      }

      particles[entry%nPara].push_back(fastjet::PseudoJet(Px, Py, Pz, E));
    }
  
    if(!doCalo){
      for(unsigned int rI = 0; rI < rho_p->size(); ++rI){rho_p->at(rI) = 0.0;}

      for(unsigned int pI = 0; pI < particles[entry%nPara].size(); ++pI){
	Int_t etaPos = -1;

	for(unsigned int eI = 0; eI < etaBins_p->size(); ++eI){
	  if(particles[entry%nPara][pI].eta() >= etaBins_p->at(eI) && particles[entry%nPara][pI].eta() < etaBins_p->at(eI+1)){
	    etaPos = eI;
	    break;
	  }
	}

	if(etaPos == -1){
	  if(particles[entry%nPara][pI].eta() <= etaBins_p->at(etaBins_p->size()-1)) etaPos = rho_p->size()-1;
	  else std::cout << "WARNING: Particle eta \'" << particles[entry%nPara][pI].eta() << "\' not find in etaBins" << std::endl;
	}

	rho_p->at(etaPos) += particles[entry%nPara][pI].E();
      }

      //      std::cout << "RHO CHECK: " << std::endl;
      
      for(unsigned int rI = 0; rI < rho_p->size(); ++rI){
	//	std::cout << " " << rI << "/" << rho_p->size() << ": " << rho_p->at(rI) << ", ";
	rho_p->at(rI) /= (deltaEta*2.*TMath::Pi());
	//	std::cout << rho_p->at(rI) << "." << std::endl;
      }
    }
    else{
      for(unsigned int rI = 0; rI < rho_p->size(); ++rI){
	rho_p->at(rI) /= 1000.;
      }
    }

    //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    runVect.push_back(run_);
    lumiVect.push_back(lumi_);
    evtVect.push_back(evt_);
    fcalA_etVect.push_back(fcalA_et_);
    fcalC_etVect.push_back(fcalC_et_);
    
    cent_ = centTable.getCent(fcalA_et_ + fcalC_et_);
    centVect.push_back(cent_);
    rhoVect.push_back((*rho_p));

    if(inATLASFileName.size() != 0){
      int entry2 = entry;
      if(!sameFileATLAS){
	std::string runStr = std::to_string(run_);
	std::string lumiStr = std::to_string(lumi_);
	std::string evtStr = std::to_string(evt_);
	
	while(runStr.size() < 10){runStr = "0" + runStr;}
	while(lumiStr.size() < 10){lumiStr = "0" + lumiStr;}
	while(evtStr.size() < 10){evtStr = "0" + evtStr;}
	
        entry2 = -1;
	if(runLumiEvtStrToEntry.count(runStr + lumiStr + evtStr) != 0) entry2 = runLumiEvtStrToEntry[runStr + lumiStr + evtStr];
      
	if(entry < 0){
	  std::cout << "NO ATLAS EVENT MATCH for (run, lumi, evt): " << run_ << ", " << lumi_ << ", " << evt_ << std::endl;
	  continue;
	}
      }
      
      atlasTree_p->GetEntry(entry2);
      
      for(unsigned int aI = 0; aI < akt4hi_em_xcalib_jet_pt_p->size(); ++aI){
	atlasPt[entry%nPara].push_back(akt4hi_em_xcalib_jet_pt_p->at(aI));
	atlasEta[entry%nPara].push_back(akt4hi_em_xcalib_jet_eta_p->at(aI));
	atlasPhi[entry%nPara].push_back(akt4hi_em_xcalib_jet_phi_p->at(aI));
      }

      if(truthChgStr.size() != 0){
	std::vector<fastjet::PseudoJet> tempParticles;

	for(unsigned int tI = 0; tI < truth_pt_p->size(); ++tI){
	  if(TMath::Abs(truth_chg_p->at(tI)) < 0.01) continue;
	  if(truth_pt_p->at(tI) < 0.1) continue;
	  if(TMath::Abs(truth_eta_p->at(tI)) > maxJtAbsEta + rParam) continue;
	  tL.SetPtEtaPhiM(truth_pt_p->at(tI), truth_eta_p->at(tI), truth_phi_p->at(tI), 0.0);
	  if(tL.E() < 0.1) continue;

	  tempParticles.push_back(fastjet::PseudoJet(tL.Px(), tL.Py(), tL.Pz(), tL.E()));
	}

	fastjet::ClusterSequence cs(tempParticles, jet_def);
	std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets(5.));
	for(unsigned int aI = 0; aI < tempJets.size(); ++aI){
	  truthPt[entry%nPara].push_back(tempJets.at(aI).pt());
	  truthEta[entry%nPara].push_back(tempJets.at(aI).eta());
	  truthPhi[entry%nPara].push_back(tempJets.at(aI).phi_std());
	}		
      }
      else{
	for(unsigned int aI = 0; aI < truth_pt_p->size(); ++aI){
	  truthPt[entry%nPara].push_back(truth_pt_p->at(aI));
	  truthEta[entry%nPara].push_back(truth_eta_p->at(aI));
	  truthPhi[entry%nPara].push_back(truth_phi_p->at(aI));
	}			
      }
    }
    
    if((entry + 1)%nPara == 0 || entry == nEntries-1){
      preCluster.stop();
      cluster1.start();

  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

#pragma omp parallel
      {
#pragma omp for
	for(Int_t pI = 0; pI < nPara; ++pI){	  
	  getJetsFromParticles(rhoVect[pI], (*etaBins_p), particles[pI], &(jets[pI]));
	}      
      }
    
      for(Int_t pI = 0; pI < nPara; ++pI){
	for(Int_t jI = 0; jI < nJtAlgo; ++jI){
	  njt_[jI] = 0;
	
	  for(fastjet::PseudoJet& ijet : jets[pI][jtAlgos[jI]] ) {
	    if(ijet.pt() < minJtPt) continue;
	    if(TMath::Abs(ijet.eta()) >= maxJtAbsEta) continue;
	  
	    jtpt_[jI][njt_[jI]] = ijet.pt();
	    jtphi_[jI][njt_[jI]] = ijet.phi_std();
	    jteta_[jI][njt_[jI]] = ijet.eta();
	    atlasmatchpos_[jI][njt_[jI]] = -1;
	    truthmatchpos_[jI][njt_[jI]] = -1;
	    
	    ++njt_[jI];
	  }
	}	
            
	if(doATLASFile){
       	  njtATLAS_ = 0;
	  for(unsigned int aI = 0; aI < atlasPt[pI].size(); ++aI){
	    jtptATLAS_[njtATLAS_] = atlasPt[pI][aI];
	    jtetaATLAS_[njtATLAS_] = atlasEta[pI][aI];
	    jtphiATLAS_[njtATLAS_] = atlasPhi[pI][aI];
	    ++njtATLAS_;
	  }

	  if(doTruth){
	    njtTruth_ = 0;
	    for(unsigned int aI = 0; aI < truthPt[pI].size(); ++aI){
	      jtptTruth_[njtTruth_] = truthPt[pI][aI];
	      jtetaTruth_[njtTruth_] = truthEta[pI][aI];
	      jtphiTruth_[njtTruth_] = truthPhi[pI][aI];
	      ++njtTruth_;
	    }
	  }

	  for(Int_t jI = 0; jI < nJtAlgo; ++jI){
	    std::vector<bool> isMatched;
	    for(unsigned int aI = 0; aI < atlasPt[pI].size(); ++aI){
	      isMatched.push_back(false);
	    }

	    for(Int_t jI2 = 0; jI2 < njt_[jI]; ++jI2){
	      for(unsigned int aI = 0; aI < atlasPt[pI].size(); ++aI){
		if(isMatched[aI]) continue;
		if(getDR(atlasEta[pI][aI], atlasPhi[pI][aI], jteta_[jI][jI2], jtphi_[jI][jI2]) < rParam){
		  atlasmatchpos_[jI][jI2] = aI;
		  isMatched[aI] = true;
		  break;
		}
	      }
	    }

	    if(doTruth){
	      isMatched.clear();
	      for(unsigned int aI = 0; aI < truthPt[pI].size(); ++aI){
		isMatched.push_back(false);
	      }
	      
	      for(Int_t jI2 = 0; jI2 < njt_[jI]; ++jI2){
		for(unsigned int aI = 0; aI < truthPt[pI].size(); ++aI){
		  if(isMatched[aI]) continue;
		  if(getDR(truthEta[pI][aI], truthPhi[pI][aI], jteta_[jI][jI2], jtphi_[jI][jI2]) < rParam){
		    truthmatchpos_[jI][jI2] = aI;
		    isMatched[aI] = true;
		    break;
		  }
		}
	      }
	    }	   

	  }

	  atlasPt[pI].clear();
	  atlasEta[pI].clear();
	  atlasPhi[pI].clear();

	  truthPt[pI].clear();
	  truthEta[pI].clear();
	  truthPhi[pI].clear();
	}
      
	run_ = runVect[pI];
	lumi_ = lumiVect[pI];
	evt_ = evtVect[pI];
	fcalA_et_ = fcalA_etVect[pI];
	fcalC_et_ = fcalC_etVect[pI];
	cent_ = centVect[pI];
	
	clusterJetsCS_p->Fill();
      }

      runVect.clear();
      lumiVect.clear();
      evtVect.clear();
      fcalA_etVect.clear();
      fcalC_etVect.clear();
      
      centVect.clear();

      rhoVect.clear();

      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      for(Int_t pI = 0; pI < nPara; ++pI){
	particles[pI].clear();
      }
    }
  }  

  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  postLoop.start();
  
  inFile_p->Close();
  delete inFile_p;
  
  if(doATLASFile && !sameFileATLAS){
    inATLASFile_p->Close();
    delete inATLASFile_p;
  }
  
  outFile_p->cd();

  clusterJetsCS_p->Write("", TObject::kOverwrite);
  delete clusterJetsCS_p;

  TDirectoryFile* paramDir_p = (TDirectoryFile*)outFile_p->mkdir("paramDir");
  paramDir_p->cd();

  std::map<std::string, std::string> paramMap;

  paramMap["nEvents"] = std::to_string(nEntries);
  paramMap["nJtAlgo"] = std::to_string(nJtAlgo);
  std::string jtAlgoStr = "";
  for(unsigned int jI = 0; jI < jtAlgos.size(); ++jI){
    jtAlgoStr = jtAlgoStr + jtAlgos[jI] + ",";
  }
  paramMap["jtAlgos"] = jtAlgoStr;

  paramMap["caloTrackStr"] = caloTrackStr;
  paramMap["doATLAS"] = std::to_string(doATLASFile);
  paramMap["doTruth"] = std::to_string(doTruth);
  
  for(auto const& iter : paramMap){
    TNamed tempName(iter.first.c_str(), iter.second.c_str());
    tempName.Write("", TObject::kOverwrite);
  }

  paramDir_p->Close();
  delete paramDir_p;
  
  outFile_p->Close();
  delete outFile_p;

  postLoop.stop();
  total.stop();
  
  std::cout << "Timing: " << std::endl;
  std::cout << " TOTAL (WALL): " << total.totalWall() << std::endl;
  std::cout << " PRELOOP: " << preLoop.totalWall() << "/" << total.totalWall() << "=" << preLoop.totalWall()/total.totalWall()<< std::endl;
  std::cout << " CLUSTER1: " << cluster1.totalWall() << "/" << total.totalWall() << "=" << cluster1.totalWall()/total.totalWall()<< std::endl;
  std::cout << " CLUSTER2SUB1: " << cluster2Sub1.totalWall() << "/" << total.totalWall() << "=" << cluster2Sub1.totalWall()/total.totalWall()<< std::endl;
  std::cout << " CLUSTER2SUB2: " << cluster2Sub2.totalWall() << "/" << total.totalWall() << "=" << cluster2Sub2.totalWall()/total.totalWall()<< std::endl;
  std::cout << " CLUSTER2SUB3: " << cluster2Sub3.totalWall() << "/" << total.totalWall() << "=" << cluster2Sub3.totalWall()/total.totalWall()<< std::endl;
  std::cout << " CLUSTER2SUB4: " << cluster2Sub4.totalWall() << "/" << total.totalWall() << "=" << cluster2Sub4.totalWall()/total.totalWall()<< std::endl;
  std::cout << " POSTCLUSTER: " << postCluster.totalWall() << "/" << total.totalWall() << "=" << postCluster.totalWall()/total.totalWall()<< std::endl;
  std::cout << " PRECLUSTER: " << preCluster.totalWall() << "/" << total.totalWall() << "=" << preCluster.totalWall()/total.totalWall()<< std::endl;
  std::cout << " POSTLOOP: " << postLoop.totalWall() << "/" << total.totalWall() << "=" << postLoop.totalWall()/total.totalWall() << std::endl;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 4){
    std::cout << "Usage: ./bin/clusterToCS.exe <inFileName> <inATLASFileName-default=\'\'> <caloTrackStr-default=\'calo\'>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += clusterToCS(argv[1]);
  else if(argc == 3) retVal += clusterToCS(argv[1], argv[2]);
  else if(argc == 4) retVal += clusterToCS(argv[1], argv[2], argv[3]);
  return retVal;
}
