//cpp
#include <iostream>
#include <string>
#include <vector>

//OMP
#include <omp.h>
#include <thread>

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
#include "include/plotUtilities.h"
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
const Float_t minRhoJtPt = 15.;
const Float_t minJtPt = 10.;
const Float_t maxJtAbsEta = 3.;
const Float_t maxTrkJtAbsEta = 2.4;
const Int_t nMaxJets = 500;
const Float_t deltaEta = 0.1;

const std::vector<std::string> baseCS = {"CSJetByJet", "CSGlobal", "CSGlobalIter"};
const std::vector<int> alphaParams = {1};

int ghostPos(std::vector<float> bins_, double ghostVal)
{
  int ghostPos = -1;
  if(ghostVal<=bins_.at(0)) ghostPos = 0;
  else if(ghostVal>=bins_.at(bins_.size()-1)) ghostPos = bins_.size()-2;//-2 because etabins are 1 greater than rhobins
  else{
    for(unsigned int ie = 0; ie < bins_.size()-1; ++ie){
      if(ghostVal>=bins_.at(ie) && ghostVal<bins_.at(ie+1)){
	ghostPos = ie;
	break;
      }
    }
  }
  
  return ghostPos;
}

int ghostEtaPos(std::vector<float> etaBins_, fastjet::PseudoJet ghost){return ghostPos(etaBins_, ghost.eta());}

int ghostPhiPos(std::vector<float> phiBins_, fastjet::PseudoJet ghost){return ghostPos(phiBins_, ghost.phi_std());}

void rescaleGhosts(std::vector<float> rho_, std::vector<float> etaBins_, std::vector<fastjet::PseudoJet>* ghosts)
{
  for(fastjet::PseudoJet& ighost : (*ghosts)){
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

void getJetsFromParticles(std::vector<float> rho_, std::vector<float> etaBins_, std::vector<fastjet::PseudoJet> particles, std::map<std::string, std::vector<fastjet::PseudoJet> >* jets, std::vector<cppWatch*> cpp)
{
  cpp[0]->start();

  //  std::cout << "STARTING PARTICLES: " << std::endl;
  //  for(unsigned int pI = 0; pI < particles.size(); ++pI){
  //    std::cout << " " << pI << "/" << particles.size() << ": " << particles[pI].pt() << ", " << particles[pI].eta() << ", " << particles[pI].phi_std() << ", " << particles[pI].m() << std::endl;
  //  }
  
  if(jets->count("NoSub") == 0){
    std::cout << "NoSub is missing from input jets collection. return" << std::endl;
    return;
  }

  for(auto const & iter : (*jets)){
    (*jets)[iter.first].clear();
  }

  std::vector<fastjet::PseudoJet> particles4GeV;
  for(unsigned int pI = 0; pI < particles.size(); ++pI){
    if(particles[pI].pt() < 4.) continue;

    particles4GeV.push_back(particles[pI]);
  }

  fastjet::ClusterSequenceArea csA(particles, jet_def, area_def);
  ((*jets)["NoSub"]) = fastjet::sorted_by_pt(csA.inclusive_jets(minJtPt));

  fastjet::ClusterSequence cs4(particles4GeV, jet_def);  
  ((*jets)["4GeVCut"]) = fastjet::sorted_by_pt(cs4.inclusive_jets(minJtPt));

  std::vector<fastjet::PseudoJet> globalGhosts, globalGhostsIter, subtracted_particles, subtracted_particles_clean;
  subtracted_particles.reserve(particles.size());
  subtracted_particles_clean.reserve(particles.size());

  for(fastjet::PseudoJet& ijet : ((*jets)["NoSub"]) ) {
    std::vector<fastjet::PseudoJet> realConst, ghosts;
    fastjet::SelectorIsPureGhost().sift(ijet.constituents(), ghosts, realConst);

    rescaleGhosts(rho_, etaBins_, &ghosts);
    /*
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
    */
    globalGhosts.insert(std::end(globalGhosts), std::begin(ghosts), std::end(ghosts));

    const Int_t nRealConst = realConst.size();
    if(nRealConst == 0) continue;
 
    //    std::cout << "PURE PARTICLES (pt, eta, phi, m): " << std::endl;
    //    for(unsigned int aI = 0; aI < realConst.size(); ++aI){
      //      std::cout << " " << aI << "/" << realConst.size() << ": " << realConst[aI].pt() << ", " << realConst[aI].eta() << ", " << realConst[aI].phi_std() << ", " << realConst[aI].m() << std::endl;
    //      realConst[aI].set_user_index(aI);
    //    }

    for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
      fastjet::contrib::ConstituentSubtractor subtractor;
      subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
      subtractor.set_max_distance(rParam);
      subtractor.set_alpha(alphaParams[aI]);
      subtractor.set_remove_all_zero_pt_particles(true);
      subtractor.set_max_eta(maxGlobalAbsEta);
      subtractor.set_keep_original_masses();
      subtracted_particles = subtractor.do_subtraction(realConst, ghosts);
   
      //      std::cout << "SUBTRACTED PARTICLES (pt, eta, phi, m): " << std::endl;

      for(unsigned int pI = 0; pI < subtracted_particles.size(); ++pI){
	//	std::cout << " " << subtracted_particles[pI].user_index() << ": " << subtracted_particles[pI].pt() << ", " << subtracted_particles[pI].eta() << ", " << subtracted_particles[pI].phi_std() << ", " << subtracted_particles[pI].m() << std::endl;

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

  cpp[0]->stop();
  cpp[1]->start();

  for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
    fastjet::contrib::ConstituentSubtractor subtractor;
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
    subtractor.set_max_distance(rParam);
    subtractor.set_alpha(alphaParams[aI]);
    subtractor.set_max_eta(maxGlobalAbsEta);
    subtractor.set_remove_all_zero_pt_particles(true);
    subtractor.set_keep_original_masses();
    subtracted_particles = subtractor.do_subtraction(particles, globalGhosts, &globalGhostsIter);

    for(unsigned int pI = 0; pI < subtracted_particles.size(); ++pI){
      if(subtracted_particles[pI].pt() < 0.1) continue;

      subtracted_particles_clean.push_back(subtracted_particles[pI]);
    }
    fastjet::ClusterSequence cs(subtracted_particles_clean, jet_def);

    std::string jtStr = baseCS[1] + "Alpha" + std::to_string(alphaParams[aI]);
    ((*jets)[jtStr]) = fastjet::sorted_by_pt(cs.inclusive_jets(minJtPt));
    subtracted_particles_clean.clear();

    //Hacked iterative subtractor - use remaining ghosts to recalc rho
    for(unsigned int rI = 0; rI < rho_.size(); ++rI){
      rho_[rI] = 0.0;
    }

    //Recalc rho based on the remaining ghosts
    for(fastjet::PseudoJet& ighost : globalGhostsIter){
      int ghostPos = ghostEtaPos(etaBins_, ighost);
      rho_[ghostPos] += ighost.E();
    }

    for(unsigned int rI = 0; rI < rho_.size(); ++rI){
      rho_[rI] /= 2.*TMath::Pi()*(etaBins_[rI+1] - etaBins_[rI]);
    }

    rescaleGhosts(rho_, etaBins_, &globalGhosts);        
    subtracted_particles = subtractor.do_subtraction(subtracted_particles, globalGhosts);
    for(unsigned int pI = 0; pI < subtracted_particles.size(); ++pI){
      if(subtracted_particles[pI].pt() < 0.1) continue;

      subtracted_particles_clean.push_back(subtracted_particles[pI]);
    }
    
    fastjet::ClusterSequence csIter(subtracted_particles_clean, jet_def);

    jtStr = baseCS[2] + "Alpha" + std::to_string(alphaParams[aI]);
    ((*jets)[jtStr]) = fastjet::sorted_by_pt(csIter.inclusive_jets(minJtPt));
  }

  cpp[1]->stop();
  cpp[2]->start();

  //We will use the above 'hacked' version for constituent subtraction for now - below murders on timing, unclear why
  /*  
  for(unsigned int aI = 0; aI < alphaParams.size(); ++aI){
    cpp[3]->start();

    fastjet::contrib::IterativeConstituentSubtractor subtractor;
    subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
    std::vector<double> max_distances = {0.2, 0.4};
    std::vector<double> alphas = {(double)alphaParams[aI], (double)alphaParams[aI]};
    subtractor.set_parameters(max_distances, alphas);
    subtractor.set_max_eta(maxGlobalAbsEta);
    //    subtractor.set_remove_all_zero_pt_particles(true);
    subtractor.set_ghost_removal(true);
    
    cpp[3]->stop();
    cpp[4]->start();

    subtracted_particles = subtractor.do_subtraction(particles, globalGhosts);

    cpp[4]->stop();
    cpp[5]->start();

    for(unsigned int pI = 0; pI < subtracted_particles.size(); ++pI){
      if(subtracted_particles[pI].pt() < 0.1) continue;    
      subtracted_particles_clean.push_back(subtracted_particles[pI]);
    }
    fastjet::ClusterSequence cs(subtracted_particles_clean, jet_def);


    std::string jtStr = baseCS[2] + "Alpha" + std::to_string(alphaParams[aI]);
    ((*jets)[jtStr]) = fastjet::sorted_by_pt(cs.inclusive_jets(5.));
    subtracted_particles_clean.clear();

    cpp[5]->stop();
  }
  */
  cpp[2]->stop();

  return;
}


int clusterToCS(std::string inFileName, std::string inATLASFileName = "", std::string caloTrackStr = "calo", std::string jzStr = "")
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

  cppWatch total, preLoop, preCluster, postCluster, postLoop;
  total.start();
  preLoop.start();

  int nthreads = std::thread::hardware_concurrency();
  const Int_t nParaMax = 4;
  const Int_t nPara = TMath::Min(nParaMax, TMath::Max(nthreads/2, 1));
  cppWatch inCluster1[nParaMax];
  cppWatch inCluster2[nParaMax];
  cppWatch inCluster3[nParaMax];
  cppWatch inCluster4[nParaMax];
  cppWatch inCluster5[nParaMax];
  cppWatch inCluster6[nParaMax];
  
  std::string outFileName = inFileName.substr(0, inFileName.rfind(".root"));
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName = "output/" + dateStr + "/" + outFileName + "_NPara" + std::to_string(nPara) + "_CS_" + dateStr + ".root";

  Int_t run_, evt_;
  UInt_t lumi_;

  Int_t jzVal_ = -1;
  if(jzStr.find("JZ1") != std::string::npos) jzVal_ = 1;
  else if(jzStr.find("JZ2") != std::string::npos) jzVal_ = 2;
  else if(jzStr.find("JZ3") != std::string::npos) jzVal_ = 3;
  else if(jzStr.find("JZ4") != std::string::npos) jzVal_ = 4;
  else if(jzStr.find("JZ5") != std::string::npos) jzVal_ = 5;

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

  std::vector<float>* etaBinsOut_p=new std::vector<float>;
  std::vector<float>* rhoOut_p=new std::vector<float>;
  std::vector<float>* rhoCorrOut_p=new std::vector<float>;  

  const Int_t nJtAlgo = 5;
  std::vector<std::string> jtAlgos = {"NoSub", "4GeVCut"};
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
  Int_t truth4GeVmatchpos_[nJtAlgo][nMaxJets];

  Int_t njtATLAS_;
  Float_t jtptATLAS_[nMaxJets];
  Float_t jtetaATLAS_[nMaxJets];
  Float_t jtphiATLAS_[nMaxJets];

  Int_t njtTruth_;
  Float_t jtptTruth_[nMaxJets];
  Float_t jtetaTruth_[nMaxJets];
  Float_t jtphiTruth_[nMaxJets];
  Float_t jtchgptTruth_[nMaxJets];
  Float_t jtchgetaTruth_[nMaxJets];
  Float_t jtchgphiTruth_[nMaxJets];

  Int_t njtTruth4GeV_;
  Float_t jtptTruth4GeV_[nMaxJets];
  Float_t jtetaTruth4GeV_[nMaxJets];
  Float_t jtphiTruth4GeV_[nMaxJets];
  Float_t jtchgptTruth4GeV_[nMaxJets];
  Float_t jtchgetaTruth4GeV_[nMaxJets];
  Float_t jtchgphiTruth4GeV_[nMaxJets];

  std::vector<float>* etaBins_p=nullptr;
  std::vector<float>* phiBins_p=nullptr;
  std::vector<float>* rho_p=nullptr;
  std::vector<float>* rho2_p=nullptr;
  std::vector<std::vector<TLorentzVector> > etaPhiTower;
  std::vector<std::vector<bool> > etaPhiGoodTower;
  std::vector<int> etaNGoodTowers;
  
  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  if(!doCalo){
    etaBins_p = new std::vector<float>;
    phiBins_p = new std::vector<float>;
    rho_p = new std::vector<float>;
    rho2_p = new std::vector<float>;
    
    Int_t nEtaBins = 2*(maxJtAbsEta + rParam)/deltaEta + 2;
    for(Int_t eI = 0; eI < nEtaBins; ++eI){
      Double_t etaVal = -(maxJtAbsEta + rParam) + eI*deltaEta;
      if(TMath::Abs(etaVal) < deltaEta/2.) etaVal = 0.0;
      etaBins_p->push_back(etaVal);
      if(eI != 0){
	rho_p->push_back(0);
	rho2_p->push_back(0);
      }
    }

    Int_t nPhiBins = 72;
    for(Int_t pI = 0; pI < nPhiBins+1; ++pI){
      Double_t phiVal = -TMath::Pi() + pI*2.*TMath::Pi()/((Double_t)nPhiBins);
      phiBins_p->push_back(phiVal);
    }

    for(Int_t eI = 0; eI < nEtaBins; ++eI){
      etaPhiTower.push_back({});
      etaPhiGoodTower.push_back({});
      etaNGoodTowers.push_back(nPhiBins);

      for(Int_t pI = 0; pI < nPhiBins; ++pI){
	etaPhiTower[eI].push_back(TLorentzVector(0,0,0,0));
	etaPhiGoodTower[eI].push_back(true);
      }
    }
    
    std::cout << "PhiBins: " << std::endl;
    for(unsigned int pI = 0; pI < phiBins_p->size(); ++pI){
      std::cout << " " << pI << "/" << phiBins_p->size() << ": " << phiBins_p->at(pI) << std::endl;
    }

    std::cout << "EtaBins: " << std::endl;
    for(unsigned int pI = 0; pI < etaBins_p->size(); ++pI){
      std::cout << " " << pI << "/" << etaBins_p->size() << ": " << etaBins_p->at(pI) << std::endl;
    }
  }
  
  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* clusterJetsCS_p = new TTree("clusterJetsCS", "");
  clusterJetsCS_p->Branch("run", &run_, "run/I");
  clusterJetsCS_p->Branch("lumi", &lumi_, "lumi/i");
  clusterJetsCS_p->Branch("evt", &evt_, "evt/I");

  clusterJetsCS_p->Branch("jzVal", &jzVal_, "jzVal/I");

  clusterJetsCS_p->Branch("fcalA_et", &fcalA_et_, "fcalA_et/F");
  clusterJetsCS_p->Branch("fcalC_et", &fcalC_et_, "fcalC_et/F");

  clusterJetsCS_p->Branch("cent", &cent_, "cent/F");
  clusterJetsCS_p->Branch("etaBins", &etaBinsOut_p);
  clusterJetsCS_p->Branch("rho", &rhoOut_p);
  clusterJetsCS_p->Branch("rhoCorr", &rhoCorrOut_p);

  for(Int_t jI = 0; jI < nJtAlgo; ++jI){
    clusterJetsCS_p->Branch(("njt" + jtAlgos[jI]).c_str(), &(njt_[jI]), ("njt" + jtAlgos[jI] + "/I").c_str());
    clusterJetsCS_p->Branch(("jtpt" + jtAlgos[jI]).c_str(), jtpt_[jI], ("jtpt" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    clusterJetsCS_p->Branch(("jteta" + jtAlgos[jI]).c_str(), jteta_[jI], ("jteta" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    clusterJetsCS_p->Branch(("jtphi" + jtAlgos[jI]).c_str(), jtphi_[jI], ("jtphi" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/F").c_str());
    if(doATLASFile){
      clusterJetsCS_p->Branch(("atlasmatchpos" + jtAlgos[jI]).c_str(), atlasmatchpos_[jI], ("atlasmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());
      if(doTruth){
	clusterJetsCS_p->Branch(("truthmatchpos" + jtAlgos[jI]).c_str(), truthmatchpos_[jI], ("truthmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());
	if(!doCalo) clusterJetsCS_p->Branch(("truth4GeVmatchpos" + jtAlgos[jI]).c_str(), truth4GeVmatchpos_[jI], ("truth4GeVmatchpos" + jtAlgos[jI] + "[njt" + jtAlgos[jI] + "]/I").c_str());
      }

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
      clusterJetsCS_p->Branch("jtchgptTruth", jtchgptTruth_, "jtchgptTruth[njtTruth]/F");
      clusterJetsCS_p->Branch("jtchgetaTruth", jtchgetaTruth_, "jtchgetaTruth[njtTruth]/F");
      clusterJetsCS_p->Branch("jtchgphiTruth", jtchgphiTruth_, "jtchgphiTruth[njtTruth]/F");

      if(!doCalo){
	clusterJetsCS_p->Branch("njtTruth4GeV", &njtTruth4GeV_, "njtTruth4GeV/I");
	clusterJetsCS_p->Branch("jtptTruth4GeV", jtptTruth4GeV_, "jtptTruth4GeV[njtTruth4GeV]/F");
	clusterJetsCS_p->Branch("jtetaTruth4GeV", jtetaTruth4GeV_, "jtetaTruth4GeV[njtTruth4GeV]/F");
	clusterJetsCS_p->Branch("jtphiTruth4GeV", jtphiTruth4GeV_, "jtphiTruth4GeV[njtTruth4GeV]/F");     
	clusterJetsCS_p->Branch("jtchgptTruth4GeV", jtchgptTruth4GeV_, "jtchgptTruth4GeV[njtTruth4GeV]/F");
	clusterJetsCS_p->Branch("jtchgphiTruth4GeV", jtchgphiTruth4GeV_, "jtchgphiTruth4GeV[njtTruth4GeV]/F");
	clusterJetsCS_p->Branch("jtchgetaTruth4GeV", jtchgetaTruth4GeV_, "jtchgetaTruth4GeV[njtTruth4GeV]/F");
      }
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

  const Int_t nEntries = TMath::Min(50000, (Int_t)clusterTree_p->GetEntries());
  const Int_t nDiv = TMath::Max(1, nEntries/400);

  /*
  clusterJetsCS_p->Branch("run", &run_, "run/I");
  clusterJetsCS_p->Branch("lumi", &lumi_, "lumi/i");
  clusterJetsCS_p->Branch("evt", &evt_, "evt/I");

  clusterJetsCS_p->Branch("fcalA_et", &fcalA_et_, "fcalA_et/F");
  clusterJetsCS_p->Branch("fcalC_et", &fcalC_et_, "fcalC_et/F");

  clusterJetsCS_p->Branch("cent", &cent_, "cent/F");
  */
  std::vector<Int_t> runVect;
  std::vector<UInt_t> lumiVect;
  std::vector<Int_t> evtVect;
  std::vector<float> fcalA_etVect;
  std::vector<float> fcalC_etVect;
  std::vector<float> centVect;
  std::vector<std::vector<float> > rhoVect;
  std::vector<std::vector<float> > rhoCorrVect;
  std::vector<std::vector<fastjet::PseudoJet> > particles;
  std::vector<std::map<std::string, std::vector<fastjet::PseudoJet> > > jets;
  
  std::vector<std::vector<float> > atlasPt, atlasPhi, atlasEta;
  std::vector<std::vector<float> > truthPt, truthPhi, truthEta, truthChgPt, truthChgPhi, truthChgEta;
  std::vector<std::vector<float> > truth4GeVPt, truth4GeVPhi, truth4GeVEta, truth4GeVChgPt, truth4GeVChgPhi, truth4GeVChgEta;
  
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
    truthChgPt.push_back({});
    truthChgEta.push_back({});
    truthChgPhi.push_back({});

    truthPt[pI].reserve(500);
    truthPhi[pI].reserve(500);
    truthEta[pI].reserve(500);
    truthChgPt[pI].reserve(500);
    truthChgEta[pI].reserve(500);
    truthChgPhi[pI].reserve(500);

    truth4GeVPt.push_back({});
    truth4GeVPhi.push_back({});
    truth4GeVEta.push_back({});
    truth4GeVChgPt.push_back({});
    truth4GeVChgPhi.push_back({});
    truth4GeVChgEta.push_back({});

    truth4GeVPt[pI].reserve(500);
    truth4GeVPhi[pI].reserve(500);
    truth4GeVEta[pI].reserve(500);
    truth4GeVChgPt[pI].reserve(500);
    truth4GeVChgEta[pI].reserve(500);
    truth4GeVChgPhi[pI].reserve(500);

    jets.push_back({});
    for(Int_t jI = 0; jI < nJtAlgo; ++jI){
      jets[pI][jtAlgos[jI]] = {};
      jets[pI][jtAlgos[jI]].reserve(500);
    }
  }

  //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  preLoop.stop();  
  TLorentzVector tL;
  
  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    preCluster.start();
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    clusterTree_p->GetEntry(entry);

    if(rhoOut_p->size() == 0){
      for(unsigned int eI = 0; eI < etaBins_p->size(); ++eI){
	etaBinsOut_p->push_back(etaBins_p->at(eI));
      }

      for(unsigned int rI = 0; rI < rho_p->size(); ++rI){
	rhoOut_p->push_back(0.0);
	rhoCorrOut_p->push_back(0.0);
      }
    }
    else{
      for(unsigned int rI = 0; rI < rhoOut_p->size(); ++rI){
	(*rhoOut_p)[rI] = 0.0;
	(*rhoCorrOut_p)[rI] = 0.0;
      }
    }

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
	E = tL.E();
	Px = tL.Px();
	Py = tL.Py();
	Pz = tL.Pz();
      }

      particles[entry%nPara].push_back(fastjet::PseudoJet(Px, Py, Pz, E));
    }

    //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(!doCalo){
      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      for(unsigned int rI = 0; rI < rho_p->size(); ++rI){
	rho_p->at(rI) = 0.0;
	rho2_p->at(rI) = 0.0;

	etaNGoodTowers[rI] = etaPhiTower[rI].size();
	for(unsigned int pI = 0; pI < etaPhiTower[rI].size(); ++pI){
	  etaPhiTower[rI][pI] = TLorentzVector(0,0,0,0);
	  etaPhiGoodTower[rI][pI] = true;
	}
      }

      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
     
      for(unsigned int pI = 0; pI < particles[entry%nPara].size(); ++pI){
	Int_t etaPos = ghostEtaPos(*etaBins_p, particles[entry%nPara][pI]);
	Int_t phiPos = ghostPhiPos(*phiBins_p, particles[entry%nPara][pI]);

	Float_t pX = particles[entry%nPara][pI].px();
	Float_t pY = particles[entry%nPara][pI].py();
	Float_t pZ = particles[entry%nPara][pI].pz();
	Float_t E = particles[entry%nPara][pI].E();
	
	rho_p->at(etaPos) += E;
	rho2_p->at(etaPos) += E*E;
	TLorentzVector tL(pX, pY, pZ, E);
	etaPhiTower[etaPos][phiPos] += tL;
      }

      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      std::vector<float> mean, sigma;
      std::vector<fastjet::PseudoJet> towerParticles;
      for(unsigned int rI = 0; rI < rho_p->size(); ++rI){
	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
	mean.push_back(rho_p->at(rI)/(Double_t)etaPhiTower[rI].size());
	sigma.push_back(TMath::Sqrt(rho2_p->at(rI)/(Double_t)etaPhiTower[rI].size() - mean[rI]*mean[rI]));

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	Double_t etaVal = (etaBins_p->at(rI) + etaBins_p->at(rI+1))/2.;

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

	for(unsigned int pI = 0; pI < etaPhiTower[rI].size(); ++pI){
	  if(etaPhiTower[rI][pI].E() > mean[rI] + sigma[rI]){
	    Double_t phiVal = (phiBins_p->at(pI) + phiBins_p->at(pI+1))/2.;
	    Double_t newE = etaPhiTower[rI][pI].E() - (mean[rI] + sigma[rI]);
	    Double_t newEt = newE/std::cosh(etaVal);
	    Double_t Px = newEt*std::cos(phiVal);
	    Double_t Py = newEt*std::sin(phiVal);
	    Double_t Pz = newEt*std::sinh(etaVal);

	    towerParticles.push_back(fastjet::PseudoJet(Px, Py, Pz, newE));
	  }
	}

	//std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      }

      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      fastjet::ClusterSequence cs(towerParticles, jet_def);
      std::vector<fastjet::PseudoJet> towerJets = fastjet::sorted_by_pt(cs.inclusive_jets(minRhoJtPt));

      for(unsigned int rI = 0; rI < rho_p->size(); ++rI){
	rhoOut_p->at(rI) = rho_p->at(rI)/(2.*TMath::Pi()*deltaEta);
	rho_p->at(rI) = 0.0;

	Double_t etaVal = (etaBins_p->at(rI) + etaBins_p->at(rI+1))/2.;
	for(unsigned int pI = 0; pI < etaPhiTower[rI].size(); ++pI){
	  Double_t phiVal = (phiBins_p->at(pI) + phiBins_p->at(pI+1))/2.;

	  for(unsigned int tI = 0; tI < towerJets.size(); ++tI){
	    if(getDR(etaVal, phiVal, towerJets[tI].eta(), towerJets[tI].phi_std()) < rParam){
	      etaPhiGoodTower[rI][pI] = false;
	      --(etaNGoodTowers[rI]);
	      break;
	    }
	  }
	}
      }      

      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
      
      for(unsigned int pI = 0; pI < particles[entry%nPara].size(); ++pI){
	Float_t etaVal = particles[entry%nPara][pI].eta();
	Float_t phiVal = particles[entry%nPara][pI].phi_std();
	Int_t etaPos = ghostPos(*etaBins_p, etaVal);
	Int_t phiPos = ghostPos(*phiBins_p, phiVal);
	Float_t E = particles[entry%nPara][pI].E();

	if(etaPhiGoodTower[etaPos][phiPos]) rho_p->at(etaPos) += E;
      }

      
      for(unsigned int rI = 0; rI < rho_p->size(); ++rI){
	Double_t areaFactor = ((Double_t)etaNGoodTowers[rI])/(Double_t)etaPhiGoodTower[rI].size();
	rho_p->at(rI) /= (deltaEta*2.*TMath::Pi()*areaFactor);
	rhoCorrOut_p->at(rI) = rho_p->at(rI);
      }
    }
    else{
      for(unsigned int rI = 0; rI < rho_p->size(); ++rI){
	rho_p->at(rI) /= 1000.;
	rhoOut_p->at(rI) = rho_p->at(rI);
	rhoCorrOut_p->at(rI) = rho_p->at(rI);
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
    rhoVect.push_back((*rhoOut_p));
    rhoCorrVect.push_back((*rhoCorrOut_p));

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
	std::vector<fastjet::PseudoJet> tempParticles, tempParticles4GeV;
	std::vector<fastjet::PseudoJet> tempChgParticles, tempChgParticles2;

	for(unsigned int tI = 0; tI < truth_pt_p->size(); ++tI){
	  if(truth_pt_p->at(tI) < 0.1) continue;
	  if(TMath::Abs(truth_eta_p->at(tI)) > maxJtAbsEta + rParam) continue;
	  tL.SetPtEtaPhiM(truth_pt_p->at(tI), truth_eta_p->at(tI), truth_phi_p->at(tI), 0.0);
	  if(tL.E() < 0.1) continue;

	  int userIndex = 1;
	  if(TMath::Abs(truth_chg_p->at(tI)) < 0.01) userIndex = 0;
	  
	  tempParticles.push_back(fastjet::PseudoJet(tL.Px(), tL.Py(), tL.Pz(), tL.E()));
	  tempParticles[tempParticles.size()-1].set_user_index(userIndex);

	  if(truth_pt_p->at(tI) > 4) tempParticles4GeV.push_back(tempParticles[tempParticles.size()-1]);
	}
	
	fastjet::ClusterSequence cs(tempParticles, jet_def);
	std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets(minJtPt));
	for(unsigned int aI = 0; aI < tempJets.size(); ++aI){
	  truthPt[entry%nPara].push_back(tempJets.at(aI).pt());

	  tempChgParticles = tempJets.at(aI).constituents();
	  for(unsigned int cI = 0; cI < tempChgParticles.size(); ++cI){
	    if(tempChgParticles[cI].user_index() == 1) tempChgParticles2.push_back(tempChgParticles[cI]);
	  }

	  fastjet::PseudoJet chgJet = join(tempChgParticles2);
	  truthChgPt[entry%nPara].push_back(chgJet.pt());
	  truthChgPhi[entry%nPara].push_back(chgJet.phi_std());
	  truthChgEta[entry%nPara].push_back(chgJet.eta());
	  tempChgParticles2.clear();
	  truthEta[entry%nPara].push_back(tempJets.at(aI).eta());
	  truthPhi[entry%nPara].push_back(tempJets.at(aI).phi_std());
	}		

	fastjet::ClusterSequence cs4GeV(tempParticles4GeV, jet_def);
	tempJets = fastjet::sorted_by_pt(cs4GeV.inclusive_jets(minJtPt));
	for(unsigned int aI = 0; aI < tempJets.size(); ++aI){
	  truth4GeVPt[entry%nPara].push_back(tempJets.at(aI).pt());

	  tempChgParticles = tempJets.at(aI).constituents();
	  for(unsigned int cI = 0; cI < tempChgParticles.size(); ++cI){
	    if(tempChgParticles[cI].user_index() == 1) tempChgParticles2.push_back(tempChgParticles[cI]);
	  }

	  fastjet::PseudoJet chgJet = join(tempChgParticles2);
	  truth4GeVChgPt[entry%nPara].push_back(chgJet.pt());
	  truth4GeVChgEta[entry%nPara].push_back(chgJet.eta());
	  truth4GeVChgPhi[entry%nPara].push_back(chgJet.phi_std());
	  tempChgParticles2.clear();
	  truth4GeVEta[entry%nPara].push_back(tempJets.at(aI).eta());
	  truth4GeVPhi[entry%nPara].push_back(tempJets.at(aI).phi_std());
	}		
      }
      else{
	for(unsigned int aI = 0; aI < truth_pt_p->size(); ++aI){
	  truthPt[entry%nPara].push_back(truth_pt_p->at(aI));
	  truthEta[entry%nPara].push_back(truth_eta_p->at(aI));
	  truthPhi[entry%nPara].push_back(truth_phi_p->at(aI));

	  truthChgPt[entry%nPara].push_back(-1);
	  truthChgEta[entry%nPara].push_back(-1);
	  truthChgPhi[entry%nPara].push_back(-1);
	}			
      }
    }
    
    preCluster.stop();
    if((entry + 1)%nPara == 0 || entry == nEntries-1){
#pragma omp parallel
      {
#pragma omp for
	for(Int_t pI = 0; pI < nPara; ++pI){	  
	  getJetsFromParticles(rhoVect[pI], (*etaBins_p), particles[pI], &(jets[pI]), {&(inCluster1[pI]), &(inCluster2[pI]), &(inCluster3[pI]), &(inCluster4[pI]), &(inCluster5[pI]), &(inCluster6[pI])});
	}      
      }
        
      postCluster.start();

      for(Int_t pI = 0; pI < nPara; ++pI){
	for(Int_t jI = 0; jI < nJtAlgo; ++jI){
	  njt_[jI] = 0;
	
	  for(fastjet::PseudoJet& ijet : jets[pI][jtAlgos[jI]] ) {
	    if(ijet.pt() < minJtPt) continue;
	    if(TMath::Abs(ijet.eta()) >= maxJtAbsEta) continue;
	    if(!doCalo){
	      if(TMath::Abs(ijet.eta()) >= maxTrkJtAbsEta) continue;
	    }
	  
	    jtpt_[jI][njt_[jI]] = ijet.pt();
	    jtphi_[jI][njt_[jI]] = ijet.phi_std();
	    jteta_[jI][njt_[jI]] = ijet.eta();
	    atlasmatchpos_[jI][njt_[jI]] = -1;
	    truthmatchpos_[jI][njt_[jI]] = -1;
	    truth4GeVmatchpos_[jI][njt_[jI]] = -1;
	    
	    ++njt_[jI];
	  }
	}	
            
	if(doATLASFile){
       	  njtATLAS_ = 0;
	  for(unsigned int aI = 0; aI < atlasPt[pI].size(); ++aI){
	    if(atlasPt[pI][aI] < minJtPt) continue;
            if(TMath::Abs(atlasEta[pI][aI]) >= maxJtAbsEta) continue;
            if(!doCalo){
              if(TMath::Abs(atlasEta[pI][aI]) >= maxTrkJtAbsEta) continue;
            }

	    jtptATLAS_[njtATLAS_] = atlasPt[pI][aI];
	    jtetaATLAS_[njtATLAS_] = atlasEta[pI][aI];
	    jtphiATLAS_[njtATLAS_] = atlasPhi[pI][aI];
	    ++njtATLAS_;
	  }

	  if(doTruth){
	    njtTruth_ = 0;
	    for(unsigned int aI = 0; aI < truthPt[pI].size(); ++aI){
	      if(truthPt[pI][aI] < minJtPt) continue;
	      if(TMath::Abs(truthEta[pI][aI]) >= maxJtAbsEta) continue;
	      if(!doCalo){
		if(TMath::Abs(truthEta[pI][aI]) >= maxTrkJtAbsEta) continue;
	      }

	      jtptTruth_[njtTruth_] = truthPt[pI][aI];
	      jtetaTruth_[njtTruth_] = truthEta[pI][aI];
	      jtphiTruth_[njtTruth_] = truthPhi[pI][aI];
	      jtchgptTruth_[njtTruth_] = truthChgPt[pI][aI];
	      jtchgphiTruth_[njtTruth_] = truthChgPhi[pI][aI];
	      jtchgetaTruth_[njtTruth_] = truthChgEta[pI][aI];
	      ++njtTruth_;
	    }

	    if(!doCalo){
	      njtTruth4GeV_ = 0;
	      for(unsigned int aI = 0; aI < truth4GeVPt[pI].size(); ++aI){
		if(truth4GeVPt[pI][aI] < minJtPt) continue;
		if(TMath::Abs(truth4GeVEta[pI][aI]) >= maxJtAbsEta) continue;
		if(!doCalo){
		  if(TMath::Abs(truth4GeVEta[pI][aI]) >= maxTrkJtAbsEta) continue;
		}

		jtptTruth4GeV_[njtTruth4GeV_] = truth4GeVPt[pI][aI];
		jtetaTruth4GeV_[njtTruth4GeV_] = truth4GeVEta[pI][aI];
		jtphiTruth4GeV_[njtTruth4GeV_] = truth4GeVPhi[pI][aI];
		jtchgptTruth4GeV_[njtTruth4GeV_] = truth4GeVChgPt[pI][aI];
		jtchgphiTruth4GeV_[njtTruth4GeV_] = truth4GeVChgPhi[pI][aI];
		jtchgetaTruth4GeV_[njtTruth4GeV_] = truth4GeVChgEta[pI][aI];
		++njtTruth4GeV_;
	      }
	    }
	  }

	  for(Int_t jI = 0; jI < nJtAlgo; ++jI){
	    std::vector<bool> isMatched;
	    for(int aI = 0; aI < njtATLAS_; ++aI){
	      isMatched.push_back(false);
	    }
	  
	    for(Int_t jI2 = 0; jI2 < njt_[jI]; ++jI2){
	      for(int aI = 0; aI < njtATLAS_; ++aI){
		if(isMatched[aI]) continue;
		if(getDR(jtetaATLAS_[aI], jtphiATLAS_[aI], jteta_[jI][jI2], jtphi_[jI][jI2]) < rParam){
		  atlasmatchpos_[jI][jI2] = aI;
		  isMatched[aI] = true;
		  break;
		}
	      }
	    }

	    if(doTruth){
	      isMatched.clear();
	      for(int aI = 0; aI < njtTruth_; ++aI){
		isMatched.push_back(false);
	      }
	      
	      for(Int_t jI2 = 0; jI2 < njt_[jI]; ++jI2){
		for(int aI = 0; aI < njtTruth_; ++aI){
		  if(isMatched[aI]) continue;
		  if(getDR(jtetaTruth_[aI], jtphiTruth_[aI], jteta_[jI][jI2], jtphi_[jI][jI2]) < rParam){
		    truthmatchpos_[jI][jI2] = aI;
		    isMatched[aI] = true;
		    break;
		  }
		}
	      }

	      if(!doCalo){
		isMatched.clear();
		for(int aI = 0; aI < njtTruth4GeV_; ++aI){
		  isMatched.push_back(false);
		}
		
		for(Int_t jI2 = 0; jI2 < njt_[jI]; ++jI2){
		  for(int aI = 0; aI < njtTruth4GeV_; ++aI){
		    if(isMatched[aI]) continue;
		    if(getDR(jtetaTruth4GeV_[aI], jtphiTruth4GeV_[aI], jteta_[jI][jI2], jtphi_[jI][jI2]) < rParam){
		      truth4GeVmatchpos_[jI][jI2] = aI;
		      isMatched[aI] = true;
		      break;
		    }
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
	  truthChgPt[pI].clear();
	  truthChgPhi[pI].clear();
	  truthChgEta[pI].clear();

	  if(!doCalo){
	    truth4GeVPt[pI].clear();
	    truth4GeVEta[pI].clear();
	    truth4GeVPhi[pI].clear();
	    truth4GeVChgPt[pI].clear();
	    truth4GeVChgPhi[pI].clear();
	    truth4GeVChgEta[pI].clear();
	  }
	}
      
	run_ = runVect[pI];
	lumi_ = lumiVect[pI];
	evt_ = evtVect[pI];
	fcalA_et_ = fcalA_etVect[pI];
	fcalC_et_ = fcalC_etVect[pI];
	cent_ = centVect[pI];

	(*rhoOut_p) = rhoVect[pI];
	(*rhoCorrOut_p) = rhoCorrVect[pI];
	
	clusterJetsCS_p->Fill();
      }
    
      runVect.clear();
      lumiVect.clear();
      evtVect.clear();
      fcalA_etVect.clear();
      fcalC_etVect.clear();
      
      centVect.clear();

      rhoVect.clear();
      rhoCorrVect.clear();

      //std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      for(Int_t pI = 0; pI < nPara; ++pI){
	particles[pI].clear();
      }

      postCluster.stop();
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

  paramMap["minJtPt"] = prettyString(minJtPt, 2, false);
  paramMap["minRhoJtPt"] = prettyString(minRhoJtPt, 2, false);

  if(doCalo) paramMap["maxJtAbsEta"] = prettyString(maxJtAbsEta, 2, false);
  else paramMap["maxJtAbsEta"] = prettyString(maxTrkJtAbsEta, 2, false);
  
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
  std::cout << " PRELOOP: " << preLoop.totalWall() << "/" << total.totalWall() << "=" << preLoop.totalWall()/total.totalWall() << std::endl;
  std::cout << " PRECLUSTER: " << preCluster.totalWall() << "/" << total.totalWall() << "=" << preCluster.totalWall()/total.totalWall()<< std::endl;
  std::cout << " INCLUSTER1: " << inCluster1[0].totalWall() << "/" << total.totalWall() << "=" << inCluster1[0].totalWall()/total.totalWall()<< std::endl;
  std::cout << " INCLUSTER2: " << inCluster2[0].totalWall() << "/" << total.totalWall() << "=" << inCluster2[0].totalWall()/total.totalWall()<< std::endl;
  std::cout << " INCLUSTER3: " << inCluster3[0].totalWall() << "/" << total.totalWall() << "=" << inCluster3[0].totalWall()/total.totalWall()<< std::endl;
  std::cout << "  SUBCLUSTER3_1: " << inCluster4[0].totalWall() << "/" << total.totalWall() << "=" << inCluster4[0].totalWall()/total.totalWall()<< std::endl;
  std::cout << "  SUBCLUSTER3_2: " << inCluster5[0].totalWall() << "/" << total.totalWall() << "=" << inCluster5[0].totalWall()/total.totalWall()<< std::endl;
  std::cout << "  SUBCLUSTER3_3: " << inCluster6[0].totalWall() << "/" << total.totalWall() << "=" << inCluster6[0].totalWall()/total.totalWall()<< std::endl;
  std::cout << " POSTCLUSTER: " << postCluster.totalWall() << "/" << total.totalWall() << "=" << postCluster.totalWall()/total.totalWall()<< std::endl;
  std::cout << " POSTLOOP: " << postLoop.totalWall() << "/" << total.totalWall() << "=" << postLoop.totalWall()/total.totalWall() << std::endl;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 5){
    std::cout << "Usage: ./bin/clusterToCS.exe <inFileName> <inATLASFileName-default=\'\'> <caloTrackStr-default=\'calo\'> <jzStr-default=\'\'>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += clusterToCS(argv[1]);
  else if(argc == 3) retVal += clusterToCS(argv[1], argv[2]);
  else if(argc == 4) retVal += clusterToCS(argv[1], argv[2], argv[3]);
  else if(argc == 5) retVal += clusterToCS(argv[1], argv[2], argv[3], argv[4]);
  return retVal;
}
