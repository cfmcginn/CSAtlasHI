//AUTHOR Chris McGinn
//cpp
#include <iostream>
#include <string>
#include <vector>

//FASTJET
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh" 
//PYTHIA
#include "Pythia8/Pythia.h"

//ROOT
#include "TDatime.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/etaPhiFunc.h"
#include "include/pdgToMassInGeV.h"

int constituentTest(const std::string hydFileName, const ULong64_t nEntries = 10000)
{
  if(!checkFile(hydFileName) || hydFileName.find(".root") == std::string::npos){
    std::cout << "Input hydFileName \'" << hydFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  const std::string outFileName = "output/" + dateStr + "/constTest_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p = new TTree("outTree", "");

  const Int_t nMaxPart = 50000;
  Int_t nGen_;
  Float_t genPt_[nMaxPart];
  Float_t genPhi_[nMaxPart];
  Float_t genEta_[nMaxPart];
  Int_t genID_[nMaxPart];

  Float_t rho_;
  
  const Int_t nMaxJets = 500;
  Int_t nJt_ = 0;
  Float_t jtptNoSub_[nMaxJets];
  Float_t jtphiNoSub_[nMaxJets];
  Float_t jtetaNoSub_[nMaxJets];
  Float_t jtptAreaSub_[nMaxJets];
  Float_t jtphiAreaSub_[nMaxJets];
  Float_t jtetaAreaSub_[nMaxJets];
  Float_t jtptConstSub_[nMaxJets];
  Float_t jtphiConstSub_[nMaxJets];
  Float_t jtetaConstSub_[nMaxJets];
  Float_t refpt_[nMaxJets];
  Float_t refeta_[nMaxJets];
  Float_t refphi_[nMaxJets];

  outTree_p->Branch("rho", &rho_, "rho/F");
  outTree_p->Branch("nJt", &nJt_, "nJt/I");
  outTree_p->Branch("jtptNoSub", jtptNoSub_, "jtptNoSub[nJt]/F");
  outTree_p->Branch("jtphiNoSub", jtphiNoSub_, "jtphiNoSub[nJt]/F");
  outTree_p->Branch("jtetaNoSub", jtetaNoSub_, "jtetaNoSub[nJt]/F");
  outTree_p->Branch("jtptAreaSub", jtptAreaSub_, "jtptAreaSub[nJt]/F");
  outTree_p->Branch("jtphiAreaSub", jtphiAreaSub_, "jtphiAreaSub[nJt]/F");
  outTree_p->Branch("jtetaAreaSub", jtetaAreaSub_, "jtetaAreaSub[nJt]/F");
  outTree_p->Branch("jtptConstSub", jtptConstSub_, "jtptConstSub[nJt]/F");
  outTree_p->Branch("jtphiConstSub", jtphiConstSub_, "jtphiConstSub[nJt]/F");
  outTree_p->Branch("jtetaConstSub", jtetaConstSub_, "jtetaConstSub[nJt]/F");
  outTree_p->Branch("refpt", refpt_, "refpt[nJt]/F");
  outTree_p->Branch("refphi", refphi_, "refphi[nJt]/F");
  outTree_p->Branch("refeta", refeta_, "refeta[nJt]/F");

  //GRAB HYDJET
  TFile* hydFile_p = new TFile(hydFileName.c_str(), "READ");
  TTree* hydTree_p = (TTree*)hydFile_p->Get("genTree");

  Int_t currHYD = 0;
  
  hydTree_p->SetBranchStatus("*", 0);
  hydTree_p->SetBranchStatus("nGen", 1);
  hydTree_p->SetBranchStatus("genPt", 1);
  hydTree_p->SetBranchStatus("genPhi", 1);
  hydTree_p->SetBranchStatus("genEta", 1);
  hydTree_p->SetBranchStatus("genID", 1);

  hydTree_p->SetBranchAddress("nGen", &nGen_);
  hydTree_p->SetBranchAddress("genPt", genPt_);
  hydTree_p->SetBranchAddress("genPhi", genPhi_);
  hydTree_p->SetBranchAddress("genEta", genEta_);
  hydTree_p->SetBranchAddress("genID", genID_);

  //PYTHIA
  Pythia8::Pythia pythia;
  pythia.readString("Beams:eCM = 5020.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  pythia.readString("PhaseSpace:pTHatMin = 50.");
  pythia.init();
  
  std::cout << "Processing " << nEntries << " events..." << std::endl;

  const double rParam = 0.4;
  const Double_t minPartPt = 0.5;
  const Double_t maxPartAbsEta = 5.1;
  const Double_t minJtPt = 50.;
  pdgToMassInGeV pdgToM;
  
  fastjet::JetDefinition jet_defE(fastjet::antikt_algorithm, rParam, fastjet::E_scheme);
  fastjet::JetDefinition jet_defEBkg(fastjet::kt_algorithm, rParam, fastjet::E_scheme);

  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);
  double ghost_maxrap = 5.0; // e.g. if particles go up to y=5
  fastjet::AreaDefinition area_def(fastjet::active_area, fastjet::GhostedAreaSpec(ghost_maxrap));
  fastjet::JetMedianBackgroundEstimator bge(fastjet::SelectorAbsRapMax(4.5), jet_defEBkg, area_def);:
  fastjet::JetMedianBackgroundEstimator bge2(fastjet::SelectorAbsRapMax(4.5), jet_defEBkg, area_def);
  fastjet::Subtractor subtractor(&bge);

  fastjet::BackgroundJetScalarPtDensity *scalarPtDensity = new fastjet::BackgroundJetScalarPtDensity();
  bge2.set_jet_density_class(scalarPtDensity);  
  fastjet::contrib::ConstituentSubtractor subtractorCS(&bge2);

  subtractor.set_use_rho_m(true);
  subtractor.set_safe_mass(true);
  
  while(nEntries > (ULong64_t)outTree_p->GetEntries()){
    if(!pythia.next()) continue;

    std::vector<fastjet::PseudoJet> particles;
    for(Int_t pI = 0; pI < pythia.event.size(); ++pI){
      if(!pythia.event[pI].isFinal()) continue;
      if(pythia.event[pI].pT() < minPartPt) continue;
      if(TMath::Abs(pythia.event[pI].eta()) > maxPartAbsEta) continue;

      if(TMath::Abs(pythia.event[pI].id()) == 12) continue;
      if(TMath::Abs(pythia.event[pI].id()) == 14) continue;
      if(TMath::Abs(pythia.event[pI].id()) == 16) continue;

      particles.push_back(fastjet::PseudoJet(pythia.event[pI].px(), pythia.event[pI].py(), pythia.event[pI].pz(), pythia.event[pI].e()));
    }

    fastjet::ClusterSequenceArea csE(particles, jet_defE, area_def);
    std::vector<fastjet::PseudoJet> jetsE = fastjet::sorted_by_pt(csE.inclusive_jets());
    
    nJt_ = 0;
    //    std::cout << "JETS SIZE: " << jetsE.size() << std::endl;
    for(unsigned int jI = 0; jI < jetsE.size(); ++jI){
      if(jetsE[jI].pt() < minJtPt) continue;

      refpt_[nJt_] = jetsE[jI].pt();
      refphi_[nJt_] = jetsE[jI].phi_std();
      refeta_[nJt_] = jetsE[jI].eta();

      jtptNoSub_[nJt_] = -999;
      jtphiNoSub_[nJt_] = -999;
      jtetaNoSub_[nJt_] = -999;

      jtptAreaSub_[nJt_] = -999;
      jtphiAreaSub_[nJt_] = -999;
      jtetaAreaSub_[nJt_] = -999;

      jtptConstSub_[nJt_] = -999;
      jtphiConstSub_[nJt_] = -999;
      jtetaConstSub_[nJt_] = -999;

      ++nJt_;
    }
    
    if(nJt_ == 0) continue;

    hydTree_p->GetEntry(currHYD);
    ++currHYD;
    
    if(outTree_p->GetEntries()%nDiv == 0) std::cout << " Entry " << outTree_p->GetEntries() << std::endl;

    TLorentzVector tL;
    for(Int_t pI = 0; pI < nGen_; ++pI){
      if(genPt_[pI] < minPartPt) continue;
      if(TMath::Abs(genEta_[pI]) > maxPartAbsEta) continue;

      if(TMath::Abs(genID_[pI]) == 12) continue;
      if(TMath::Abs(genID_[pI]) == 14) continue;
      if(TMath::Abs(genID_[pI]) == 16) continue;

      tL.SetPtEtaPhiM(genPt_[pI], genEta_[pI], genPhi_[pI], pdgToM.getMassFromID(TMath::Abs(genID_[pI])));
      particles.push_back(fastjet::PseudoJet(tL.Px(), tL.Py(), tL.Pz(), tL.E()));
    }

    fastjet::ClusterSequenceArea csE2(particles, jet_defE, area_def);
    std::vector<fastjet::PseudoJet> jetsE2 = fastjet::sorted_by_pt(csE2.inclusive_jets());
    bge.set_particles(particles);
    bge2.set_particles(particles);
    rho_ = bge.rho();
    std::vector<fastjet::PseudoJet> subtracted_jets = fastjet::sorted_by_pt(subtractor(jetsE2));

    std::vector<bool> isUsed;
    for(unsigned int jI = 0; jI < jetsE2.size(); ++jI){isUsed.push_back(false);}

    for(Int_t jI = 0; jI < nJt_; ++jI){
      for(unsigned int jI2 = 0; jI2 < jetsE2.size(); ++jI2){
	if(isUsed[jI2]) continue;

	if(getDR(refeta_[jI], refphi_[jI], jetsE2[jI2].eta(), jetsE2[jI2].phi_std()) < 0.2){
	  jtptNoSub_[jI] = jetsE2[jI2].pt();
	  jtphiNoSub_[jI] = jetsE2[jI2].phi_std();
	  jtetaNoSub_[jI] = jetsE2[jI2].eta();	  

	  fastjet::PseudoJet subtracted_jet = subtractorCS(jetsE2[jI2]);

	  jtptConstSub_[jI] = subtracted_jet.pt();
	  jtphiConstSub_[jI] = subtracted_jet.phi_std();
	  jtetaConstSub_[jI] = subtracted_jet.eta();	  

	  
	  isUsed[jI2] = true;
	}     
      }
    }    


    isUsed.clear();
    for(unsigned int jI = 0; jI < subtracted_jets.size(); ++jI){isUsed.push_back(false);}

    for(Int_t jI = 0; jI < nJt_; ++jI){
      for(unsigned int jI2 = 0; jI2 < subtracted_jets.size(); ++jI2){
	if(isUsed[jI2]) continue;

	if(getDR(refeta_[jI], refphi_[jI], subtracted_jets[jI2].eta(), subtracted_jets[jI2].phi_std()) < 0.2){
	  jtptAreaSub_[jI] = subtracted_jets[jI2].pt();
	  jtphiAreaSub_[jI] = subtracted_jets[jI2].phi_std();
	  jtetaAreaSub_[jI] = subtracted_jets[jI2].eta();

	  isUsed[jI2] = true;
	}     
      }
    }    

    
    particles.clear();
    outTree_p->Fill();
  }

  hydFile_p->Close();
  delete hydFile_p;
  
  outFile_p->cd();

  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;
  
  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 3){
    std::cout << "Usage: ./bin/constituentTest.exe <inHydFileName> <nEntries>. return 1" << std::endl;
    return 1;
  }
  
  int retVal;
  if(argc == 2) retVal += constituentTest(argv[1]);
  else if(argc == 3) retVal += constituentTest(argv[1], std::stol(argv[2]));
  return retVal;
}
