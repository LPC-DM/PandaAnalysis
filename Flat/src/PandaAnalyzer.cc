#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

using namespace panda;
using namespace std;

PandaAnalyzer::PandaAnalyzer(int debug_/*=0*/) {
  DEBUG = debug_;

  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Calling constructor");
  gt = new GeneralTree();
  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Built GeneralTree");
  ibetas = gt->get_ibetas();
  Ns = gt->get_Ns();
  orders = gt->get_orders();
  flags["fatjet"]         = true;
  flags["puppi"]          = true;
  flags["monohiggs"]      = false;
  flags["vbf"]            = false;
  flags["firstGen"]       = true;
  flags["applyJSON"]      = true;
  flags["genOnly"]        = false;
  flags["pfCands"]        = false;
  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Called constructor");
}


PandaAnalyzer::~PandaAnalyzer() {
  if (DEBUG) PDebug("PandaAnalyzer::~PandaAnalyzer","Calling destructor");
}


void PandaAnalyzer::ResetBranches() {
  genObjects.clear();
  matchPhos.clear();
  matchEles.clear();
  matchLeps.clear();
  gt->Reset();
  if (DEBUG) PDebug("PandaAnalyzer::ResetBranches","Reset");
}


void PandaAnalyzer::SetOutputFile(TString fOutName) {
  fOut = new TFile(fOutName,"RECREATE");
  tOut = new TTree("events","events");

  fOut->WriteTObject(hDTotalMCWeight);    

  gt->monohiggs = flags["monohiggs"];
  gt->vbf       = flags["vbf"];
  gt->fatjet    = flags["fatjet"];

  // fill the signal weights
  for (auto& id : wIDs) 
    gt->signal_weights[id] = 1;

  // Build the input tree here 
  gt->WriteTree(tOut);

  if (DEBUG) PDebug("PandaAnalyzer::SetOutputFile","Created output in "+fOutName);
}


int PandaAnalyzer::Init(TTree *t, TH1D *hweights, TTree *weightNames)
{
  if (DEBUG) PDebug("PandaAnalyzer::Init","Starting initialization");
  if (!t || !hweights) {
    PError("PandaAnalyzer::Init","Malformed input!");
    return 0;
  }
  tIn = t;

  event.setStatus(*t, {"!*"}); // turn everything off first

  TString jetname = (flags["puppi"]) ? "puppi" : "chs";
  panda::utils::BranchList readlist({"runNumber", "lumiNumber", "eventNumber", "rho", 
                                     "isData", "npv", "npvTrue", "weight", "chsAK4Jets", 
                                     "electrons", "muons", "taus", "photons", 
                                     "pfMet", "caloMet", "puppiMet", "rawMet", 
                                     "recoil","metFilters","genMet",});
  readlist.setVerbosity(0);

  if (flags["fatjet"])
   readlist += {jetname+"CA15Jets", "subjets", jetname+"CA15Subjets","Subjets"};
  
  if (flags["pfCands"])
    readlist.push_back("pfCandidates");

  if (isData) {
    readlist.push_back("triggers");
  } else {
   readlist.push_back("genParticles");
   readlist.push_back("genReweight");
  }


  event.setAddress(*t, readlist); // pass the readlist so only the relevant branches are turned on
  if (DEBUG) PDebug("PandaAnalyzer::Init","Set addresses");

  hDTotalMCWeight = new TH1F("hDTotalMCWeight","hDTotalMCWeight",1,0,2);
  hDTotalMCWeight->SetBinContent(1,hweights->GetBinContent(1));

  if (weightNames) {
    if (weightNames->GetEntries()!=377 && weightNames->GetEntries()!=22) {
      PError("PandaAnalyzer::Init",
          TString::Format("Reweighting failed because only found %u weights!",
                          unsigned(weightNames->GetEntries())));
      return 1;
    }
    TString *id = new TString();
    weightNames->SetBranchAddress("id",&id);
    unsigned nW = weightNames->GetEntriesFast();
    for (unsigned iW=0; iW!=nW; ++iW) {
      weightNames->GetEntry(iW);
      wIDs.push_back(*id);
    }
  } else if (processType==kSignal) {
    PError("PandaAnalyzer::Init","This is a signal file, but the weights are missing!");
    return 2;
  }


  // manipulate the output tree
  if (isData) {
    std::vector<TString> droppable = {"mcWeight","scale","scaleUp",
                                      "scaleDown","pdf.*","gen.*","sf_.*"};
    gt->RemoveBranches(droppable,{"sf_phoPurity"});
  }
  if (flags["genOnly"]) {
    std::vector<TString> keepable = {"mcWeight","scale","scaleUp",
                                     "scaleDown","pdf*","gen*","fj1*",
                                     "nFatjet","sf_tt*","sf_qcdTT*"};
    gt->RemoveBranches({".*"},keepable);
  }

  if (!flags["fatjet"]) {
    gt->RemoveBranches({"fj1.*"});
  } else if (flags["pfCands"]) {
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    double radius = 1.5;
    double sdZcut = 0.15;
    double sdBeta = 1.;
    activeArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
    areaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*activeArea);
    jetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,radius);
    softDrop = new fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);
  } else { 
    std::vector<TString> droppable = {"fj1NConst","fj1NSDConst","fj1EFrac100","fj1SDEFrac100"};
    gt->RemoveBranches(droppable);
  }

  if (DEBUG) PDebug("PandaAnalyzer::Init","Finished configuration");

  return 0;
}


panda::GenParticle const *PandaAnalyzer::MatchToGen(double eta, double phi, double radius, int pdgid) {
  panda::GenParticle const* found=NULL;
  double r2 = radius*radius;
  pdgid = abs(pdgid);

  unsigned int counter=0;
  for (map<panda::GenParticle const*,float>::iterator iG=genObjects.begin();
      iG!=genObjects.end(); ++iG) {
    if (found!=NULL)
      break;
    if (pdgid!=0 && abs(iG->first->pdgid)!=pdgid)
      continue;
    if (DeltaR2(eta,phi,iG->first->eta(),iG->first->phi())<r2)
      found = iG->first;
  }

  return found;
}


void PandaAnalyzer::Terminate() {
  fOut->WriteTObject(tOut);
  fOut->Close();

  for (auto *f : fCorrs)
    if (f)
      f->Close();
  for (auto *h : h1Corrs)
    delete h;
  for (auto *h : h2Corrs)
    delete h;

  delete btagCalib;
  delete sj_btagCalib;
  for (auto *reader : btagReaders )
    delete reader;

  for (auto& iter : ak8UncReader)
    delete iter.second;

  delete ak8JERReader;

  for (auto& iter : ak4UncReader)
    delete iter.second;

  for (auto& iter : ak4ScaleReader)
    delete iter.second;

  delete ak4JERReader;

  delete activeArea;
  delete areaDef;
  delete jetDef;
  delete softDrop;

  delete hDTotalMCWeight;
  if (DEBUG) PDebug("PandaAnalyzer::Terminate","Finished with output");
}

void PandaAnalyzer::OpenCorrection(CorrectionType ct, TString fpath, TString hname, int dim) {
  fCorrs[ct] = TFile::Open(fpath);
  if (dim==1) 
    h1Corrs[ct] = new THCorr1((TH1D*)fCorrs[ct]->Get(hname));
  else
    h2Corrs[ct] = new THCorr2((TH2D*)fCorrs[ct]->Get(hname));
}

double PandaAnalyzer::GetCorr(CorrectionType ct, double x, double y) {
  if (h1Corrs[ct]!=0) {
    return h1Corrs[ct]->Eval(x); 
  } else if (h2Corrs[ct]!=0) {
    return h2Corrs[ct]->Eval(x,y);
  } else {
    PError("PandaAnalyzer::GetCorr",
       TString::Format("No correction is defined for CorrectionType=%u",ct));
    return 1;
  }
}

void PandaAnalyzer::SetDataDir(const char *s) {
  TString dirPath(s);
  dirPath += "/";

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Starting loading of data");

  // pileup
  OpenCorrection(cNPV,dirPath+"moriond17/normalized_npv.root","data_npv_Wmn",1);
  OpenCorrection(cPU,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeights",1);

  // electrons
  OpenCorrection(cEleVeto,dirPath+"moriond17/scaleFactor_electron_summer16.root",
                 "scaleFactor_electron_vetoid_RooCMSShape_pu_0_100",2);
  OpenCorrection(cEleTight,dirPath+"moriond17/scaleFactor_electron_summer16.root",
                 "scaleFactor_electron_tightid_RooCMSShape_pu_0_100",2);
  OpenCorrection(cEleReco,dirPath+"moriond17/scaleFactor_electron_reco_summer16.root",
                 "scaleFactor_electron_reco_RooCMSShape_pu_0_100",2);

  // muons
  OpenCorrection(cMuLooseID,dirPath+"moriond17/muon_scalefactors_37ifb.root",
                 "scalefactors_MuonLooseId_Muon",2);
  OpenCorrection(cMuLooseIso,dirPath+"moriond17/muon_scalefactors_37ifb.root",
                 "scalefactors_Iso_MuonLooseId",2);
  OpenCorrection(cMuTightID,dirPath+"moriond17/muon_scalefactors_37ifb.root",
                 "scalefactors_TightId_Muon",2);
  OpenCorrection(cMuTightIso,dirPath+"moriond17/muon_scalefactors_37ifb.root",
                 "scalefactors_Iso_MuonTightId",2);
  OpenCorrection(cMuReco,dirPath+"moriond17/Tracking_12p9.root","htrack2",1);

  // photons
  OpenCorrection(cPho,dirPath+"moriond17/scalefactors_80x_medium_photon_37ifb.root",
                 "EGamma_SF2D",2);

  // triggers
  OpenCorrection(cTrigMET,dirPath+"moriond17/metTriggerEfficiency_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);
  OpenCorrection(cTrigEle,dirPath+"moriond17/eleTrig.root","hEffEtaPt",2);
  OpenCorrection(cTrigPho,dirPath+"moriond17/photonTriggerEfficiency_photon_TH1F.root",
                 "hden_photonpt_clone_passed",1);
  OpenCorrection(cTrigMETZmm,dirPath+"moriond17/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded scale factors");

  // kfactors
  TFile *fKFactor = new TFile(dirPath+"kfactors.root"); 
  fCorrs[cZNLO] = fKFactor; // just for garbage collection

  TH1D *hZLO    = (TH1D*)fKFactor->Get("ZJets_LO/inv_pt");
  TH1D *hWLO    = (TH1D*)fKFactor->Get("WJets_LO/inv_pt");
  TH1D *hALO    = (TH1D*)fKFactor->Get("GJets_LO/inv_pt_G");

  h1Corrs[cZNLO] = new THCorr1((TH1D*)fKFactor->Get("ZJets_012j_NLO/nominal"));
  h1Corrs[cWNLO] = new THCorr1((TH1D*)fKFactor->Get("WJets_012j_NLO/nominal"));
  h1Corrs[cANLO] = new THCorr1((TH1D*)fKFactor->Get("GJets_1j_NLO/nominal_G"));

  h1Corrs[cZEWK] = new THCorr1((TH1D*)fKFactor->Get("EWKcorr/Z"));
  h1Corrs[cWEWK] = new THCorr1((TH1D*)fKFactor->Get("EWKcorr/W"));
  h1Corrs[cAEWK] = new THCorr1((TH1D*)fKFactor->Get("EWKcorr/photon"));

  h1Corrs[cZEWK]->GetHist()->Divide(h1Corrs[cZNLO]->GetHist());     
  h1Corrs[cWEWK]->GetHist()->Divide(h1Corrs[cWNLO]->GetHist());     
  h1Corrs[cAEWK]->GetHist()->Divide(h1Corrs[cANLO]->GetHist());

  h1Corrs[cZNLO]->GetHist()->Divide(hZLO);    
  h1Corrs[cWNLO]->GetHist()->Divide(hWLO);    
  h1Corrs[cANLO]->GetHist()->Divide(hALO);

  OpenCorrection(cANLO2j,dirPath+"moriond17/histo_photons_2jet.root","Func",1);

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded k factors");

  TFile *fKFactor_VBFZ = new TFile(dirPath+"vbf16/kqcd/kfactor_VBF_zjets_v2.root");
  h1Corrs[cVBF_ZNLO] = new THCorr1((TH1D*)fKFactor_VBFZ->Get("bosonPt_NLO_vbf_relaxed"));
  h1Corrs[cVBF_ZNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFZ->Get("bosonPt_LO_vbf_relaxed"));
  h1Corrs[cVBF_ZNLO]->GetHist()->Multiply((TH1D*)fKFactor_VBFZ->Get("bosonPt_LO_monojet"));
  h1Corrs[cVBF_ZNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFZ->Get("bosonPt_NLO_monojet"));
  h1Corrs[cVBFTight_ZNLO] = new THCorr1((TH1D*)fKFactor_VBFZ->Get("bosonPt_NLO_vbf"));
  h1Corrs[cVBFTight_ZNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFZ->Get("bosonPt_LO_vbf"));
  h1Corrs[cVBFTight_ZNLO]->GetHist()->Multiply((TH1D*)fKFactor_VBFZ->Get("bosonPt_LO_monojet"));
  h1Corrs[cVBFTight_ZNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFZ->Get("bosonPt_NLO_monojet"));

  TFile *fKFactor_VBFW = new TFile(dirPath+"vbf16/kqcd/kfactor_VBF_wjets_v2.root");
  h1Corrs[cVBF_WNLO] = new THCorr1((TH1D*)fKFactor_VBFW->Get("bosonPt_NLO_vbf_relaxed"));
  h1Corrs[cVBF_WNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFW->Get("bosonPt_LO_vbf_relaxed"));
  h1Corrs[cVBF_WNLO]->GetHist()->Multiply((TH1D*)fKFactor_VBFW->Get("bosonPt_LO_monojet"));
  h1Corrs[cVBF_WNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFW->Get("bosonPt_NLO_monojet"));
  h1Corrs[cVBFTight_WNLO] = new THCorr1((TH1D*)fKFactor_VBFW->Get("bosonPt_NLO_vbf"));
  h1Corrs[cVBFTight_WNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFW->Get("bosonPt_LO_vbf"));
  h1Corrs[cVBFTight_WNLO]->GetHist()->Multiply((TH1D*)fKFactor_VBFW->Get("bosonPt_LO_monojet"));
  h1Corrs[cVBFTight_WNLO]->GetHist()->Divide((TH1D*)fKFactor_VBFW->Get("bosonPt_NLO_monojet"));

  OpenCorrection(cVBF_EWKZ,dirPath+"vbf16/kewk/kFactor_ZToNuNu_pT_Mjj.root",
                 "TH2F_kFactor",2);
  OpenCorrection(cVBF_EWKW,dirPath+"vbf16/kewk/kFactor_WToLNu_pT_Mjj.root",
                 "TH2F_kFactor",2);

  OpenCorrection(cVBF_TrigMET,dirPath+"vbf16/trig/metTriggerEfficiency_mjj_vbf.root",
                 "h_eff",2);
  OpenCorrection(cVBF_TrigMETZmm,dirPath+"vbf16/trig/metTriggerEfficiency_mjj_vbf_zmm.root",
                 "h_eff",2);

  OpenCorrection(cBadECALJets,dirPath+"vbf16/hotjets-runBCDEFGH.root",
                 "h2jet",2);

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded VBF k factors");

  // btag SFs
  btagCalib = new BTagCalibration("csvv2",(dirPath+"moriond17/CSVv2_Moriond17_B_H.csv").Data());
  btagReaders[bJetL] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
  btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
  btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
  btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

  sj_btagCalib = new BTagCalibration("csvv2",(dirPath+"moriond17/subjet_CSVv2_Moriond17_B_H.csv").Data());
  btagReaders[bSubJetL] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
  btagReaders[bSubJetL]->load(*sj_btagCalib,BTagEntry::FLAV_B,"lt");
  btagReaders[bSubJetL]->load(*sj_btagCalib,BTagEntry::FLAV_C,"lt");
  btagReaders[bSubJetL]->load(*sj_btagCalib,BTagEntry::FLAV_UDSG,"incl");

  //the Medium working
  btagReaders[bJetM] = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"});
  btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
  btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
  btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded btag SFs");

  // mSD corr
  MSDcorr = new TFile(dirPath+"/puppiCorr.root");
  puppisd_corrGEN = (TF1*)MSDcorr->Get("puppiJECcorr_gen");;
  puppisd_corrRECO_cen = (TF1*)MSDcorr->Get("puppiJECcorr_reco_0eta1v3");
  puppisd_corrRECO_for = (TF1*)MSDcorr->Get("puppiJECcorr_reco_1v3eta2v5");

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded mSD correction");

  TString jecV = "V4", jecReco = "23Sep2016"; 
  TString jecVFull = jecReco+jecV;
  ak8UncReader["MC"] = new JetCorrectionUncertainty(
     (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK8PFPuppi.txt").Data()
    );
  std::vector<TString> eraGroups = {"BCD","EF","G","H"};
  for (auto e : eraGroups) {
    ak8UncReader["data"+e] = new JetCorrectionUncertainty(
       (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK8PFPuppi.txt").Data()
      );
  }

  ak8JERReader = new JERReader(dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt",
                               dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt");


  ak4UncReader["MC"] = new JetCorrectionUncertainty(
     (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK4PFPuppi.txt").Data()
    );
  for (auto e : eraGroups) {
    ak4UncReader["data"+e] = new JetCorrectionUncertainty(
       (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK4PFPuppi.txt").Data()
      );
  }

  ak4JERReader = new JERReader(dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt",
                               dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt");

  std::vector<JetCorrectorParameters> params = {
    JetCorrectorParameters(
      (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L1FastJet_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2Relative_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L3Absolute_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2L3Residual_AK4PFPuppi.txt").Data())
  };
  ak4ScaleReader["MC"] = new FactorizedJetCorrector(params);
  if (DEBUG>1) PDebug("PandaAnalyzer::SetDataDir","Loaded JES for AK4 MC");
  for (auto e : eraGroups) {
    params = {
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L1FastJet_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2Relative_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L3Absolute_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2L3Residual_AK4PFPuppi.txt").Data())
    };
    ak4ScaleReader["data"+e] = new FactorizedJetCorrector(params);
    if (DEBUG>1) PDebug("PandaAnalyzer::SetDataDir","Loaded JES for AK4 "+e);
  }

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded JES/R");

}


void PandaAnalyzer::AddGoodLumiRange(int run, int l0, int l1) {
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) { // don't know about this run yet
    std::vector<LumiRange> newLumiList;
    newLumiList.emplace_back(l0,l1);
    goodLumis[run] = newLumiList;
  } else {
    run_->second.emplace_back(l0,l1);
  }
}


bool PandaAnalyzer::PassGoodLumis(int run, int lumi) {
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) {
    // matched no run
    if (DEBUG) 
      PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
    return false;
  }

  // found the run, now look for a lumi range
  for (auto &range : run_->second) {
    if (range.Contains(lumi)) {
      if (DEBUG) 
        PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
      return true;
    }
  }

  // matched no lumi range
  if (DEBUG) 
    PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
  return false;
}


bool PandaAnalyzer::PassPreselection() {
  // TODO: refactor this function
  // was originally written this way to handle more complex conditions
  // like triggers, but could probably clean it up with a Condition class
  
  if (preselBits==0)
    return true;
  bool isGood=false;

  if (preselBits & kFatjet) {
    if (gt->fj1Pt>250)
      isGood = true;
  }

  float max_puppi = std::max({gt->puppimet, gt->puppiUZmag, gt->puppiUWmag, gt->puppiUAmag});
  float max_pf = std::max({gt->pfmet, gt->pfUZmag, gt->pfUWmag, gt->pfUAmag, gt->pfUWWmag});
  float max_pfUp = std::max({gt->pfmetUp, gt->pfUZmagUp, gt->pfUWmagUp, gt->pfUAmagUp, gt->pfUWWmagUp});
  float max_pfDown = std::max({gt->pfmetDown, gt->pfUZmagDown, gt->pfUWmagDown, gt->pfUAmagDown, gt->pfUWWmagDown});

  if (preselBits & kRecoil) {
    if ( max_pfDown>200 || max_pf>200 || max_pfUp>200 || max_puppi>200 ) {
      isGood = true;
    }
  }
  if (preselBits & kMonotop) {
    if (gt->nFatjet>=1 && gt->fj1Pt>200) {
      if ( max_pf>200 || max_puppi>200) {
        isGood = true;
      }
    }
  }
  if (preselBits & kMonojet) {
    if (gt->nJet>=1 && gt->jet1Pt>100) {
      if ( max_pfDown>250 || max_pf>250 || max_pfUp>250 || max_puppi>250 ) {
        isGood = true;
      }
    }
  }
  if (preselBits & kMonohiggs) {
    if (gt->nFatjet>=1 && gt->fj1Pt>150) {
      if ( max_pf>250 || max_puppi>250) {
        isGood = true;
      }
    }
  }

  return isGood;
}


void PandaAnalyzer::CalcBJetSFs(BTagType bt, int flavor,
                double eta, double pt, double eff, double uncFactor,
                double &sf, double &sfUp, double &sfDown) 
{
  if (flavor==5) {
    sf     = btagReaders[bt]->eval_auto_bounds("central",BTagEntry::FLAV_B,eta,pt);
    sfUp   = btagReaders[bt]->eval_auto_bounds("up",BTagEntry::FLAV_B,eta,pt);
    sfDown = btagReaders[bt]->eval_auto_bounds("down",BTagEntry::FLAV_B,eta,pt);
  } else if (flavor==4) {
    sf     = btagReaders[bt]->eval_auto_bounds("central",BTagEntry::FLAV_C,eta,pt);
    sfUp   = btagReaders[bt]->eval_auto_bounds("up",BTagEntry::FLAV_C,eta,pt);
    sfDown = btagReaders[bt]->eval_auto_bounds("down",BTagEntry::FLAV_C,eta,pt);
  } else {
    sf     = btagReaders[bt]->eval_auto_bounds("central",BTagEntry::FLAV_UDSG,eta,pt);
    sfUp   = btagReaders[bt]->eval_auto_bounds("up",BTagEntry::FLAV_UDSG,eta,pt);
    sfDown = btagReaders[bt]->eval_auto_bounds("down",BTagEntry::FLAV_UDSG,eta,pt);
  }

  sfUp = uncFactor*(sfUp-sf)+sf;
  sfDown = uncFactor*(sfDown-sf)+sf;
  return;
}

void PandaAnalyzer::EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
			       GeneralTree::BTagShift shift,GeneralTree::BTagJet jettype, bool do2) 
{
  float sf0 = 1, sf1 = 1, sfGT0 = 1, sf2=1;
  float prob_mc0=1, prob_data0=1;
  float prob_mc1=0, prob_data1=0;
  unsigned int nC = cands.size();

  for (unsigned int iC=0; iC!=nC; ++iC) {
    double sf_i = sfs[iC];
    double eff_i = cands[iC].eff;
    prob_mc0 *= (1-eff_i); // denominator for 0 b-tag
    prob_data0 *= (1-sf_i*eff_i); // nominator for 0 b-tag
    float tmp_mc1=1, tmp_data1=1;
    for (unsigned int jC=0; jC!=nC; ++jC) {
      if (iC==jC) continue;
      double sf_j = sfs[jC];
      double eff_j = cands[jC].eff;
      tmp_mc1 *= (1-eff_j); 
      tmp_data1 *= (1-eff_j*sf_j);
    }
    prob_mc1 += eff_i * tmp_mc1; // denominator for 1 b-tag
    prob_data1 += eff_i * sf_i * tmp_data1; // nominator for 1 b-tag
  }
  
  if (nC>0) {
    sf0 = prob_data0/prob_mc0;
    sf1 = prob_data1/prob_mc1;
    sfGT0 = (1-prob_data0)/(1-prob_mc0); // all b-tag (negating 0 b-tag)
  }

  GeneralTree::BTagParams p;
  p.shift = shift;
  p.jet = jettype;
  p.tag=GeneralTree::b0; gt->sf_btags[p] = sf0;
  p.tag=GeneralTree::b1; gt->sf_btags[p] = sf1;
  p.tag=GeneralTree::bGT0; gt->sf_btags[p] = sfGT0;

  if (do2) {
    float prob_mc2=0, prob_data2=0;
    unsigned int nC = cands.size();

    for (unsigned int iC=0; iC!=nC; ++iC) {
      double sf_i = sfs[iC], eff_i = cands[iC].eff;
      for (unsigned int jC=iC+1; jC!=nC; ++jC) {
        double sf_j = sfs[jC], eff_j = cands[jC].eff;
        float tmp_mc2=1, tmp_data2=1;
        for (unsigned int kC=0; kC!=nC; ++kC) {
          if (kC==iC || kC==jC) continue;
          double sf_k = sfs[kC], eff_k = cands[kC].eff;
          tmp_mc2 *= (1-eff_k); 
          tmp_data2 *= (1-eff_k*sf_k);
        }
        prob_mc2 += eff_i * eff_j * tmp_mc2; // denominator for 2 b-tag
        prob_data2 += eff_i * sf_i * eff_j * sf_j * tmp_data2; // nominator for 2 b-tag
      }
    }

    if (nC>1) {
      sf2 = prob_data2/prob_mc2;
    }

    p.tag=GeneralTree::b2; gt->sf_btags[p] = sf2;
  }

}

float PandaAnalyzer::GetMSDCorr(Float_t puppipt, Float_t puppieta) {

  float genCorr   = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  genCorr = puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta) <= 1.3 ){
    recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  }
  else{
    recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  }
  totalWeight = genCorr * recoCorr;

  return totalWeight;
}

void PandaAnalyzer::RegisterTrigger(TString path, std::vector<unsigned> &idxs) {
  unsigned idx = event.registerTrigger(path);
  if (DEBUG>1) PDebug("PandaAnalyzer::RegisterTrigger",
            TString::Format("At %u found trigger=%s",idx,path.Data()));
  idxs.push_back(idx);
}

// run
void PandaAnalyzer::Run() {

  // INITIALIZE --------------------------------------------------------------------------

  unsigned int nEvents = tIn->GetEntries();
  unsigned int nZero = 0;
  if (lastEvent>=0 && lastEvent<(int)nEvents)
    nEvents = lastEvent;
  if (firstEvent>=0)
    nZero = firstEvent;

  if (!fOut || !tIn) {
    PError("PandaAnalyzer::Run","NOT SETUP CORRECTLY");
    exit(1);
  }

  // get bounds
  float genBosonPtMin=150, genBosonPtMax=1000;
  if (!isData) {
    genBosonPtMin = h1Corrs[cZNLO]->GetHist()->GetBinCenter(1);
    genBosonPtMax = h1Corrs[cZNLO]->GetHist()->GetBinCenter(h1Corrs[cZNLO]->GetHist()->GetNbinsX());
  }

  panda::FatJetCollection* fatjets(0);
  if (flags["fatjet"]) {
   if (flags["puppi"])
    fatjets = &event.puppiCA15Jets;
   else
    fatjets = &event.chsCA15Jets;
  }

  panda::JetCollection* jets(0);
  // seems like now we always use chs? - yeah this was overridden to be consistent with PF MET
  jets = &event.chsAK4Jets;

  // these are bins of b-tagging eff in pT and eta, derived in 8024 TT MC
  // TODO: don't hardcode these 
  std::vector<double> vbtagpt {20.0,50.0,80.0,120.0,200.0,300.0,400.0,500.0,700.0,1000.0};
  std::vector<double> vbtageta {0.0,0.5,1.5,2.5};
  std::vector<std::vector<double>> lfeff  = {{0.081,0.065,0.060,0.063,0.072,0.085,0.104,0.127,0.162},
                       {0.116,0.097,0.092,0.099,0.112,0.138,0.166,0.185,0.222},
                       {0.173,0.145,0.149,0.175,0.195,0.225,0.229,0.233,0.250}};
  std::vector<std::vector<double>> ceff = {{0.377,0.389,0.391,0.390,0.391,0.375,0.372,0.392,0.435},
                      {0.398,0.407,0.416,0.424,0.424,0.428,0.448,0.466,0.500},
                      {0.375,0.389,0.400,0.425,0.437,0.459,0.481,0.534,0.488}};
  std::vector<std::vector<double>> beff = {{0.791,0.815,0.825,0.835,0.821,0.799,0.784,0.767,0.760},
                      {0.794,0.816,0.829,0.836,0.823,0.804,0.798,0.792,0.789},
                      {0.739,0.767,0.780,0.789,0.776,0.771,0.779,0.787,0.806}};
  Binner btagpt(vbtagpt);
  Binner btageta(vbtageta);

  JetCorrectionUncertainty *uncReader=0;
  JetCorrectionUncertainty *uncReaderAK4=0;
  FactorizedJetCorrector *scaleReaderAK4=0;

  std::vector<unsigned int> metTriggers;
  std::vector<unsigned int> eleTriggers;
  std::vector<unsigned int> phoTriggers;
  std::vector<unsigned int> muTriggers;

  if (isData) {
    std::vector<TString> metTriggerPaths = {
          "HLT_PFMET170_NoiseCleaned",
          "HLT_PFMET170_HBHECleaned",
          "HLT_PFMET170_JetIdCleaned",
          "HLT_PFMET170_NotCleaned",
          "HLT_PFMET170_HBHE_BeamHaloCleaned",
          "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight",
          "HLT_PFMETNoMu110_NoiseCleaned_PFMHTNoMu110_IDTight",
          "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight",
          "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight",
          "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight",
          "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"
    };
    std::vector<TString> eleTriggerPaths = {
//                    "HLT_Ele25_eta2p1_WPTight_Gsf",
//                     "HLT_Ele27_eta2p1_WPLoose_Gsf",

//                    "HLT_Ele23_CaloIdL_TrackIdL_IsoVL",
//                    "HLT_Ele22_eta2p1_WP75_Gsf",
//                    "HLT_Ele23_WPLoose_Gsf",
          "HLT_Ele27_WP85_Gsf",
          "HLT_Ele27_WPLoose_Gsf",
          "HLT_Ele105_CaloIdVT_GsfTrkIdT",
          "HLT_Ele27_WPTight_Gsf",
          "HLT_Ele30_WPTight_Gsf",
          "HLT_Ele27_eta2p1_WPTight_Gsf",
          "HLT_Ele32_eta2p1_WPTight_Gsf",
          "HLT_Ele35_WPLoose_Gsf",
//                    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",

          "HLT_ECALHT800"
    };

    std::vector<TString> phoTriggerPaths = {
          "HLT_Photon175",
          "HLT_Photon165_HE10",
          "HLT_Photon36_R9Id90_HE10_IsoM",
          "HLT_Photon50_R9Id90_HE10_IsoM",
          "HLT_Photon75_R9Id90_HE10_IsoM",
          "HLT_Photon90_R9Id90_HE10_IsoM",
          "HLT_Photon120_R9Id90_HE10_IsoM",
          "HLT_Photon165_R9Id90_HE10_IsoM",
          "HLT_Photon300_NoHE",

          "HLT_ECALHT800"
    };

    if (DEBUG>1) PDebug("PandaAnalyzer::Run","Loading MET triggers");
    for (auto path : metTriggerPaths) {
      RegisterTrigger(path,metTriggers);
    }
    if (DEBUG>1) PDebug("PandaAnalyzer::Run","Loading SingleElectron triggers");
    for (auto path : eleTriggerPaths) {
      RegisterTrigger(path,eleTriggers);
    }
    if (DEBUG>1) PDebug("PandaAnalyzer::Run","Loading SinglePhoton triggers");
    for (auto path : phoTriggerPaths) {
      RegisterTrigger(path,phoTriggers);
    }

  }

  float EGMSCALE = isData ? 1 : 1;

  // set up reporters
  unsigned int iE=0;
  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,10);
  TimeReporter tr("PandaAnalyzer::Run",DEBUG);

  bool applyJSON = flags["applyJSON"];
  bool doMonoH = flags["monohiggs"];
  bool doVBF = flags["vbf"];
  bool doFatjet = flags["fatjet"];
  

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();
    ResetBranches();
    event.getEntry(*tIn,iE);

    tr.TriggerEvent(TString::Format("GetEntry %u",iE));
    if (DEBUG>2) {
      PDebug("PandaAnalyzer::Run::Dump","");
      event.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.photons.print(std::cout, 2);       // photon branch
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.muons.print(std::cout, 2);         // muons branch
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.electrons.print(std::cout, 2);     // electron branch
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.chsAK4Jets.print(std::cout, 2);    // chsAK4Jets branch
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.pfMet.print(std::cout, 2);         // pfMet branch
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.metMuOnlyFix.print(std::cout, 2);  // metMuOnlyFix branch
      std::cout << std::endl;
    }

    if ( (preselBits&kMonotop) || (preselBits&kMonohiggs) || 
      (preselBits&kMonojet) || (preselBits&kRecoil) ) 
    {
     if (event.recoil.max<175)
       continue;
    } 

    // event info
    //gt->mcWeight = (event.weight>0) ? 1 : -1;
    gt->mcWeight = event.weight;
    gt->runNumber = event.runNumber;
    gt->lumiNumber = event.lumiNumber;
    gt->eventNumber = event.eventNumber;
    gt->npv = event.npv;
    gt->pu = event.npvTrue;
    gt->metFilter = (event.metFilters.pass()) ? 1 : 0;
    // these two are not need since we use muon-fixed MET
    // gt->metFilter = (gt->metFilter==1 && !event.metFilters.badMuons) ? 1 : 0;
    // gt->metFilter = (gt->metFilter==1 && !event.metFilters.duplicateMuons) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badPFMuons) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badChargedHadrons) ? 1 : 0;

    if (isData) {
      // check the json
      if (applyJSON && !PassGoodLumis(gt->runNumber,gt->lumiNumber))
        continue;

      // save triggers
      for (auto iT : metTriggers) {
       if (event.triggerFired(iT)) {
        gt->trigger |= kMETTrig;
        break;
       }
      }
      for (auto iT : eleTriggers) {
       if (event.triggerFired(iT)) {
        gt->trigger |= kSingleEleTrig;
        break;
       }
      }
      for (auto iT : phoTriggers) {
       if (event.triggerFired(iT)) {
        gt->trigger |= kSinglePhoTrig;
        break;
       }
      }
    } else {
      gt->sf_npv = GetCorr(cNPV,gt->npv);
      gt->sf_pu = GetCorr(cPU,gt->pu);
    }

    if (uncReader==0) {
      if (isData) {
        TString thisEra = eras.getEra(gt->runNumber);
        for (auto &iter : ak8UncReader) {
          if (! iter.first.Contains("data"))
            continue;
          if (iter.first.Contains(thisEra)) {
            uncReader = iter.second;
            uncReaderAK4 = ak4UncReader[iter.first];
            scaleReaderAK4 = ak4ScaleReader[iter.first];
            break;
          }
        }
      } else {
        uncReader = ak8UncReader["MC"];
        uncReaderAK4 = ak4UncReader["MC"];
        scaleReaderAK4 = ak4ScaleReader["MC"];
      }
    }

    tr.TriggerEvent("initialize");

    // met
    gt->pfmetRaw = event.rawMet.pt;
    gt->pfmet = event.pfMet.pt;
    gt->pfmetphi = event.pfMet.phi;
    gt->pfmetUp = event.pfMet.ptCorrUp;
    gt->pfmetDown = event.pfMet.ptCorrDown;
    gt->calomet = event.caloMet.pt;
    gt->sumETRaw = event.pfMet.sumETRaw;
    gt->puppimet = event.puppiMet.pt;
    gt->puppimetphi = event.puppiMet.phi;
    TLorentzVector vPFMET, vPuppiMET;
    vPFMET.SetPtEtaPhiM(gt->pfmet,0,gt->pfmetphi,0);
    vPuppiMET.SetPtEtaPhiM(gt->puppimet,0,gt->puppimetphi,0);
    TVector2 vMETNoMu; vMETNoMu.SetMagPhi(gt->pfmet,gt->pfmetphi); //       for trigger eff

    tr.TriggerEvent("met");

    gt->isGS = 0;

    //electrons
    std::vector<panda::Lepton*> looseLeps, tightLeps;
    for (auto& ele : event.electrons) {
     float pt = ele.pt()*EGMSCALE; float eta = ele.eta(); float aeta = fabs(eta);
      if (pt<10 || aeta>2.5 /* || (aeta>1.4442 && aeta<1.566) */)
        continue;
      if (!ele.veto)
        continue;
      if (!ElectronIP(ele.eta(),ele.dxy,ele.dz))
        continue;
      looseLeps.push_back(&ele);
      if (doVBF) {
        matchLeps.push_back(&ele);
        matchEles.push_back(&ele);
      }
      gt->nLooseElectron++;
    }

    // muons
    for (auto& mu : event.muons) {
     float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
      if (pt<10 || aeta>2.4)
        continue;
      if (!mu.loose)
        continue;
      if (!MuonIsolation(pt,eta,mu.combIso(),panda::kLoose))
        continue;
      looseLeps.push_back(&mu);
      if (doVBF)
        matchLeps.push_back(&mu);
      gt->nLooseMuon++;
      TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
      vMETNoMu += vMu;
    }
    gt->pfmetnomu = vMETNoMu.Mod();

    // now consider all leptons
    gt->nLooseLep = looseLeps.size();
    if (gt->nLooseLep>0) {
     auto ptsort([](panda::Lepton const* l1, panda::Lepton const* l2)->bool {
       return l1->pt() > l2->pt();
      });
     int nToSort = TMath::Min(3,gt->nLooseLep);
     std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),ptsort);
    }
    int lep_counter=1;
    for (auto* lep : looseLeps) {
      if (lep_counter==1) {
       gt->looseLep1Pt = lep->pt();
       gt->looseLep1Eta = lep->eta();
       gt->looseLep1Phi = lep->phi();
      } else if (lep_counter==2) {
	gt->looseLep2Pt = lep->pt();
	gt->looseLep2Eta = lep->eta();
	gt->looseLep2Phi = lep->phi();
      } else {
        break;
      }
      // now specialize lepton types
      panda::Muon *mu = dynamic_cast<panda::Muon*>(lep);
      if (mu!=NULL) {
        bool isTight = ( mu->tight &&
			 MuonIsolation(mu->pt(),mu->eta(),mu->combIso(),panda::kTight) &&
			 mu->pt()>20 && fabs(mu->eta())<2.4 );
        if (lep_counter==1) {
          gt->looseLep1PdgId = mu->charge*-13;
          gt->looseLep1IsHLTSafe = 1;
          if (isTight) {
            gt->nTightMuon++;
            gt->looseLep1IsTight = 1;
            if (!doVBF)
              matchLeps.push_back(lep);
          }
        } else if (lep_counter==2) {
          gt->looseLep2PdgId = mu->charge*-13;
          gt->looseLep2IsHLTSafe = 1;
          if (isTight) {
            gt->nTightMuon++;
            gt->looseLep2IsTight = 1;
          }
          if (!doVBF && (isTight || gt->looseLep1IsTight))
            matchLeps.push_back(lep);
        }
      } else {
        panda::Electron *ele = dynamic_cast<panda::Electron*>(lep);
        bool isTight = ( ele->tight &&
			 /*ElectronIsolation(ele->pt,ele->eta,ele->iso,PElectron::kTight) &&*/
			 ele->pt()>40 && fabs(ele->eta())<2.5 );
        if (lep_counter==1) {
          gt->looseLep1Pt *= EGMSCALE;
          gt->looseLep1PdgId = ele->charge*-11;
          gt->looseLep1IsHLTSafe = ele->hltsafe ? 1 : 0;
          if (isTight) {
            gt->nTightElectron++;
            gt->looseLep1IsTight = 1;
            if (!doVBF) {
              matchLeps.push_back(lep);
              matchEles.push_back(lep);
            }
          }
        } else if (lep_counter==2) {
          gt->looseLep2Pt *= EGMSCALE;
          gt->looseLep2PdgId = ele->charge*-11;
          gt->looseLep2IsHLTSafe = ele->hltsafe ? 1 : 0;
          if (isTight) {
            gt->nTightElectron++;
            gt->looseLep2IsTight = 1;
          }
          if (!doVBF && (isTight || gt->looseLep1IsTight)) {
            matchLeps.push_back(lep);
            matchEles.push_back(lep);
          }
        }
      }
      ++lep_counter;
    }
    gt->nTightLep = gt->nTightElectron + gt->nTightMuon;
    if (gt->nLooseLep>0) {
      panda::Lepton* lep1 = looseLeps[0];
      gt->mT = MT(lep1->pt(),lep1->phi(),gt->pfmet,gt->pfmetphi);
    }
    if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId==0) {
      TLorentzVector v1,v2;
      panda::Lepton *lep1=looseLeps[0], *lep2=looseLeps[1];
      v1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());
      v2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());
      gt->diLepMass = (v1+v2).M();
    } else {
      gt->diLepMass = -1;
    }
    
    tr.TriggerEvent("leptons");
    
    // photons
    std::vector<panda::Photon*> loosePhos;
    for (auto& pho : event.photons) {
      if (!pho.loose || !pho.csafeVeto)
        continue;
      float pt = pho.pt() * EGMSCALE;
      if (pt<1) continue;
      float eta = pho.eta(), phi = pho.phi();
      if (pt<15 || fabs(eta)>2.5)
        continue;
      /*
	if (IsMatched(&matchEles,0.16,eta,phi))
        continue;
      */
      loosePhos.push_back(&pho);
      gt->nLoosePhoton++;
      if (gt->nLoosePhoton==1) {
        gt->loosePho1Pt = pt;
        gt->loosePho1Eta = eta;
        gt->loosePho1Phi = phi;
      }
      if ( pho.medium &&
	   pt>175 /*&& fabs(eta)<1.4442*/ ) { // apply eta cut offline
        if (gt->nLoosePhoton==1)
          gt->loosePho1IsTight=1;
        gt->nTightPhoton++;
        matchPhos.push_back(&pho);
      }
    }
    
    // TODO - store in a THCorr
    if (isData && gt->nLoosePhoton>0) {
      if (gt->loosePho1Pt>=175 && gt->loosePho1Pt<200)
        gt->sf_phoPurity = 0.04802;
      else if (gt->loosePho1Pt>=200 && gt->loosePho1Pt<250)
        gt->sf_phoPurity = 0.04241;
      else if (gt->loosePho1Pt>=250 && gt->loosePho1Pt<300)
        gt->sf_phoPurity = 0.03641;
      else if (gt->loosePho1Pt>=300 && gt->loosePho1Pt<350)
        gt->sf_phoPurity = 0.0333;
      else if (gt->loosePho1Pt>=350)
        gt->sf_phoPurity = 0.02544;
    }
    
    tr.TriggerEvent("photons");
    
    // trigger efficiencies
    gt->sf_eleTrig=1; gt->sf_metTrig=1; gt->sf_phoTrig=1; gt->sf_metTrigZmm=1;
    gt->sf_metTrigVBF=1; gt->sf_metTrigZmmVBF=1;
    if (!isData) {
      gt->sf_metTrig = GetCorr(cTrigMET,gt->pfmetnomu);
      gt->sf_metTrigZmm = GetCorr(cTrigMETZmm,gt->pfmetnomu);
      
      if (gt->nLooseElectron>0 && abs(gt->looseLep1PdgId)==11
          && gt->looseLep1IsTight==1) {
        float eff1=0, eff2=0;
        eff1 = GetCorr(cTrigEle,gt->looseLep1Eta,gt->looseLep1Pt);
        if (gt->nLooseElectron>1 && abs(gt->looseLep2PdgId)==11) {
          eff2 = GetCorr(cTrigEle,gt->looseLep2Eta,gt->looseLep2Pt);
        }
        gt->sf_eleTrig = 1 - (1-eff1)*(1-eff2);
      } // done with ele trig SF

      if (gt->nLoosePhoton>0 && gt->loosePho1IsTight)
        gt->sf_phoTrig = GetCorr(cTrigPho,gt->loosePho1Pt);
    }

    tr.TriggerEvent("triggers");

    // recoil!
    TLorentzVector vpfUp; vpfUp.SetPtEtaPhiM(gt->pfmetUp,0,gt->pfmetphi,0);
    TLorentzVector vpfDown; vpfDown.SetPtEtaPhiM(gt->pfmetDown,0,gt->pfmetphi,0);
    TLorentzVector vObj1, vObj2;
    TLorentzVector vpuppiUW, vpuppiUZ, vpuppiUA;
    TLorentzVector vpfUW, vpfUZ, vpfUA, vpfUWW;
    TLorentzVector vpuppiU, vpfU;
    int whichRecoil = 0; // -1=photon, 0=MET, 1,2=nLep, 3 = WW
    if (gt->nLooseLep>0) {
      panda::Lepton *lep1 = looseLeps.at(0);
      vObj1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());

      // one lep => W
      vpuppiUW = vPuppiMET+vObj1; gt->puppiUWmag=vpuppiUW.Pt(); gt->puppiUWphi=vpuppiUW.Phi();
      vpfUW = vPFMET+vObj1; gt->pfUWmag=vpfUW.Pt(); gt->pfUWphi=vpfUW.Phi();
      
      TLorentzVector vpfUWUp = vpfUp+vObj1; gt->pfUWmagUp = vpfUWUp.Pt();
      TLorentzVector vpfUWDown = vpfDown+vObj1; gt->pfUWmagDown = vpfUWDown.Pt();

      if (gt->nLooseLep>1 && gt->looseLep1PdgId+gt->looseLep2PdgId==0) {
        // two OS lep => Z
        panda::Lepton *lep2 = looseLeps.at(1);
        vObj2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());

        vpuppiUZ=vpuppiUW+vObj2; gt->puppiUZmag=vpuppiUZ.Pt(); gt->puppiUZphi=vpuppiUZ.Phi();
        vpfUZ=vpfUW+vObj2; gt->pfUZmag=vpfUZ.Pt(); gt->pfUZphi=vpfUZ.Phi();

        TLorentzVector vpfUZUp = vpfUWUp+vObj2; gt->pfUZmagUp = vpfUZUp.Pt();
        TLorentzVector vpfUZDown = vpfUWDown+vObj2; gt->pfUZmagDown = vpfUZDown.Pt();

        vpuppiU = vpuppiUZ; vpfU = vpfUZ;
        whichRecoil = 2;
      } else if ( gt->nLooseLep>1 && ( (gt->looseLep1PdgId==11 && gt->looseLep2PdgId==-13) 
				       || ( gt->looseLep1PdgId==-11 && gt->looseLep2PdgId==13 ) 
				       || ( gt->looseLep1PdgId==13 && gt->looseLep2PdgId==-11 ) 
				       || ( gt->looseLep1PdgId==-13 && gt->looseLep2PdgId==11 )  
				       ) ){
	// two W with opposite flavor
	panda::Lepton *lep2 = looseLeps.at(1);
	vObj2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());

        vpfUWW=vpfUW+vObj2; gt->pfUWWmag=vpfUWW.Pt(); gt->pfUWWphi=vpfUWW.Phi();

        TLorentzVector vpfUWWUp = vpfUWUp+vObj2; gt->pfUWWmagUp = vpfUWWUp.Pt();
        TLorentzVector vpfUWWDown = vpfUWDown+vObj2; gt->pfUWWmagDown = vpfUWWDown.Pt();
	vpfU = vpfUWW;
	whichRecoil = 3;
      }
      else {
        vpuppiU = vpuppiUW; vpfU = vpfUW;
        whichRecoil = 1;
      }
    }
    if (gt->nLoosePhoton>0) {
      panda::Photon *pho = loosePhos.at(0);
      vObj1.SetPtEtaPhiM(pho->pt(),pho->eta(),pho->phi(),0.);

      vpuppiUA=vPuppiMET+vObj1; gt->puppiUAmag=vpuppiUA.Pt(); gt->puppiUAphi=vpuppiUA.Phi();
      vpfUA=vPFMET+vObj1; gt->pfUAmag=vpfUA.Pt(); gt->pfUAphi=vpfUA.Phi();

      TLorentzVector vpfUAUp = vpfUp+vObj1; gt->pfUAmagUp = vpfUAUp.Pt();
      TLorentzVector vpfUADown = vpfDown+vObj1; gt->pfUAmagDown = vpfUADown.Pt();

      if (gt->nLooseLep==0) {
        vpuppiU = vpuppiUA; vpfU = vpfUA;
        whichRecoil = -1;
      }
    }
    if (gt->nLooseLep==0 && gt->nLoosePhoton==0) {
      vpuppiU = vPuppiMET;
      vpfU = vPFMET;
      whichRecoil = 0;
    }
    gt->puppiUmag = vpuppiU.Pt();
    gt->puppiUphi = vpuppiU.Phi();
    gt->pfUmag = vpfU.Pt();
    gt->pfUphi = vpfU.Phi();

    tr.TriggerEvent("recoils");
    /*
    std::vector<int> LepIds = {11,13,15,-11,-13,-15};
    std::vector<int> NeuIds = {12,14,16,-12,-14,-16};
    std::vector<int> HadIds = {1,2,3,4,5,-1,-2,-3,-4,-5};
    int nGen = event.genParticles.size();
    for (int iG=0; iG!=nGen; ++iG) {
      auto& part(event.genParticles.at(iG));
      int pdgid = part.pdgid;
      unsigned int abspdgid = abs(pdgid);
      bool isLastCopy=true;
      for (int jG=0; jG!=nGen; ++jG) {
	if (event.genParticles.at(jG).parent.get() == &part) {
	  isLastCopy=false;
	  break;
	}
      }
      if (!isLastCopy)
	continue;
    } //looking for targets     
    
    tr.TriggerEvent("GenParticle analysis");
    */
    panda::FatJet *fj1=0;
    gt->nFatjet=0;
    if (doFatjet) {
      int fatjet_counter=-1;
      for (auto& fj : *fatjets) {
        ++fatjet_counter;
        float pt = fj.pt();
        float rawpt = fj.rawPt;
        float eta = fj.eta();
        float mass = fj.m();
        float ptcut = 150;
        if (doMonoH)
          ptcut = 150;
	
        if (pt<ptcut || fabs(eta)>2.4 || !fj.monojet)
          continue;
	
        float phi = fj.phi();
        if (IsMatched(&matchLeps,2.25,eta,phi) || IsMatched(&matchPhos,2.25,eta,phi)) {
          continue;
        }
	
        gt->nFatjet++;
        if (gt->nFatjet==1) {
          fj1 = &fj;
          if (fatjet_counter==0)
            gt->fj1IsClean = 1;
          else
            gt->fj1IsClean = 0;
          gt->fj1Pt = pt;
          gt->fj1Eta = eta;
          gt->fj1Phi = phi;
          gt->fj1M = mass;
          gt->fj1MSD = fj.mSD;
          gt->fj1RawPt = rawpt;

          // do a bit of jet energy scaling
          // uncReader->setJetEta(eta); uncReader->setJetPt(pt);
          // double scaleUnc = uncReader->getUncertainty(true);
          double scaleUnc = (fj.ptCorrUp - gt->fj1Pt) / gt->fj1Pt; 
          gt->fj1PtScaleUp    = gt->fj1Pt  * (1 + 2*scaleUnc);
          gt->fj1PtScaleDown  = gt->fj1Pt  * (1 - 2*scaleUnc);
          gt->fj1MSDScaleUp   = gt->fj1MSD * (1 + 2*scaleUnc);
          gt->fj1MSDScaleDown = gt->fj1MSD * (1 - 2*scaleUnc);

          // do some jet energy smearing
          if (isData) {
            gt->fj1PtSmeared = gt->fj1Pt;
            gt->fj1PtSmearedUp = gt->fj1Pt;
            gt->fj1PtSmearedDown = gt->fj1Pt;
            gt->fj1MSDSmeared = gt->fj1MSD;
            gt->fj1MSDSmearedUp = gt->fj1MSD;
            gt->fj1MSDSmearedDown = gt->fj1MSD;
          } else {
            double smear=1, smearUp=1, smearDown=1;
            ak8JERReader->getStochasticSmear(pt,eta,event.rho,smear,smearUp,smearDown);

            gt->fj1PtSmeared = smear*gt->fj1Pt;
            gt->fj1PtSmearedUp = smearUp*gt->fj1Pt;
            gt->fj1PtSmearedDown = smearDown*gt->fj1Pt;

            gt->fj1MSDSmeared = smear*gt->fj1MSD;
            gt->fj1MSDSmearedUp = smearUp*gt->fj1MSD;
            gt->fj1MSDSmearedDown = smearDown*gt->fj1MSD;
          }

          // now have to do this mess with the subjets...
          TLorentzVector sjSum, sjSumUp, sjSumDown, sjSumSmear;
          for (unsigned int iSJ=0; iSJ!=fj.subjets.size(); ++iSJ) {
            auto& subjet = fj.subjets.objAt(iSJ);
            // now correct...
            double factor=1;
            if (fabs(subjet.eta())<5.191) {
              scaleReaderAK4->setJetPt(subjet.pt());
              scaleReaderAK4->setJetEta(subjet.eta());
              scaleReaderAK4->setJetPhi(subjet.phi());
              scaleReaderAK4->setJetE(subjet.e());
              scaleReaderAK4->setRho(event.rho);
              scaleReaderAK4->setJetA(0);
              scaleReaderAK4->setJetEMF(-99.0);
              factor = scaleReaderAK4->getCorrection();
            }
            TLorentzVector vCorr = factor * subjet.p4();
            sjSum += vCorr;
            double corr_pt = vCorr.Pt();

            // now vary
            uncReaderAK4->setJetEta(subjet.eta()); uncReaderAK4->setJetPt(corr_pt);
            double scaleUnc = uncReaderAK4->getUncertainty(true);
            sjSumUp += (1 + 2*scaleUnc) * vCorr;
            sjSumDown += (1 - 2*scaleUnc) * vCorr;

            // now smear...
            double smear=1, smearUp=1, smearDown=1;
            ak4JERReader->getStochasticSmear(corr_pt,subjet.eta(),event.rho,smear,smearUp,smearDown);
            sjSumSmear += smear * vCorr;
          }
          gt->fj1PtScaleUp_sj = gt->fj1Pt * (sjSumUp.Pt()/sjSum.Pt());
          gt->fj1PtScaleDown_sj = gt->fj1Pt * (sjSumDown.Pt()/sjSum.Pt());
          gt->fj1PtSmeared_sj = gt->fj1Pt * (sjSumSmear.Pt()/sjSum.Pt());
          gt->fj1MSDScaleUp_sj = gt->fj1MSD * (sjSumUp.Pt()/sjSum.Pt());
          gt->fj1MSDScaleDown_sj = gt->fj1MSD * (sjSumDown.Pt()/sjSum.Pt());
          gt->fj1MSDSmeared_sj = gt->fj1MSD * (sjSumSmear.Pt()/sjSum.Pt());


          // mSD correction
          float corrweight=1.;
          corrweight = GetMSDCorr(pt,eta);
          gt->fj1MSD_corr = corrweight*gt->fj1MSD;

          // now we do substructure
          gt->fj1Tau32 = clean(fj.tau3/fj.tau2);
          gt->fj1Tau32SD = clean(fj.tau3SD/fj.tau2SD);
          gt->fj1Tau21 = clean(fj.tau2/fj.tau1);
          gt->fj1Tau21SD = clean(fj.tau2SD/fj.tau1SD);

          for (auto ibeta : ibetas) {
            for (auto N : Ns) {
              for (auto order : orders) {
                GeneralTree::ECFParams p;
                p.order = order; p.N = N; p.ibeta = ibeta;
                if (gt->fj1IsClean || true)
                  gt->fj1ECFNs[p] = fj.get_ecf(order,N,ibeta);
                else
                  gt->fj1ECFNs[p] = fj.get_ecf(order,N,ibeta);
              }
            }
          } //loop over betas
          gt->fj1HTTMass = fj.htt_mass;
          gt->fj1HTTFRec = fj.htt_frec;

          std::vector<panda::MicroJet const*> subjets;
          for (unsigned iS(0); iS != fj.subjets.size(); ++iS)
           subjets.push_back(&fj.subjets.objAt(iS));

          auto csvsort([](panda::MicroJet const* j1, panda::MicroJet const* j2)->bool {
            return j1->csv > j2->csv;
           });

          std::sort(subjets.begin(),subjets.end(),csvsort);
          if (subjets.size()>0) {
            gt->fj1MaxCSV = subjets.at(0)->csv;
            gt->fj1MinCSV = subjets.back()->csv;
            if (subjets.size()>1) {
              gt->fj1SubMaxCSV = subjets.at(1)->csv;
            }
          }
          gt->fj1DoubleCSV = fj.double_sub;

          if (doMonoH) {
            for (unsigned int iSJ=0; iSJ!=fj.subjets.size(); ++iSJ) {
              auto& subjet = fj.subjets.objAt(iSJ);
              gt->fj1sjPt[iSJ]=subjet.pt();
              gt->fj1sjEta[iSJ]=subjet.eta();
              gt->fj1sjPhi[iSJ]=subjet.phi();
              gt->fj1sjM[iSJ]=subjet.m();
              gt->fj1sjCSV[iSJ]=subjet.csv;
              gt->fj1sjQGL[iSJ]=subjet.qgl;
            }
          }
        }
      }
      tr.TriggerSubEvent("fatjet basics");

      if (flags["pfCands"] && fj1) {
        VPseudoJet particles = ConvertPFCands(event.pfCandidates,flags["puppi"],0);
        fastjet::ClusterSequenceArea seq(particles,*jetDef,*areaDef);
        VPseudoJet allJets(seq.inclusive_jets(0.));
        fastjet::PseudoJet *pj1=0;
        double minDR2 = 999;
        for (auto &jet : allJets) {
          double dr2 = DeltaR2(jet.eta(),jet.phi_std(),fj1->eta(),fj1->phi());
          if (dr2<minDR2) {
            minDR2 = dr2;
            pj1 = &jet;
          }
        }
        if (pj1) {
          VPseudoJet constituents = fastjet::sorted_by_pt(pj1->constituents());

          gt->fj1NConst = constituents.size();
          double eTot=0, eTrunc=0;
          for (unsigned iC=0; iC!=gt->fj1NConst; ++iC) {
            double e = constituents.at(iC).E();
            eTot += e;
            if (iC<100)
              eTrunc += e;
          }
          gt->fj1EFrac100 = eTrunc/eTot;


          fastjet::PseudoJet sdJet = (*softDrop)(*pj1);
          VPseudoJet sdConstituents = fastjet::sorted_by_pt(sdJet.constituents());
          gt->fj1NSDConst = sdConstituents.size();
          eTot=0; eTrunc=0;
          for (unsigned iC=0; iC!=gt->fj1NSDConst; ++iC) {
            double e = sdConstituents.at(iC).E();
            eTot += e;
            if (iC<100)
              eTrunc += e;
          }
          gt->fj1SDEFrac100 = eTrunc/eTot;

        }
        tr.TriggerSubEvent("fatjet reclustering");
      }
    }

    tr.TriggerEvent("fatjet");

    // first identify interesting jets
    vector<panda::Jet*> cleanedJets, isoJets, btaggedJets, centralJets, forwardJets;
    vector<int> btagindices;
    TLorentzVector vJet;
    panda::Jet *jet1=0, *jet2=0, *jet3=0, *jet4=0, *jet5=0;
    panda::Jet *jot1=0, *jot2=0;
    panda::Jet *jotUp1=0, *jotUp2=0;
    panda::Jet *jotDown1=0, *jotDown2=0;
    gt->dphipuppimet=999; gt->dphipfmet=999;
    gt->dphipuppiUW=999; gt->dphipfUW=999;
    gt->dphipuppiUZ=999; gt->dphipfUZ=999;
    gt->dphipuppiUA=999; gt->dphipfUA=999;
    gt->dphipfUWW=999;
    float maxJetEta = (doVBF) ? 4.7 : 4.5;
    float maxIsoEta = (doMonoH) ? maxJetEta : 2.5;
    unsigned nJetDPhi = (doVBF) ? 4 : 4;

    gt->badECALFilter = 1;
    for (auto& jet : *jets) {
     
     // only do eta-phi checks here
     if (abs(jet.eta()) > maxJetEta)
        continue;
     // NOTE:
     // For VBF we require nTightLep>0, but in monotop looseLep1IsTight
     // No good reason to do that, should switch to former
     // Should update jet cleaning accordingly (just check all loose objects)
     if (IsMatched(&matchLeps,0.16,jet.eta(),jet.phi()) ||
         IsMatched(&matchPhos,0.16,jet.eta(),jet.phi()))
        continue;
     if (doVBF && !jet.loose)
       continue;

     if (doVBF && jet.pt()>20 && fabs(jet.eta())<2.4 && jet.csv>0.8484) {
        ++(gt->jetNMBtags);
     }


     if (jet.pt()>30) { // nominal jets
      cleanedJets.push_back(&jet);
      if (cleanedJets.size()==1) {
        jot1 = &jet;
        gt->jot1Pt = jet.pt();
        gt->jot1Eta = jet.eta();
        gt->jot1Phi = jet.phi();
        if (doVBF && fabs(gt->jot1Eta)<2.4) { // if it's a central jet, must jot1 ID requirements
          gt->jot1VBFID = jet.monojet;
        } else { // if leading jet is not central, leave the event be
          gt->jot1VBFID = 1;
        }
      } else if (cleanedJets.size()==2) {
        jot2 = &jet;
        gt->jot2Pt = jet.pt();
        gt->jot2Eta = jet.eta();
        gt->jot2Phi = jet.phi();
      }

      if (cleanedJets.size()<3) {
        bool isBad = GetCorr(cBadECALJets,jet.eta(),jet.phi()) > 0;
        if (isBad)
          gt->badECALFilter = 0;
      }

      float csv = (fabs(jet.eta())<2.5) ? jet.csv : -1;
      if (fabs(jet.eta())<2.4) {
        centralJets.push_back(&jet);
        if (centralJets.size()==1) {
          jet1 = &jet;
          gt->jet1Pt = jet.pt();
          gt->jet1Eta = jet.eta();
          gt->jet1Phi = jet.phi();
          gt->jet1CSV = csv;
          gt->jet1IsTight = jet.monojet ? 1 : 0;
        } else if (centralJets.size()==2) {
          jet2 = &jet;
          gt->jet2Pt = jet.pt();
          gt->jet2Eta = jet.eta();
          gt->jet2Phi = jet.phi();
          gt->jet2CSV = csv;
        } else if (centralJets.size()==3) {
	  jet3 = &jet;
          gt->jet3Pt = jet.pt();
          gt->jet3Eta = jet.eta();
          gt->jet3Phi = jet.phi();
          gt->jet3CSV = csv;
	} else if (centralJets.size()==4) {
	  jet4 = &jet;
          gt->jet4Pt = jet.pt();
          gt->jet4Eta = jet.eta();
          gt->jet4Phi = jet.phi();
          gt->jet4CSV = csv;
	} else if (centralJets.size()==5) {
	  jet5 = &jet;
          gt->jet5Pt = jet.pt();
          gt->jet5Eta = jet.eta();
          gt->jet5Phi = jet.phi();
          gt->jet5CSV = csv;
	}
	
      }
      else{
	forwardJets.push_back(&jet);
      }

      if (doMonoH) {
        gt->jetPt[cleanedJets.size()-1]=jet.pt();
        gt->jetEta[cleanedJets.size()-1]=jet.eta();
        gt->jetPhi[cleanedJets.size()-1]=jet.phi();
        gt->jetE[cleanedJets.size()-1]=jet.m();
        gt->jetCSV[cleanedJets.size()-1]=csv;
        gt->jetQGL[cleanedJets.size()-1]=jet.qgl;
      }

      // compute dphi wrt mets
      if (cleanedJets.size() <= nJetDPhi) {
        vJet.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());
        gt->dphipuppimet = std::min(fabs(vJet.DeltaPhi(vPuppiMET)),(double)gt->dphipuppimet);
        gt->dphipfmet = std::min(fabs(vJet.DeltaPhi(vPFMET)),(double)gt->dphipfmet);
        gt->dphipuppiUA = std::min(fabs(vJet.DeltaPhi(vpuppiUA)),(double)gt->dphipuppiUA);
        gt->dphipuppiUW = std::min(fabs(vJet.DeltaPhi(vpuppiUW)),(double)gt->dphipuppiUW);
        gt->dphipuppiUZ = std::min(fabs(vJet.DeltaPhi(vpuppiUZ)),(double)gt->dphipuppiUZ);
        gt->dphipfUA = std::min(fabs(vJet.DeltaPhi(vpfUA)),(double)gt->dphipfUA);
        gt->dphipfUW = std::min(fabs(vJet.DeltaPhi(vpfUW)),(double)gt->dphipfUW);
        gt->dphipfUZ = std::min(fabs(vJet.DeltaPhi(vpfUZ)),(double)gt->dphipfUZ);
	gt->dphipfUWW = std::min(fabs(vJet.DeltaPhi(vpfUWW)),(double)gt->dphipfUWW);
      }
      // btags
      if (csv>0.5426) {
        ++(gt->jetNBtags);
        if (doMonoH) {
          btaggedJets.push_back(&jet);
          btagindices.push_back(cleanedJets.size()-1);
        }
        if (!doVBF && csv>0.8484) 
          ++(gt->jetNMBtags);
      }

      bool isIsoJet = ( (gt->nFatjet==0) || 
               (fabs(jet.eta())<maxIsoEta 
               && DeltaR2(gt->fj1Eta,gt->fj1Phi,jet.eta(),jet.phi())>2.25) ); 

      if (isIsoJet) {
        isoJets.push_back(&jet);
        if (csv>0.5426)
          ++gt->isojetNBtags;
        if (isoJets.size()==1) {
         gt->isojet1Pt = jet.pt();
         gt->isojet1CSV = jet.csv;
        } else if (isoJets.size()==2) {
         gt->isojet2Pt = jet.pt();
         gt->isojet2CSV = jet.csv;
        }
        if (doMonoH)
          gt->jetIso[cleanedJets.size()-1]=1;
      } else {
        if (doMonoH)
          gt->jetIso[cleanedJets.size()-1]=0;
      }
     }

     // do jes variation OUTSIDE of pt>30 check
     if (jet.ptCorrUp>30) {
      if (jet.ptCorrUp > gt->jot1PtUp) {
        if (jotUp1) {
          jotUp2 = jotUp1;
          gt->jot2PtUp = gt->jot1PtUp;
          gt->jot2EtaUp = gt->jot1EtaUp;
        }
        jotUp1 = &jet;
        gt->jot1PtUp = jet.ptCorrUp;
        gt->jot1EtaUp = jet.eta();
      } else if (jet.ptCorrUp > gt->jot2PtUp) {
        jotUp2 = &jet;
        gt->jot2PtUp = jet.ptCorrUp;
        gt->jot2EtaUp = jet.eta();
      }
     }
     if (jet.ptCorrDown>30) {
      if (jet.ptCorrDown > gt->jot1PtDown) {
        if (jotDown1) {
          jotDown2 = jotDown1;
          gt->jot2PtDown = gt->jot1PtDown;
          gt->jot2EtaDown = gt->jot1EtaDown;
        }
        jotDown1 = &jet;
        gt->jot1PtDown = jet.ptCorrDown;
        gt->jot1EtaDown = jet.eta();
      } else if (jet.ptCorrDown > gt->jot2PtDown) {
        jotDown2 = &jet;
        gt->jot2PtDown = jet.ptCorrDown;
        gt->jot2EtaDown = jet.eta();
      }
     }

    } // VJet loop

    switch (whichRecoil) {
      case -1: // photon
        gt->dphipuppiU = gt->dphipuppiUA;
        gt->dphipfU = gt->dphipfUA;
        break;
      case 0: // MET
        gt->dphipuppiU = gt->dphipuppimet;
        gt->dphipfU = gt->dphipfmet;
        break;
      case 1:
        gt->dphipuppiU = gt->dphipuppiUW;
        gt->dphipfU = gt->dphipfUW;
        break;
      case 2:
        gt->dphipuppiU = gt->dphipuppiUZ;
        gt->dphipfU = gt->dphipfUZ;
        break;
      case 3:
	//gt->dphipuppiU = gt->dphipuppiUZ;
        gt->dphipfU = gt->dphipfUWW;
	break;
      default: // impossible
        break;
    }

    gt->nJet = centralJets.size();
    gt->nJot = cleanedJets.size();
    if (gt->nJot>1 && doVBF) {
     TLorentzVector vj1, vj2;
     vj1.SetPtEtaPhiM(jot1->pt(),jot1->eta(),jot1->phi(),jot1->m());
     vj2.SetPtEtaPhiM(jot2->pt(),jot2->eta(),jot2->phi(),jot2->m());
     gt->jot12Mass = (vj1+vj2).M();
     gt->jot12DPhi = vj1.DeltaPhi(vj2);
     gt->jot12DEta = fabs(jot1->eta()-jot2->eta());

     if (jotUp1 && jotUp2) {
       vj1.SetPtEtaPhiM(jotUp1->ptCorrUp,jotUp1->eta(),jotUp1->phi(),jotUp1->m());
       vj2.SetPtEtaPhiM(jotUp2->ptCorrUp,jotUp2->eta(),jotUp2->phi(),jotUp2->m());
       gt->jot12MassUp = (vj1+vj2).M();
       gt->jot12DPhiUp = vj1.DeltaPhi(vj2);
       gt->jot12DEtaUp = fabs(jotUp1->eta()-jotUp2->eta());
     }
     
     if (jotDown1 && jotDown2) {
       vj1.SetPtEtaPhiM(jotDown1->ptCorrDown,jotDown1->eta(),jotDown1->phi(),jotDown1->m());
       vj2.SetPtEtaPhiM(jotDown2->ptCorrDown,jotDown2->eta(),jotDown2->phi(),jotDown2->m());
       gt->jot12MassDown = (vj1+vj2).M();
       gt->jot12DPhiDown = vj1.DeltaPhi(vj2);
       gt->jot12DEtaDown = fabs(jotDown1->eta()-jotDown2->eta());
     }
    }

    tr.TriggerEvent("jets");


    if (doMonoH) {
      // Higgs reconstruction for resolved analysis - highest pt pair of b jets
      float tmp_hbbpt=-99;
      float tmp_hbbeta=-99;
      float tmp_hbbphi=-99;
      float tmp_hbbm=-99;
      int tmp_hbbjtidx1=-1;
      int tmp_hbbjtidx2=-1;
      for (unsigned int i = 0;i<btaggedJets.size();i++){
        panda::Jet *jet_1 = btaggedJets.at(i);
        TLorentzVector hbbdaughter1;
        hbbdaughter1.SetPtEtaPhiM(jet_1->pt(),jet_1->eta(),jet_1->phi(),jet_1->m());
        for (unsigned int j = i+1;j<btaggedJets.size();j++){
          panda::Jet *jet_2 = btaggedJets.at(j);
          TLorentzVector hbbdaughter2;
          hbbdaughter2.SetPtEtaPhiM(jet_2->pt(),jet_2->eta(),jet_2->phi(),jet_2->m());
          TLorentzVector hbbsystem = hbbdaughter1 + hbbdaughter2;
          if (hbbsystem.Pt()>tmp_hbbpt){
            tmp_hbbpt = hbbsystem.Pt();
            tmp_hbbeta = hbbsystem.Eta();
            tmp_hbbphi = hbbsystem.Phi();
            tmp_hbbm = hbbsystem.M();
            tmp_hbbjtidx1 = btagindices.at(i);
            tmp_hbbjtidx2 = btagindices.at(j);
          }
        }
      }
      gt->hbbpt = tmp_hbbpt;
      gt->hbbeta = tmp_hbbeta;
      gt->hbbphi = tmp_hbbphi;
      gt->hbbm = tmp_hbbm;
      gt->hbbjtidx[0] = tmp_hbbjtidx1;
      gt->hbbjtidx[1] = tmp_hbbjtidx2;

      tr.TriggerEvent("monohiggs");
    }

    for (auto& tau : event.taus) {
      if (doVBF) {
        if (!tau.decayMode || !tau.decayModeNew)
          continue;
        if (!tau.looseIsoMVAOld)
          continue;
      } else {
        if (!tau.decayMode || !tau.decayModeNew)
          continue;
        if (!tau.looseIsoMVA)
          continue;
      }
      /*
      if (tau.isoDeltaBetaCorr>5)
        continue;
      */
      if (tau.pt()<18 || fabs(tau.eta())>2.3)
        continue;
      if (IsMatched(&matchLeps,0.16,tau.eta(),tau.phi()))
        continue;
      gt->nTau++;
    }

    tr.TriggerEvent("taus");

    if (!PassPreselection())
      continue;

    std::cout<<"Pass Preselection"<<std::endl;

    tr.TriggerEvent("presel");

    // identify interesting gen particles for fatjet matching
    unsigned int pdgidTarget=0;
    if (!isData && processType>=kTT) {
      switch(processType) {
        case kTop:
        case kTT:
        case kSignal:
          pdgidTarget=6;
          break;
        case kV:
          pdgidTarget=24;
          break;
        case kH:
          pdgidTarget=25;
          break;
        default:
          // processType>=kTT means we should never get here
          PError("PandaAnalyzer::Run","Reached an unknown process type");
      }

      std::vector<int> targets;

      int nGen = event.genParticles.size();
      for (int iG=0; iG!=nGen; ++iG) {
        auto& part(event.genParticles.at(iG));
        int pdgid = part.pdgid;
        unsigned int abspdgid = abs(pdgid);
        if (abspdgid == pdgidTarget)
          targets.push_back(iG);
      } //looking for targets

      for (int iG : targets) {
        auto& part(event.genParticles.at(iG));

        // check there is no further copy:
        bool isLastCopy=true;
        for (int jG : targets) {
          if (event.genParticles.at(jG).parent.get() == &part) {
            isLastCopy=false;
            break;
          }
        }
        if (!isLastCopy)
          continue;

        // (a) check it is a hadronic decay and if so, (b) calculate the size
        if (processType==kTop||processType==kTT) {

          // first look for a W whose parent is the top at iG, or a W further down the chain
          panda::GenParticle const* lastW(0);
          for (int jG=0; jG!=nGen; ++jG) {
            GenParticle const& partW(event.genParticles.at(jG));
            if (TMath::Abs(partW.pdgid)==24 && partW.pdgid*part.pdgid>0) {
              // it's a W and has the same sign as the top
              if (!lastW && partW.parent.get() == &part) {
                lastW = &partW;
              } else if (lastW && partW.parent.get() == lastW) {
                lastW = &partW;
              }
            }
          } // looking for W
          if (!lastW) {// ???
            continue;
          }
          auto& partW(*lastW);

          // now look for b or W->qq
          int iB=-1, iQ1=-1, iQ2=-1;
          double size=0, sizeW=0;
          for (int jG=0; jG!=nGen; ++jG) {
            auto& partQ(event.genParticles.at(jG));
            int pdgidQ = partQ.pdgid;
            unsigned int abspdgidQ = TMath::Abs(pdgidQ);
            if (abspdgidQ>5)
              continue;
            if (abspdgidQ==5 && iB<0 && partQ.parent.get() == &part) {
              // only keep first copy
              iB = jG;
              size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),size);
            } else if (abspdgidQ<5 && partQ.parent.get() == &partW) {
              if (iQ1<0) {
                iQ1 = jG;
                size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                         size);
                sizeW = TMath::Max(DeltaR2(partW.eta(),partW.phi(),partQ.eta(),partQ.phi()),
                         sizeW);
              } else if (iQ2<0) {
                iQ2 = jG;
                size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                         size);
                sizeW = TMath::Max(DeltaR2(partW.eta(),partW.phi(),partQ.eta(),partQ.phi()),
                         sizeW);
              }
            }
            if (iB>=0 && iQ1>=0 && iQ2>=0)
              break;
          } // looking for quarks


          bool isHadronic = (iB>=0 && iQ1>=0 && iQ2>=0); // all 3 quarks were found
          if (isHadronic)
            genObjects[&part] = size;

          bool isHadronicW = (iQ1>=0 && iQ2>=0);
          if (isHadronicW)
            genObjects[&partW] = sizeW;

        } else { // these are W,Z,H - 2 prong decays

          int iQ1=-1, iQ2=-1;
          double size=0;
          for (int jG=0; jG!=nGen; ++jG) {
            auto& partQ(event.genParticles.at(jG));
            int pdgidQ = partQ.pdgid;
            unsigned int abspdgidQ = TMath::Abs(pdgidQ);
            if (abspdgidQ>5)
              continue;
            if (partQ.parent.get() == &part) {
              if (iQ1<0) {
                iQ1=jG;
                size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                         size);
              } else if (iQ2<0) {
                iQ2=jG;
                size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                         size);
              }
            }
            if (iQ1>=0 && iQ2>=0)
              break;
          } // looking for quarks

          bool isHadronic = (iQ1>=0 && iQ2>=0); // both quarks were found

          // add to collection
          if (isHadronic)
            genObjects[&part] = size;
        }

      } // loop over targets
    } // process is interesting

    tr.TriggerEvent("gen matching");

    if (!isData && gt->nFatjet>0) {
      // first see if jet is matched
      auto* matched = MatchToGen(fj1->eta(),fj1->phi(),1.5,pdgidTarget);
      if (matched!=NULL) {
        gt->fj1IsMatched = 1;
        gt->fj1GenPt = matched->pt();
        gt->fj1GenSize = genObjects[matched];
      } else {
        gt->fj1IsMatched = 0;
      }
      if (pdgidTarget==6) { // matched to top; try for W
        auto* matchedW = MatchToGen(fj1->eta(),fj1->phi(),1.5,24);
        if (matchedW!=NULL) {
          gt->fj1IsWMatched = 1;
          gt->fj1GenWPt = matchedW->pt();
          gt->fj1GenWSize = genObjects[matchedW];
        } else {
          gt->fj1IsWMatched = 0;
        }
      }

      bool found_b_from_g=false;
      int bs_inside_cone=0;
      int has_gluon_splitting=0;
      panda::GenParticle const* first_b_mo(0);
      // now get the highest pT gen particle inside the jet cone
      for (auto& gen : event.genParticles) {
        float pt = gen.pt();
        int pdgid = gen.pdgid;
        if (pt>(gt->fj1HighestPtGenPt)
          && DeltaR2(gen.eta(),gen.phi(),fj1->eta(),fj1->phi())<2.25) {
          gt->fj1HighestPtGenPt = pt;
          gt->fj1HighestPtGen = pdgid;
        }

        if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
          continue;

        //count bs and cs
        int apdgid = abs(pdgid);
        if (apdgid!=5 && apdgid!=4) 
          continue;

        if (gen.pt()>5) {
          gt->nHF++;
          if (apdgid==5)
            gt->nB++;
        }

        if (DeltaR2(gen.eta(),gen.phi(),fj1->eta(),fj1->phi())<2.25) {
          gt->fj1NHF++;
          if (apdgid==5) {
            if (gen.parent.isValid() && gen.parent->pdgid==21 && gen.parent->pt()>20) {
              if (!found_b_from_g) {
                found_b_from_g=true;
                first_b_mo=gen.parent.get();
                bs_inside_cone+=1;
              } else if (gen.parent.get()==first_b_mo) {
                bs_inside_cone+=1;
                has_gluon_splitting=1;
              } else {
                bs_inside_cone+=1;
              }
            } else {
              bs_inside_cone+=1;
            }
          }
        }
      }

      gt->fj1Nbs=bs_inside_cone;
      gt->fj1gbb=has_gluon_splitting;
    
      // now get the subjet btag SFs
      vector<btagcand> sj_btagcands;
      vector<double> sj_sf_cent, sj_sf_bUp, sj_sf_bDown, sj_sf_mUp, sj_sf_mDown;
      unsigned int nSJ = fj1->subjets.size();
      for (unsigned int iSJ=0; iSJ!=nSJ; ++iSJ) {
        auto& subjet = fj1->subjets.objAt(iSJ);
        int flavor=0;
        for (auto& gen : event.genParticles) {
          int apdgid = abs(gen.pdgid);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(subjet.eta(),subjet.phi(),gen.eta(),gen.phi());
          if (dr2<0.09) {
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid;
              break;
            } else {
              flavor=0;
            }
          }
        } // finding the subjet flavor

        float pt = subjet.pt();
        float btagUncFactor = 1;
        float eta = subjet.eta();
        double eff(1),sf(1),sfUp(1),sfDown(1);
        unsigned int binpt = btagpt.bin(pt);
        unsigned int bineta = btageta.bin(fabs(eta));
        if (flavor==5) {
          eff = beff[bineta][binpt];
        } else if (flavor==4) {
          eff = ceff[bineta][binpt];
        } else {
          eff = lfeff[bineta][binpt];
        }
        CalcBJetSFs(bSubJetL,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
        sj_btagcands.push_back(btagcand(iSJ,flavor,eff,sf,sfUp,sfDown));
        sj_sf_cent.push_back(sf);
        if (flavor>0) {
          sj_sf_bUp.push_back(sfUp); sj_sf_bDown.push_back(sfDown);
          sj_sf_mUp.push_back(sf); sj_sf_mDown.push_back(sf);
        } else {
          sj_sf_bUp.push_back(sf); sj_sf_bDown.push_back(sf);
          sj_sf_mUp.push_back(sfUp); sj_sf_mDown.push_back(sfDown);
        }

      } // loop over subjets

      EvalBTagSF(sj_btagcands,sj_sf_cent,GeneralTree::bCent,GeneralTree::bSubJet);
      EvalBTagSF(sj_btagcands,sj_sf_bUp,GeneralTree::bBUp,GeneralTree::bSubJet);
      EvalBTagSF(sj_btagcands,sj_sf_bDown,GeneralTree::bBDown,GeneralTree::bSubJet);
      EvalBTagSF(sj_btagcands,sj_sf_mUp,GeneralTree::bMUp,GeneralTree::bSubJet);
      EvalBTagSF(sj_btagcands,sj_sf_mDown,GeneralTree::bMDown,GeneralTree::bSubJet);

    }

    tr.TriggerEvent("fatjet gen-matching");
    
    if (!isData) {
      // now get the jet btag SFs
      vector<btagcand> btagcands;
      vector<btagcand> btagcands_alt;
      vector<double> sf_cent, sf_bUp, sf_bDown, sf_mUp, sf_mDown;
      vector<double> sf_cent_alt, sf_bUp_alt, sf_bDown_alt, sf_mUp_alt, sf_mDown_alt;

      unsigned int nJ = centralJets.size();
      for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
        panda::Jet *jet = centralJets.at(iJ);
        bool isIsoJet=false;
        if (std::find(isoJets.begin(), isoJets.end(), jet) != isoJets.end())
          isIsoJet = true;
        int flavor=0;
        float genpt=0;
	float PX=0;
	float PY=0;
	float PZ=0;
	float EN=0;
	int ID=0;
	int status=0;
	float MASS=0;
        for (auto& gen : event.genParticles) {
          int apdgid = abs(gen.pdgid);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(jet->eta(),jet->phi(),gen.eta(),gen.phi());
          if (dr2<0.09) {
            genpt = gen.pt();
	    PX = gen.px();
	    PY = gen.py();
	    PZ = gen.pz();
	    EN = gen.e();
	    //status = gen.status();
	    ID = gen.pdgid;
	    MASS = gen.m();
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid; //flavor variable strictly for btag usage, to get pdgid of genparticle, refer to XXGenFlav
              break;
            } else {
              flavor=0;
            }
          }
        } // finding the jet flavor
        float pt = jet->pt();
        float btagUncFactor = 1;
        float eta = jet->eta();
        double eff(1),sf(1),sfUp(1),sfDown(1);
        unsigned int binpt = btagpt.bin(pt);
        unsigned int bineta = btageta.bin(fabs(eta));
        if (flavor==5)
          eff = beff[bineta][binpt];
        else if (flavor==4)
          eff = ceff[bineta][binpt];
        else
          eff = lfeff[bineta][binpt];
        if (jet==centralJets.at(0)) {
          gt->jet1Flav = flavor;
          gt->jet1GenPt = genpt;
	  //Gen level of jets
	  gt->jet1GenPx = PX;
	  gt->jet1GenPy = PY;
	  gt->jet1GenPz = PZ;
	  gt->jet1GenEn = EN;
	  //gt->genjet1status = status;
	  gt->jet1GenFlav = ID;
	  gt->jet1GenMass = MASS;
        } else if (jet==centralJets.at(1)) {
          gt->jet2Flav = flavor;
          gt->jet2GenPt = genpt;
	  gt->jet2GenPx = PX;
          gt->jet2GenPy = PY;
          gt->jet2GenPz = PZ;
          gt->jet2GenEn = EN;
          //gt->genjet2status = status;
          gt->jet2GenFlav = ID;
	  gt->jet2GenMass = MASS;

        }else if (jet==centralJets.at(2)) {
	  gt->jet3Flav = flavor;
          gt->jet3GenPt = genpt;
	  gt->jet3GenPx = PX;
          gt->jet3GenPy = PY;
          gt->jet3GenPz = PZ;
          gt->jet3GenEn = EN;
          //gt->genjet3status = status;
          gt->jet3GenFlav = ID;
	  gt->jet3GenMass = MASS;
	}else if (jet==centralJets.at(3)) {
	  gt->jet4Flav = flavor;
          gt->jet4GenPt = genpt;
	  gt->jet4GenPx = PX;
          gt->jet4GenPy = PY;
          gt->jet4GenPz = PZ;
          gt->jet4GenEn = EN;
	  gt->jet4GenFlav = ID;
	  gt->jet4GenMass = MASS;
	}else if (jet==centralJets.at(4)) {
	  gt->jet5Flav = flavor;
          gt->jet5GenPt = genpt;
	  gt->jet5GenPx = PX;
          gt->jet5GenPy = PY;
          gt->jet5GenPz = PZ;
          gt->jet5GenEn = EN;
	  gt->jet5GenFlav = ID;
	  gt->jet5GenMass = MASS;
	}


        if (isIsoJet) {
          if (jet==isoJets.at(0))
            gt->isojet1Flav = flavor;
          else if (jet==isoJets.at(1))
            gt->isojet2Flav = flavor;

          CalcBJetSFs(bJetL,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
          btagcands.push_back(btagcand(iJ,flavor,eff,sf,sfUp,sfDown));
          sf_cent.push_back(sf);

          if (flavor>0) {
            sf_bUp.push_back(sfUp); sf_bDown.push_back(sfDown);
            sf_mUp.push_back(sf); sf_mDown.push_back(sf);
          } else {
            sf_bUp.push_back(sf); sf_bDown.push_back(sf);
            sf_mUp.push_back(sfUp); sf_mDown.push_back(sfDown);
          }

        }

        // I'm killing this temporarily as we don't use it -SN
	// I AM ENABLING THIS PART
        if (doMonoH){
          // alternate stuff for inclusive jet collection (also different b tagging WP)
	  //std::cout<<"B-tagging Medium SF evaluating"<<std::endl;
          double sf_alt(1),sfUp_alt(1),sfDown_alt(1);

	  //Calculating medium working point SF
          CalcBJetSFs(bJetM,flavor,eta,pt,eff,btagUncFactor,sf_alt,sfUp_alt,sfDown_alt); //
          btagcands_alt.push_back(btagcand(iJ,flavor,eff,sf_alt,sfUp_alt,sfDown_alt));
          sf_cent_alt.push_back(sf_alt);
	  ////
          if (flavor>0) {
            sf_bUp_alt.push_back(sfUp_alt); sf_bDown_alt.push_back(sfDown_alt);
            sf_mUp_alt.push_back(sf_alt); sf_mDown_alt.push_back(sf_alt);
          } else {
            sf_bUp_alt.push_back(sf_alt); sf_bDown_alt.push_back(sf_alt);
            sf_mUp_alt.push_back(sfUp_alt); sf_mDown_alt.push_back(sfDown_alt);
	    }
	    
	}
        
      } // loop over jets
      /*
      // Loop over forward jet, eta> 2.5
      unsigned int nFJ = forwardJets.size();
      if ( nFJ > 2 ){
	for (unsigned int iFJ=0; iFJ!=nFJ; ++iFJ) {
	  panda::Jet *iforwardjet = centralJets.at(iFJ);

	  TLorentzVector forwardjet1(0.,0.,0.,0.);
	  forwardjet1.SetPxPyPzE(iforwardjet->px(),iforwardjet->py(),iforwardjet->pz(),iforwardjet->e());
	  for (unsigned int jFJ=1; jFJ!=nFJ; ++jFJ) {
	    if (iFJ==jFJ) continue;
	    panda::Jet *jforwardjet = centralJets.at(jFJ);
	    TLorentzVector forwardjet2(0.,0.,0.,0.);
	    forwardjet2.SetPxPyPzE(jforwardjet->px(),jforwardjet->py(),jforwardjet->pz(),jforwardjet->e());
	  
	    double dijetmass = (forwardjet1+forwardjet2).M();
	    double deltaeta = iforwardjet->eta() - jforwardjet->eta();
	    if (deltaeta > 4.2 && dijetmass>800 ){
	      gt->forwjet1Pt = iforwardjet->pt();
	      gt->forwjet1Phi = iforwardjet->phi();
	      gt->forwjet1Eta = iforwardjet->eta();
	      gt->forwjet2Pt = jforwardjet->pt();
	      gt->forwjet2Phi = jforwardjet->phi();
	      gt->forwjet2Eta = jforwardjet->eta();
	      break;
	    }
	  }
	} // loop over forward jet
      }
      */
      //ENABLE 2-btag SF 
      // isojet with Low working points
      EvalBTagSF(btagcands,sf_cent,GeneralTree::bCent,GeneralTree::bJet,false);
      EvalBTagSF(btagcands,sf_bUp,GeneralTree::bBUp,GeneralTree::bJet,false);
      EvalBTagSF(btagcands,sf_bDown,GeneralTree::bBDown,GeneralTree::bJet,false);
      EvalBTagSF(btagcands,sf_mUp,GeneralTree::bMUp,GeneralTree::bJet,false);
      EvalBTagSF(btagcands,sf_mDown,GeneralTree::bMDown,GeneralTree::bJet,false);

      // see above, also needs to use the new functions -SN
      
      if (flags["monohiggs"]) {

	EvalBTagSF(btagcands_alt,sf_cent_alt,GeneralTree::bCent,GeneralTree::bMedJet,true);
	EvalBTagSF(btagcands_alt,sf_bUp_alt, GeneralTree::bBUp, GeneralTree::bMedJet,true);
	EvalBTagSF(btagcands_alt,sf_bDown_alt, GeneralTree::bBDown, GeneralTree::bMedJet,true);
	EvalBTagSF(btagcands_alt,sf_mUp_alt, GeneralTree::bMUp, GeneralTree::bMedJet,true);
	EvalBTagSF(btagcands_alt,sf_mDown_alt, GeneralTree::bMDown, GeneralTree::bMedJet,true);
	
        //EvalBTagSF(btagcands_alt,sf_cent_alt, gt->sf_btag0_alt,gt->sf_btag1_alt,gt->sf_btag2_alt,gt->sf_btagGT0_alt);
	//EvalBTagSF(btagcands_alt,sf_bUp_alt, gt->sf_btag0BUp_alt,gt->sf_btag1BUp_alt,gt->sf_btag2BUp_alt,gt->sf_btagGT0BUp_alt);
	//EvalBTagSF(btagcands_alt,sf_bDown_alt, gt->sf_btag0BDown_alt,gt->sf_btag1BDown_alt,gt->sf_btag2BDown_alt,gt->sf_btagGT0BDown_alt);
	//EvalBTagSF(btagcands_alt,sf_mUp_alt, gt->sf_btag0MUp_alt,gt->sf_btag1MUp_alt,gt->sf_btag2MUp_alt,gt->sf_btagGT0MUp_alt);
	//EvalBTagSF(btagcands_alt,sf_mDown_alt, gt->sf_btag0MDown_alt,gt->sf_btag1MDown_alt,gt->sf_btag2MDown_alt,gt->sf_btagGT0MDown_alt);
	
      }
      
    }

    tr.TriggerEvent("ak4 gen-matching");

    // ttbar pT weight
    gt->sf_tt = 1; gt->sf_tt_ext = 1; gt->sf_tt_bound = 1;
    gt->sf_tt8TeV = 1; gt->sf_tt8TeV_ext = 1; gt->sf_tt8TeV_bound = 1;
    gt->sf_qcdTT = 1;
    if (!isData && processType==kTT) {
      gt->genWPlusPt = -1; gt->genWMinusPt = -1;
      for (auto& gen : event.genParticles) {
        if (abs(gen.pdgid)!=24)
          continue;
        if (flags["firstGen"]) {
          if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
            continue; // must be first copy
        }
        if (gen.pdgid>0) {
         gt->genWPlusPt = gen.pt();
         gt->genWPlusEta = gen.eta();
        } else {
         gt->genWMinusPt = gen.pt();
         gt->genWMinusEta = gen.eta();
        }
        if (flags["firstGen"]) {
          if (gt->genWPlusPt>0 && gt->genWMinusPt>0)
            break;
        }
      }
      TLorentzVector vT,vTbar;
      float pt_t=0, pt_tbar=0;
      for (auto& gen : event.genParticles) {
        if (abs(gen.pdgid)!=6)
          continue;
        if (flags["firstGen"]) {
          if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
            continue; // must be first copy
        }
        if (gen.pdgid>0) {
         pt_t = gen.pt();
         gt->genTopPt = gen.pt();
         gt->genTopEta = gen.eta();
         vT.SetPtEtaPhiM(gen.pt(),gen.eta(),gen.phi(),gen.m());
        } else {
         pt_tbar = gen.pt();
         gt->genAntiTopPt = gen.pt();
         gt->genAntiTopEta = gen.eta();
         vTbar.SetPtEtaPhiM(gen.pt(),gen.eta(),gen.phi(),gen.m());
        }
        if (flags["firstGen"]) {
          if (pt_t>0 && pt_tbar>0)
            break;
        }
      }
      if (pt_t>0 && pt_tbar>0) {
        TLorentzVector vTT = vT+vTbar;
        gt->genTTPt = vTT.Pt(); gt->genTTEta = vTT.Eta();
        gt->sf_tt8TeV       = TMath::Sqrt(TMath::Exp(0.156-0.00137*TMath::Min((float)400.,pt_t)) *
                         TMath::Exp(0.156-0.00137*TMath::Min((float)400.,pt_tbar)));
        gt->sf_tt           = TMath::Sqrt(TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_t)) *
                         TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_tbar)));
        gt->sf_tt8TeV_ext   = TMath::Sqrt(TMath::Exp(0.156-0.00137*pt_t) *
                         TMath::Exp(0.156-0.00137*pt_tbar));
        gt->sf_tt_ext       = TMath::Sqrt(TMath::Exp(0.0615-0.0005*pt_t) *
                         TMath::Exp(0.0615-0.0005*pt_tbar));
        gt->sf_tt8TeV_bound = TMath::Sqrt(((pt_t>400) ? 1 : TMath::Exp(0.156-0.00137*pt_t)) *
                         ((pt_tbar>400) ? 1 : TMath::Exp(0.156-0.00137*pt_tbar)));
        gt->sf_tt_bound     = TMath::Sqrt(((pt_t>400) ? 1 : TMath::Exp(0.0615-0.0005*pt_t)) *
                         ((pt_tbar>400) ? 1 : TMath::Exp(0.0615-0.0005*pt_tbar)));
      }

      if (pt_t>0)
        gt->sf_qcdTT *= TTNLOToNNLO(pt_t);
      if (pt_tbar>0) 
        gt->sf_qcdTT *= TTNLOToNNLO(pt_tbar);
      gt->sf_qcdTT = TMath::Sqrt(gt->sf_qcdTT);

    }

    tr.TriggerEvent("tt SFs");

    // derive ewk/qcd weights
    gt->sf_qcdV=1; gt->sf_ewkV=1;
    gt->sf_qcdV_VBF=1;
    if (!isData) {
      bool found = processType!=kA && processType!=kZ && processType!=kW
                     && processType!=kZEWK && processType!=kWEWK;
      //gen info of boson stored according to pdgid
      int target=24;
      if (processType==kZ || processType==kZEWK) target=23;
      if (processType==kA) target=22;

      for (auto& gen : event.genParticles) {
        if (found) break;
        int apdgid = abs(gen.pdgid);
	float PX = gen.px();
	float PY = gen.py();
	float PZ = gen.pz();
	float EN = gen.e();
	int status = 0;
	int ID = gen.pdgid;
	float MASS = gen.m();
        if (apdgid==target)     {
          bool foundChild = false;
          for (auto& child : event.genParticles) {

            if (abs(child.pdgid) != target)
              continue;
            if (child.parent.isValid() && child.parent.get() == &(gen)) {
              foundChild = true; 
              break;
            }
          }
          if (foundChild)
            continue;
          if (processType==kZ) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
	    gt->genBosonGenPx= PX;
	    gt->genBosonGenPy= PY;
	    gt->genBosonGenPz= PZ;
	    gt->genBosonGenEn= EN;
	    gt->genBosonGenFlav = ID;
	    gt->genBosonGenMass = MASS;
            gt->sf_qcdV = GetCorr(cZNLO,gt->genBosonPt);
            gt->sf_ewkV = GetCorr(cZEWK,gt->genBosonPt);
            gt->sf_qcdV_VBF = GetCorr(cVBF_ZNLO,gt->genBosonPt) * gt->sf_qcdV;
            gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_ZNLO,gt->genBosonPt) * gt->sf_qcdV;
            found=true;
          } else if (processType==kW) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
	    gt->genBosonGenPx= PX;
            gt->genBosonGenPy= PY;
            gt->genBosonGenPz= PZ;
            gt->genBosonGenEn= EN;
            gt->genBosonGenFlav = ID;
	    gt->genBosonGenMass = MASS;
            gt->sf_qcdV = GetCorr(cWNLO,gt->genBosonPt);
            gt->sf_ewkV = GetCorr(cWEWK,gt->genBosonPt);
            gt->sf_qcdV_VBF = GetCorr(cVBF_WNLO,gt->genBosonPt) * gt->sf_qcdV;
            gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_WNLO,gt->genBosonPt) * gt->sf_qcdV;
            found=true;
          } else if (processType==kZEWK) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
	    gt->genBosonGenPx= PX;
            gt->genBosonGenPy= PY;
            gt->genBosonGenPz= PZ;
            gt->genBosonGenEn= EN;
            gt->genBosonGenFlav = ID;
	    gt->genBosonGenMass = MASS;
            gt->sf_qcdV_VBF = GetCorr(cVBF_EWKZ,gt->genBosonPt,gt->jot12Mass);
            gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
          } else if (processType==kWEWK) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
	    gt->genBosonGenPx= PX;
            gt->genBosonGenPy= PY;
            gt->genBosonGenPz= PZ;
            //gt->genBosonGenEn= EN;
            gt->genBosonGenFlav = ID;
	    gt->genBosonGenMass = MASS;
            gt->sf_qcdV_VBF = GetCorr(cVBF_EWKW,gt->genBosonPt,gt->jot12Mass);
            gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
          } else if (processType==kA) {
            // take the highest pT
            if (gen.pt() > gt->trueGenBosonPt) {
              gt->trueGenBosonPt = gen.pt();
              gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
	      gt->genBosonGenPx= PX;
	      gt->genBosonGenPy= PY;
	      gt->genBosonGenPz= PZ;
	      gt->genBosonGenEn= EN;
	      gt->genBosonGenFlav = ID;
	      gt->genBosonGenMass = MASS;
              gt->sf_qcdV = GetCorr(cANLO,gt->genBosonPt);
              gt->sf_ewkV = GetCorr(cAEWK,gt->genBosonPt);
              gt->sf_qcdV2j = GetCorr(cANLO2j,gt->genBosonPt);
            }
          }
        } // target matches
      }
    }

    tr.TriggerEvent("qcd/ewk SFs");

    gt->sf_metTrigVBF=1; gt->sf_metTrigZmmVBF=1;
    if (!isData) {
      gt->sf_metTrigVBF = GetCorr(cVBF_TrigMET,gt->pfmetnomu,gt->jot12Mass);
      gt->sf_metTrigZmmVBF = GetCorr(cVBF_TrigMETZmm,gt->pfmetnomu,gt->jot12Mass);
    }

    if (!isData && processType==kSignal) {
      bool found=false, foundbar=false;
      TLorentzVector vMediator(0,0,0,0);
      for (auto& gen : event.genParticles) {
        if (found && foundbar)
          break;
        if (abs(gen.pdgid) != 18)
          continue;
        if (gen.parent.isValid() && gen.parent->pdgid == gen.pdgid)
          continue;
        if (gen.pdgid == 18 && !found) {
          found = true;
          vMediator += gen.p4();
        } else if (gen.pdgid == -18 && !foundbar) {
          foundbar = true;
          vMediator += gen.p4();
        }
      }
      if (found && foundbar) {
        gt->trueGenBosonPt = vMediator.Pt();
        gt->genBosonPt = bound(gt->trueGenBosonPt,175,1200);
      }
      // gt->trueGenBosonPt = event.genMet.pt;
      // gt->genBosonPt = bound(gt->trueGenBosonPt,175,1200);
      // tr.TriggerEvent("signal gen kinematics");
    }

    //lepton SFs
    gt->sf_lepID=1; gt->sf_lepIso=1; gt->sf_lepTrack=1;
    if (!isData) {
      for (unsigned int iL=0; iL!=TMath::Min(gt->nLooseLep,2); ++iL) {
        auto* lep = looseLeps.at(iL);
        float pt = lep->pt(), eta = lep->eta(), aeta = TMath::Abs(eta);
        bool isTight = (iL==0 && gt->looseLep1IsTight) || (iL==1 && gt->looseLep2IsTight);
        auto* mu = dynamic_cast<panda::Muon*>(lep);
        if (mu!=NULL) {
          if (isTight) {
            gt->sf_lepID *= GetCorr(cMuTightID,aeta,pt);
            gt->sf_lepIso *= GetCorr(cMuTightIso,aeta,pt);
          } else {
            gt->sf_lepID *= GetCorr(cMuLooseID,aeta,pt);
            gt->sf_lepIso *= GetCorr(cMuLooseIso,aeta,pt);
          }
          gt->sf_lepTrack *= GetCorr(cMuReco,gt->npv);
        } else {
          if (isTight) {
            gt->sf_lepID *= GetCorr(cEleTight,eta,pt);
          } else {
            gt->sf_lepID *= GetCorr(cEleVeto,eta,pt);
          }
          gt->sf_lepTrack *= GetCorr(cEleReco,eta,pt);
        }
      }
    }

    tr.TriggerEvent("lepton SFs");

    //photon SF
    gt->sf_pho=1;
    if (!isData && gt->nLoosePhoton>0) {
      float pt = gt->loosePho1Pt, eta = gt->loosePho1Eta;
      if (gt->loosePho1IsTight)
        gt->sf_pho = GetCorr(cPho,eta,pt);
    }

    tr.TriggerEvent("photon SFs");

    // scale and PDF weights, if they exist
    gt->scaleUp = 1; gt->scaleDown = 1;
    gt->pdfUp = 1; gt->pdfDown = 1;
    if (!isData) {
      gt->pdfUp = 1 + event.genReweight.pdfDW;
      gt->pdfDown = 1 - event.genReweight.pdfDW;
      auto &genReweight = event.genReweight;
      for (unsigned iS=0; iS!=6; ++iS) {
        float s=1;
        switch (iS) {
          case 0:
            s = genReweight.r1f2DW; break;
          case 1:
            s = genReweight.r1f5DW; break;
          case 2:
            s = genReweight.r2f1DW; break;
          case 3:
            s = genReweight.r2f2DW; break;
          case 4:
            s = genReweight.r5f1DW; break;
          case 5:
            s = genReweight.r5f5DW; break;
          default:
            break;
        }
        gt->scale[iS] = s; 
        gt->scaleUp = max(float(gt->scaleUp),float(s));
        gt->scaleDown = min(float(gt->scaleDown),float(s));
      }
      tr.TriggerEvent("qcd uncertainties");

      unsigned nW = wIDs.size();
      if (nW) {
        for (unsigned iW=0; iW!=nW; ++iW) {
          gt->signal_weights[wIDs[iW]] = event.genReweight.genParam[iW];
        }
      }
    }

    gt->Fill();
    
  } // entry loop

  if (DEBUG) { PDebug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()

