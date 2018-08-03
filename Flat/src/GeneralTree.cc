#include "../interface/GeneralTree.h"
#include <iostream>

#define NJET 20
#define NSUBJET 2

GeneralTree::GeneralTree() {
//STARTCUSTOMCONST
  for (unsigned iS=0; iS!=6; ++iS) {
    scale[iS] = 1;
  }
  
  for (auto ibeta : ibetas) {
    for (auto N : Ns) {
      for (auto order : orders) {
        ECFParams p;
        p.ibeta = ibeta;
        p.N = N;
        p.order = order;
        ecfParams.push_back(p);
        fj1ECFNs[p] = -1;
      }
    }
  }

  for (unsigned iShift=0; iShift!=bNShift; ++iShift) {
    for (unsigned iJet=0; iJet!=bNJet; ++iJet) {
      for (unsigned iTags=0; iTags!=bNTags; ++iTags) {
        BTagParams p;
        p.jet = (BTagJet)iJet;
        p.tag = (BTagTags)iTags;
        p.shift = (BTagShift)iShift;
        btagParams.push_back(p);
        sf_btags[p] = 1;
      }
    }
  }
  
  for (unsigned iShift=0; iShift<nCsvShifts; iShift++) {
    csvShift shift = csvShifts[iShift];
    sf_csvWeights[shift] = 1;
  }
  for (unsigned int iSJ=0; iSJ!=NSUBJET; ++iSJ) {
    fj1sjPt[iSJ] = -1;
    fj1sjEta[iSJ] = -1;
    fj1sjPhi[iSJ] = -1;
    fj1sjM[iSJ] = -1;
    fj1sjCSV[iSJ] = -1;
    fj1sjQGL[iSJ] = -1;
  }
  for (unsigned int iJ=0; iJ!=NJET; ++iJ) {
    jetPt[iJ] = -1;
    jetPtUp[iJ] = -1;
    jetPtDown[iJ] = -1;
    jetEta[iJ] = -1;
    jetPhi[iJ] = -1;
    jetM[iJ] = -1;
    jetE[iJ] = -1;
    jetCSV[iJ] = -1;
    jetCMVA[iJ] = -1;
    jetIso[iJ] = -1;
    jetQGL[iJ] = -1;
    jetLeadingLepPt[iJ] = -1;
    jetLeadingLepPtRel[iJ] = -1;
    jetLeadingLepDeltaR[iJ] = -1;
    jetLeadingTrkPt[iJ] = -1;
    jetvtxPt[iJ] = -1;
    jetvtxMass[iJ] = -1;
    jetvtx3Dval[iJ] = -1;
    jetvtx3Derr[iJ] = -1;
    jetvtxNtrk[iJ] = -1;
    jetNLep[iJ] = -1;
    jetEMFrac[iJ] = -1;
    jetHadFrac[iJ] = -1;
    jetGenFlavor[iJ] = -1;
  }

//ENDCUSTOMCONST
}

GeneralTree::~GeneralTree() {
//STARTCUSTOMDEST
//ENDCUSTOMDEST
}

void GeneralTree::SetAuxTree(TTree *t) {
//STARTCUSTOMAUX
  for (auto p : ecfParams) { 
    TString ecfn(makeECFString(p));
    t->Branch("fj1"+ecfn,&(fj1ECFNs[p]),"fj1"+ecfn+"/F");
  }
  t->Branch("fj1Tau32SD",&(fj1Tau32SD),"fj1Tau32SD/F");
  t->Branch("fj1HTTFRec",&(fj1HTTFRec),"fj1HTTFRec/F");
//ENDCUSTOMAUX
}

void GeneralTree::Reset() {
//STARTCUSTOMRESET
  for (unsigned iS=0; iS!=6; ++iS) {
    scale[iS] = 1;
  }

  for (auto p : ecfParams) { 
    fj1ECFNs[p] = -1;
  }

  for (auto p : btagParams) { 
    sf_btags[p] = 1;
  }
  for (unsigned iShift=0; iShift<nCsvShifts; iShift++) {
    csvShift shift = csvShifts[iShift];
    sf_csvWeights[shift] = -1;
  }
  for (unsigned int iSJ=0; iSJ!=NSUBJET; ++iSJ) {
    fj1sjPt[iSJ] = -99;
    fj1sjEta[iSJ] = -99;
    fj1sjPhi[iSJ] = -99;
    fj1sjM[iSJ] = -99;
    fj1sjCSV[iSJ] = -99;
    fj1sjQGL[iSJ] = -99;
  }
  for (unsigned int iJ=0; iJ!=NJET; ++iJ) {
    jetPt[iJ] = -99;
    jetPtUp[iJ] = -99;
    jetPtDown[iJ] = -99;
    jetEta[iJ] = -99;
    jetPhi[iJ] = -99;
    jetM[iJ] = -99;
    jetE[iJ] = -99;
    jetCSV[iJ] = -99;
    jetCMVA[iJ] = -99;
    jetIso[iJ] = -99;
    jetQGL[iJ] = -99;
    jetLeadingLepPt[iJ] = -1;
    jetLeadingLepPtRel[iJ] = -1;
    jetLeadingLepDeltaR[iJ] = -1;
    jetLeadingTrkPt[iJ] = -1;
    jetvtxPt[iJ] = -1;
    jetvtxMass[iJ] = -1;
    jetvtx3Dval[iJ] = -1;
    jetvtx3Derr[iJ] = -1;
    jetvtxNtrk[iJ] = -1;
    jetEMFrac[iJ] = -1;
    jetHadFrac[iJ] = -1;
    jetNLep[iJ] = -99;
    jetGenPt[iJ] = -99;
    jetGenFlavor[iJ] = -99;
  }
  for (unsigned int iL=0; iL!=NLEP; ++iL) {
    muonPt[iL] = -99;
    muonEta[iL] = -99;
    muonPhi[iL] = -99;
    muonD0[iL] = -99;
    muonDZ[iL] = -99;
    muonSfLoose[iL] = -99;
    muonSfMedium[iL] = -99;
    muonSfTight[iL] = -99;
    muonSfUnc[iL] = -99;
    muonSfReco[iL] = -99;
    muonSelBit[iL] = -99;
    muonPdgId[iL] = -99;
    muonIsSoftMuon[iL] = -99;
    muonCombIso[iL] = -99;
    //muonIsGlobalMuon[iL] = -99;
    //muonIsTrackerMuon[iL] = -99;
    //muonNValidMuon[iL] = -99;
    //muonNValidPixel[iL] = -99;
    //muonTrkLayersWithMmt[iL] = -99;
    //muonPixLayersWithMmt[iL] = -99;
    //muonNMatched[iL] = -99;
    //muonChi2LocalPosition[iL] = -99;
    //muonTrkKink[iL] = -99;
    //muonValidFraction[iL] = -99;
    //muonNormChi2[iL] = -99;
    //muonSegmentCompatibility[iL] = -99;

    electronPt[iL] = -99;
    electronEta[iL] = -99;
    electronPhi[iL] = -99;
    electronD0[iL] = -99;
    electronDZ[iL] = -99;
    electronSfLoose[iL] = -99;
    electronSfMedium[iL] = -99;
    electronSfTight[iL] = -99;
    electronSfMvaWP90[iL] = -99;
    electronSfMvaWP80[iL] = -99;
    electronSfUnc[iL] = -99;
    electronSfReco[iL] = -99;
    electronSelBit[iL] = -99;
    electronPdgId[iL] = -99;
    //electronChIsoPh[iL] = -99;
    //electronNhIsoPh[iL] = -99;
    //electronPhIsoPh[iL] = -99;
    //electronEcalIso[iL] = -99;
    //electronHcalIso[iL] = -99;
    //electronTrackIso[iL] = -99;
    //electronIsoPUOffset[iL] = -99;
    //electronSieie[iL] = -99;
    //electronSipip[iL] = -99;
    //electronDEtaInSeed[iL] = -99;
    //electronDPhiIn[iL] = -99;
    //electronEseed[iL] = -99;
    //electronHOverE[iL] = -99;
    //electronEcalE[iL] = -99;
    //electronTrackP[iL] = -99;
    electronNMissingHits[iL] = -99;
    electronTripleCharge[iL] = -99;
    electronCombIso[iL] = -99;
  }
  for (auto iter=signal_weights.begin(); iter!=signal_weights.end(); ++iter) {
    signal_weights[iter->first] = 1; 
  }

//ENDCUSTOMRESET
    jot1PhiUp = -1;
    jot1PhiDown = -1;
    jot2PhiUp = -1;
    jot2PhiDown = -1;
    loosePho1SelBit = 0;
    looseGenPho1PdgId = 0;
    genFatJetNProngs = 0;
    genFatJetPt = -1;
    fj1NBPartons = 0;
    fj1NCPartons = 0;
    fj1Rho2 = -1;
    fj1RawRho2 = -1;
    fj1Rho = -1;
    fj1RawRho = -1;
    fj1NPartons = 0;
    fj1PartonM = -1;
    fj1PartonPt = -1;
    fj1PartonEta = -1;
    trkmetphi = -1;
    sf_zzUnc = 1;
    sf_zz = 1;
    sf_wz = 1;
    sf_vh = 1;
    sf_vhUp = 1;
    sf_vhDown = 1;
    genLep1Pt = -1;
    genLep1Eta = -1;
    genLep1Phi = -1;
    genLep1PdgId = 0;
    genLep2Pt = -1;
    genLep2Eta = -1;
    genLep2Phi = -1;
    genLep2PdgId = 0;
    genLep3Pt = -1;
    genLep3Eta = -1;
    genLep3Phi = -1;
    genLep3PdgId = 0;
    genLep4Pt = -1;
    genLep4Eta = -1;
    genLep4Phi = -1;
    genLep4PdgId = 0;
    looseGenLep1PdgId = 0;
    looseGenLep2PdgId = 0;
    looseGenLep3PdgId = 0;
    looseGenLep4PdgId = 0;
    whichRecoil = 0;
    genJet1Pt = -1;
    genJet2Pt = -1;
    genJet1Eta = -1;
    genJet2Eta = -1;
    genMjj = -1;
    sf_qcdV_VBF2lTight = 1;
    sf_qcdV_VBF2l = 1;
    barrelHTMiss = -1;
    barrelJet12Pt = -1;
    barrelJet1Pt = -1;
    barrelJet1Eta = -1;
    barrelHT = -1;
    sf_puUp = 1;
    sf_puDown = 1;
    genMuonPt = -1;
    genMuonEta = -1;
    genElectronPt = -1;
    genElectronEta = -1;
    genTauPt = -1;
    genTauEta = -1;
    badECALFilter = 0;
    sf_qcdV_VBFTight = 1;
    sf_metTrigVBF = 1;
    sf_metTrigZmmVBF = 1;
    sumETRaw = -1;
    jot1VBFID = 0;
    sf_metTrigZmm = 1;
    sf_qcdV_VBF = 1;
    jetNMBtags = 0;
    pfmetRaw = -1;
    nAK8jet = 0;
    ak81Pt = -1;
    ak81Eta = -1;
    ak81Phi = -1;
    ak81MaxCSV = -1;
    nB = 0;
    nBGenJets = 0;
    fj1MSDScaleUp_sj = -1;
    fj1MSDScaleDown_sj = -1;
    fj1MSDSmeared_sj = -1;
    fj1MSDSmearedUp_sj = -1;
    fj1MSDSmearedDown_sj = -1;
    fj1PtScaleUp_sj = -1;
    fj1PtScaleDown_sj = -1;
    fj1PtSmeared_sj = -1;
    fj1PtSmearedUp_sj = -1;
    fj1PtSmearedDown_sj = -1;
    jot2EtaUp = -1;
    jot2EtaDown = -1;
    jot1EtaUp = -1;
    jot1EtaDown = -1;
    jot1PtUp = -1;
    jot1PtDown = -1;
    jot2PtUp = -1;
    jot2PtDown = -1;
    jot12MassUp = -1;
    jot12DEtaUp = -1;
    jot12DPhiUp = -1;
    jot12MassDown = -1;
    jot12DEtaDown = -1;
    jot12DPhiDown = -1;
    pfmetUp = -1;
    pfmetDown = -1;
    pfUWmagUp = -1;
    pfUZmagUp = -1;
    pfUWWmagUp = -1;
    pfUAmagUp = -1;
    pfUmagUp = -1;
    pfUWmagDown = -1;
    pfUZmagDown = -1;
    pfUWWmagDown = -1;
    pfUAmagDown = -1;
    pfUmagDown = -1;
    nJot = 0;
    nJot_jesUp = 0;
    nJot_jesDown = 0;
    jot1Phi = -1;
    jot1Pt = -1;
    jot1GenPt = -1;
    jot1Eta = -1;
    jot2Phi = -1;
    jot2Pt = -1;
    jot2GenPt = -1;
    jot2Eta = -1;
    jot12DPhi = -1;
    isGS = 0;
    fj1SubMaxCSV = -1;
    runNumber = 0;
    lumiNumber = 0;
    eventNumber = 0;
    npv = 0;
    pu = 0;
    mcWeight = -1;
    trigger = 0;
    metFilter = 0;
    egmFilter = 0;
    filter_maxRecoil = -1;
    filter_whichRecoil = -1;
    sf_ewkV = 1;
    sf_qcdV = 1;
    sf_ewkV2j = 1;
    sf_qcdV2j = 1;
    sf_qcdTT = 1;
    sf_lepID = 1;
    sf_lepIso = 1;
    sf_lepTrack = 1;
    sf_pho = 1;
    sf_eleTrig = 1;
    sf_muTrig = 1;
    sf_phoTrig = 1;
    sf_metTrig = 1;
    sf_pu = 1;
    sf_npv = 1;
    sf_tt = 1;
    sf_tt_ext = 1;
    sf_tt_bound = 1;
    sf_tt8TeV = 1;
    sf_tt8TeV_ext = 1;
    sf_tt8TeV_bound = 1;
    sf_phoPurity = 1;
    pfmet = -1;
    pfmetphi = -1;
    pfmetnomu = -1;
    pfmetsig = -1;
    puppimet = -1;
    puppimetphi = -1;
    puppimetsig = -1;
    calomet = -1;
    calometphi = -1;
    pfcalobalance = -1;
    sumET = -1;
    trkmet = -1;
    puppiUWmag = -1;
    puppiUWphi = -1;
    puppiUZmag = -1;
    puppiUZphi = -1;
    puppiUAmag = -1;
    puppiUAphi = -1;
    puppiUperp = -1;
    puppiUpara = -1;
    puppiUmag = -1;
    puppiUphi = -1;
    pfUWmag = -1;
    pfUWphi = -1;
    pfUZmag = -1;
    pfUZphi = -1;
    pfUWWmag = -1;
    pfUWWphi = -1;
    pfUAmag = -1;
    pfUAphi = -1;
    pfUperp = -1;
    pfUpara = -1;
    pfUmag = -1;
    pfUphi = -1;
    dphipfmet = -1;
    dphipuppimet = -1;
    dphipuppiUW = -1;
    dphipuppiUZ = -1;
    dphipuppiUA = -1;
    dphipfUW = -1;
    dphipfUZ = -1;
    dphipfUWW = -1;
    dphipfUA = -1;
    dphipuppiU = -1;
    dphipfU = -1;
    trueGenBosonPt = -1;
    genBosonPt = -1;
    genBosonEta = -1;
    genBosonMass = -1;
    genBosonPhi = -1;
    genWPlusPt = -1;
    genWMinusPt = -1;
    genWPlusEta = -1;
    genWMinusEta = -1;
    genTopPt = -1;
    genTopIsHad = 0;
    genTopEta = -1;
    genAntiTopPt = -1;
    genAntiTopIsHad = 0;
    genAntiTopEta = -1;
    genTTPt = -1;
    genTTEta = -1;
    nJet = 0;
    nJet_jesUp = 0;
    nJet_jesDown = 0;
    nIsoJet = 0;
    jet1Flav = 0;
    jet1Phi = -1;
    jet1Pt = -1;
    jet1GenPt = -1;
    jet1Eta = -1;
    jet1CSV = -1;
    jet1CMVA = -1;
    jet1IsTight = 0;
    jet2Flav = 0;
    jet2Phi = -1;
    jet2Pt = -1;
    jet2GenPt = -1;
    jet2Eta = -1;
    jet2CSV = -1;
    jet2CMVA = -1;
    jet2EtaUp = -1;
    jet2EtaDown = -1;
    jet1EtaUp = -1;
    jet1EtaDown = -1;
    jet1PtUp = -1;
    jet1PtDown = -1;
    jet2PtUp = -1;
    jet2PtDown = -1;
    jet3Flav = 0;
    jet3Phi = -1;
    jet3Pt = -1;
    jet3GenPt = -1;
    jet3Eta = -1;
    jet3CSV = -1;
    isojet1Pt = -1;
    isojet1CSV = -1;
    isojet1Flav = 0;
    isojet2Pt = -1;
    isojet2CSV = -1;
    isojet2Flav = 0;
    jot12Mass = -1;
    jot12DEta = -1;
    jetNBtags = 0;
    isojetNBtags = 0;
    nFatjet = 0;
    fj1Tau32 = -1;
    fj1Tau21 = -1;
    fj1Tau32SD = -1;
    fj1Tau21SD = -1;
    fj1MSD = -1;
    fj1MSDScaleUp = -1;
    fj1MSDScaleDown = -1;
    fj1MSDSmeared = -1;
    fj1MSDSmearedUp = -1;
    fj1MSDSmearedDown = -1;
    fj1MSD_corr = -1;
    fj1Pt = -1;
    fj1PtScaleUp = -1;
    fj1PtScaleDown = -1;
    fj1PtSmeared = -1;
    fj1PtSmearedUp = -1;
    fj1PtSmearedDown = -1;
    fj1Phi = -1;
    fj1Eta = -1;
    fj1M = -1;
    fj1MaxCSV = -1;
    fj1MinCSV = -1;
    fj1DoubleCSV = -1;
    fj1Nbs = 0;
    fj1gbb = 0;
    fj1GenPt = -1;
    fj1GenSize = -1;
    fj1IsMatched = 0;
    fj1GenWPt = -1;
    fj1GenWSize = -1;
    fj1IsWMatched = 0;
    fj1HighestPtGen = 0;
    fj1HighestPtGenPt = -1;
    fj1IsTight = 0;
    fj1IsLoose = 0;
    fj1RawPt = -1;
    fj1NHF = 0;
    fj1HTTMass = -1;
    fj1HTTFRec = -1;
    fj1IsClean = 0;
    fj1NConst = 0;
    fj1NSDConst = 0;
    fj1EFrac100 = -1;
    fj1SDEFrac100 = -1;
    nHF = 0;
    nLoosePhoton = 0;
    nTightPhoton = 0;
    loosePho1IsTight = 0;
    loosePho1Pt = -1;
    loosePho1Eta = -1;
    loosePho1Phi = -1;
    nLooseLep = 0;
    nLooseElectron = 0;
    nLooseMuon = 0;
    nTightLep = 0;
    nTightElectron = 0;
    nTightMuon = 0;
    diLepMass = -1;
    nTau = 0;
    mT = -1;
    bosonpt = -1;
    bosoneta = -1;
    bosonphi = -1;
    bosonm = -1;
    bosoness = -1;
    bosondr = -1;
    bosondphibrf = -1;
    bosonm_reg = -1;
    bosonpt_reg = -1;
    bosonpt_jesUp = -1;
    bosoneta_jesUp = -1;
    bosonphi_jesUp = -1;
    bosonm_jesUp = -1;
    bosonm_reg_jesUp = -1;
    bosonpt_reg_jesUp = -1;
    bosonpt_jesDown = -1;
    bosoneta_jesDown = -1;
    bosonphi_jesDown = -1;
    bosonm_jesDown = -1;
    bosonm_reg_jesDown = -1;
    bosonpt_reg_jesDown = -1;
    bosonCosThetaJJ = -1;
    bosonCosThetaCSJ1 = -1;
    topMassLep1Met = -1;
    topMassLep1Met_jesUp = -1;
    topMassLep1Met_jesDown = -1;
    topWBosonCosThetaCS = -1;
    topWBosonPt = -1;
    topWBosonEta = -1;
    topWBosonPhi = -1;
    sumEtSoft1 = -1;
    nSoft2 = -1;
    nSoft5 = -1;
    nSoft10 = -1;
    scaleUp = 1;
    scaleDown = 1;
    pdfUp = 1;
    pdfDown = 1;
}

void GeneralTree::WriteTree(TTree *t) {
  treePtr = t;
//STARTCUSTOMWRITE
  for (auto iter=signal_weights.begin(); iter!=signal_weights.end(); ++iter) {
    Book("rw_"+iter->first,&(signal_weights[iter->first]),"rw_"+iter->first+"/F");
  }

  Book("nJet",&nJet,"nJet/I");
  Book("nJot",&nJot,"nJot/I");
  if (leptonic) {
    Book("nJot_jesUp",&nJot_jesUp,"nJot_jesUp/I");
    Book("nJot_jesDown",&nJot_jesDown,"nJot_jesDown/I");
    Book("jot1PhiUp",&jot1PhiUp,"jot1PhiUp/F");
    Book("jot1PhiDown",&jot1PhiDown,"jot1PhiDown/F");
    Book("jot2PhiUp",&jot2PhiUp,"jot2PhiUp/F");
    Book("jot2PhiDown",&jot2PhiDown,"jot2PhiDown/F");
  }
  Book("nLooseLep",&nLooseLep,"nLooseLep/I");
  Book("nLooseElectron",&nLooseElectron,"nLooseElectron/I");
  Book("nLooseMuon",&nLooseMuon,"nLooseMuon/I");
  Book("nTightLep",&nTightLep,"nTightLep/I");
  Book("nTightElectron",&nTightElectron,"nTightElectron/I");
  Book("nTightMuon",&nTightMuon,"nTightMuon/I");
  Book("muonPt",muonPt,"muonPt[nLooseMuon]/F");
  Book("muonEta",muonEta,"muonEta[nLooseMuon]/F");
  Book("muonPhi",muonPhi,"muonPhi[nLooseMuon]/F");
  Book("muonSelBit",muonSelBit,"muonSelBit[nLooseMuon]/I");
  Book("muonPdgId",muonPdgId,"muonPdgId[nLooseMuon]/I");
  Book("electronPt",electronPt,"electronPt[nLooseElectron]/F");
  Book("electronEta",electronEta,"electronEta[nLooseElectron]/F");
  Book("electronPhi",electronPhi,"electronPhi[nLooseElectron]/F");
  Book("electronSelBit",electronSelBit,"electronSelBit[nLooseElectron]/I");
  Book("electronPdgId",electronPdgId,"electronPdgId[nLooseElectron]/I");
  Book("nJet_jesUp",&nJet_jesUp,"nJet_jesUp/I");
  Book("nJet_jesDown",&nJet_jesDown,"nJet_jesDown/I");
  Book("pfmetUp",&pfmetUp,"pfmetUp/F");
  Book("pfmetDown",&pfmetDown,"pfmetDown/F");
  Book("jetPt",jetPt,"jetPt[nJot]/F");
  Book("jetPtUp",jetPtUp,"jetPtUp[nJot]/F");
  Book("jetPtDown",jetPtDown,"jetPtDown[nJot]/F");
  Book("jetEta",jetEta,"jetEta[nJot]/F");
  Book("jetPhi",jetPhi,"jetPhi[nJot]/F");
  Book("jetM",jetM,"jetM[nJot]/F");
  Book("jetE",jetE,"jetE[nJot]/F");
  Book("jetCSV",jetCSV,"jetCSV[nJot]/F");
  Book("jetCMVA",jetCMVA,"jetCMVA[nJot]/F");
  Book("jetIso",jetIso,"jetIso[nJot]/F");
  Book("jetQGL",jetQGL,"jetQGL[nJot]/F");
  Book("jetLeadingLepPt",jetLeadingLepPt,"jetLeadingLepPt[nJot]/F");
  Book("jetEMFrac",jetEMFrac,"jetEMFrac[nJot]/F");
  Book("jetHadFrac",jetHadFrac,"jetHadFrac[nJot]/F");
  Book("jetNLep",jetNLep,"jetNLep[nJot]/I");
  Book("jetGenPt",jetGenPt,"jetGenPt[nJot]/F");
  Book("jetGenFlavor",jetGenFlavor,"jetGenFlavor[nJot]/I");
  Book("jet1EtaUp",&jet1EtaUp,"jet1EtaUp/F");
  Book("jet1EtaDown",&jet1EtaDown,"jet1EtaDown/F");
  Book("jet1PtUp",&jet1PtUp,"jet1PtUp/F");
  Book("jet1PtDown",&jet1PtDown,"jet1PtDown/F");
  Book("jet2EtaUp",&jet2EtaUp,"jet2EtaUp/F");
  Book("jet2EtaDown",&jet2EtaDown,"jet2EtaDown/F");
  Book("jet2PtUp",&jet2PtUp,"jet2PtUp/F");
  Book("jet2PtDown",&jet2PtDown,"jet2PtDown/F");
/*
    Book("jetLeadingLepPtRel",jetLeadingLepPtRel,"jetLeadingLepPtRel[nJot]/F");
    Book("jetLeadingLepDeltaR",jetLeadingLepDeltaR,"jetLeadingLepDeltaR[nJot]/F");
    Book("jetLeadingTrkPt",jetLeadingTrkPt,"jetLeadingTrkPt[nJot]/F");
    Book("jetvtxPt",jetvtxPt,"jetvtxPt[nJot]/F");
    Book("jetvtxMass",jetvtxMass,"jetvtxMass[nJot]/F");
    Book("jetvtx3Dval",jetvtx3Dval,"jetvtx3Dval[nJot]/F");
    Book("jetvtx3Derr",jetvtx3Derr,"jetvtx3Derr[nJot]/F");
    Book("jetvtxNtrk",jetvtxNtrk,"jetvtxNtrk[nJot]/I");
    Book("fj1sjPt",fj1sjPt,"fj1sjPt[2]/F");
    Book("fj1sjPhi",fj1sjPhi,"fj1sjPhi[2]/F");
    Book("fj1sjEta",fj1sjEta,"fj1sjEta[2]/F");
    Book("fj1sjM",fj1sjM,"fj1sjM[2]/F");
    Book("fj1sjCSV",fj1sjCSV,"fj1sjCSV[2]/F");
    Book("fj1sjQGL",fj1sjQGL,"fj1sjQGL[2]/F");
    Book("fj1Nbs",&fj1Nbs,"fj1Nbs/I");
    Book("fj1gbb",&fj1gbb,"fj1gbb/I");
*/
  if (resolved) {
    Book("bosonpt",&bosonpt,"bosonpt/F");
    Book("bosoneta",&bosoneta,"bosoneta/F");
    Book("bosonphi",&bosonphi,"bosonphi/F");
    Book("bosonm",&bosonm,"bosonm/F");
    Book("bosoness",&bosoness,"bosoness/F");
    Book("bosondr",&bosondr,"bosondr/F");
    Book("bosondphibrf",&bosondphibrf,"bosondphibrf/F");
    Book("bosonm_reg",&bosonm_reg,"bosonm_reg/F");
    Book("bosonpt_reg",&bosonpt_reg,"bosonpt_reg/F");
    Book("bosonjtidx",bosonjtidx,"bosonjtidx[2]/I");
    Book("bosonpt_jesUp",&bosonpt_jesUp,"bosonpt_jesUp/F");
    Book("bosoneta_jesUp",&bosoneta_jesUp,"bosoneta_jesUp/F");
    Book("bosonphi_jesUp",&bosonphi_jesUp,"bosonphi_jesUp/F");
    Book("bosonm_jesUp",&bosonm_jesUp,"bosonm_jesUp/F");
    Book("bosonm_reg_jesUp",&bosonm_reg_jesUp,"bosonm_reg_jesUp/F");
    Book("bosonpt_reg_jesUp",&bosonpt_reg_jesUp,"bosonpt_reg_jesUp/F");
    Book("bosonpt_jesDown",&bosonpt_jesDown,"bosonpt_jesDown/F");
    Book("bosoneta_jesDown",&bosoneta_jesDown,"bosoneta_jesDown/F");
    Book("bosonphi_jesDown",&bosonphi_jesDown,"bosonphi_jesDown/F");
    Book("bosonm_jesDown",&bosonm_jesDown,"bosonm_jesDown/F");
    Book("bosonm_reg_jesDown",&bosonm_reg_jesDown,"bosonm_reg_jesDown/F");
    Book("bosonpt_reg_jesDown",&bosonpt_reg_jesDown,"bosonpt_reg_jesDown/F");
    Book("bosonjtidx",bosonjtidx,"bosonjtidx[2]/I");
    Book("bosonCosThetaJJ",&bosonCosThetaJJ,"bosonCosThetaJJ/F");
    Book("bosonCosThetaCSJ1",&bosonCosThetaCSJ1,"bosonCosThetaCSJ1/F");
  }
/*
    Book("topMassLep1Met",&topMassLep1Met,"topMassLep1Met/F");
    Book("topMassLep1Met_jesUp",&topMassLep1Met_jesUp,"topMassLep1Met_jesUp/F");
    Book("topMassLep1Met_jesDown",&topMassLep1Met_jesDown,"topMassLep1Met_jesDown/F");
    Book("topWBosonCosThetaCS",&topWBosonCosThetaCS,"topWBosonCosThetaCS/F");
    Book("topWBosonPt",&topWBosonPt,"topWBosonPt/F");
    Book("topWBosonEta",&topWBosonEta,"topWBosonEta/F");
    Book("topWBosonPhi",&topWBosonPhi,"topWBosonPhi/F");
    Book("sumEtSoft1",&sumEtSoft1,"sumEtSoft1/F");
    Book("nSoft2",&nSoft2,"nSoft2/I");
    Book("nSoft5",&nSoft5,"nSoft5/I");
    Book("nSoft10",&nSoft10,"nSoft10/I");
    Book("jetRegFac",jetRegFac,"jetRegFac[2]/F");
*/

  if (vbf || leptonic) {
    Book("jot1Phi",&jot1Phi,"jot1Phi/F");
    Book("jot1Eta",&jot1Eta,"jot1Eta/F");
    Book("jot1Pt",&jot1Pt,"jot1Pt/F");
    Book("jot2Phi",&jot2Phi,"jot2Phi/F");
    Book("jot2Eta",&jot2Eta,"jot2Eta/F");
    Book("jot2Pt",&jot2Pt,"jot2Pt/F");
    Book("jot1EtaUp",&jot1EtaUp,"jot1EtaUp/F");
    Book("jot1EtaDown",&jot1EtaDown,"jot1EtaDown/F");
    Book("jot1PtUp",&jot1PtUp,"jot1PtUp/F");
    Book("jot1PtDown",&jot1PtDown,"jot1PtDown/F");
    Book("jot2EtaUp",&jot2EtaUp,"jot2EtaUp/F");
    Book("jot2EtaDown",&jot2EtaDown,"jot2EtaDown/F");
    Book("jot2PtUp",&jot2PtUp,"jot2PtUp/F");
    Book("jot2PtDown",&jot2PtDown,"jot2PtDown/F");
  }
  if (vbf) { 
    Book("jot1GenPt",&jot1GenPt,"jot1GenPt/F");
    Book("jot2GenPt",&jot2GenPt,"jot2GenPt/F");
    Book("jot12DPhi",&jot12DPhi,"jot12DPhi/F");
    Book("jot12Mass",&jot12Mass,"jot12Mass/F");
    Book("jot12DEta",&jot12DEta,"jot12DEta/F");
    Book("pfmetUp",&pfmetUp,"pfmetUp/F");
    Book("pfmetDown",&pfmetDown,"pfmetDown/F");
    Book("pfUWmagUp",&pfUWmagUp,"pfUWmagUp/F");
    Book("pfUZmagUp",&pfUZmagUp,"pfUZmagUp/F");
    Book("pfUWWmagUp",&pfUWWmagUp,"pfUWWmagUp/F");
    Book("pfUAmagUp",&pfUAmagUp,"pfUAmagUp/F");
    Book("pfUmagUp",&pfUmagUp,"pfUmagUp/F");
    Book("pfUWmagDown",&pfUWmagDown,"pfUWmagDown/F");
    Book("pfUZmagDown",&pfUZmagDown,"pfUZmagDown/F");
    Book("pfUWWmagDown",&pfUWWmagDown,"pfUWWmagDown/F");
    Book("pfUAmagDown",&pfUAmagDown,"pfUAmagDown/F");
    Book("pfUmagDown",&pfUmagDown,"pfUmagDown/F");
    Book("jot12MassUp",&jot12MassUp,"jot12MassUp/F");
    Book("jot12DEtaUp",&jot12DEtaUp,"jot12DEtaUp/F");
    Book("jot12DPhiUp",&jot12DPhiUp,"jot12DPhiUp/F");
    Book("jot12MassDown",&jot12MassDown,"jot12MassDown/F");
    Book("jot12DEtaDown",&jot12DEtaDown,"jot12DEtaDown/F");
    Book("jot12DPhiDown",&jot12DPhiDown,"jot12DPhiDown/F");
    Book("jot1VBFID",&jot1VBFID,"jot1VBFID/I");
    Book("sf_qcdV_VBF2lTight",&sf_qcdV_VBF2lTight,"sf_qcdV_VBF2lTight/F");
    Book("sf_qcdV_VBF2l",&sf_qcdV_VBF2l,"sf_qcdV_VBF2l/F");
    Book("barrelHTMiss",&barrelHTMiss,"barrelHTMiss/F");
    Book("barrelJet12Pt",&barrelJet12Pt,"barrelJet12Pt/F");
    Book("barrelJet1Pt",&barrelJet1Pt,"barrelJet1Pt/F");
    Book("barrelJet1Eta",&barrelJet1Eta,"barrelJet1Eta/F");
    Book("barrelHT",&barrelHT,"barrelHT/F");
    Book("sf_puUp",&sf_puUp,"sf_puUp/F");
    Book("sf_puDown",&sf_puDown,"sf_puDown/F");
    Book("genMuonPt",&genMuonPt,"genMuonPt/F");
    Book("genMuonEta",&genMuonEta,"genMuonEta/F");
    Book("genElectronPt",&genElectronPt,"genElectronPt/F");
    Book("genElectronEta",&genElectronEta,"genElectronEta/F");
    Book("genTauPt",&genTauPt,"genTauPt/F");
    Book("genTauEta",&genTauEta,"genTauEta/F");
    Book("sf_qcdV_VBFTight",&sf_qcdV_VBFTight,"sf_qcdV_VBFTight/F");
    Book("sf_metTrigVBF",&sf_metTrigVBF,"sf_metTrigVBF/F");
    Book("sf_metTrigZmmVBF",&sf_metTrigZmmVBF,"sf_metTrigZmmVBF/F");
    Book("sumETRaw",&sumETRaw,"sumETRaw/F");
    Book("sf_metTrigZmm",&sf_metTrigZmm,"sf_metTrigZmm/F");
    Book("sf_qcdV_VBF",&sf_qcdV_VBF,"sf_qcdV_VBF/F");
  }
  Book("scale",scale,"scale[6]/F");

  if (boosted) {
  for (auto p : ecfParams) { 
    TString ecfn(makeECFString(p));
    Book("fj1"+ecfn,&(fj1ECFNs[p]),"fj1"+ecfn+"/F");
    }
  }
  for (auto p : btagParams) {
    TString btagn(makeBTagSFString(p));
    Book(btagn,&(sf_btags[p]),btagn+"/F");
  }
  if (leptonic) {
    // Per-leg lepton scale factors, only needed for the nice analyses :-)
    Book("muonSfLoose",muonSfLoose,"muonSfLoose[nLooseMuon]/F");
    Book("muonSfMedium",muonSfMedium,"muonSfMedium[nLooseMuon]/F");
    Book("muonSfTight",muonSfTight,"muonSfTight[nLooseMuon]/F");
    Book("muonSfUnc",muonSfUnc,"muonSfUnc[nLooseMuon]/F");
    Book("muonSfReco",muonSfReco,"muonSfReco[nLooseMuon]/F");
    Book("electronSfLoose",electronSfLoose,"electronSfLoose[nLooseElectron]/F");
    Book("electronSfMedium",electronSfMedium,"electronSfMedium[nLooseElectron]/F");
    Book("electronSfTight",electronSfTight,"electronSfTight[nLooseElectron]/F");
    Book("electronSfMvaWP90",electronSfMvaWP90,"electronSfMvaWP90[nLooseElectron]/F");
    Book("electronSfMvaWP80",electronSfMvaWP80,"electronSfMvaWP80[nLooseElectron]/F");
    Book("electronSfUnc",electronSfUnc,"electronSfUnc[nLooseElectron]/F");
    Book("electronSfReco",electronSfReco,"electronSfReco[nLooseElectron]/F");
    // Advanced muon properties for nerds
    Book("muonD0",muonD0,"muonD0[nLooseMuon]/F");
    Book("muonDZ",muonDZ,"muonDZ[nLooseMuon]/F");
    Book("muonIsSoftMuon",muonIsSoftMuon,"muonIsSoftMuon[nLooseMuon]/I");
    Book("muonCombIso",muonCombIso,"muonCombIso[nLooseMuon]/F");
    // Advanced electron properties for nerds
    Book("electronD0",electronD0,"electronD0[nLooseElectron]/F");
    Book("electronDZ",electronDZ,"electronDZ[nLooseElectron]/F");
    Book("electronNMissingHits",electronNMissingHits,"electronNMissingHits[nLooseElectron]/I");
    Book("electronTripleCharge",electronTripleCharge,"electronTripleCharge[nLooseElectron]/I");
    Book("electronCombIso",electronCombIso,"electronCombIso[nLooseElectron]/F");
    // Gen study
    Book("sf_zz",&sf_zz,"sf_zz/F");
    Book("sf_zzUnc",&sf_zzUnc,"sf_zzUnc/F");
    Book("sf_wz",&sf_wz,"sf_wz/F");
    Book("sf_vh",&sf_vh,"sf_vh/F");
    Book("sf_vhUp",&sf_vhUp,"sf_vhUp/F");
    Book("sf_vhDown",&sf_vhDown,"sf_vhDown/F");
    Book("genLep1Pt",&genLep1Pt,"genLep1Pt/F");
    Book("genLep1Eta",&genLep1Eta,"genLep1Eta/F");
    Book("genLep1Phi",&genLep1Phi,"genLep1Phi/F");
    Book("genLep1PdgId",&genLep1PdgId,"genLep1PdgId/I");
    Book("genLep2Pt",&genLep2Pt,"genLep2Pt/F");
    Book("genLep2Eta",&genLep2Eta,"genLep2Eta/F");
    Book("genLep2Phi",&genLep2Phi,"genLep2Phi/F");
    Book("genLep2PdgId",&genLep2PdgId,"genLep2PdgId/I");
    Book("genLep3Pt",&genLep3Pt,"genLep3Pt/F");
    Book("genLep3Eta",&genLep3Eta,"genLep3Eta/F");
    Book("genLep3Phi",&genLep3Phi,"genLep3Phi/F");
    Book("genLep3PdgId",&genLep3PdgId,"genLep3PdgId/I");
    Book("genLep4Pt",&genLep4Pt,"genLep4Pt/F");
    Book("genLep4Eta",&genLep4Eta,"genLep4Eta/F");
    Book("genLep4Phi",&genLep4Phi,"genLep4Phi/F");
    Book("genLep4PdgId",&genLep4PdgId,"genLep4PdgId/I");
    Book("looseGenLep1PdgId",&looseGenLep1PdgId,"looseGenLep1PdgId/I");
    Book("looseGenLep2PdgId",&looseGenLep2PdgId,"looseGenLep2PdgId/I");
    Book("looseGenLep3PdgId",&looseGenLep3PdgId,"looseGenLep3PdgId/I");
    Book("looseGenLep4PdgId",&looseGenLep4PdgId,"looseGenLep4PdgId/I");
  } else {
    Book("sf_lepID",&sf_lepID,"sf_lepID/F");
    Book("sf_lepIso",&sf_lepIso,"sf_lepIso/F");
    Book("sf_lepTrack",&sf_lepTrack,"sf_lepTrack/F");
  }
  if (fatjet) {
    Book("fj1Pt",&fj1Pt,"fj1Pt/F");
    Book("fj1Phi",&fj1Phi,"fj1Phi/F");
    Book("fj1Eta",&fj1Eta,"fj1Eta/F");
    Book("fj1M",&fj1M,"fj1M/F");
    Book("nFatjet",&nFatjet,"nFatjet/I");
    Book("nIsoJet",&nIsoJet,"nIsoJet/I");
    Book("nAK8jet",&nAK8jet,"nAK8jet/I");
    Book("ak81Pt",&ak81Pt,"ak81Pt/F");
    Book("ak81Eta",&ak81Eta,"ak81Eta/F");
    Book("ak81Phi",&ak81Phi,"ak81Phi/F");
  }
  if (boosted) {
    Book("fj1Nbs",&fj1Nbs,"fj1Nbs/I");
    Book("fj1gbb",&fj1gbb,"fj1gbb/I");
    Book("fj1Tau32",&fj1Tau32,"fj1Tau32/F");
    Book("fj1Tau21",&fj1Tau21,"fj1Tau21/F");
    Book("fj1Tau32SD",&fj1Tau32SD,"fj1Tau32SD/F");
    Book("fj1Tau21SD",&fj1Tau21SD,"fj1Tau21SD/F");
    Book("fj1MSD",&fj1MSD,"fj1MSD/F");
    Book("fj1MSDScaleUp",&fj1MSDScaleUp,"fj1MSDScaleUp/F");
    Book("fj1MSDScaleDown",&fj1MSDScaleDown,"fj1MSDScaleDown/F");
    Book("fj1MSDSmeared",&fj1MSDSmeared,"fj1MSDSmeared/F");
    Book("fj1MSDSmearedUp",&fj1MSDSmearedUp,"fj1MSDSmearedUp/F");
    Book("fj1MSDSmearedDown",&fj1MSDSmearedDown,"fj1MSDSmearedDown/F");
    Book("fj1MSD_corr",&fj1MSD_corr,"fj1MSD_corr/F");
    Book("fj1PtScaleUp",&fj1PtScaleUp,"fj1PtScaleUp/F");
    Book("fj1PtScaleDown",&fj1PtScaleDown,"fj1PtScaleDown/F");
    Book("fj1PtSmeared",&fj1PtSmeared,"fj1PtSmeared/F");
    Book("fj1PtSmearedUp",&fj1PtSmearedUp,"fj1PtSmearedUp/F");
    Book("fj1PtSmearedDown",&fj1PtSmearedDown,"fj1PtSmearedDown/F");
    Book("fj1MaxCSV",&fj1MaxCSV,"fj1MaxCSV/F");
    Book("fj1MinCSV",&fj1MinCSV,"fj1MinCSV/F");
    Book("fj1DoubleCSV",&fj1DoubleCSV,"fj1DoubleCSV/F");
    Book("fj1GenPt",&fj1GenPt,"fj1GenPt/F");
    Book("fj1GenSize",&fj1GenSize,"fj1GenSize/F");
    Book("fj1IsMatched",&fj1IsMatched,"fj1IsMatched/I");
    Book("fj1GenWPt",&fj1GenWPt,"fj1GenWPt/F");
    Book("fj1GenWSize",&fj1GenWSize,"fj1GenWSize/F");
    Book("fj1IsWMatched",&fj1IsWMatched,"fj1IsWMatched/I");
    Book("fj1HighestPtGen",&fj1HighestPtGen,"fj1HighestPtGen/I");
    Book("fj1HighestPtGenPt",&fj1HighestPtGenPt,"fj1HighestPtGenPt/F");
    Book("fj1IsTight",&fj1IsTight,"fj1IsTight/I");
    Book("fj1IsLoose",&fj1IsLoose,"fj1IsLoose/I");
    Book("fj1RawPt",&fj1RawPt,"fj1RawPt/F");
    Book("fj1NHF",&fj1NHF,"fj1NHF/I");
  //  Book("fj1HTTMass",&fj1HTTMass,"fj1HTTMass/F");
  //  Book("fj1HTTFRec",&fj1HTTFRec,"fj1HTTFRec/F");
    Book("fj1IsClean",&fj1IsClean,"fj1IsClean/I");
    Book("fj1NConst",&fj1NConst,"fj1NConst/I");
    Book("fj1NSDConst",&fj1NSDConst,"fj1NSDConst/I");
    Book("fj1EFrac100",&fj1EFrac100,"fj1EFrac100/F");
    Book("fj1SDEFrac100",&fj1SDEFrac100,"fj1SDEFrac100/F");
/*
    Book("fj1MSDScaleUp_sj",&fj1MSDScaleUp_sj,"fj1MSDScaleUp_sj/F");
    Book("fj1MSDScaleDown_sj",&fj1MSDScaleDown_sj,"fj1MSDScaleDown_sj/F");
    Book("fj1MSDSmeared_sj",&fj1MSDSmeared_sj,"fj1MSDSmeared_sj/F");
    Book("fj1MSDSmearedUp_sj",&fj1MSDSmearedUp_sj,"fj1MSDSmearedUp_sj/F");
    Book("fj1MSDSmearedDown_sj",&fj1MSDSmearedDown_sj,"fj1MSDSmearedDown_sj/F");
    Book("fj1PtScaleUp_sj",&fj1PtScaleUp_sj,"fj1PtScaleUp_sj/F");
    Book("fj1PtScaleDown_sj",&fj1PtScaleDown_sj,"fj1PtScaleDown_sj/F");
    Book("fj1PtSmeared_sj",&fj1PtSmeared_sj,"fj1PtSmeared_sj/F");
    Book("fj1PtSmearedUp_sj",&fj1PtSmearedUp_sj,"fj1PtSmearedUp_sj/F");
    Book("fj1PtSmearedDown_sj",&fj1PtSmearedDown_sj,"fj1PtSmearedDown_sj/F");
*/
    Book("fj1SubMaxCSV",&fj1SubMaxCSV,"fj1SubMaxCSV/F");
    Book("isojet1Pt",&isojet1Pt,"isojet1Pt/F");
    Book("isojet1CSV",&isojet1CSV,"isojet1CSV/F");
    Book("isojet1Flav",&isojet1Flav,"isojet1Flav/I");
    Book("isojet2Pt",&isojet2Pt,"isojet2Pt/F");
    Book("isojet2CSV",&isojet2CSV,"isojet2CSV/F");
    Book("isojet2Flav",&isojet2Flav,"isojet2Flav/I");
    Book("isojetNBtags",&isojetNBtags,"isojetNBtags/I");
    Book("fj1NPartons",&fj1NPartons,"fj1NPartons/I");
    Book("fj1PartonM",&fj1PartonM,"fj1PartonM/F");
    Book("fj1PartonPt",&fj1PartonPt,"fj1PartonPt/F");
    Book("fj1PartonEta",&fj1PartonEta,"fj1PartonEta/F");
    Book("ak81MaxCSV",&ak81MaxCSV,"ak81MaxCSV/F");
    Book("genFatJetPt",&genFatJetPt,"genFatJetPt/F");
    Book("fj1NBPartons",&fj1NBPartons,"fj1NBPartons/I");
    Book("fj1NCPartons",&fj1NCPartons,"fj1NCPartons/I");
    Book("genFatJetNProngs",&genFatJetNProngs,"genFatJetNProngs/I");
  }
  if (hfCounting) {
    Book("nHF",&nHF,"nHF/I");
    Book("nB",&nB,"nB/I");
    Book("nBGenJets",&nBGenJets,"nBGenJets/I");
  }
  if (btagWeights) { 
    for (unsigned iShift=0; iShift<nCsvShifts; iShift++) {
      csvShift shift = csvShifts[iShift];
      TString csvWeightString = makeCsvWeightString(shift, useCMVA);
      Book(csvWeightString, &(sf_csvWeights[shift]), csvWeightString+"/F");
    }
  }
  if (photonic) {
    Book("loosePho1SelBit",&loosePho1SelBit,"loosePho1SelBit/I");
    Book("looseGenPho1PdgId",&looseGenPho1PdgId,"looseGenPho1PdgId/I");
  }
//ENDCUSTOMWRITE

 // Book("genFatJetNProngs",&genFatJetNProngs,"genFatJetNProngs/I");
  Book("genJet1Pt",&genJet1Pt,"genJet1Pt/F");
  Book("genJet2Pt",&genJet2Pt,"genJet2Pt/F");
  Book("genJet1Eta",&genJet1Eta,"genJet1Eta/F");
  Book("genJet2Eta",&genJet2Eta,"genJet2Eta/F");
  Book("genMjj",&genMjj,"genMjj/F");
  Book("jet1Flav",&jet1Flav,"jet1Flav/I");
  Book("jet1Phi",&jet1Phi,"jet1Phi/F");
  Book("jet1Pt",&jet1Pt,"jet1Pt/F");
  Book("jet1GenPt",&jet1GenPt,"jet1GenPt/F");
  Book("jet1Eta",&jet1Eta,"jet1Eta/F");
  Book("jet1CSV",&jet1CSV,"jet1CSV/F");
  Book("jet1CMVA",&jet1CMVA,"jet1CMVA/F");
  Book("jet1IsTight",&jet1IsTight,"jet1IsTight/I");
  Book("jet2Flav",&jet2Flav,"jet2Flav/I");
  Book("jet2Phi",&jet2Phi,"jet2Phi/F");
  Book("jet2Pt",&jet2Pt,"jet2Pt/F");
  Book("jet2GenPt",&jet2GenPt,"jet2GenPt/F");
  Book("jet2Eta",&jet2Eta,"jet2Eta/F");
  Book("jet2CSV",&jet2CSV,"jet2CSV/F");
  Book("jet2CMVA",&jet2CMVA,"jet2CMVA/F");
/*
  Book("fj1Rho2",&fj1Rho2,"fj1Rho2/F");
  Book("fj1RawRho2",&fj1RawRho2,"fj1RawRho2/F");
  Book("fj1Rho",&fj1Rho,"fj1Rho/F");
  Book("fj1RawRho",&fj1RawRho,"fj1RawRho/F");
*/
  Book("puppiUWmag",&puppiUWmag,"puppiUWmag/F");
  Book("puppiUWphi",&puppiUWphi,"puppiUWphi/F");
  Book("puppiUZmag",&puppiUZmag,"puppiUZmag/F");
  Book("puppiUZphi",&puppiUZphi,"puppiUZphi/F");
  Book("puppiUAmag",&puppiUAmag,"puppiUAmag/F");
  Book("puppiUAphi",&puppiUAphi,"puppiUAphi/F");
  Book("puppiUperp",&puppiUperp,"puppiUperp/F");
  Book("puppiUpara",&puppiUpara,"puppiUpara/F");
  Book("puppiUmag",&puppiUmag,"puppiUmag/F");
  Book("puppiUphi",&puppiUphi,"puppiUphi/F");
  Book("pfUWmag",&pfUWmag,"pfUWmag/F");
  Book("pfUWWmag",&pfUWWmag,"pfUWWmag/F");
  Book("pfUWphi",&pfUWphi,"pfUWphi/F");
  Book("pfUZmag",&pfUZmag,"pfUZmag/F");
  Book("pfUZphi",&pfUZphi,"pfUZphi/F");
  Book("pfUAmag",&pfUAmag,"pfUAmag/F");
  Book("pfUAphi",&pfUAphi,"pfUAphi/F");
  Book("pfUperp",&pfUperp,"pfUperp/F");
  Book("pfUpara",&pfUpara,"pfUpara/F");
  Book("pfUmag",&pfUmag,"pfUmag/F");
  Book("pfUphi",&pfUphi,"pfUphi/F");
  Book("dphipuppiUW",&dphipuppiUW,"dphipuppiUW/F");
  Book("dphipuppiUZ",&dphipuppiUZ,"dphipuppiUZ/F");
  Book("dphipuppiUA",&dphipuppiUA,"dphipuppiUA/F");
  Book("dphipfUW",&dphipfUW,"dphipfUW/F");
  Book("dphipfUWW",&dphipfUWW,"dphipfUWW/F");
  Book("dphipfUZ",&dphipfUZ,"dphipfUZ/F");
  Book("dphipfUA",&dphipfUA,"dphipfUA/F");
  Book("dphipuppiU",&dphipuppiU,"dphipuppiU/F");
  Book("dphipfU",&dphipfU,"dphipfU/F");
//  Book("isGS",&isGS,"isGS/I");
//  Book("trkmetphi",&trkmetphi,"trkmetphi/F");
  Book("whichRecoil",&whichRecoil,"whichRecoil/I");
  Book("badECALFilter",&badECALFilter,"badECALFilter/I");
  Book("jetNMBtags",&jetNMBtags,"jetNMBtags/I");
  Book("pfmetRaw",&pfmetRaw,"pfmetRaw/F");
  Book("runNumber",&runNumber,"runNumber/I");
  Book("lumiNumber",&lumiNumber,"lumiNumber/I");
  Book("eventNumber",&eventNumber,"eventNumber/l");
  Book("npv",&npv,"npv/I");
  Book("pu",&pu,"pu/I");
  Book("mcWeight",&mcWeight,"mcWeight/F");
  Book("trigger",&trigger,"trigger/I");
  Book("metFilter",&metFilter,"metFilter/I");
  Book("egmFilter",&egmFilter,"egmFilter/I");
  Book("filter_maxRecoil",&filter_maxRecoil,"filter_maxRecoil/F");
  Book("filter_whichRecoil",&filter_whichRecoil,"filter_whichRecoil/F");
  Book("sf_ewkV",&sf_ewkV,"sf_ewkV/F");
  Book("sf_qcdV",&sf_qcdV,"sf_qcdV/F");
  Book("sf_ewkV2j",&sf_ewkV2j,"sf_ewkV2j/F");
  Book("sf_qcdV2j",&sf_qcdV2j,"sf_qcdV2j/F");
  Book("sf_qcdTT",&sf_qcdTT,"sf_qcdTT/F");
  Book("sf_pho",&sf_pho,"sf_pho/F");
  Book("sf_eleTrig",&sf_eleTrig,"sf_eleTrig/F");
  Book("sf_muTrig",&sf_muTrig,"sf_muTrig/F");
  Book("sf_phoTrig",&sf_phoTrig,"sf_phoTrig/F");
  Book("sf_metTrig",&sf_metTrig,"sf_metTrig/F");
  Book("sf_pu",&sf_pu,"sf_pu/F");
  Book("sf_npv",&sf_npv,"sf_npv/F");
  Book("sf_tt",&sf_tt,"sf_tt/F");
  Book("sf_tt_ext",&sf_tt_ext,"sf_tt_ext/F");
  Book("sf_tt_bound",&sf_tt_bound,"sf_tt_bound/F");
  Book("sf_tt8TeV",&sf_tt8TeV,"sf_tt8TeV/F");
  Book("sf_tt8TeV_ext",&sf_tt8TeV_ext,"sf_tt8TeV_ext/F");
  Book("sf_tt8TeV_bound",&sf_tt8TeV_bound,"sf_tt8TeV_bound/F");
  Book("sf_phoPurity",&sf_phoPurity,"sf_phoPurity/F");
  Book("pfmet",&pfmet,"pfmet/F");
  Book("pfmetphi",&pfmetphi,"pfmetphi/F");
  Book("pfmetnomu",&pfmetnomu,"pfmetnomu/F");
  Book("pfmetsig",&pfmetsig,"pfmetsig/F");
  Book("puppimet",&puppimet,"puppimet/F");
  Book("puppimetphi",&puppimetphi,"puppimetphi/F");
  Book("puppimetsig",&puppimetsig,"puppimetsig/F");
  Book("calomet",&calomet,"calomet/F");
  Book("calometphi",&calometphi,"calometphi/F");
//  Book("pfcalobalance",&pfcalobalance,"pfcalobalance/F");
  Book("sumET",&sumET,"sumET/F");
//  Book("trkmet",&trkmet,"trkmet/F");
  Book("dphipfmet",&dphipfmet,"dphipfmet/F");
  Book("dphipuppimet",&dphipuppimet,"dphipuppimet/F");
  Book("trueGenBosonPt",&trueGenBosonPt,"trueGenBosonPt/F");
  Book("genBosonPt",&genBosonPt,"genBosonPt/F");
  Book("genBosonEta",&genBosonEta,"genBosonEta/F");
  Book("genBosonMass",&genBosonMass,"genBosonMass/F");
  Book("genBosonPhi",&genBosonPhi,"genBosonPhi/F");
  Book("genWPlusPt",&genWPlusPt,"genWPlusPt/F");
  Book("genWMinusPt",&genWMinusPt,"genWMinusPt/F");
  Book("genWPlusEta",&genWPlusEta,"genWPlusEta/F");
  Book("genWMinusEta",&genWMinusEta,"genWMinusEta/F");
  Book("genTopPt",&genTopPt,"genTopPt/F");
  Book("genTopIsHad",&genTopIsHad,"genTopIsHad/I");
  Book("genTopEta",&genTopEta,"genTopEta/F");
  Book("genAntiTopPt",&genAntiTopPt,"genAntiTopPt/F");
  Book("genAntiTopIsHad",&genAntiTopIsHad,"genAntiTopIsHad/I");
  Book("genAntiTopEta",&genAntiTopEta,"genAntiTopEta/F");
  Book("genTTPt",&genTTPt,"genTTPt/F");
  Book("genTTEta",&genTTEta,"genTTEta/F");
  Book("jetNBtags",&jetNBtags,"jetNBtags/I");
  Book("nLoosePhoton",&nLoosePhoton,"nLoosePhoton/I");
  Book("nTightPhoton",&nTightPhoton,"nTightPhoton/I");
  Book("loosePho1IsTight",&loosePho1IsTight,"loosePho1IsTight/I");
  Book("loosePho1Pt",&loosePho1Pt,"loosePho1Pt/F");
  Book("loosePho1Eta",&loosePho1Eta,"loosePho1Eta/F");
  Book("loosePho1Phi",&loosePho1Phi,"loosePho1Phi/F");
  Book("diLepMass",&diLepMass,"diLepMass/F");
  Book("nTau",&nTau,"nTau/I");
  Book("mT",&mT,"mT/F");
  Book("scaleUp",&scaleUp,"scaleUp/F");
  Book("scaleDown",&scaleDown,"scaleDown/F");
  Book("pdfUp",&pdfUp,"pdfUp/F");
  Book("pdfDown",&pdfDown,"pdfDown/F");
}
