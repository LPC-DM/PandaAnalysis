#include "../interface/GeneralTree.h"

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
    jetEta[iJ] = -1;
    jetPhi[iJ] = -1;
    jetE[iJ] = -1;
    jetCSV[iJ] = -1;
    jetIso[iJ] = -1;
    jetQGL[iJ] = -1;
  }

//ENDCUSTOMCONST
}

GeneralTree::~GeneralTree() {
//STARTCUSTOMDEST
//ENDCUSTOMDEST
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
    jetEta[iJ] = -99;
    jetPhi[iJ] = -99;
    jetE[iJ] = -99;
    jetCSV[iJ] = -99;
    jetIso[iJ] = -99;
    jetQGL[iJ] = -99;
  }

  for (auto iter=signal_weights.begin(); iter!=signal_weights.end(); ++iter) {
    signal_weights[iter->first] = 1; // does pair::second return a reference?
  }

//ENDCUSTOMRESET
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
    nB = 0;
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
    looseLep1IsHLTSafe = 0;
    looseLep2IsHLTSafe = 0;
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
    puppimet = -1;
    puppimetphi = -1;
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
    nIsoJet = 0;
    jet1Flav = 0;
    jet1Phi = -1;
    jet1Pt = -1;
    jet1GenPt = -1;
    jet1Eta = -1;
    jet1CSV = -1;
    jet1IsTight = 0;
    jet2Flav = 0;
    jet2Phi = -1;
    jet2Pt = -1;
    jet2GenPt = -1;
    jet2Eta = -1;
    jet2CSV = -1;
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
    looseLep1PdgId = 0;
    looseLep2PdgId = 0;
    looseLep1IsTight = 0;
    looseLep2IsTight = 0;
    looseLep1Pt = -1;
    looseLep1Eta = -1;
    looseLep1Phi = -1;
    looseLep2Pt = -1;
    looseLep2Eta = -1;
    looseLep2Phi = -1;
    diLepMass = -1;
    nTau = 0;
    mT = -1;
    hbbpt = -1;
    hbbeta = -1;
    hbbphi = -1;
    hbbm = -1;
    scaleUp = -1;
    scaleDown = -1;
    pdfUp = -1;
    pdfDown = -1;
}

void GeneralTree::WriteTree(TTree *t) {
  treePtr = t;
//STARTCUSTOMWRITE
  for (auto iter=signal_weights.begin(); iter!=signal_weights.end(); ++iter) {
    Book("rw_"+iter->first,&(signal_weights[iter->first]),"rw_"+iter->first+"/F");
  }

  Book("nJet",&nJet,"nJet/I");
  if (monohiggs) {
    Book("jetPt",jetPt,"jetPt[nJet]/F");
    Book("jetEta",jetEta,"jetEta[nJet]/F");
    Book("jetPhi",jetPhi,"jetPhi[nJet]/F");
    Book("jetE",jetE,"jetE[nJet]/F");
    Book("jetCSV",jetCSV,"jetCSV[nJet]/F");
    Book("jetIso",jetIso,"jetIso[nJet]/F");
    Book("jetQGL",jetQGL,"jetQGL[nJet]/F");
    Book("fj1sjPt",fj1sjPt,"fj1sjPt[2]/F");
    Book("fj1sjPhi",fj1sjPhi,"fj1sjPhi[2]/F");
    Book("fj1sjEta",fj1sjEta,"fj1sjEta[2]/F");
    Book("fj1sjM",fj1sjM,"fj1sjM[2]/F");
    Book("fj1sjCSV",fj1sjCSV,"fj1sjCSV[2]/F");
    Book("fj1sjQGL",fj1sjQGL,"fj1sjQGL[2]/F");
    Book("fj1Nbs",&fj1Nbs,"fj1Nbs/I");
    Book("fj1gbb",&fj1gbb,"fj1gbb/I");
    Book("hbbpt",&hbbpt,"hbbpt/F");
    Book("hbbeta",&hbbeta,"hbbeta/F");
    Book("hbbphi",&hbbphi,"hbbphi/F");
    Book("hbbm",&hbbm,"hbbm/F");
    Book("hbbjtidx",hbbjtidx,"hbbjtidx[2]/I");
  }

  if (vbf) { 
    Book("nJot",&nJot,"nJot/I");
    Book("jot1Phi",&jot1Phi,"jot1Phi/F");
    Book("jot1Pt",&jot1Pt,"jot1Pt/F");
    Book("jot1GenPt",&jot1GenPt,"jot1GenPt/F");
    Book("jot1Eta",&jot1Eta,"jot1Eta/F");
    Book("jot2Phi",&jot2Phi,"jot2Phi/F");
    Book("jot2Pt",&jot2Pt,"jot2Pt/F");
    Book("jot2GenPt",&jot2GenPt,"jot2GenPt/F");
    Book("jot2Eta",&jot2Eta,"jot2Eta/F");
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
    Book("jot1EtaUp",&jot1EtaUp,"jot1EtaUp/F");
    Book("jot1EtaDown",&jot1EtaDown,"jot1EtaDown/F");
    Book("jot1PtUp",&jot1PtUp,"jot1PtUp/F");
    Book("jot1PtDown",&jot1PtDown,"jot1PtDown/F");
    Book("jot2PtUp",&jot2PtUp,"jot2PtUp/F");
    Book("jot2PtDown",&jot2PtDown,"jot2PtDown/F");
    Book("jot12MassUp",&jot12MassUp,"jot12MassUp/F");
    Book("jot12DEtaUp",&jot12DEtaUp,"jot12DEtaUp/F");
    Book("jot12DPhiUp",&jot12DPhiUp,"jot12DPhiUp/F");
    Book("jot12MassDown",&jot12MassDown,"jot12MassDown/F");
    Book("jot12DEtaDown",&jot12DEtaDown,"jot12DEtaDown/F");
    Book("jot12DPhiDown",&jot12DPhiDown,"jot12DPhiDown/F");
    Book("jot2EtaUp",&jot2EtaUp,"jot2EtaUp/F");
    Book("jot2EtaDown",&jot2EtaDown,"jot2EtaDown/F");
    Book("jot1VBFID",&jot1VBFID,"jot1VBFID/I");
  }
  Book("scale",scale,"scale[6]/F");

  for (auto p : ecfParams) { 
    TString ecfn(makeECFString(p));
    Book("fj1"+ecfn,&(fj1ECFNs[p]),"fj1"+ecfn+"/F");
  }

  for (auto p : btagParams) {
    TString btagn(makeBTagSFString(p));
    Book(btagn,&(sf_btags[p]),btagn+"/F");
  }
//ENDCUSTOMWRITE
    Book("badECALFilter",&badECALFilter,"badECALFilter/I");
    Book("sf_qcdV_VBFTight",&sf_qcdV_VBFTight,"sf_qcdV_VBFTight/F");
    Book("sf_metTrigVBF",&sf_metTrigVBF,"sf_metTrigVBF/F");
    Book("sf_metTrigZmmVBF",&sf_metTrigZmmVBF,"sf_metTrigZmmVBF/F");
    Book("sumETRaw",&sumETRaw,"sumETRaw/F");
    Book("sf_metTrigZmm",&sf_metTrigZmm,"sf_metTrigZmm/F");
    Book("sf_qcdV_VBF",&sf_qcdV_VBF,"sf_qcdV_VBF/F");
    Book("jetNMBtags",&jetNMBtags,"jetNMBtags/I");
    Book("pfmetRaw",&pfmetRaw,"pfmetRaw/F");
    Book("nB",&nB,"nB/I");
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
    Book("isGS",&isGS,"isGS/I");
    Book("fj1SubMaxCSV",&fj1SubMaxCSV,"fj1SubMaxCSV/F");
    Book("looseLep1IsHLTSafe",&looseLep1IsHLTSafe,"looseLep1IsHLTSafe/I");
    Book("looseLep2IsHLTSafe",&looseLep2IsHLTSafe,"looseLep2IsHLTSafe/I");
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
    Book("sf_lepID",&sf_lepID,"sf_lepID/F");
    Book("sf_lepIso",&sf_lepIso,"sf_lepIso/F");
    Book("sf_lepTrack",&sf_lepTrack,"sf_lepTrack/F");
    Book("sf_pho",&sf_pho,"sf_pho/F");
    Book("sf_eleTrig",&sf_eleTrig,"sf_eleTrig/F");
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
    Book("puppimet",&puppimet,"puppimet/F");
    Book("puppimetphi",&puppimetphi,"puppimetphi/F");
    Book("calomet",&calomet,"calomet/F");
    Book("calometphi",&calometphi,"calometphi/F");
    Book("pfcalobalance",&pfcalobalance,"pfcalobalance/F");
    Book("sumET",&sumET,"sumET/F");
    Book("trkmet",&trkmet,"trkmet/F");
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
    Book("pfUWphi",&pfUWphi,"pfUWphi/F");
    Book("pfUZmag",&pfUZmag,"pfUZmag/F");
    Book("pfUZphi",&pfUZphi,"pfUZphi/F");
    Book("pfUWWmag",&pfUWWmag,"pfUWWmag/F");
    Book("pfUWWphi",&pfUWWphi,"pfUWWphi/F");
    Book("pfUAmag",&pfUAmag,"pfUAmag/F");
    Book("pfUAphi",&pfUAphi,"pfUAphi/F");
    Book("pfUperp",&pfUperp,"pfUperp/F");
    Book("pfUpara",&pfUpara,"pfUpara/F");
    Book("pfUmag",&pfUmag,"pfUmag/F");
    Book("pfUphi",&pfUphi,"pfUphi/F");
    Book("dphipfmet",&dphipfmet,"dphipfmet/F");
    Book("dphipuppimet",&dphipuppimet,"dphipuppimet/F");
    Book("dphipuppiUW",&dphipuppiUW,"dphipuppiUW/F");
    Book("dphipuppiUZ",&dphipuppiUZ,"dphipuppiUZ/F");
    Book("dphipuppiUA",&dphipuppiUA,"dphipuppiUA/F");
    Book("dphipfUW",&dphipfUW,"dphipfUW/F");
    Book("dphipfUWW",&dphipfUWW,"dphipfUWW/F");
    Book("dphipfUZ",&dphipfUZ,"dphipfUZ/F");
    Book("dphipfUA",&dphipfUA,"dphipfUA/F");
    Book("dphipuppiU",&dphipuppiU,"dphipuppiU/F");
    Book("dphipfU",&dphipfU,"dphipfU/F");
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
    Book("nIsoJet",&nIsoJet,"nIsoJet/I");
    Book("jet1Flav",&jet1Flav,"jet1Flav/I");
    Book("jet1Phi",&jet1Phi,"jet1Phi/F");
    Book("jet1Pt",&jet1Pt,"jet1Pt/F");
    Book("jet1GenPt",&jet1GenPt,"jet1GenPt/F");
    Book("jet1Eta",&jet1Eta,"jet1Eta/F");
    Book("jet1CSV",&jet1CSV,"jet1CSV/F");
    Book("jet1IsTight",&jet1IsTight,"jet1IsTight/I");
    Book("jet2Flav",&jet2Flav,"jet2Flav/I");
    Book("jet2Phi",&jet2Phi,"jet2Phi/F");
    Book("jet2Pt",&jet2Pt,"jet2Pt/F");
    Book("jet2GenPt",&jet2GenPt,"jet2GenPt/F");
    Book("jet2Eta",&jet2Eta,"jet2Eta/F");
    Book("jet2CSV",&jet2CSV,"jet2CSV/F");
    Book("jet3Flav",&jet3Flav,"jet3Flav/I");
    Book("jet3Phi",&jet3Phi,"jet3Phi/F");
    Book("jet3Pt",&jet3Pt,"jet3Pt/F");
    Book("jet3GenPt",&jet3GenPt,"jet3GenPt/F");
    Book("jet3Eta",&jet3Eta,"jet3Eta/F");
    Book("jet3CSV",&jet3CSV,"jet3CSV/F");
    Book("isojet1Pt",&isojet1Pt,"isojet1Pt/F");
    Book("isojet1CSV",&isojet1CSV,"isojet1CSV/F");
    Book("isojet1Flav",&isojet1Flav,"isojet1Flav/I");
    Book("isojet2Pt",&isojet2Pt,"isojet2Pt/F");
    Book("isojet2CSV",&isojet2CSV,"isojet2CSV/F");
    Book("isojet2Flav",&isojet2Flav,"isojet2Flav/I");
    Book("jetNBtags",&jetNBtags,"jetNBtags/I");
    Book("isojetNBtags",&isojetNBtags,"isojetNBtags/I");
    Book("nFatjet",&nFatjet,"nFatjet/I");
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
    Book("fj1Pt",&fj1Pt,"fj1Pt/F");
    Book("fj1PtScaleUp",&fj1PtScaleUp,"fj1PtScaleUp/F");
    Book("fj1PtScaleDown",&fj1PtScaleDown,"fj1PtScaleDown/F");
    Book("fj1PtSmeared",&fj1PtSmeared,"fj1PtSmeared/F");
    Book("fj1PtSmearedUp",&fj1PtSmearedUp,"fj1PtSmearedUp/F");
    Book("fj1PtSmearedDown",&fj1PtSmearedDown,"fj1PtSmearedDown/F");
    Book("fj1Phi",&fj1Phi,"fj1Phi/F");
    Book("fj1Eta",&fj1Eta,"fj1Eta/F");
    Book("fj1M",&fj1M,"fj1M/F");
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
    Book("fj1HTTMass",&fj1HTTMass,"fj1HTTMass/F");
    Book("fj1HTTFRec",&fj1HTTFRec,"fj1HTTFRec/F");
    Book("fj1IsClean",&fj1IsClean,"fj1IsClean/I");
    Book("fj1NConst",&fj1NConst,"fj1NConst/I");
    Book("fj1NSDConst",&fj1NSDConst,"fj1NSDConst/I");
    Book("fj1EFrac100",&fj1EFrac100,"fj1EFrac100/F");
    Book("fj1SDEFrac100",&fj1SDEFrac100,"fj1SDEFrac100/F");
    Book("nHF",&nHF,"nHF/I");
    Book("nLoosePhoton",&nLoosePhoton,"nLoosePhoton/I");
    Book("nTightPhoton",&nTightPhoton,"nTightPhoton/I");
    Book("loosePho1IsTight",&loosePho1IsTight,"loosePho1IsTight/I");
    Book("loosePho1Pt",&loosePho1Pt,"loosePho1Pt/F");
    Book("loosePho1Eta",&loosePho1Eta,"loosePho1Eta/F");
    Book("loosePho1Phi",&loosePho1Phi,"loosePho1Phi/F");
    Book("nLooseLep",&nLooseLep,"nLooseLep/I");
    Book("nLooseElectron",&nLooseElectron,"nLooseElectron/I");
    Book("nLooseMuon",&nLooseMuon,"nLooseMuon/I");
    Book("nTightLep",&nTightLep,"nTightLep/I");
    Book("nTightElectron",&nTightElectron,"nTightElectron/I");
    Book("nTightMuon",&nTightMuon,"nTightMuon/I");
    Book("looseLep1PdgId",&looseLep1PdgId,"looseLep1PdgId/I");
    Book("looseLep2PdgId",&looseLep2PdgId,"looseLep2PdgId/I");
    Book("looseLep1IsTight",&looseLep1IsTight,"looseLep1IsTight/I");
    Book("looseLep2IsTight",&looseLep2IsTight,"looseLep2IsTight/I");
    Book("looseLep1Pt",&looseLep1Pt,"looseLep1Pt/F");
    Book("looseLep1Eta",&looseLep1Eta,"looseLep1Eta/F");
    Book("looseLep1Phi",&looseLep1Phi,"looseLep1Phi/F");
    Book("looseLep2Pt",&looseLep2Pt,"looseLep2Pt/F");
    Book("looseLep2Eta",&looseLep2Eta,"looseLep2Eta/F");
    Book("looseLep2Phi",&looseLep2Phi,"looseLep2Phi/F");
    Book("diLepMass",&diLepMass,"diLepMass/F");
    Book("nTau",&nTau,"nTau/I");
    Book("mT",&mT,"mT/F");
    Book("scaleUp",&scaleUp,"scaleUp/F");
    Book("scaleDown",&scaleDown,"scaleDown/F");
    Book("pdfUp",&pdfUp,"pdfUp/F");
    Book("pdfDown",&pdfDown,"pdfDown/F");
}

