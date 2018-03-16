#ifndef PandaAnalyzer_h
#define PandaAnalyzer_h

// STL
#include "vector"
#include <unordered_set>
#include "map"
#include <string>
#include <cmath>

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

#include "AnalyzerUtilities.h"
#include "GeneralTree.h"

// btag
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "BTagCalibrationStandalone.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// Utils
#include "PandaAnalysis/Utilities/interface/RoccoR.h"
#include "PandaAnalysis/Utilities/interface/CSVHelper.h"
#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"

// TMVA
#include "TMVA/Reader.h"

/////////////////////////////////////////////////////////////////////////////
// some misc definitions

#define NMAXPF 100
#define NMAXSV 10


/////////////////////////////////////////////////////////////////////////////
// PandaAnalyzer definition
class PandaAnalyzer {
public :
    // configuration enums
    enum PreselectionBit {
     kMonotop    =(1<<0),
     kMonohiggs  =(1<<2),
     kMonojet    =(1<<3),
     kPassTrig   =(1<<4),
     kVBF        =(1<<5),
     kRecoil     =(1<<6),
     kFatjet     =(1<<7),
     kFatjet450  =(1<<8),
     kRecoil50   =(1<<9),
     kGenBosonPt =(1<<10),
     kGenFatJet  =(1<<11),
     kVHBB       =(1<<12),
     kLepton     =(1<<13),
     kLeptonFake =(1<<14)
    };
    
    enum LepSelectionBit {
     kLoose     =(1<<0),
     kFake      =(1<<1),
     kMedium    =(1<<2),
     kTight     =(1<<3),
     kDxyz      =(1<<4),
     kEleMvaWP90=(1<<5),
     kEleMvaWP80=(1<<6)
    };

    enum PhoSelectionBit {
     pMedium    =(1<<0),
     pTight     =(1<<1),
     pHighPt    =(1<<2),
     pCsafeVeto =(1<<3),
     pPixelVeto =(1<<4),
     pTrkVeto   =(1<<5)
    };

    enum TriggerBits {
        kMETTrig       = 0,
        kSingleEleTrig,
        kSinglePhoTrig,
        kSingleMuTrig,
        kDoubleMuTrig,
        kDoubleEleTrig,
        kEMuTrig,
        kJetHTTrig,
        kMuFakeTrig,
        kEleFakeTrig,
        kNTrig
    };


    //////////////////////////////////////////////////////////////////////////////////////

    PandaAnalyzer(int debug_=0);
    ~PandaAnalyzer();
    int Init(TTree *tree, TH1D *hweights, TTree *weightNames=0);
    void SetOutputFile(TString fOutName);
    void ResetBranches();
    void Run();
    void Terminate();
    void SetDataDir(const char *s);
    void SetPreselectionBit(PreselectionBit b,bool on=true) {
        if (on) 
            preselBits |= b;
        else 
            preselBits &= ~b;
    }
    void AddGoodLumiRange(int run, int l0, int l1);

    // public configuration
    void SetAnalysis(Analysis *a) { analysis = a; }
    bool isData=false;              // to do gen matching, etc
    int firstEvent=-1;
    int lastEvent=-1;               // max events to process; -1=>all

private:
    enum CorrectionType { //!< enum listing relevant corrections applied to MC
        cNPV=0,       //!< npv weight
        cPU,          //!< true pu weight
        cPUUp,        //!< true pu weight
        cPUDown,      //!< true pu weight
        cEleVeto,     //!< monojet SF, Veto ID for e
        cEleLoose,    //!< monojet SF, Tight ID for e
        cEleMedium,   //!< monojet SF, Tight ID for e
        cEleTight,    //!< monojet SF, Tight ID for e
        cEleMvaWP90,
        cEleMvaWP80,
        cEleReco,     //!< monojet SF, tracking for e
        cWmHEwkCorr,     //!< W(l-V)H Ewk Corr weight  
        cWmHEwkCorrUp,   //!< W(l-V)H Ewk Corr weight Up  
        cWmHEwkCorrDown, //!< W(l-V)H Ewk Corr weight Down  
        cWpHEwkCorr,     //!< W(l+v)H Ewk Corr weight  
        cWpHEwkCorrUp,   //!< W(l+v)H Ewk Corr weight Up  
        cWpHEwkCorrDown, //!< W(l+v)H Ewk Corr weight Down  
        cZnnHEwkCorr,     //!< Z(vv)H Ewk Corr weight  
        cZnnHEwkCorrUp,   //!< Z(vv)H Ewk Corr weight Up  
        cZnnHEwkCorrDown, //!< Z(vv)H Ewk Corr weight Down  
        cZllHEwkCorr,     //!< Z(ll)H Ewk Corr weight  
        cZllHEwkCorrUp,   //!< Z(ll)H Ewk Corr weight Up  
        cZllHEwkCorrDown, //!< Z(ll)H Ewk Corr weight Down  
        cWZEwkCorr,
        cqqZZQcdCorr,
        cMuLooseID,   //!< MUO POG SF, Loose ID for mu 
        cMuMediumID,  //!< MUO POG SF, Tight ID for mu 
        cMuTightID,   //!< MUO POG SF, Tight ID for mu 
        cMuLooseIso,  //!< MUO POG SF, Loose Iso for mu 
        cMuMediumIso, //!< MUO POG SF, Loose Iso for mu 
        cMuTightIso,  //!< MUO POG SF, Tight Iso for mu 
        cMuReco,      //!< MUO POG SF, tracking for mu
        cPho,         //!< EGM POG SF, contains ID for gamma
        cTrigMET,     //!< MET trigger eff        
        cTrigMETZmm,  //!< Zmumu MET trigger eff
        cTrigEle,     //!< Ele trigger eff        
        cTrigMu,      //!< Mu trigger eff        
        cTrigPho,     //!< Pho trigger eff        
        cZNLO,        //!< NLO weights for QCD Z,W,A,A+2j
        cWNLO,
        cANLO,
        cANLO2j,
        cZEWK,        //!< EWK weights for QCD Z,W,A,A+2j
        cWEWK,
        cAEWK,
        cVBF_ZNLO,    //!< NLO weights for QCD Z,W in VBF phase space
        cVBF_ZllNLO,  
        cVBF_WNLO,
        cVBFTight_ZNLO,    //!< NLO weights for QCD Z,W in tight VBF phase space
        cVBFTight_ZllNLO,  
        cVBFTight_WNLO,
        cVBF_EWKZ,    //!< k-factors for EWK Z,W in VBF phase space
        cVBF_EWKW,
        cVBF_TrigMET, //!< MET trigger eff as a f'n of mjj/met 
        cVBF_TrigMETZmm,
        cBadECALJets,  //!< bad ECAL clusters to filter jets
        cN
    };

    enum BTagType {
        bJetL=0,
        bSubJetL,
        bJetM,
        bSubJetM,
        bN
    };

    class btagcand {
        public:
            btagcand(unsigned int i, int f,double e,double cent,double up,double down) {
                idx = i;
                flav = f;
                eff = e;
                sf = cent;
                sfup = up;
                sfdown = down;
            }
            ~btagcand() { }
            int flav, idx;
            double eff, sf, sfup, sfdown;
    };

    struct GenJetInfo {
    public:
      float pt=-1, eta=-1, phi=-1, m=-1;
      float msd=-1;
      float tau3=-1, tau2=-1, tau1=-1;
      float tau3sd=-1, tau2sd=-1, tau1sd=-1;
      int nprongs=-1;
      float partonpt=-1, partonm=-1;
      std::vector<std::vector<float>> particles;
      std::vector<std::vector<std::vector<float>>> ecfs; // uh
      void reset() {
        pt=-1; eta=-1; phi=-1; m=-1;
        msd=-1;
        tau3=-1; tau2=-1; tau1=-1;
        tau3sd=-1; tau2sd=-1; tau1sd=-1;
        nprongs=-1;
        partonpt=-1; partonm=-1;
        for (auto& v : particles) {
          for (auto& vv : v) {
            vv = 0;
          }
        }
        for (auto& v : ecfs) {
          for (auto& vv : v) {
            for (auto& vvv : vv) {
              vvv = -1;
            }
          }
        }
      }

    };

    struct JetHistory {
      int user_idx;
      int child_idx;
    };

    //////////////////////////////////////////////////////////////////////////////////////

    bool PassGoodLumis(int run, int lumi);
    bool PassPreselection();
    void OpenCorrection(CorrectionType,TString,TString,int);
    double GetCorr(CorrectionType ct,double x, double y=0);
    double GetError(CorrectionType ct,double x, double y=0);
    void RegisterTriggers(); 
    void GetMETSignificance(); 

    // these are functions used for analysis-specific tasks inside Run.
    // ideally the return type is void (e.g. they are stateful functions),
    // but that is not always possible (e.g. RecoilPresel)
    void CalcBJetSFs(BTagType bt, int flavor, double eta, double pt, 
                     double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
    void ComplicatedLeptons();
    void ComplicatedPhotons();
    void EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
                    GeneralTree::BTagShift shift,GeneralTree::BTagJet jettype, bool do2=false);
    void IncrementAuxFile(bool close=false);
    void IncrementGenAuxFile(bool close=false);
    void FatjetBasics();
    void FatjetMatching();
    void FatjetPartons();
    void FatjetRecluster();
    void FillPFTree();
    void GenFatJet();
    void GenJetsNu();
    void GenStudyEWK();
    float GetMSDCorr(float, float); 
    void HeavyFlavorCounting();
    void IsoJet(panda::Jet&);
    void JetBRegressionInfo(panda::Jet&);
    void JetBasics();
    void JetBtagSFs();
    void JetCMVAWeights();
    void JetHbbBasics(panda::Jet&);
    void JetHbbReco();
    void JetHbbSoftActivity();
    void JetVBFBasics(panda::Jet&);
    void JetVBFSystem();
    void JetVaryJES(panda::Jet&);
    void LeptonSFs();
    bool PFChargedPhotonMatch(const panda::Photon& photon);
    void PhotonSFs();
    void QCDUncs();
    void Recoil();
    bool RecoilPresel();
    void SaveGenLeptons();
    void SetupJES();
    void SignalInfo();
    void SignalReweights();
    void SimpleLeptons();
    void SimplePhotons();
    void Taus();
    void TopPTReweight();
    void TriggerEffs();
    void VJetsReweight();
    double WeightEWKCorr(float pt, int type);
    double WeightZHEWKCorr(float baseCorr);

    // templated functions
    template <typename T> void CountGenPartons(std::unordered_set<const T*>&, const panda::Collection<T>&);
    template <typename T> void FillGenTree(panda::Collection<T>& genParticles);
    template <typename T> void MatchGenJets(T& genJets);

    //////////////////////////////////////////////////////////////////////////////////////

    int DEBUG = 0; //!< debug verbosity level
    Analysis *analysis = 0; //!< configure what to run
    TimeReporter *tr = 0; //!< profile time usage
    float FATJETMATCHDR2 = 2.25;

    //////////////////////////////////////////////////////////////////////////////////////

    // stuff for matching objects
    std::map<panda::GenParticle const*,float> genObjects; 
        //!< particles we want to match the jets to, and the 'size' of the daughters
    panda::GenParticle const* MatchToGen(double eta, double phi, double r2, int pdgid=0);   
        //!< private function to match a jet; returns NULL if not found
    std::map<int,std::vector<LumiRange>> goodLumis;
    std::vector<panda::Particle*> matchPhos, matchEles, matchLeps;
    std::map<int, int> pdgToQ; 
    
    // fastjet reclustering
    fastjet::JetDefinition *jetDef=0;
    fastjet::JetDefinition *jetDefKt=0;
    fastjet::contrib::SoftDrop *softDrop=0;
    fastjet::contrib::Njettiness *tauN=0;
    fastjet::AreaDefinition *areaDef=0;
    fastjet::GhostedAreaSpec *activeArea=0;
    fastjet::JetDefinition *jetDefGen=0;
    fastjet::JetDefinition *softTrackJetDefinition=0;

    //////////////////////////////////////////////////////////////////////////////////////

    // CMSSW-provided utilities
    BTagCalibration *btagCalib=0;
    BTagCalibration *sj_btagCalib=0;
    std::vector<BTagCalibrationReader*> btagReaders = std::vector<BTagCalibrationReader*>(bN,0); 
      //!< maps BTagType to a reader 

    Binner btagpt = Binner({});
    Binner btageta = Binner({});
    std::vector<std::vector<double>> lfeff, ceff, beff;
    TMVA::Reader *bjetregReader=0; 

    std::map<TString,JetCorrectionUncertainty*> ak8UncReader; //!< calculate JES unc on the fly
    JERReader *ak8JERReader{0}; //!< fatjet jet energy resolution reader
    std::map<TString,JetCorrectionUncertainty*> ak4UncReader; //!< calculate JES unc on the fly
    std::map<TString,FactorizedJetCorrector*> ak4ScaleReader; //!< calculate JES on the fly
    JERReader *ak4JERReader{0}; //!< fatjet jet energy resolution reader
    JetCorrectionUncertainty *uncReader=0;           
    JetCorrectionUncertainty *uncReaderAK4=0;        
    FactorizedJetCorrector *scaleReaderAK4=0;        

    EraHandler eras = EraHandler(2016); //!< determining data-taking era, to be used for era-dependent JEC
    ParticleGridder *grid = 0;
    pandaecf::ECFNManager *ecfnMan = 0;

    //////////////////////////////////////////////////////////////////////////////////////

    // files and histograms containing weights
    std::vector<TFile*>   fCorrs  = std::vector<TFile*>  (cN,0); //!< files containing corrections
    std::vector<THCorr1*> h1Corrs = std::vector<THCorr1*>(cN,0); //!< histograms for binned corrections
    std::vector<THCorr2*> h2Corrs = std::vector<THCorr2*>(cN,0); //!< histograms for binned corrections
    std::vector<TF1Corr*> f1Corrs = std::vector<TF1Corr*>(cN,0); //!< TF1s for continuous corrections
    TFile *MSDcorr=0;
    TF1* puppisd_corrGEN=0;
    TF1* puppisd_corrRECO_cen=0;
    TF1* puppisd_corrRECO_for=0;
    RoccoR *rochesterCorrection=0;
    TRandom3 rng;
    CSVHelper *csvReweighter=0, *cmvaReweighter=0;

    //////////////////////////////////////////////////////////////////////////////////////

    // IO for the analyzer
    TString fOutPath;
    TFile *fOut=0;     // output file is owned by PandaAnalyzer
    TTree *tOut=0;
    GeneralTree *gt=0; // essentially a wrapper around tOut
    TString auxFilePath="";
    unsigned auxCounter=0;
    TFile *fAux=0; // auxillary file
    TTree *tAux=0;
    TH1F *hDTotalMCWeight=0;
    TTree *tIn=0;    // input tree to read
    unsigned int preselBits=0;
    panda::EventAnalysis event;

    //////////////////////////////////////////////////////////////////////////////////////

    // configuration read from output tree
    std::vector<int> ibetas;
    std::vector<int> Ns; 
    std::vector<int> orders;

    //////////////////////////////////////////////////////////////////////////////////////

    // stuff that gets passed between modules
    //
    // NB: ensure that any global vectors/maps that are per-event
    // are reset properly in ResetBranches(), or you can really
    // mess up behavior
    std::vector<TString> wIDs;
    std::vector<TriggerHandler> triggerHandlers = std::vector<TriggerHandler>(kNTrig);

    std::vector<panda::Lepton*> looseLeps, tightLeps;
    std::vector<panda::Photon*> loosePhos;
    int looseLep1PdgId, looseLep2PdgId, looseLep3PdgId, looseLep4PdgId;

    TLorentzVector vPFMET, vPuppiMET;
    TVector2 vMETNoMu;
    TLorentzVector vpfUW, vpfUZ, vpfUWW, vpfUA, vpfU;
    TLorentzVector vpuppiUW, vpuppiUZ, vpuppiUWW, vpuppiUA, vpuppiU;
    
    panda::FatJet *fj1 = 0;
    panda::FatJetCollection *fatjets = 0;
    std::vector<panda::Jet*> cleanedJets, isoJets, centralJets, bCandJets;
    TLorentzVector vJet, vBarrelJets;
    panda::JetCollection *jets = 0;
    panda::Jet *jot1 = 0, *jot2 = 0;
    panda::Jet *jotUp1 = 0, *jotUp2 = 0;
    panda::Jet *jotDown1 = 0, *jotDown2 = 0;
    panda::Jet *jetUp1 = 0, *jetUp2 = 0;
    panda::Jet *jetDown1 = 0, *jetDown2 = 0;
    float jetPtThreshold=30;
    float bJetPtThreshold=30;
    std::map<panda::Jet*,int> bCandJetGenFlavor;
    std::map<panda::Jet*,float> bCandJetGenPt;

    std::vector<panda::GenJet> genJetsNu;
    float genBosonPtMin, genBosonPtMax;

    float *bjetreg_vars = 0;

    std::vector<std::vector<float>> pfInfo;
    std::vector<std::vector<float>> svInfo; 
    float fjmsd, fjpt, fjrawpt, fjeta, fjphi;
    int NPFPROPS = 9, NSVPROPS = 13;

    GenJetInfo genJetInfo;
    int NGENPROPS = 8; 
    
    float minSoftTrackPt=0.3; // 300 MeV
};


/** templated functions **/

#include "TemplatedPandaAnalyzer.h"

#endif

