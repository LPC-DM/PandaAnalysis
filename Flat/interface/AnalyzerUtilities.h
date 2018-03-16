#ifndef ANALYZERUTILS_H
#define ANALYZERUTILS_H

// PandaProd Objects
#include "PandaTree/Objects/interface/EventAnalysis.h"

// PANDACore
#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/DataTools.h"
#include "PandaCore/Tools/interface/JERReader.h"

// fastjet
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/MeasureDefinition.hh"

// root
#include "TRotation.h"

////////////////////////////////////////////////////////////////////////////////////

class JetTree {
  public:
    JetTree(fastjet::PseudoJet& root): _root(root) { }
    ~JetTree() { }
    std::vector<int> GetTerminals() {
      std::vector<int> v;
      _root.GetTerminals(v);
      return v;
    }

    struct compare : public std::unary_function<fastjet::PseudoJet, bool> {
      explicit compare(int x_) : _x(x_) {}
      bool operator() (const fastjet::PseudoJet& y_) { return _x == y_.user_index(); }
      int _x;
    };

  private:
    class Node {
      public:
        Node(fastjet::PseudoJet& pj_);
        ~Node() { delete l; delete r; }
        void GetTerminals(std::vector<int>&);
        fastjet::PseudoJet _pj;
        Node *l=0, *r=0;   
    };

    Node _root;
};

////////////////////////////////////////////////////////////////////////////////////

class ParticleGridder {
  public:
    ParticleGridder(unsigned etaN, unsigned phiN, float etaMax=5);
    ~ParticleGridder() { clear(); delete hEta_; delete hPhi_; }
    void clear();
    void add(panda::Particle& p);
    std::vector<TLorentzVector>& get();
  private:
    float etaMax_, phiMax_;
    TH1F *hEta_=0, *hPhi_=0;
    std::vector<std::vector<std::vector<TLorentzVector*>>> collections_;
    std::vector<TLorentzVector> particles_;
    std::vector<std::pair<int,int>> nonEmpty_;
    std::vector<TLorentzVector> gridded_; 
};

////////////////////////////////////////////////////////////////////////////////////

class JetRotation {
  public:
    JetRotation(float x1, float y1, float z1,
                float x2, float y2, float z2);
    void Rotate(float& x, float& y, float& z);
  private:
    TRotation r_toz;
    TRotation r_inxy;
};

////////////////////////////////////////////////////////////////////////////////////

typedef std::vector<fastjet::PseudoJet> VPseudoJet;
VPseudoJet ConvertPFCands(std::vector<const panda::PFCand*> &incoll, bool puppi, double minPt=0.001);
VPseudoJet ConvertPFCands(panda::RefVector<panda::PFCand> &incoll, bool puppi, double minPt=0.001);
VPseudoJet ConvertPFCands(panda::PFCandCollection &incoll, bool puppi, double minPt=0.001);

////////////////////////////////////////////////////////////////////////////////////

enum ProcessType { 
    kNoProcess,
    kZ,
    kW,
    kA,
    kZEWK,
    kWEWK,
    kTT,
    kTop, // used for non-ttbar top
    kV, // used for non V+jets W or Z
    kH,
    kSignal,
};

class Analysis {
public:
  Analysis(TString name_ = "") { name = name_; }
  ~Analysis() {}
  TString name;
  ProcessType processType=kNoProcess;
  bool ak8 = false;
  bool applyMCTriggers = false;
  bool bjetRegression = false;
  bool btagSFs = true;
  bool btagWeights = false;
  bool complicatedLeptons = false;
  bool complicatedPhotons = false;
  bool deep = false;
  bool deepAntiKtSort = false;
  bool deepGen = false;
  bool deepGenGrid = false;
  bool deepKtSort = false;
  bool deepSVs = false;
  bool deepTracks = false;
  bool fatjet = true;
  bool firstGen = true;
  bool genOnly = false;
  bool hbb = false;
  bool hfCounting = false;
  bool jetFlavorPartons = true;
  bool jetFlavorJets = false;
 // bool monoh = false;
  bool boosted = false;
  bool puppi_jets = true;
  bool recluster = false;
  bool reclusterGen = false;
  bool recoil = true;
  bool rerunJES = false;
  bool useCMVA = false;
  bool varyJES = false;
  bool vbf = false;
};

////////////////////////////////////////////////////////////////////////////////////

class LumiRange {
public:
    LumiRange(int l0_,int l1_):
        l0(l0_),
        l1(l1_)
     { }
    ~LumiRange() {}
    bool Contains(int l) {
        return l0<=l && l<=l1;
    }
private:
    int l0, l1;
};

////////////////////////////////////////////////////////////////////////////////////

class TriggerHandler {  
public:
  TriggerHandler() {};
  ~TriggerHandler() {};
  void addTriggers(std::vector<TString> paths_) { 
    for (auto &path : paths_) {
      paths.push_back(path); 
      indices.push_back(-1); 
    }
  }
  void registerTrigger(unsigned my_idx, int panda_idx) { indices[my_idx] = panda_idx; }
  std::vector<int> indices;
  std::vector<TString> paths;
};

////////////////////////////////////////////////////////////////////////////////////

template <typename T>
class TCorr {
public:
  TCorr(T*) {} 
  virtual ~TCorr() {}
  virtual double Eval(double x)=0; 
protected:
  T *h=0;
};


class TF1Corr : public TCorr<TF1> {
public:
  TF1Corr(TF1 *f_):
    TCorr(f_) 
  {
//    f_->SetDirectory(0);
    h = f_;
  }
  ~TF1Corr() {} 
  double Eval(double x) {
    return h->Eval(x);
  }

  TF1 *GetFunc() { return h; }
};

template <typename T>
class THCorr : public TCorr<T> {
public:
  // wrapper around TH* to do corrections
  THCorr(T *h_):
    TCorr<T>(h_)
  {
    this->h = h_;
    dim = this->h->GetDimension();
    TAxis *thurn = this->h->GetXaxis(); 
    lo1 = thurn->GetBinCenter(1);
    hi1 = thurn->GetBinCenter(thurn->GetNbins());
    if (dim>1) {
      TAxis *taxis = this->h->GetYaxis();
      lo2 = taxis->GetBinCenter(1);
      hi2 = taxis->GetBinCenter(taxis->GetNbins());
    }
  }
  ~THCorr() {} // does not own histogram!
  double Eval(double x) {
    if (dim!=1) {
      PError("THCorr1::Eval",
        TString::Format("Trying to access a non-1D histogram (%s)!",this->h->GetName()));
      return -1;
    }
    return getVal(this->h,bound(x,lo1,hi1));
  }

  double Eval(double x, double y) {
    if (dim!=2) {
      PError("THCorr1::Error",
       TString::Format("Trying to access a non-2D histogram (%s)!",this->h->GetName()));
      return -1;
    }
    return getVal(this->h,bound(x,lo1,hi1),bound(y,lo2,hi2));
  }

  double Error(double x) {
    if (dim!=1) {
      PError("THCorr1::Error",
        TString::Format("Trying to access a non-1D histogram (%s)!",this->h->GetName()));
      return -1;
    }
    return getError(this->h,bound(x,lo1,hi1));
  }

  double Error(double x, double y) {
    if (dim!=2) {
      PError("THCorr1::Eval",
       TString::Format("Trying to access a non-2D histogram (%s)!",this->h->GetName()));
      return -1;
    }
    return getError(this->h,bound(x,lo1,hi1),bound(y,lo2,hi2));
  }

  T *GetHist() { return this->h; }

private:
  int dim;
  double lo1, lo2, hi1, hi2;
};

typedef THCorr<TH1D> THCorr1;
typedef THCorr<TH2D> THCorr2;

////////////////////////////////////////////////////////////////////////////////////

bool ElectronIP(double eta, double dxy, double dz); 
bool MuonIP(double dxy, double dz); 
double TTNLOToNNLO(double pt);
bool IsMatched(std::vector<panda::Particle*>*objects,
               double deltaR2, double eta, double phi);

////////////////////////////////////////////////////////////////////////////////////

#endif
