#include "PandaAnalysis/Flat/interface/AnalyzerUtilities.h"
#include <cassert>

using namespace fastjet;
using namespace std;

////////////////////////////////////////////////////////////////////////////////////

JetTree::Node::Node(PseudoJet& pj_):
  _pj(pj_)
{
  PseudoJet dau1, dau2;
  if (_pj.has_parents(dau1, dau2)) {
    l = new Node(dau1);
    r = new Node(dau2);
  }
}

void JetTree::Node::GetTerminals(vector<int>& terminals_)
{
  if (l && r) {
    l->GetTerminals(terminals_);
    r->GetTerminals(terminals_);
  } else {
    terminals_.push_back(_pj.user_index());
  }
}

////////////////////////////////////////////////////////////////////////////////////

ParticleGridder::ParticleGridder(unsigned etaN, unsigned phiN, float etaMax) 
{
  float phiMax = TMath::Pi();
  hEta_ = new TH1F("eta","eta",etaN*2,-etaMax,etaMax);
  hPhi_ = new TH1F("phi","phi",phiN*2,-phiMax,phiMax);
  etaMax_ = etaMax; phiMax_ = phiMax - 0.000001;
  collections_.resize(etaN*2);
  for (auto &v : collections_)
    v.resize(phiN*2);
}

void ParticleGridder::clear() 
{
  particles_.clear();
  for (auto &v : collections_) {
    for (auto &vv : v) {
      vv.clear();
    }
  }
  gridded_.clear();
  nonEmpty_.clear();
}

void ParticleGridder::add(panda::Particle& p) 
{
  particles_.push_back(p.p4());
}


std::vector<TLorentzVector>& ParticleGridder::get() 
{
  for (auto &p : particles_) {
    float eta = p.Eta();
    float phi = p.Phi();
    if (fabs(eta) >= etaMax_)
      continue;
    phi = bound(phi, -phiMax_, phiMax_);
    int etabin = hEta_->FindBin(eta)-1;
    int phibin = hPhi_->FindBin(phi)-1;
    collections_.at(etabin).at(phibin).push_back(&p);
    std::pair<int,int>bin(etabin,phibin);
    if (std::find(nonEmpty_.begin(),nonEmpty_.end(),bin) == nonEmpty_.end())
      nonEmpty_.push_back(bin);
  }

  // grid is filled
  gridded_.clear();
  for (auto &bin : nonEmpty_) {
    int iEta = bin.first, iPhi = bin.second;
    std::vector<TLorentzVector*> inBin = collections_.at(iEta).at(iPhi);
    if (inBin.size() > 0) {
      TLorentzVector vSum;
      for (auto *p : inBin) {
        vSum += *p;
      }
      float eta = hEta_->GetBinCenter(iEta);
      float phi = hPhi_->GetBinCenter(iPhi);
      vSum.SetPtEtaPhiM(vSum.Pt(), eta, phi, vSum.M());
      gridded_.push_back(vSum);
    }
  }
  return gridded_;
}

////////////////////////////////////////////////////////////////////////////////////

JetRotation::JetRotation(float x1, float y1, float z1, float x2, float y2, float z2)
{
  TVector3 axis1(x1, y1, z1); // this axis gets rotated onto the z-axis
  TVector3 axis2(x2, y2, z2); // this axis will get rotated into the x-z plane 
  TVector3 axisz(0, 0, 1);
  TVector3 axisx(1, 0, 0);

  r_toz.Rotate(axis1.Angle(axisz), axis1.Cross(axisz)); 
  assert((r_toz*axis1).Angle(axisz) < 0.0001); // allow some rounding

  axis2 = r_toz * axis2;      // first rotate it as before 
  axis2.SetZ(0);              // zero-out the z-component 
  r_inxy.Rotate(axis2.Angle(axisx), axis2.Cross(axisx));
  assert((r_inxy*axis2).Angle(axisx) < 0.0001);
}

void JetRotation::Rotate(float& x, float& y, float& z) 
{
  TVector3 v(x, y, z);
  v = r_toz * v;
  v = r_inxy * v;
  x = v.x(); y = v.y(); z = v.z();
}

////////////////////////////////////////////////////////////////////////////////////

VPseudoJet ConvertPFCands(std::vector<const panda::PFCand*> &incoll, bool puppi, double minPt) 
{
  VPseudoJet vpj;
  vpj.reserve(incoll.size());
  int idx = -1;
  for (auto *incand : incoll) {
    double factor = puppi ? incand->puppiW() : 1;
    idx++;
    if (factor*incand->pt()<minPt)
      continue;
    vpj.emplace_back(factor*incand->px(),factor*incand->py(),
                     factor*incand->pz(),factor*incand->e());
    vpj.back().set_user_index(idx);
  }
  return vpj;
}

VPseudoJet ConvertPFCands(panda::RefVector<panda::PFCand> &incoll, bool puppi, double minPt) 
{
  std::vector<const panda::PFCand*> outcoll;
  outcoll.reserve(incoll.size());
  for (auto incand : incoll)
    outcoll.push_back(incand.get());

  return ConvertPFCands(outcoll, puppi, minPt);
}

VPseudoJet ConvertPFCands(panda::PFCandCollection &incoll, bool puppi, double minPt) 
{
  std::vector<const panda::PFCand*> outcoll;
  outcoll.reserve(incoll.size());
  for (auto &incand : incoll)
    outcoll.push_back(&incand);

  return ConvertPFCands(outcoll, puppi, minPt);
}

////////////////////////////////////////////////////////////////////////////////////

double TTNLOToNNLO(double pt) 
{
    double a = 0.1102;
    double b = 0.1566;
    double c = -3.685e-4;
    double d = 1.098;

    return TMath::Min(1.25,
                        a*TMath::Exp(-b*pow(pt,2)+1) + c*pt + d);
}

bool ElectronIP(double eta, double dxy, double dz) 
{
  double aeta = fabs(eta);
  if (aeta<1.4442) {
    return (dxy < 0.05 && dz < 0.10) ;
  } else {
    return (dxy < 0.10 && dz < 0.20);
  }
}

bool MuonIP(double dxy, double dz) 
{
  return (dxy < 0.02 && dz < 0.10);
}

bool IsMatched(std::vector<panda::Particle*>*objects,
               double deltaR2, double eta, double phi) 
{
  for (auto *x : *objects) {
    if (x->pt()>0) {
      if ( DeltaR2(x->eta(),x->phi(),eta,phi) < deltaR2 )
        return true;
    }
  }
  return false;
}

