#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1

using namespace panda;
using namespace std;


void PandaAnalyzer::TopPTReweight()
{
      if (analysis->processType != kTT)
        return;

      gt->genWPlusPt = -1; gt->genWMinusPt = -1;
      for (auto& gen : event.genParticles) {
        if (abs(gen.pdgid)!=24)
          continue;
        if (analysis->firstGen) {
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
        if (analysis->firstGen) {
          if (gt->genWPlusPt>0 && gt->genWMinusPt>0)
            break;
        }
      }
      TLorentzVector vT,vTbar;
      float pt_t=0, pt_tbar=0;
      for (auto& gen : event.genParticles) {
        if (abs(gen.pdgid)!=6)
          continue;
        if (analysis->firstGen) {
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
        if (analysis->firstGen) {
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

    tr->TriggerEvent("tt SFs");
}

void PandaAnalyzer::VJetsReweight() 
{
      // calculate the mjj 
      TLorentzVector vGenJet;
      if (analysis->vbf) {
        // first find high pT leptons
        std::vector<GenParticle*> genLeptons;
        for (auto &gp : event.genParticles) {
          if (!gp.finalState)
            continue;
          unsigned id = abs(gp.pdgid);
          if ((id == 11 || id == 13) &&
              (gp.pt() > 20 && fabs(gp.eta()) < 4.7))
            genLeptons.push_back(&gp);
        }
        unsigned nGenJet = 0;
        TLorentzVector v;
        for (auto &gj : event.ak4GenJets) {
          bool matchesLep = false;
          for (auto *gl : genLeptons) {
            if (DeltaR2(gj.eta(), gj.phi(), gl->eta(), gl->phi()) < 0.16) {
              matchesLep = true;
              break;
            }
          }
          if (matchesLep)
            continue;

          v.SetPtEtaPhiM(gj.pt(), gj.eta(), gj.phi(), gj.m());
          if (nGenJet == 0) { gt->genJet1Pt = gj.pt(); gt->genJet1Eta = gj.eta(); }
          else if (nGenJet == 1) { gt->genJet2Pt = gj.pt(); gt->genJet2Eta = gj.eta(); }

          vGenJet += v;
          nGenJet++;
          if (nGenJet == 2)
            break;
        }
      }
      gt->genMjj = vGenJet.M();

      bool found = analysis->processType!=kA 
                   && analysis->processType!=kZ 
                   && analysis->processType!=kW
                   && analysis->processType!=kZEWK 
                   && analysis->processType!=kWEWK;
      if (found)
        return;

      int target=24;
      if (analysis->processType==kZ || analysis->processType==kZEWK) target=23;
      if (analysis->processType==kA) target=22;

      for (auto& gen : event.genParticles) {
        if (found) break;
        int apdgid = abs(gen.pdgid);
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
          if (analysis->processType==kZ) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonMass = gen.m();
            gt->genBosonEta = gen.eta();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
            gt->sf_qcdV = GetCorr(cZNLO,gt->genBosonPt);
            gt->sf_ewkV = GetCorr(cZEWK,gt->genBosonPt);
            if (analysis->vbf) {
              gt->sf_qcdV_VBF = GetCorr(cVBF_ZNLO,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBF2l = GetCorr(cVBF_ZllNLO,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_ZNLO,gt->genBosonPt);
              gt->sf_qcdV_VBF2lTight = GetCorr(cVBFTight_ZllNLO,gt->genBosonPt);
            }
            found=true;
          } else if (analysis->processType==kW) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonMass = gen.m();
            gt->genBosonEta = gen.eta();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
            gt->sf_qcdV = GetCorr(cWNLO,gt->genBosonPt);
            gt->sf_ewkV = GetCorr(cWEWK,gt->genBosonPt);
            if (analysis->vbf) {
              gt->sf_qcdV_VBF = GetCorr(cVBF_WNLO,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_WNLO,gt->genBosonPt);
            }
            found=true;
          } else if (analysis->processType==kZEWK) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonMass = gen.m();
            gt->genBosonEta = gen.eta();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
            if (analysis->vbf) {
              gt->sf_qcdV_VBF = GetCorr(cVBF_EWKZ,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
            }
          } else if (analysis->processType==kWEWK) {
            gt->trueGenBosonPt = gen.pt();
            gt->genBosonMass = gen.m();
            gt->genBosonEta = gen.eta();
            gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
            if (analysis->vbf) {
              gt->sf_qcdV_VBF = GetCorr(cVBF_EWKW,gt->genBosonPt,gt->genMjj);
              gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
            }
          } else if (analysis->processType==kA) {
            // take the highest pT
            if (gen.pt() > gt->trueGenBosonPt) {
              gt->trueGenBosonPt = gen.pt();
              gt->genBosonMass = gen.m();
              gt->genBosonEta = gen.eta();
              gt->genBosonPt = bound(gen.pt(),genBosonPtMin,genBosonPtMax);
              gt->sf_qcdV = GetCorr(cANLO,gt->genBosonPt);
              gt->sf_ewkV = GetCorr(cAEWK,gt->genBosonPt);
              gt->sf_qcdV2j = GetCorr(cANLO2j,gt->genBosonPt);
            }
          }
        } // target matches
      } // gen particle loop ends

      //now for the cases where we did not find a gen boson
      if (gt->genBosonPt < 0) {

        TLorentzVector vpt(0,0,0,0);

        for (auto& part : event.genParticles) {
          int pdgid = part.pdgid;
          unsigned int abspdgid = abs(pdgid);

          if ((abspdgid == 11 || abspdgid == 13) &&
              (part.statusFlags == GenParticle::kIsPrompt || 
               part.statusFlags == GenParticle::kIsTauDecayProduct || 
               part.statusFlags == GenParticle::kIsPromptTauDecayProduct || 
               part.statusFlags == GenParticle::kIsDirectTauDecayProduct || 
               part.statusFlags == GenParticle::kIsDirectPromptTauDecayProduct )) {

            //ideally you want to have dressed leptons (lepton + photon), 
            //but we have in any ways have a photon veto in the analysis
            if (IsMatched(&matchLeps,0.01,part.eta(),part.phi()))
              vpt += part.p4();
          }
          
          if ((abspdgid == 12 || abspdgid == 14 || abspdgid == 16) && part.finalState==1) {
            vpt += part.p4();
          }
        }
        
        gt->genBosonPt = bound(vpt.Pt(),genBosonPtMin,genBosonPtMax);
        gt->trueGenBosonPt = vpt.Pt();
        gt->genBosonMass = vpt.M();
        gt->genBosonEta = vpt.Eta();

        if (analysis->processType==kZ) {
          gt->sf_qcdV = GetCorr(cZNLO,gt->genBosonPt);
          gt->sf_ewkV = GetCorr(cZEWK,gt->genBosonPt);
          if (analysis->vbf) {
            gt->sf_qcdV_VBF = GetCorr(cVBF_ZNLO,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBF2l = GetCorr(cVBF_ZllNLO,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_ZNLO,gt->genBosonPt);
            gt->sf_qcdV_VBF2lTight = GetCorr(cVBFTight_ZllNLO,gt->genBosonPt);
          }
        } 
        else if (analysis->processType==kW) {
          gt->sf_qcdV = GetCorr(cWNLO,gt->genBosonPt);
          gt->sf_ewkV = GetCorr(cWEWK,gt->genBosonPt);
          if (analysis->vbf) {
            gt->sf_qcdV_VBF = GetCorr(cVBF_WNLO,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBFTight = GetCorr(cVBFTight_WNLO,gt->genBosonPt);
          }
        } 
        else if (analysis->processType==kZEWK) {
          if (analysis->vbf) {
            gt->sf_qcdV_VBF = GetCorr(cVBF_EWKZ,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
          }
        } 
        else if (analysis->processType==kWEWK) {
          if (analysis->vbf) {
            gt->sf_qcdV_VBF = GetCorr(cVBF_EWKW,gt->genBosonPt,gt->genMjj);
            gt->sf_qcdV_VBFTight = gt->sf_qcdV_VBF; // for consistency
          }
        }  
      }

    tr->TriggerEvent("qcd/ewk SFs");
}


void PandaAnalyzer::SignalInfo()
{
      if (analysis->processType != kSignal)
        return;

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
}

void PandaAnalyzer::QCDUncs()
{
      gt->pdfUp = 1 + event.genReweight.pdfDW;
      gt->pdfDown = 1 - event.genReweight.pdfDW;
      auto &genReweight = event.genReweight;
      for (unsigned iS=0; iS!=6; ++iS) {
        float s = 0;
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
        s += 1;
        gt->scale[iS] = s; 
        gt->scaleUp = max(float(gt->scaleUp),float(s));
        gt->scaleDown = min(float(gt->scaleDown),float(s));
      }
      tr->TriggerEvent("qcd uncertainties");
}

void PandaAnalyzer::SignalReweights()
{
      unsigned nW = wIDs.size();
      if (nW) {
        for (unsigned iW=0; iW!=nW; ++iW) {
          gt->signal_weights[wIDs[iW]] = event.genReweight.genParam[iW];
        }
      }
}
double PandaAnalyzer::WeightEWKCorr(float pt, int type) {
  double parWZ08[2] = { 2.85714,-0.05714};
  double parZZ08[2] = {-4.57143,-0.06857};
  double parWZ14[3] = {3.69800,-0.0726117,0.0000318044};
  double parZZ14[3] = {-0.586985000,-0.099845900,0.0000445083};
  double corrA = 0.0;
  double corrB = 0.0;
  if     (type == 0) { // WZ13
    corrA = (parWZ08[0]+parWZ08[1]*pt)/100.;
    corrB = (parWZ14[0]+parWZ14[1]*pt+parWZ14[2]*pt*pt)/100.;
  }
  else if (type == 1) { // ZZ13
    corrA = (parZZ08[0]+parZZ08[1]*pt)/100.;
    corrB = (parZZ14[0]+parZZ14[1]*pt+parZZ14[2]*pt*pt)/100.;
  }
  double corr = corrB - (corrB-corrA)/6.;

  if (corr >= 0.0) return 1.0;
  return (1.0+corr);
}

double PandaAnalyzer::WeightZHEWKCorr(float baseCorr) {
  return (baseCorr+0.31+0.11)/((1-0.053)+0.31+0.11);
}
void PandaAnalyzer::GenStudyEWK() {
  gt->genLep1Pt = 0;
  gt->genLep1Eta = -1;
  gt->genLep1Phi = -1;
  gt->genLep1PdgId = 0;
  gt->genLep2Pt = 0;
  gt->genLep2Eta = -1;
  gt->genLep2Phi = -1;
  gt->genLep2PdgId = 0;
  gt->genLep3Pt = 0;
  gt->genLep3Eta = -1;
  gt->genLep3Phi = -1;
  gt->genLep3PdgId = 0;
  gt->genLep4Pt = 0;
  gt->genLep4Eta = -1;
  gt->genLep4Phi = -1;
  gt->genLep4PdgId = 0;
  gt->looseGenLep1PdgId = 0;
  gt->looseGenLep2PdgId = 0;
  gt->looseGenLep3PdgId = 0;
  gt->looseGenLep4PdgId = 0;
  gt->looseGenPho1PdgId = 0;
  if (isData) return;
  TLorentzVector v1,v2,v3,v4,p1;
  if (gt->nLooseLep>=1) {
    panda::Lepton *lep1=looseLeps[0];
    v1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());
  }
  if (gt->nLooseLep>=2) {
    panda::Lepton *lep2=looseLeps[1];
    v2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());
  }
  if (gt->nLooseLep>=3) {
    panda::Lepton *lep3=looseLeps[2];
    v3.SetPtEtaPhiM(lep3->pt(),lep3->eta(),lep3->phi(),lep3->m());
  }
  if (gt->nLooseLep>=4) {
    panda::Lepton *lep4=looseLeps[3];
    v4.SetPtEtaPhiM(lep4->pt(),lep4->eta(),lep4->phi(),lep4->m());
  }
  if (gt->nLoosePhoton>=1) {
    panda::Photon *pho = loosePhos[0];
    p1.SetPtEtaPhiM(pho->pt(),pho->eta(),pho->phi(),0.);
  }
  // gen lepton matching
  std::vector<int> targetsLepton;
  std::vector<int> targetsPhoton;
  std::vector<int> targetsV;
  std::vector<int> targetsH;
  std::vector<int> targetsTop;
  std::vector<int> targetsN;

  int nGen = event.genParticles.size();
  for (int iG=0; iG!=nGen; ++iG) {
    auto& part(event.genParticles.at(iG));
    int pdgid = part.pdgid;
    unsigned int abspdgid = abs(pdgid);
    
     if ((abspdgid == 11 || abspdgid == 13) 
         && part.finalState 
         && (part.testFlag(GenParticle::kIsPrompt) 
             || part.statusFlags == GenParticle::kIsPrompt 
             || part.testFlag(GenParticle::kIsTauDecayProduct) 
             || part.testFlag(GenParticle::kIsPromptTauDecayProduct) 
             || part.testFlag(GenParticle::kIsDirectTauDecayProduct) 
             || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) 
             || (part.parent.isValid() && abs(part.parent->pdgid) == 15)
            )
        )
    {
      targetsLepton.push_back(iG);
    }

    if (abspdgid == 22 && part.finalState)
      targetsPhoton.push_back(iG);

    if (abspdgid == 23 || abspdgid == 24)
      targetsV.push_back(iG);

    if (abspdgid == 25)
      targetsH.push_back(iG);

    if (abspdgid == 6)
      targetsTop.push_back(iG);

    if (abspdgid == 12 || abspdgid == 14 || abspdgid == 16)
      targetsN.push_back(iG);

  } //looking for targets

  tr->TriggerSubEvent("check gen infos");

  for (int jG : targetsPhoton) {
    auto& partph(event.genParticles.at(jG));
    
    if (partph.pt() <= 10) continue; // ignore low pt photons

    if (p1.Pt() > 0 && DeltaR2(partph.eta(),partph.phi(),p1.Eta(),p1.Phi()) < 0.01) {
        gt->looseGenPho1PdgId = 3;
    }    
  }

  TLorentzVector rhoP4(0,0,0,0);
  double bosonPtMin = 1000000000;
  for (int iG : targetsLepton) {
    auto& part(event.genParticles.at(iG));
    TLorentzVector dressedLepton;
    dressedLepton.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());
  
    rhoP4 = rhoP4 + dressedLepton;
    for (int jG : targetsPhoton) {
      auto& partj(event.genParticles.at(jG));

      if (DeltaR2(part.eta(),part.phi(),partj.eta(),partj.phi()) < 0.01) {
        TLorentzVector photonV;
        photonV.SetPtEtaPhiM(partj.pt(),partj.eta(),partj.phi(),partj.m());
        dressedLepton += photonV;
      }
    }

    int pdgId = part.pdgid;
    if (part.testFlag(GenParticle::kIsTauDecayProduct) 
        || part.testFlag(GenParticle::kIsPromptTauDecayProduct) 
        || part.testFlag(GenParticle::kIsDirectTauDecayProduct) 
        || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) 
        ||(part.parent.isValid() && abs(part.parent->pdgid) == 15)) 
    {
      pdgId = 15 * part.pdgid/abs(part.pdgid);
    }
    if (dressedLepton.Pt() > gt->genLep1Pt) {
      gt->genLep4Pt    = gt->genLep3Pt; 
      gt->genLep4Eta   = gt->genLep3Eta;
      gt->genLep4Phi   = gt->genLep3Phi;
      gt->genLep4PdgId = gt->genLep3PdgId; 
      gt->genLep3Pt    = gt->genLep2Pt; 
      gt->genLep3Eta   = gt->genLep2Eta;
      gt->genLep3Phi   = gt->genLep2Phi;
      gt->genLep3PdgId = gt->genLep2PdgId; 
      gt->genLep2Pt    = gt->genLep1Pt; 
      gt->genLep2Eta   = gt->genLep1Eta;
      gt->genLep2Phi   = gt->genLep1Phi;
      gt->genLep2PdgId = gt->genLep1PdgId; 
      gt->genLep1Pt    = dressedLepton.Pt();
      gt->genLep1Eta   = dressedLepton.Eta();
      gt->genLep1Phi   = dressedLepton.Phi();
      gt->genLep1PdgId = pdgId;
    } 
    else if (dressedLepton.Pt() > gt->genLep2Pt) {
      gt->genLep4Pt    = gt->genLep3Pt; 
      gt->genLep4Eta   = gt->genLep3Eta;
      gt->genLep4Phi   = gt->genLep3Phi;
      gt->genLep4PdgId = gt->genLep3PdgId; 
      gt->genLep3Pt    = gt->genLep2Pt; 
      gt->genLep3Eta   = gt->genLep2Eta;
      gt->genLep3Phi   = gt->genLep2Phi;
      gt->genLep3PdgId = gt->genLep2PdgId; 
      gt->genLep2Pt    = dressedLepton.Pt();
      gt->genLep2Eta   = dressedLepton.Eta();
      gt->genLep2Phi   = dressedLepton.Phi();
      gt->genLep2PdgId = pdgId; 
    }
    else if (dressedLepton.Pt() > gt->genLep3Pt) {
      gt->genLep4Pt    = gt->genLep3Pt; 
      gt->genLep4Eta   = gt->genLep3Eta;
      gt->genLep4Phi   = gt->genLep3Phi;
      gt->genLep4PdgId = gt->genLep3PdgId; 
      gt->genLep3Pt    = dressedLepton.Pt();
      gt->genLep3Eta   = dressedLepton.Eta();
      gt->genLep3Phi   = dressedLepton.Phi();
      gt->genLep3PdgId = pdgId; 
    }
    else if (dressedLepton.Pt() > gt->genLep4Pt) {
      gt->genLep4Pt    = dressedLepton.Pt();
      gt->genLep4Eta   = dressedLepton.Eta();
      gt->genLep4Phi   = dressedLepton.Phi();
      gt->genLep4PdgId = pdgId; 
    }
    panda::Muon *mu; panda::Electron *ele;
    if (v1.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v1.Eta(),v1.Phi()) < 0.01) {
      if (part.testFlag(GenParticle::kIsTauDecayProduct) 
          || part.testFlag(GenParticle::kIsPromptTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) 
          || (part.parent.isValid() && abs(part.parent->pdgid) == 15)) 
      {
        gt->looseGenLep1PdgId = 2;
      } else if (part.testFlag(GenParticle::kIsPrompt) 
                 || part.statusFlags == GenParticle::kIsPrompt) 
      {
        gt->looseGenLep1PdgId = 1;
      }
      if (part.pdgid != looseLep1PdgId) 
        gt->looseGenLep1PdgId = -1 * gt->looseGenLep1PdgId;
    }

    if (v2.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v2.Eta(),v2.Phi()) < 0.01) {
      if (part.testFlag(GenParticle::kIsTauDecayProduct) 
          || part.testFlag(GenParticle::kIsPromptTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) 
          || (part.parent.isValid() && abs(part.parent->pdgid) == 15)) 
      {
        gt->looseGenLep2PdgId = 2;
      } else if (part.testFlag(GenParticle::kIsPrompt) 
                 || part.statusFlags == GenParticle::kIsPrompt) 
      {
        gt->looseGenLep2PdgId = 1;
      }
      if (part.pdgid != looseLep2PdgId) 
        gt->looseGenLep2PdgId = -1 * gt->looseGenLep2PdgId;
    }
    if (v3.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v3.Eta(),v3.Phi()) < 0.01) {
      if (part.testFlag(GenParticle::kIsTauDecayProduct) 
          || part.testFlag(GenParticle::kIsPromptTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) 
          || (part.parent.isValid() && abs(part.parent->pdgid) == 15)) 
      {
        gt->looseGenLep3PdgId = 2;
      } else if (part.testFlag(GenParticle::kIsPrompt) 
                 || part.statusFlags == GenParticle::kIsPrompt) 
      {
        gt->looseGenLep3PdgId = 1;
      }
      if (part.pdgid != looseLep3PdgId) 
        gt->looseGenLep3PdgId = -1 * gt->looseGenLep3PdgId;
    }
    if (v4.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v4.Eta(),v4.Phi()) < 0.01) {
      if (part.testFlag(GenParticle::kIsTauDecayProduct) 
          || part.testFlag(GenParticle::kIsPromptTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) 
          || (part.parent.isValid() && abs(part.parent->pdgid) == 15)) 
      {
        gt->looseGenLep4PdgId = 2;
      } else if (part.testFlag(GenParticle::kIsPrompt) 
                 || part.statusFlags == GenParticle::kIsPrompt) 
      {
        gt->looseGenLep4PdgId = 1;
      }
      if (part.pdgid != looseLep4PdgId) 
        gt->looseGenLep4PdgId = -1 * gt->looseGenLep4PdgId;
    }
    if (p1.Pt() > 0 && DeltaR2(part.eta(),part.phi(),p1.Eta(),p1.Phi()) < 0.01) {
      if (part.testFlag(GenParticle::kIsTauDecayProduct) 
          || part.testFlag(GenParticle::kIsPromptTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectTauDecayProduct) 
          || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) 
          || (part.parent.isValid() && abs(part.parent->pdgid) == 15)) 
      {
        gt->looseGenPho1PdgId = 2;
      } else if (part.testFlag(GenParticle::kIsPrompt) 
                 || part.statusFlags == GenParticle::kIsPrompt) 
      {
        gt->looseGenPho1PdgId = 1;
      }
    }
  }
  
  tr->TriggerSubEvent("gen leptons 2");
  
  unsigned char nNeutrinos=0;
  for (int iG : targetsN) {
    auto& part(event.genParticles.at(iG));
    TLorentzVector neutrino;
    neutrino.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());
    // check there is no further copy:
    bool isLastCopy=true;
    for (int kG : targetsN) {
      if (event.genParticles.at(kG).parent.isValid() && event.genParticles.at(kG).parent.get() == &part) {
        isLastCopy=false;
        break;
      }
    }
    if (!isLastCopy)
      continue;
    rhoP4 = rhoP4 + neutrino;
    nNeutrinos++;
  }

  tr->TriggerSubEvent("gen neutrinos");
  
  TLorentzVector higgsBosons(0,0,0,0);
  int nHBosons = 0;
  for (int iG : targetsH) {
    auto& part(event.genParticles.at(iG));
    TLorentzVector boson;
    boson.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());

    // check there is no further copy:
    bool isLastCopy=true;
    for (int kG : targetsH) {
      if (event.genParticles.at(kG).parent.isValid() 
          && event.genParticles.at(kG).parent.get() == &part) {
        isLastCopy=false;
        break;
      }
    }
    if (!isLastCopy) continue;
  
    if (abs(part.pdgid) == 25) {higgsBosons = higgsBosons + boson; nHBosons++;}
  }
  tr->TriggerSubEvent("gen higgs");

  TLorentzVector zBosons(0,0,0,0);
  TLorentzVector wBosons(0,0,0,0);
  int nZBosons = 0; int nWBosons = 0; vector<bool> wBosonQ; wBosonQ.reserve(8);
  for (int iG : targetsV) {
    auto& part(event.genParticles.at(iG));
    TLorentzVector boson;
    boson.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());

    // check there is no further copy:
    bool isLastCopy=true;
    for (int kG : targetsV) {
      if (event.genParticles.at(kG).parent.isValid() 
          && event.genParticles.at(kG).parent.get() == &part) {
        isLastCopy=false;
        break;
      }
    }
    if (!isLastCopy) continue;
  
    if (boson.Pt() < bosonPtMin) bosonPtMin = boson.Pt();
    if (abs(part.pdgid) == 23) {zBosons = zBosons + boson; nZBosons++;}
    if (abs(part.pdgid) == 24) {wBosons = wBosons + boson; nWBosons++; wBosonQ.push_back(part.pdgid>0); }
  }
  if (nZBosons+nWBosons == 0) bosonPtMin = 0;

  if (nZBosons >= 2) {
    double rho = 0.0; if (rhoP4.P() > 0) rho = rhoP4.Pt()/rhoP4.P();
    double ZZCorr[2] {1,1};
    ZZCorr[0] = WeightEWKCorr(bosonPtMin,1);
    float GENmZZ = zBosons.M();
    ZZCorr[1] = GetCorr(cqqZZQcdCorr,2,GENmZZ); // final state = 2 is fixed
    gt->sf_zz = ZZCorr[0]*ZZCorr[1];
    if (rho <= 0.3) gt->sf_zzUnc = (1.0+TMath::Abs((ZZCorr[0]-1)*(15.99/9.89-1)));
    else            gt->sf_zzUnc = (1.0+TMath::Abs((ZZCorr[0]-1)               ));
  } else {
    gt->sf_zz    = 1.0;
    gt->sf_zzUnc = 1.0;
  }

  if (nWBosons == 1 && nZBosons == 1) {
    TLorentzVector WZBoson = wBosons + zBosons;
    gt->sf_wz = GetCorr(cWZEwkCorr,WZBoson.M());
  } else {
    gt->sf_wz = 1.0;
  }
  if (analysis->processType==kH) { 
    if (nZBosons == 1 && gt->genLep2PdgId!=0) {
      gt->sf_vh     = WeightZHEWKCorr(GetCorr(cZllHEwkCorr,bound(zBosons.Pt(),0,499.999)));
      gt->sf_vhUp   = WeightZHEWKCorr(GetCorr(cZllHEwkCorrUp,bound(zBosons.Pt(),0,499.999)));
      gt->sf_vhDown = WeightZHEWKCorr(GetCorr(cZllHEwkCorrDown,bound(zBosons.Pt(),0,499.999)));
      gt->genBosonMass = zBosons.M();
      gt->genBosonEta  = zBosons.Eta();
      gt->genBosonPt   = zBosons.Pt();
      gt->trueGenBosonPt = higgsBosons.Pt();
    } else if (nZBosons==1 && nNeutrinos>0) {
      gt->sf_vh     = WeightZHEWKCorr(GetCorr(cZnnHEwkCorr,bound(zBosons.Pt(),0,499.999)));
      gt->sf_vhUp   = WeightZHEWKCorr(GetCorr(cZnnHEwkCorrUp,bound(zBosons.Pt(),0,499.999)));
      gt->sf_vhDown = WeightZHEWKCorr(GetCorr(cZnnHEwkCorrDown,bound(zBosons.Pt(),0,499.999)));
      gt->genBosonMass = zBosons.M();
      gt->genBosonEta  = zBosons.Eta();
      gt->genBosonPt   = zBosons.Pt();
      gt->trueGenBosonPt = higgsBosons.Pt();
    } else if (nWBosons==1 && gt->genLep1PdgId!=0 && nNeutrinos>0) {
      if (wBosonQ[0]>0) {
        gt->sf_vh     = WeightZHEWKCorr(GetCorr(cWpHEwkCorr,bound(wBosons.Pt(),0,499.999)));
        gt->sf_vhUp   = WeightZHEWKCorr(GetCorr(cWpHEwkCorrUp,bound(wBosons.Pt(),0,499.999)));
        gt->sf_vhDown = WeightZHEWKCorr(GetCorr(cWpHEwkCorrDown,bound(wBosons.Pt(),0,499.999)));
        gt->genWPlusPt = wBosons.Pt();
        gt->genWPlusEta = wBosons.Eta();
      } else {
        gt->sf_vh     = WeightZHEWKCorr(GetCorr(cWmHEwkCorr,bound(wBosons.Pt(),0,499.999)));
        gt->sf_vhUp   = WeightZHEWKCorr(GetCorr(cWmHEwkCorrUp,bound(wBosons.Pt(),0,499.999)));
        gt->sf_vhDown = WeightZHEWKCorr(GetCorr(cWmHEwkCorrDown,bound(wBosons.Pt(),0,499.999)));
        gt->genWMinusPt = wBosons.Pt();
        gt->genWMinusEta = wBosons.Eta();
      }
      gt->genBosonMass = higgsBosons.M();
      gt->genBosonEta  = higgsBosons.Eta();
      gt->genBosonPt   = higgsBosons.Pt();
    } else {
      gt->sf_vh     = 1.0;
      gt->sf_vhUp   = 1.0;
      gt->sf_vhDown = 1.0;
    }
  }

  tr->TriggerSubEvent("boson corrections");

  tr->TriggerEvent("EWK gen study");
}
