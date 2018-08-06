#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>
#include "PandaAnalysis/Utilities/src/NeutrinoSolver.cc"
#include "PandaAnalysis/Utilities/interface/Helicity.h"
#include "TLorentzRotation.h"
#define EGMSCALE 1
#define TWOPI 6.28318531

using namespace panda;
using namespace std;


void PandaAnalyzer::SetupJES()
{
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
}

void PandaAnalyzer::JetBasics() 
{
  gt->barrelJet12Pt = 0;
  gt->barrelHT = 0;
  unsigned nBarrelJets = 0;
  panda::Jet *jet1=0, *jet2=0;
  jot1=0; jot2=0;
  jotUp1=0; jotUp2=0;
  jotDown1=0; jotDown2=0;
  jetUp1=0; jetUp2=0;
  jetDown1=0; jetDown2=0;
  gt->dphipuppimet=999; gt->dphipfmet=999;
  gt->dphipuppiUW=999; gt->dphipfUW=999;
  gt->dphipuppiUZ=999; gt->dphipfUZ=999;
  gt->dphipuppiUA=999; gt->dphipfUA=999;
  float maxJetEta = (analysis->vbf) ? 4.7 : 4.5;
  unsigned nJetDPhi = (analysis->vbf) ? 4 : 5;

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
    if (analysis->vbf && !jet.loose)
      continue;

    if (jet.pt()>jetPtThreshold || jet.pt()>bJetPtThreshold) { // nominal or b jets

      // get jet flavor
      int flavor=0;
      float genpt=0;
      if (analysis->jetFlavorPartons) {
        // here we try to match to hard partons
        for (auto& gen : event.genParticles) {
          int apdgid = abs(gen.pdgid);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(jet.eta(),jet.phi(),gen.eta(),gen.phi());
          if (dr2<0.09) {
            genpt = gen.pt();
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid;
              break;
            } else {
              flavor=0;
            }
          }
        } 
      } else if (analysis->jetFlavorJets) {
        // or we can match to gen jets (probably better)
        for (auto &gen : event.ak4GenJets) {
          if (DeltaR2(gen.eta(), gen.phi(), jet.eta(), jet.phi()) < 0.09) {
            int apdgid = abs(gen.pdgid);
            genpt = gen.pt();
            if (apdgid == 4 || apdgid == 5) {
              flavor = apdgid;
              break;
            } else {
              flavor = 0;
            }
          }
        }
      }

      if (jet.pt()>bJetPtThreshold && fabs(jet.eta())<2.4) { // b jets

        if (jet.csv > 0.5426) {
          // loose WP
          ++(gt->jetNBtags);
          if (jet.csv > 0.8484) {
            // medum WP
            if (!analysis->vbf || 
                (analysis->vbf && fabs(jet.eta())<2.4))
            {
              ++(gt->jetNMBtags);
            }
          }
        }

        bCandJets.push_back(&jet);
        bCandJetGenFlavor[&jet] = flavor;
        bCandJetGenPt[&jet] = genpt;
      }

      if (jet.pt()>jetPtThreshold) { // nominal jets
        if ((analysis->monojet || analysis->boson || analysis->boosted || analysis->lepmonotop) && cleanedJets.size() >= NJET) 
          continue;
        cleanedJets.push_back(&jet);
        // Set jetGenPt, jetGenFlavor for these jets
        // This will be overwritten later if reclusterGen is turned on
        if (analysis->monojet || analysis->boson || analysis->boosted || analysis->lepmonotop) {
          gt->jetGenFlavor[cleanedJets.size()-1] = flavor;
          gt->jetGenPt    [cleanedJets.size()-1] = genpt ;
        }

        if (cleanedJets.size()<3) {
          bool isBad = GetCorr(cBadECALJets,jet.eta(),jet.phi()) > 0;
          if (isBad)
            gt->badECALFilter = 0;
        }

        //if (analysis->fatjet)
        if (analysis->boosted)
          IsoJet(jet);

        float csv = (fabs(jet.eta())<2.5) ? jet.csv : -1;
        float cmva = (fabs(jet.eta())<2.5) ? jet.cmva : -1;
        if (fabs(jet.eta())<2.4) {
          centralJets.push_back(&jet  );
          if (centralJets.size()==1) {
            jet1 = &jet;
            gt->jet1Pt = jet.pt();
            gt->jet1Eta = jet.eta();
            gt->jet1Phi = jet.phi();
            gt->jet1CSV = csv;
            gt->jet1CMVA = cmva;
            gt->jet1IsTight = jet.monojet ? 1 : 0;
            gt->jet1Flav = flavor;
            gt->jet1GenPt = genpt;
          } else if (centralJets.size()==2) {
            jet2 = &jet;
            gt->jet2Pt = jet.pt();
            gt->jet2Eta = jet.eta();
            gt->jet2Phi = jet.phi();
            gt->jet2CSV = csv;
            gt->jet2CMVA = cmva;
            gt->jet2Flav = flavor;
            gt->jet2GenPt = genpt;
          }
        }

        vJet.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());

        if (analysis->vbf || analysis->complicatedLeptons)
          JetVBFBasics(jet);

        if (analysis->boosted || analysis->boson || analysis->monojet || analysis->lepmonotop) {
          JetHbbBasics(jet);
          if (analysis->bjetRegression)
            JetBRegressionInfo(jet);
        }

        // compute dphi wrt mets
        if (cleanedJets.size() <= nJetDPhi) {
          gt->dphipuppimet = std::min(fabs(vJet.DeltaPhi(vPuppiMET)),(double)gt->dphipuppimet);
          gt->dphipfmet = std::min(fabs(vJet.DeltaPhi(vPFMET)),(double)gt->dphipfmet);
          if (analysis->recoil) {
            gt->dphipuppiUA = std::min(fabs(vJet.DeltaPhi(vpuppiUA)),(double)gt->dphipuppiUA);
            gt->dphipuppiUW = std::min(fabs(vJet.DeltaPhi(vpuppiUW)),(double)gt->dphipuppiUW);
            gt->dphipuppiUZ = std::min(fabs(vJet.DeltaPhi(vpuppiUZ)),(double)gt->dphipuppiUZ);
            gt->dphipfUA = std::min(fabs(vJet.DeltaPhi(vpfUA)),(double)gt->dphipfUA);
            gt->dphipfUW = std::min(fabs(vJet.DeltaPhi(vpfUW)),(double)gt->dphipfUW);
            gt->dphipfUZ = std::min(fabs(vJet.DeltaPhi(vpfUZ)),(double)gt->dphipfUZ);
          }
        }
      } // end jetPt>jetPtThreshold
    }

    if (analysis->varyJES) {
      if (fabs(jet.eta())<2.4) {
        if (jet.ptCorrUp>jetPtThreshold) 
          gt->nJet_jesUp++;
        if (jet.ptCorrDown>jetPtThreshold) 
          gt->nJet_jesDown++;
      }
      JetVaryJES(jet);

      if (jet.ptCorrUp>jetPtThreshold) 
        gt->nJot_jesUp++;
      if (jet.ptCorrDown>jetPtThreshold) 
        gt->nJot_jesDown++;
    }

  } // jet loop
  gt->barrelHTMiss = vBarrelJets.Pt();

  switch (gt->whichRecoil) {
    case 0: // MET
      gt->dphipuppiU = gt->dphipuppimet;
      gt->dphipfU = gt->dphipfmet;
      break;
    case -1: // photon
      gt->dphipuppiU = gt->dphipuppiUA;
      gt->dphipfU = gt->dphipfUA;
      break;
    case 1:
      gt->dphipuppiU = gt->dphipuppiUW;
      gt->dphipfU = gt->dphipfUW;
      break;
    case 2:
      gt->dphipuppiU = gt->dphipuppiUZ;
      gt->dphipfU = gt->dphipfUZ;
      break;
    default: // c'est impossible !
      break;
  }

  gt->nJet = centralJets.size();
  gt->nJot = cleanedJets.size();

  if (analysis->vbf) {
    JetVBFSystem();
  }

  tr->TriggerEvent("jets");

}

void PandaAnalyzer::JetHbbBasics(panda::Jet& jet)
{
  float csv = (fabs(jet.eta())<2.5) ? jet.csv : -1;
  float cmva = (fabs(jet.eta())<2.5) ? jet.cmva : -1;
  unsigned N = cleanedJets.size()-1;
  gt->jetPt[N]=jet.pt();
  gt->jetPtUp[N]=jet.ptCorrUp;
  gt->jetPtDown[N]=jet.ptCorrDown;
  gt->jetEta[N]=jet.eta();
  gt->jetPhi[N]=jet.phi();
  gt->jetM[N]=jet.m();
  gt->jetE[N]=jet.e();
  gt->jetCSV[N]=csv;
  gt->jetCMVA[N]=cmva;
  gt->jetQGL[N]=jet.qgl;

  tr->TriggerSubEvent("H->bb jet");
}

void PandaAnalyzer::JetBRegressionInfo(panda::Jet& jet)
{
  unsigned N = cleanedJets.size()-1;
  gt->jetEMFrac[N] = jet.cef + jet.nef;
  gt->jetHadFrac[N] = jet.chf + jet.nhf;
  gt->jetLeadingLepPt[N] = 0;
  gt->jetLeadingTrkPt[N] = 0;
  gt->jetNLep[N] = 0;
  for (const panda::Ref<panda::PFCand> &c_iter : jet.constituents) {
    if (!c_iter.isValid())
      continue;
    auto *pf = c_iter.get();
    if (pf->q() != 0) {
      float pt = pf->pt();
      gt->jetLeadingTrkPt[N] = max(pt, gt->jetLeadingTrkPt[N]);
      unsigned pdgid = abs(pf->pdgId());
      if (pdgid == 11 || pdgid == 13) {
        gt->jetNLep[N]++;
        if (pt > gt->jetLeadingLepPt[N]) {
          gt->jetLeadingLepPt[N] = pt;
          gt->jetLeadingLepPtRel[N] = pf->p4().Perp(jet.p4().Vect());
          gt->jetLeadingLepDeltaR[N] = sqrt(DeltaR2(pf->eta(), pf->phi(), jet.eta(), jet.phi()));
        }
      }
    }
  }

  auto& vert = jet.secondaryVertex;
  if (vert.isValid()) {
    gt->jetvtxPt[N] = vert->pt();
    gt->jetvtxMass[N] = vert->m();
    gt->jetvtx3Dval[N] = vert->vtx3DVal;
    gt->jetvtx3Derr[N] = vert->vtx3DeVal;
    gt->jetvtxNtrk[N] = vert->ntrk;
  }

  tr->TriggerSubEvent("b-jet reg info");
}

void PandaAnalyzer::JetVBFBasics(panda::Jet& jet)
{
  if (cleanedJets.size()==1) {
    jot1 = &jet;
    gt->jot1Pt = jet.pt();
    gt->jot1Eta = jet.eta();
    gt->jot1Phi = jet.phi();
    if (analysis->vbf && fabs(gt->jot1Eta)<2.4) { // if it's a central jet, must pass ID requirements
      gt->jot1VBFID = jet.monojet ? 1 : 0;
    } else { // if leading jet is not central, leave the event be
      gt->jot1VBFID = 1;
    }
  } else if (cleanedJets.size()==2) {
    jot2 = &jet;
    gt->jot2Pt = jet.pt();
    gt->jot2Eta = jet.eta();
    gt->jot2Phi = jet.phi();
  }

  if (fabs(jet.eta())<3.0) {
    gt->barrelHT += jet.pt();
    vBarrelJets += vJet;
    if (gt->barrelJet1Pt <= 0) {
      gt->barrelJet1Pt = jet.pt();
      gt->barrelJet1Eta = jet.eta();
    }
  }
  tr->TriggerSubEvent("VBF jet");

}

void PandaAnalyzer::IsoJet(panda::Jet& jet) 
{
  float maxIsoEta = (analysis->boosted) ? 4.5 : 2.5;
  bool isIsoJet = ( (gt->nFatjet==0) || 
      (fabs(jet.eta())<maxIsoEta 
       && DeltaR2(gt->fj1Eta,gt->fj1Phi,jet.eta(),jet.phi())>FATJETMATCHDR2) ); 

  if (isIsoJet) {
    isoJets.push_back(&jet);
    gt->nIsoJet++;
    float csv = (fabs(jet.eta())<2.5) ? jet.csv : -1;
    if (csv>0.5426)
      ++gt->isojetNBtags;
    if (isoJets.size()==1) {
      gt->isojet1Pt = jet.pt();
      gt->isojet1CSV = jet.csv;
    } else if (isoJets.size()==2) {
      gt->isojet2Pt = jet.pt();
      gt->isojet2CSV = jet.csv;
    }
    gt->jetIso[cleanedJets.size()-1]=1;
  } else {
    gt->jetIso[cleanedJets.size()-1]=0;
  }
  tr->TriggerSubEvent("iso jets");
}

void PandaAnalyzer::JetVaryJES(panda::Jet& jet)
{
  // do jes variation OUTSIDE of nominal jet check 
  if (jet.ptCorrUp>jetPtThreshold) { // Vary the pT up
    // forward and central jets:
    if (jet.ptCorrUp > gt->jot1PtUp) {
      if (jotUp1) {
        jotUp2 = jotUp1;
        gt->jot2PtUp = gt->jot1PtUp;
        gt->jot2EtaUp = gt->jot1EtaUp;
        gt->jot2PhiUp = gt->jot1PhiUp;
      }
      jotUp1 = &jet;
      gt->jot1PtUp = jet.ptCorrUp;
      gt->jot1EtaUp = jet.eta();
    } else if (jet.ptCorrUp > gt->jot2PtUp) {
      jotUp2 = &jet;
      gt->jot2PtUp = jet.ptCorrUp;
      gt->jot2EtaUp = jet.eta();
      gt->jot2PhiUp = jet.phi();
    }
    // central only jets:
    if (fabs(jet.eta()) < 2.4) {
      if (jet.ptCorrUp > gt->jet1PtUp) {
        if (jetUp1) {
          jetUp2 = jetUp1;
          gt->jet2PtUp = gt->jet1PtUp;
          gt->jet2EtaUp = gt->jet1EtaUp;
        }
        jetUp1 = &jet;
        gt->jet1PtUp = jet.ptCorrUp;
        gt->jet1EtaUp = jet.eta();
      } else if (jet.ptCorrUp > gt->jet2PtUp) {
        jetUp2 = &jet;
        gt->jet2PtUp = jet.ptCorrUp;
        gt->jet2EtaUp = jet.eta();
      }
    }
  }
  if (jet.ptCorrDown>jetPtThreshold) { // Vary the pT down:
    // forward and central jets:
    if (jet.ptCorrDown > gt->jot1PtDown) {
      if (jotDown1) {
        jotDown2 = jotDown1;
        gt->jot2PtDown = gt->jot1PtDown;
        gt->jot2EtaDown = gt->jot1EtaDown;
        gt->jot2PhiDown = gt->jot1PhiDown;
      }
      jotDown1 = &jet;
      gt->jot1PtDown = jet.ptCorrDown;
      gt->jot1EtaDown = jet.eta();
    } else if (jet.ptCorrDown > gt->jot2PtDown) {
      jotDown2 = &jet;
      gt->jot2PtDown = jet.ptCorrDown;
      gt->jot2EtaDown = jet.eta();
      gt->jot2PhiDown = jet.phi();
    }
    // central only jets:
    if (fabs(jet.eta()) < 2.4) {
      if (jet.ptCorrDown > gt->jet1PtDown) {
        if (jetDown1) {
          jetDown2 = jetDown1;
          gt->jet2PtDown = gt->jet1PtDown;
          gt->jet2EtaDown = gt->jet1EtaDown;
        }
        jetDown1 = &jet;
        gt->jet1PtDown = jet.ptCorrDown;
        gt->jet1EtaDown = jet.eta();
      } else if (jet.ptCorrDown > gt->jet2PtDown) {
        jetDown2 = &jet;
        gt->jet2PtDown = jet.ptCorrDown;
        gt->jet2EtaDown = jet.eta();
      }
    }
  }
  tr->TriggerSubEvent("vary jet JES");
}

void PandaAnalyzer::JetVBFSystem() 
{
  if (gt->nJot>1 && analysis->vbf) {
    TLorentzVector vj1, vj2;
    vj1.SetPtEtaPhiM(jot1->pt(),jot1->eta(),jot1->phi(),jot1->m());
    vj2.SetPtEtaPhiM(jot2->pt(),jot2->eta(),jot2->phi(),jot2->m());
    gt->jot12Mass = (vj1+vj2).M();
    gt->jot12DPhi = vj1.DeltaPhi(vj2);
    gt->jot12DEta = fabs(jot1->eta()-jot2->eta());

    if (analysis->varyJES && jotUp1 && jotUp2) {
      vj1.SetPtEtaPhiM(jotUp1->ptCorrUp,jotUp1->eta(),jotUp1->phi(),jotUp1->m());
      vj2.SetPtEtaPhiM(jotUp2->ptCorrUp,jotUp2->eta(),jotUp2->phi(),jotUp2->m());
      gt->jot12MassUp = (vj1+vj2).M();
      gt->jot12DPhiUp = vj1.DeltaPhi(vj2);
      gt->jot12DEtaUp = fabs(jotUp1->eta()-jotUp2->eta());
    }

    if (analysis->varyJES && jotDown1 && jotDown2) {
      vj1.SetPtEtaPhiM(jotDown1->ptCorrDown,jotDown1->eta(),jotDown1->phi(),jotDown1->m());
      vj2.SetPtEtaPhiM(jotDown2->ptCorrDown,jotDown2->eta(),jotDown2->phi(),jotDown2->m());
      gt->jot12MassDown = (vj1+vj2).M();
      gt->jot12DPhiDown = vj1.DeltaPhi(vj2);
      gt->jot12DEtaDown = fabs(jotDown1->eta()-jotDown2->eta());
    }
  }

  tr->TriggerSubEvent("VBF jet system");
}

void PandaAnalyzer::JetBosonReco() 
{
  float tmp_bosonpt=-99;
  float tmp_bosoneta=-99;
  float tmp_bosonphi=-99;
  float tmp_bosonm=-99;
  int tmp_bosonjtidx1=-1;
  int tmp_bosonjtidx2=-1;
  float tmp_bosoness=-999;
  float tmp_bosondr=-99;
  float tmp_bosondphibrf=-99;
  if (cleanedJets.size() > 1) {
     TLorentzVector bosondaughter1, bosondaughter2, bosonsystem,
                    bosondaughter1_jesUp, bosondaughter2_jesUp, bosonsystem_jesUp,
                    bosondaughter1_jesDown, bosondaughter2_jesDown, bosonsystem_jesDown;
     for (unsigned int i = 0;i<cleanedJets.size();i++){
         panda::Jet *jet_1 = cleanedJets.at(i);
         bosondaughter1.SetPtEtaPhiM(jet_1->pt(),jet_1->eta(),jet_1->phi(),jet_1->m());
         for (unsigned int j = i+1;j<cleanedJets.size();j++){
             panda::Jet *jet_2 = cleanedJets.at(j);
             bosondaughter2.SetPtEtaPhiM(jet_2->pt(),jet_2->eta(),jet_2->phi(),jet_2->m());
             bosonsystem = bosondaughter1 + bosondaughter2;
             double bosondr = sqrt(DeltaR2(bosondaughter1.Eta(),bosondaughter1.Phi(),bosondaughter2.Eta(),bosondaughter2.Phi()));
             double bosondphibrf = dPhiBRF(bosondaughter1, bosondaughter2, bosonsystem);
            // if (bosonsystem.Pt()>tmp_bosonpt){
             if (bosondphibrf>tmp_bosondphibrf){
                tmp_bosondphibrf = bosondphibrf;
                tmp_bosoness = (bosonsystem.Pt()*bosondr)/(2*bosonsystem.M());  
                tmp_bosondr = bosondr;
                tmp_bosonpt = bosonsystem.Pt();
                tmp_bosoneta = bosonsystem.Eta();
                tmp_bosonphi = bosonsystem.Phi();
                tmp_bosonm = bosonsystem.M();
                tmp_bosonjtidx1 = i;
                tmp_bosonjtidx2 = j;
                //tmp_bosondphibrf = dPhiBRF(bosondaughter1, bosondaughter2, bosonsystem);
             }
         }
     }
      gt->bosoness = tmp_bosoness;
      gt->bosonpt = tmp_bosonpt;
      gt->bosoneta = tmp_bosoneta;
      gt->bosonphi = tmp_bosonphi;
      gt->bosonm = tmp_bosonm;
      gt->bosondr = tmp_bosondr;
      gt->bosondphibrf = tmp_bosondphibrf;
      gt->bosonjtidx[0] = tmp_bosonjtidx1;
      gt->bosonjtidx[1] = tmp_bosonjtidx2;

    // Daughter jet energies varied Up
    bosondaughter1_jesUp.SetPtEtaPhiM(gt->jetPtUp[tmp_bosonjtidx1],gt->jetEta[tmp_bosonjtidx1],gt->jetPhi[tmp_bosonjtidx1],gt->jetM[tmp_bosonjtidx1]);
    bosondaughter2_jesUp.SetPtEtaPhiM(gt->jetPtUp[tmp_bosonjtidx2],gt->jetEta[tmp_bosonjtidx2],gt->jetPhi[tmp_bosonjtidx2],gt->jetM[tmp_bosonjtidx2]);
    //bosondaughter1_jesUp.SetPtEtaPhiM(jet_1->ptCorrUp,jet_1->eta(),jet_1->phi(),jet_1->m());
  //  bosondaughter2_jesUp.SetPtEtaPhiM(jet_2->ptCorrUp,jet_2->eta(),jet_2->phi(),jet_2->m());
    bosonsystem_jesUp = bosondaughter1_jesUp + bosondaughter2_jesUp;
    gt->bosonpt_jesUp = bosonsystem_jesUp.Pt();
    gt->bosoneta_jesUp = bosonsystem_jesUp.Eta();
    gt->bosonphi_jesUp = bosonsystem_jesUp.Phi();
    gt->bosonm_jesUp = bosonsystem_jesUp.M();
    // Daughter jet energies varied Down
    bosondaughter1_jesDown.SetPtEtaPhiM(gt->jetPtDown[tmp_bosonjtidx1],gt->jetEta[tmp_bosonjtidx1],gt->jetPhi[tmp_bosonjtidx1],gt->jetM[tmp_bosonjtidx1]);
    bosondaughter2_jesDown.SetPtEtaPhiM(gt->jetPtDown[tmp_bosonjtidx2],gt->jetEta[tmp_bosonjtidx2],gt->jetPhi[tmp_bosonjtidx2],gt->jetM[tmp_bosonjtidx2]);
    //bosondaughter1_jesDown.SetPtEtaPhiM(jet_1->ptCorrDown,jet_1->eta(),jet_1->phi(),jet_1->m());
    //bosondaughter2_jesDown.SetPtEtaPhiM(jet_2->ptCorrDown,jet_2->eta(),jet_2->phi(),jet_2->m());
    bosonsystem = bosondaughter1_jesDown + bosondaughter2_jesDown;
    gt->bosonpt_jesDown = bosonsystem_jesDown.Pt();
    gt->bosoneta_jesDown = bosonsystem_jesDown.Eta();
    gt->bosonphi_jesDown = bosonsystem_jesDown.Phi();
    gt->bosonm_jesDown = bosonsystem_jesDown.M();
    tr->TriggerSubEvent("Bare Hbb reco");
    
    
    TLorentzVector bosondaughters_corr[2], bosondaughters_corr_jesUp[2], bosondaughters_corr_jesDown[2];
    TLorentzVector bosonsystem_corr, bosonsystem_corr_jesUp, bosonsystem_corr_jesDown;
    if (analysis->bjetRegression && gt->bosonm>0.) {
      
      for (unsigned i = 0; i<2; i++) {
        // Central value for the jet energies to perform the b-jet regression
        bjetreg_vars[0] = gt->jetPt[gt->bosonjtidx[i]];
        bjetreg_vars[1] = gt->nJot;
        bjetreg_vars[2] = gt->jetEta[gt->bosonjtidx[i]];
        bjetreg_vars[3] = gt->jetE[gt->bosonjtidx[i]];
        bjetreg_vars[4] = gt->npv;
        bjetreg_vars[5] = gt->jetLeadingTrkPt[gt->bosonjtidx[i]];
        bjetreg_vars[6] = gt->jetLeadingLepPt[gt->bosonjtidx[i]];
        bjetreg_vars[7] = gt->jetNLep[gt->bosonjtidx[i]];
        bjetreg_vars[8] = gt->jetEMFrac[gt->bosonjtidx[i]];
        bjetreg_vars[9] = gt->jetHadFrac[gt->bosonjtidx[i]];
        
        // B-jet regression with jet energy varied up
        // Don't propagate the JES uncertainty to the hardest track/lepton or the EM fraction for now
        bjetreg_vars[0] = gt->jetPtUp[gt->bosonjtidx[i]];
        bjetreg_vars[3] = gt->jetE[gt->bosonjtidx[i]] * gt->jetPtUp[gt->bosonjtidx[i]] / gt->jetPt[gt->bosonjtidx[i]];
        gt->jetRegFac[i] = (bjetregReader->EvaluateRegression("BDT method"))[0];
        bosondaughters_corr_jesUp[i].SetPtEtaPhiM(
          gt->jetRegFac[i]*gt->jetPtUp[gt->bosonjtidx[i]],
          gt->jetEta[gt->bosonjtidx[i]],
          gt->jetPhi[gt->bosonjtidx[i]],
          cleanedJets.at(gt->bosonjtidx[i])->m()
        );
        // B-jet regression with jet energy varied down
        bjetreg_vars[0] = gt->jetPtDown[gt->bosonjtidx[i]];
        bjetreg_vars[3] = gt->jetE[gt->bosonjtidx[i]] * gt->jetPtDown[gt->bosonjtidx[i]] / gt->jetPt[gt->bosonjtidx[i]];
        gt->jetRegFac[i] = (bjetregReader->EvaluateRegression("BDT method"))[0];
        bosondaughters_corr_jesDown[i].SetPtEtaPhiM(
          gt->jetRegFac[i]*gt->jetPtDown[gt->bosonjtidx[i]],
          gt->jetEta[gt->bosonjtidx[i]],
          gt->jetPhi[gt->bosonjtidx[i]],
          cleanedJets.at(gt->bosonjtidx[i])->m()
        );
        // B-jet regression with central value for jet energy
        // Call this last so that the central value for jetRegFac[i] is stored in gt
        gt->jetRegFac[i] = (bjetregReader->EvaluateRegression("BDT method"))[0];
        bosondaughters_corr[i].SetPtEtaPhiM(
          gt->jetRegFac[i]*gt->jetPt[gt->bosonjtidx[i]],
          gt->jetEta[gt->bosonjtidx[i]],
          gt->jetPhi[gt->bosonjtidx[i]],
          cleanedJets.at(gt->bosonjtidx[i])->m()
        );

      }
      bosonsystem_corr = bosondaughters_corr[0] + bosondaughters_corr[1];
      gt->bosonm_reg = bosonsystem_corr.M();
      gt->bosonpt_reg = bosonsystem_corr.Pt();
      bosonsystem_corr_jesUp = bosondaughters_corr_jesUp[0] + bosondaughters_corr_jesUp[1];
      gt->bosonm_reg_jesUp = bosonsystem_corr_jesUp.M();
      gt->bosonpt_reg_jesUp = bosonsystem_corr_jesUp.Pt();
      bosonsystem_corr_jesDown = bosondaughters_corr_jesDown[0] + bosondaughters_corr_jesDown[1];
      gt->bosonm_reg_jesDown = bosonsystem_corr_jesDown.M();
      gt->bosonpt_reg_jesDown = bosonsystem_corr_jesDown.Pt();
      
      tr->TriggerSubEvent("Regr. Hbb reco");
    }
    if (gt->bosonm>0.) { 
      gt->bosonCosThetaJJ   = bosonsystem.CosTheta();
      // Collins-Soper frame calculation
      if (analysis->bjetRegression) {
        if (bosondaughters_corr[0].Pt() > bosondaughters_corr[1].Pt()) 
          gt->bosonCosThetaCSJ1 = CosThetaCollinsSoper(bosondaughters_corr[0],bosondaughters_corr[1]);
        else
          gt->bosonCosThetaCSJ1 = CosThetaCollinsSoper(bosondaughters_corr[1],bosondaughters_corr[0]);
      } else {
        if (bosondaughters_corr[0].Pt() > bosondaughters_corr[1].Pt()) 
          gt->bosonCosThetaCSJ1 = CosThetaCollinsSoper(bosondaughter1,bosondaughter2);
        else
          gt->bosonCosThetaCSJ1 = CosThetaCollinsSoper(bosondaughter2,bosondaughter1);
      }
      tr->TriggerSubEvent("Hbb spin correl.");
    }
 
    // Top mass reconstruction
    if (gt->bosonm>0. && gt->nLooseLep>0) {
      TLorentzVector leptonP4, metP4, nuP4, *jet1P4, *jet2P4, WP4, topP4;
      float dRJet1W, dRJet2W;
      bool jet1IsCloser;
      leptonP4=looseLeps[0]->p4();
      
      metP4.SetPtEtaPhiM(gt->pfmet, 0, gt->pfmetphi, 0);
      nuP4 = getNu4Momentum( leptonP4, metP4 );
      WP4 = leptonP4 + nuP4;
      // If using b-jet regression, use the regressed jets for the top mass reconstruction
      // Otherwise, use the un regressed jets
     // if (analysis->bjetRegression) {
     //   jet1P4 = &bosondaughters_corr[0]; jet2P4 = &bosondaughters_corr[1];
    //  } else {
      jet1P4 = &bosondaughter1; jet2P4 = &bosondaughter2;
    //  }
      dRJet1W=jet1P4->DeltaR(leptonP4); dRJet2W=jet2P4->DeltaR(leptonP4); jet1IsCloser = (dRJet1W < dRJet2W);
      topP4 = jet1IsCloser? (*jet1P4)+WP4 : (*jet2P4)+WP4;
      gt->topMassLep1Met = topP4.M();
      gt->topWBosonCosThetaCS = CosThetaCollinsSoper(WP4,jet1IsCloser? *jet1P4:*jet2P4);
      gt->topWBosonPt  = WP4.Pt();
      gt->topWBosonEta = WP4.Eta();
      gt->topWBosonPhi = WP4.Phi();

      if (analysis->varyJES) {
        // Top mass with Jet energy scale Up
        metP4.SetPtEtaPhiM(gt->pfmetUp, 0, gt->pfmetphi, 0);
        nuP4 = getNu4Momentum( leptonP4, metP4 );
        WP4 = leptonP4 + nuP4;
        // If using b-jet regression, use the regressed jets for the top mass reconstruction
        // Otherwise, use the un regressed jets
     //   if (analysis->bjetRegression) {
     //     jet1P4 = &bosondaughters_corr_jesUp[0]; jet2P4 = &bosondaughters_corr_jesUp[1];
      //  } else {
        jet1P4 = &bosondaughter1_jesUp; jet2P4 = &bosondaughter2_jesUp;
      //  }
        dRJet1W=jet1P4->DeltaR(leptonP4); dRJet2W=jet2P4->DeltaR(leptonP4); jet1IsCloser = (dRJet1W < dRJet2W);
        topP4 = jet1IsCloser? (*jet1P4)+WP4 : (*jet2P4)+WP4;
        gt->topMassLep1Met_jesUp = topP4.M();
        
        // Top mass with Jet energy scale Down
        metP4.SetPtEtaPhiM(gt->pfmetDown, 0, gt->pfmetphi, 0);
        nuP4 = getNu4Momentum( leptonP4, metP4 );
        WP4 = leptonP4 + nuP4;
        // If using b-jet regression, use the regressed jets for the top mass reconstruction
        // Otherwise, use the un regressed jets
      //  if (analysis->bjetRegression) {
     //     jet1P4 = &bosondaughters_corr_jesDown[0]; jet2P4 = &bosondaughters_corr_jesDown[1];
     //   } else {
        jet1P4 = &bosondaughter1_jesDown; jet2P4 = &bosondaughter2_jesDown;
     //   }
        dRJet1W=jet1P4->DeltaR(leptonP4); dRJet2W=jet2P4->DeltaR(leptonP4); jet1IsCloser = (dRJet1W < dRJet2W);
        topP4 = jet1IsCloser? (*jet1P4)+WP4 : (*jet2P4)+WP4;
        gt->topMassLep1Met_jesDown = topP4.M();
      }
      tr->TriggerSubEvent("Top(bW) reco");
    }


  }
}



void PandaAnalyzer::GenJetsNu()
{

  std::vector<fastjet::PseudoJet> finalStates;
  std::vector<panda::GenParticle*> bcs;
  for (auto &p : event.genParticles) {
    if (p.finalState && p.pt() > 0.001) {
      finalStates.emplace_back(p.px(), p.py(), p.pz(), p.e());
      continue;
    }
    unsigned apdgid = abs(p.pdgid);
    if (apdgid == 4 || apdgid == 5) {
      bcs.push_back(&p);
    }
  }

  fastjet::ClusterSequenceArea seq(finalStates, *jetDefGen, *areaDef);
  std::vector<fastjet::PseudoJet> allJets(seq.inclusive_jets(0.01));

  genJetsNu.reserve(allJets.size());
  for (auto &pj : allJets) {
    genJetsNu.emplace_back(panda::GenJet());
    genJetsNu.back().setXYZE(pj.px(), pj.py(), pj.pz(), pj.e());
    int flavor = 0;
    for (auto *bc : bcs) {
      if (DeltaR2(pj.eta(), pj.phi(), bc->eta(), bc->phi()) < 0.09) {
        flavor = abs(bc->pdgid);
        break;
      }
    }
    genJetsNu.back().pdgid = flavor;
  }
  tr->TriggerEvent("gen jets nu");
}

void PandaAnalyzer::JetHbbSoftActivity() {
  // Soft activity
  if (gt->bosonm>0.) {
    gt->sumEtSoft1=0; gt->nSoft2=0; gt->nSoft5=0; gt->nSoft10=0;
    // Define the ellipse of particles to forget about
    // https://math.stackexchange.com/questions/426150/what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate 
    // ((x-h)cos(A) + (y-k)sin(A))^2 /a^2 + ((x-h)sin(A) - (y-k)cos(A))^2 /b^2 <=1
    double ellipse_cosA, ellipse_sinA, ellipse_h, ellipse_k, ellipse_a, ellipse_b; {
      double ellipse_alpha;
      float phi1=gt->jetPhi[gt->bosonjtidx[0]], phi2=gt->jetPhi[gt->bosonjtidx[1]];
      float eta1=gt->jetEta[gt->bosonjtidx[0]], eta2=gt->jetEta[gt->bosonjtidx[1]];
      double phi1MinusPhi2 = phi1-phi2;
      double eta1MinusEta2 = eta1-eta2;
      double phi1MinusPhi2MPP = TVector2::Phi_mpi_pi(phi1MinusPhi2);
      ellipse_alpha = atan2( phi1MinusPhi2, eta1MinusEta2);
      // compute delta R using already computed qty's to save time
      ellipse_a = (sqrt(pow(eta1MinusEta2,2) + pow(TVector2::Phi_mpi_pi(phi1MinusPhi2),2)) + 1.)/2.; // Major axis 2*a = dR(b,b)+1
      ellipse_b = 1./2.; // Minor axis 2*b = 1
      ellipse_h = (eta1+eta2)/2.;
      ellipse_k = TVector2::Phi_mpi_pi(phi2 + phi1MinusPhi2MPP/2.);
      ellipse_cosA = cos(ellipse_alpha);
      ellipse_sinA = sin(ellipse_alpha);
      if (DEBUG > 10) {
        PDebug("PandaAnalyzer::JetBosonReco",Form("Calculating ellipse with (eta1,phi1)=(%.2f,%.2f), (eta2,phi2)=(%.2f,%.2f)",eta1,phi1,eta2,phi2));
        PDebug("PandaAnalyzer::JetBosonReco",Form("Found ellipse parameters (a,b,h,k,alpha)=(%.2f,%.2f,%.2f,%.2f,%.2f)",ellipse_a,ellipse_b,ellipse_h,ellipse_k,ellipse_alpha));
      }
    }

    // Find out which PF constituents to not use
    RefVector<PFCand> jet1Tracks = cleanedJets[gt->bosonjtidx[0]]->constituents,
                      jet2Tracks = cleanedJets[gt->bosonjtidx[1]]->constituents;

    // Get vector of pseudo jets for clustering
    panda::PFCandCollection &allTracks = event.pfCandidates;
    std::vector<fastjet::PseudoJet> softTracksPJ;
    softTracksPJ.reserve(allTracks.size());
    panda::PFCand *softTrack=0;
    for (auto &softTrackRef : allTracks) {
      softTrack = &softTrackRef;
      // Minimum track pT threshold (300 MeV default)
      if (softTrack->pt() < minSoftTrackPt) continue;
      // High quality track flag
      if (!softTrack->track.isValid() || !softTrack->track.get()->highPurity) continue;
      // Only consider tracks with dz < 0.2 w.r.t. the primary vertex
      if (fabs(softTrack->track.get()->dz()) > 0.2) continue;
      // Track cannot be a constituent of loose leptons or the two b-jets
      bool trackIsSpokenFor=false;
      if (!trackIsSpokenFor) for (UShort_t iJetTrack=0; iJetTrack<jet1Tracks.size(); iJetTrack++) {
        if (!jet1Tracks.at(iJetTrack).isValid()) continue;
        if (softTrack==jet1Tracks.at(iJetTrack).get()) { trackIsSpokenFor=true; break; }
      }
      if (!trackIsSpokenFor) for (UShort_t iJetTrack=0; iJetTrack<jet2Tracks.size(); iJetTrack++) {
        if (!jet2Tracks.at(iJetTrack).isValid()) continue;
        if (softTrack==jet2Tracks.at(iJetTrack).get()) { trackIsSpokenFor=true; break; }
      }
      if (!trackIsSpokenFor) for (int iLep=0; iLep<gt->nLooseLep; iLep++) {
        if (!looseLeps[iLep]->matchedPF.isValid()) continue;
        if (softTrack==looseLeps[iLep]->matchedPF.get()) { trackIsSpokenFor=true; break; }
      }
      if (trackIsSpokenFor) continue;
      // Require tracks to have the lowest |dz| with the hardest PV amongst all others
      int idxVertexWithMinAbsDz=-1; float minAbsDz=9999;
      for (int iV=0; iV!=event.vertices.size(); iV++) {
        auto& theVertex = event.vertices[iV];
        float vertexAbsDz = fabs(softTrack->dz(theVertex.position()));
        if (DEBUG > 10) PDebug("PandaAnalyzer::JetBosonReco",Form("Track has |dz| %.2f with vertex %d",vertexAbsDz,iV));
        if (vertexAbsDz >= minAbsDz) continue;
        idxVertexWithMinAbsDz = iV;
        minAbsDz = vertexAbsDz;
      }
      if (idxVertexWithMinAbsDz!=0 || minAbsDz>0.2) continue;
      if (DEBUG > 10) PDebug("PandaAnalyzer::JetBosonReco",Form("Track above 300 MeV has dz %.3f", softTrack->track.isValid()?softTrack->track.get()->dz():-1));
      // Need to add High Quality track flags :-)
      bool trackIsInHbbEllipse=false; {
        double ellipse_x = softTrack->eta();
        double ellipse_y = softTrack->phi();
        double ellipse_term1 = pow(
          (TVector2::Phi_mpi_pi(ellipse_x - ellipse_h)*ellipse_cosA + (ellipse_y - ellipse_k)*ellipse_sinA) / ellipse_a,
          2
        );
        double ellipse_term2 = pow(
          (TVector2::Phi_mpi_pi(ellipse_x - ellipse_h)*ellipse_sinA - (ellipse_y - ellipse_k)*ellipse_cosA) / ellipse_b,
          2
        );
        double ellipse_equation = (ellipse_term1 + ellipse_term2);
        trackIsInHbbEllipse = (ellipse_equation <= 1.);
      } if (trackIsInHbbEllipse) continue;
      softTracksPJ.emplace_back(softTrack->px(),softTrack->py(),softTrack->pz(),softTrack->e());
    }
    if (DEBUG > 10) PDebug("PandaAnalyzer::JetBosonReco",Form("Found %ld soft tracks that passed track quality cuts and the ellipse, jet constituency, and lepton matching vetoes",softTracksPJ.size()));
    softTrackJetDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4);
    fastjet::ClusterSequenceArea softTrackSequence(softTracksPJ, *softTrackJetDefinition, *areaDef);
    
    std::vector<fastjet::PseudoJet> softTrackJets(softTrackSequence.inclusive_jets(1.));
    if (DEBUG > 10) PDebug("PandaAnalyzer::JetBosonReco",Form("Clustered %ld jets of pT>1GeV using anti-kT algorithm (dR 0.4) from the soft tracks",softTrackJets.size()));
    for (std::vector<fastjet::PseudoJet>::size_type iSTJ=0; iSTJ<softTrackJets.size(); iSTJ++) {
      if (fabs(softTrackJets[iSTJ].eta()) > 4.7) continue;
      gt->sumEtSoft1 += softTrackJets[iSTJ].Et(); 
      if (DEBUG > 10) PDebug("PandaAnalyzer::JetBosonReco",Form("Soft jet %d has pT %.2f",(int)iSTJ,softTrackJets[iSTJ].pt()));
      if (softTrackJets[iSTJ].pt() >  2.)  gt->nSoft2++; else continue;
      if (softTrackJets[iSTJ].pt() >  5.)  gt->nSoft5++; else continue;
      if (softTrackJets[iSTJ].pt() > 10.) gt->nSoft10++; else continue;
    }
    tr->TriggerEvent("Soft activity");
  }
}
  double PandaAnalyzer::dPhiBRF(TLorentzVector v1, TLorentzVector v2, TLorentzVector v3){
   TVector3 hVelocity = v3.BoostVector();
   TLorentzRotation Boost(hVelocity);
   TLorentzRotation tosubjetRest = Boost.Inverse();
   TVector3 v1Dir = (tosubjetRest*v1).Vect();
   TVector3 v2Dir = (tosubjetRest*v2).Vect();

   float angle_mht = atan2((v1.Py()+v2.Py()),(v1.Px()+v2.Px())) + TMath::Pi();

   Float_t b1 = TMath::Pi() - acos(cos(v1Dir.Phi()-angle_mht));
   Float_t b2 = TMath::Pi() - acos(cos(v2Dir.Phi()-angle_mht));

   return TMath::Min(b1,b2);
}
