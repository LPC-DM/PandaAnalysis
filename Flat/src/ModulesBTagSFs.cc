
#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>
#include <iostream>

#define EGMSCALE 1

using namespace panda;
using namespace std;

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
                               GeneralTree::BTagShift shift,GeneralTree::BTagJet jettype, bool do2, bool do3) 
{
  float sf0 = 1, sf1 = 1, sfGT0 = 1, sf2=1, sf3=1;
  float prob_mc0=1, prob_data0=1;
  float prob_mc1=0, prob_data1=0;
  unsigned int nC = cands.size();

  for (unsigned int iC=0; iC!=nC; ++iC) {
    double sf_i = sfs[iC];
    double eff_i = cands[iC].eff;
    prob_mc0 *= (1-eff_i);
    prob_data0 *= (1-sf_i*eff_i);
    float tmp_mc1=1, tmp_data1=1;
    for (unsigned int jC=0; jC!=nC; ++jC) {
      if (iC==jC) continue;
      double sf_j = sfs[jC];
      double eff_j = cands[jC].eff;
      tmp_mc1 *= (1-eff_j);
      tmp_data1 *= (1-eff_j*sf_j);
    }
    prob_mc1 += eff_i * tmp_mc1;
    prob_data1 += eff_i * sf_i * tmp_data1;
  }
  
  if (nC>0) {
    sf0 = prob_data0/prob_mc0;
    sf1 = prob_data1/prob_mc1;
    sfGT0 = (1-prob_data0)/(1-prob_mc0);
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
        prob_mc2 += eff_i * eff_j * tmp_mc2;
        prob_data2 += eff_i * sf_i * eff_j * sf_j * tmp_data2;
      }
    }

    if (nC>1) {
      sf2 = prob_data2/prob_mc2;
    }

    p.tag=GeneralTree::b2; gt->sf_btags[p] = sf2;
  }

  if (do3) {
    float prob_mc3=0, prob_data3=0;
    unsigned int nC = cands.size();
    for (unsigned int iC=0; iC!=nC; ++iC) {
      double sf_i = sfs[iC], eff_i = cands[iC].eff;
      for (unsigned int jC=iC+1; jC!=nC; ++jC) {
        double sf_j = sfs[jC], eff_j = cands[jC].eff;
        for (unsigned int kC=iC+2; kC!=nC; ++kC) { 
          double sf_k = sfs[kC], eff_k = cands[kC].eff;
          float tmp_mc3=1, tmp_data3=1;
          for (unsigned int lC=0; lC!=nC; ++lC) {
	    if (lC==iC || lC==jC || lC==kC) continue;
	    double sf_l = sfs[lC], eff_l = cands[lC].eff;
	    tmp_mc3 *= (1-eff_l);
	    tmp_data3 *= (1-eff_l*sf_l);
          }
          prob_mc3 += eff_i * eff_j * eff_k * tmp_mc3;
          prob_data3 += eff_i * sf_i * eff_j * sf_j * eff_k * sf_k * tmp_data3;
        }
      }
    }

    if (nC>2) {
      sf3 = prob_data3/prob_mc3;
    }

    p.tag=GeneralTree::b3; gt->sf_btags[p] = sf3;
  }

}


void PandaAnalyzer::JetBtagSFs() 
{
      // now get the jet btag SFs
      vector<btagcand> btagcands;
      vector<double> sf_cent, sf_bUp, sf_bDown, sf_mUp, sf_mDown;

      unsigned int nJ = bCandJets.size();
      for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
        panda::Jet *jet = bCandJets.at(iJ);
        bool isIsoJet=false;
      
// 
        if (!analysis->boosted || // if we do not consider fatjets, everything is an isojet 
            std::find(isoJets.begin(), isoJets.end(), jet) != isoJets.end()) // otherwise, explicitly check isojet
          isIsoJet = true;
        int flavor = bCandJetGenFlavor[jet];
        // float genpt = bCandJetGenPt[jet]; // not needed right now but it's here if it becomes needed
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
        if (isIsoJet) {
          if (analysis->boosted) {
            if (jet==isoJets.at(0))
              gt->isojet1Flav = flavor;
            else if (jet==isoJets.at(1))
              gt->isojet2Flav = flavor;
          }

          BTagType wp = bJetM;
          if (analysis->boosted) wp = bJetL;          


          CalcBJetSFs(wp,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
          btagcands.emplace_back(iJ,flavor,eff,sf,sfUp,sfDown);
          sf_cent.push_back(sf);

          if (flavor>0) {
            sf_bUp.push_back(sfUp); sf_bDown.push_back(sfDown);
            sf_mUp.push_back(sf); sf_mDown.push_back(sf);
          } else {
            sf_bUp.push_back(sf); sf_bDown.push_back(sf);
            sf_mUp.push_back(sfUp); sf_mDown.push_back(sfDown);
          }

        }
      } // loop over jets

      EvalBTagSF(btagcands,sf_cent,GeneralTree::bCent,GeneralTree::bJet);
      EvalBTagSF(btagcands,sf_bUp,GeneralTree::bBUp,GeneralTree::bJet);
      EvalBTagSF(btagcands,sf_bDown,GeneralTree::bBDown,GeneralTree::bJet);
      EvalBTagSF(btagcands,sf_mUp,GeneralTree::bMUp,GeneralTree::bJet);
      EvalBTagSF(btagcands,sf_mDown,GeneralTree::bMDown,GeneralTree::bJet);

    tr->TriggerEvent("ak4 gen-matching");
}

void PandaAnalyzer::JetCMVAWeights() 
{
  for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
    GeneralTree::csvShift shift = gt->csvShifts[iShift];
    gt->sf_csvWeights[shift] = 1;
  }
  if (bCandJets.size() < 1) return;

  //get vectors of jet properties
  std::vector<double> jetPts, jetEtas, jetCSVs, jetCMVAs;
  std::vector<int> jetFlavors;
  unsigned int nJ = bCandJets.size();
  //cout<<"nJ= "<<nJ<<endl; 
  jetPts.reserve(nJ);
  jetEtas.reserve(nJ);
  jetCSVs.reserve(nJ);
  jetCMVAs.reserve(nJ);
  jetFlavors.reserve(nJ);
  for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
    panda::Jet *jet = bCandJets.at(iJ);
    jetPts.push_back(jet->pt());
    jetEtas.push_back(jet->eta());
    jetCSVs.push_back(jet->csv);
    jetCMVAs.push_back(jet->cmva);
    int flavor = bCandJetGenFlavor[jet];
    jetFlavors.push_back(flavor);
    //cout<<"jetCSV= "<<jet->csv<<endl; 
  }
  // throwaway addresses
  double csvWgtHF, csvWgtLF, csvWgtCF, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF;
  for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
    GeneralTree::csvShift shift = gt->csvShifts[iShift];
    if (analysis->useCMVA) {
      gt->sf_csvWeights[shift] = cmvaReweighter->getCSVWeight(jetPts,jetEtas,jetCMVAs,jetFlavors, shift, cmvaWgtHF, cmvaWgtLF, cmvaWgtCF);
    }
    else 
      gt->sf_csvWeights[shift] = csvReweighter->getCSVWeight(jetPts,jetEtas,jetCSVs,jetFlavors, shift, csvWgtHF, csvWgtLF, csvWgtCF);
  }

}
