#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#define EGMSCALE 1

using namespace panda;
using namespace std;

void PandaAnalyzer::IncrementGenAuxFile(bool close)
{
  if (fAux) {
    fAux->WriteTObject(tAux, "inputs", "Overwrite");
    fAux->Close();
  }
  if (close)
    return;

  TString path = TString::Format(auxFilePath.Data(),auxCounter++);
  fAux = TFile::Open(path.Data(), "RECREATE");
  if (DEBUG) PDebug("PandaAnalyzer::IncrementAuxFile", "Opening "+path);
  tAux = new TTree("inputs","inputs");
  
  genJetInfo.particles.resize(NMAXPF);
  for (unsigned i = 0; i != NMAXPF; ++i) {
    genJetInfo.particles[i].resize(NGENPROPS);
  }

  genJetInfo.ecfs.resize(3);
  for (int o = 0; o != 3; ++o) {
    genJetInfo.ecfs.at(o).resize(4);
    for (int N = 0; N != 4; ++N) {
      genJetInfo.ecfs.at(o).at(N).resize(2);
    }
  }

  tAux->Branch("eventNumber",&(gt->eventNumber),"eventNumber/l");
  tAux->Branch("nprongs",&(genJetInfo.nprongs),"nprongs/I");
  tAux->Branch("partonpt",&(genJetInfo.partonpt),"partonpt/F");
  tAux->Branch("partonm",&(genJetInfo.partonm),"partonm/F");
  tAux->Branch("pt",&(genJetInfo.pt),"pt/F");
  tAux->Branch("msd",&(genJetInfo.msd),"msd/F");
  tAux->Branch("eta",&(genJetInfo.eta),"eta/F");
  tAux->Branch("phi",&(genJetInfo.phi),"phi/F");
  tAux->Branch("m",&(genJetInfo.m),"m/F");
  tAux->Branch("tau3",&(genJetInfo.tau3),"tau3/F");
  tAux->Branch("tau2",&(genJetInfo.tau2),"tau2/F");
  tAux->Branch("tau1",&(genJetInfo.tau1),"tau1/F");
  tAux->Branch("tau3sd",&(genJetInfo.tau3sd),"tau3sd/F");
  tAux->Branch("tau2sd",&(genJetInfo.tau2sd),"tau2sd/F");
  tAux->Branch("tau1sd",&(genJetInfo.tau1sd),"tau1sd/F");
  for (int o = 1; o != 4; ++o) {
    for (int N = 1; N != 5; ++N) {
      for (int beta = 1; beta != 3; ++beta) {
        TString bname = Form("%i_%i_%i",o,N,beta);
        tAux->Branch(bname,&(genJetInfo.ecfs[o-1][N-1][beta-1]),bname+"/F");
      }
    }
  }
  tAux->Branch("kinematics",&(genJetInfo.particles));

  fOut->cd();

  if (tr)
    tr->TriggerEvent("increment aux file");
}

void PandaAnalyzer::IncrementAuxFile(bool close)
{
  if (fAux) {
    fAux->WriteTObject(tAux, "inputs", "Overwrite");
    fAux->Close();
  }
  if (close)
    return;

  TString path = TString::Format(auxFilePath.Data(),auxCounter++);
  fAux = TFile::Open(path.Data(), "RECREATE");
  if (DEBUG) PDebug("PandaAnalyzer::IncrementAuxFile", "Opening "+path);
  tAux = new TTree("inputs","inputs");
  
  pfInfo.resize(NMAXPF);
  for (unsigned i = 0; i != NMAXPF; ++i) {
    pfInfo[i].resize(NPFPROPS);
  }
  tAux->Branch("kinematics",&pfInfo);
  
  svInfo.resize(NMAXSV);
  for (unsigned i = 0; i != NMAXSV; ++i) {
    svInfo[i].resize(NSVPROPS);
  }
  tAux->Branch("svs",&svInfo);

  tAux->Branch("msd",&fjmsd,"msd/F");
  tAux->Branch("pt",&fjpt,"pt/F");
  tAux->Branch("rawpt",&fjrawpt,"rawpt/F");
  tAux->Branch("eta",&fjeta,"eta/F");
  tAux->Branch("phi",&fjphi,"phi/F");
  tAux->Branch("rho",&(gt->fj1Rho),"rho/f");
  tAux->Branch("rawrho",&(gt->fj1RawRho),"rawrho/f");
  tAux->Branch("rho2",&(gt->fj1Rho2),"rho2/f");
  tAux->Branch("rawrho2",&(gt->fj1RawRho2),"rawrho2/f");
  tAux->Branch("nPartons",&(gt->fj1NPartons),"nPartons/I");
  tAux->Branch("nBPartons",&(gt->fj1NBPartons),"nBPartons/I");
  tAux->Branch("nCPartons",&(gt->fj1NCPartons),"nCPartons/I");
  tAux->Branch("partonM",&(gt->fj1PartonM),"partonM/f");
  tAux->Branch("partonPt",&(gt->fj1PartonPt),"partonPt/f");
  tAux->Branch("partonEta",&(gt->fj1PartonEta),"partonEta/f");
  tAux->Branch("tau32",&(gt->fj1Tau32),"tau32/f");
  tAux->Branch("tau32SD",&(gt->fj1Tau32SD),"tau32SD/f");
  tAux->Branch("tau21",&(gt->fj1Tau21),"tau21/f");
  tAux->Branch("tau21SD",&(gt->fj1Tau21SD),"tau21SD/f");
  tAux->Branch("eventNumber",&(gt->eventNumber),"eventNumber/l");
  tAux->Branch("maxcsv",&(gt->fj1MaxCSV),"maxcsv/f");
  tAux->Branch("mincsv",&(gt->fj1MinCSV),"mincsv/f");
  tAux->Branch("doubleb",&(gt->fj1DoubleCSV),"doubleb/f");

  gt->SetAuxTree(tAux);

  fOut->cd();

  if (tr)
    tr->TriggerEvent("increment aux file");
}

void PandaAnalyzer::RegisterTriggers() 
{
  for (auto &th : triggerHandlers) {
    unsigned N = th.paths.size();
    for (unsigned i = 0; i != N; i++) {
      unsigned panda_idx = event.registerTrigger(th.paths.at(i));
      th.indices[i] = panda_idx;
      if (DEBUG) PDebug("PandaAnalyzer::RegisterTriggers",
        Form("Got index %d for trigger path %s", panda_idx, th.paths.at(i).Data())
      );
    }
  }
}

bool PandaAnalyzer::RecoilPresel() 
{
    if ( (preselBits&kMonotop) || (preselBits&kMonohiggs) || 
         (preselBits&kMonojet) || (preselBits&kRecoil) ) 
    {
       if (event.recoil.max<175)
         return false;
    } 
    return true;
}

void PandaAnalyzer::TriggerEffs()
{

    // trigger efficiencies
    gt->sf_metTrig = GetCorr(cTrigMET,gt->pfmetnomu);
    gt->sf_metTrigZmm = GetCorr(cTrigMETZmm,gt->pfmetnomu);

    if (gt->nLooseElectron>0) {
      panda::Electron *ele1=0, *ele2=0;
      if (gt->nLooseLep>0) ele1 = dynamic_cast<panda::Electron*>(looseLeps[0]);
      if (gt->nLooseLep>1) ele2 = dynamic_cast<panda::Electron*>(looseLeps[1]);
      float eff1=0, eff2=0;
      if (ele1 && ele1->tight) {
        eff1 = GetCorr(cTrigEle, ele1->eta(), ele1->pt());
        if (ele2 && ele2->tight)
          eff2 = GetCorr(cTrigEle, ele2->eta(), ele2->pt());
        gt->sf_eleTrig = 1 - (1-eff1)*(1-eff2);
      }
    } // done with ele trig SF
    if (gt->nLooseMuon>0) {
      panda::Muon *mu1=0, *mu2=0;
      if (gt->nLooseLep>0) mu1 = dynamic_cast<panda::Muon*>(looseLeps[0]);
      if (gt->nLooseLep>1) mu2 = dynamic_cast<panda::Muon*>(looseLeps[1]);
      float eff1=0, eff2=0;
      if (mu1 && mu1->tight) {
        eff1 = GetCorr(
          cTrigMu,
          fabs(mu1->eta()),
          TMath::Max((float)26.,TMath::Min((float)499.99,(float)mu1->pt()))
        );
        if (mu2 && mu2->tight)
          eff2 = GetCorr(
            cTrigMu,
            fabs(mu2->eta()),
            TMath::Max((float)26.,TMath::Min((float)499.99,(float)mu2->pt()))
          );
        gt->sf_muTrig = 1 - (1-eff1)*(1-eff2);
      }
    } // done with mu trig SF

    if (gt->nLoosePhoton>0 && gt->loosePho1IsTight)
      gt->sf_phoTrig = GetCorr(cTrigPho,gt->loosePho1Pt);

    if (analysis->vbf) {
      gt->sf_metTrigVBF = GetCorr(cVBF_TrigMET,gt->barrelHTMiss);
      gt->sf_metTrigZmmVBF = GetCorr(cVBF_TrigMETZmm,gt->barrelHTMiss);
    }
    tr->TriggerEvent("triggers");
}

void PandaAnalyzer::Recoil()
{
    TLorentzVector vpfUp; vpfUp.SetPtEtaPhiM(gt->pfmetUp,0,gt->pfmetphi,0);
    TLorentzVector vpfDown; vpfDown.SetPtEtaPhiM(gt->pfmetDown,0,gt->pfmetphi,0);
    TLorentzVector vpfUWUp, vpfUWDown;
    TLorentzVector vObj1, vObj2;
    gt->whichRecoil = 0; // -1=photon, 0=MET, 1,2=nLep
    if (gt->nLooseLep>0) {
      panda::Lepton *lep1 = looseLeps.at(0);
      vObj1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());

      // one lep => W
      vpuppiUW = vPuppiMET+vObj1; gt->puppiUWmag=vpuppiUW.Pt(); gt->puppiUWphi=vpuppiUW.Phi();
      vpfUW = vPFMET+vObj1; gt->pfUWmag=vpfUW.Pt(); gt->pfUWphi=vpfUW.Phi();
      
      if (analysis->varyJES) {
        vpfUWUp = vpfUp+vObj1; gt->pfUWmagUp = vpfUWUp.Pt();
        vpfUWDown = vpfDown+vObj1; gt->pfUWmagDown = vpfUWDown.Pt();
      }

      if (gt->nLooseLep>1 && looseLep1PdgId+looseLep2PdgId==0) {
	//std::cout<<"two OS lep => Z"<<std::endl;
        panda::Lepton *lep2 = looseLeps.at(1);
        vObj2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());

        vpuppiUZ=vpuppiUW+vObj2; gt->puppiUZmag=vpuppiUZ.Pt(); gt->puppiUZphi=vpuppiUZ.Phi();
        vpfUZ=vpfUW+vObj2; gt->pfUZmag=vpfUZ.Pt(); gt->pfUZphi=vpfUZ.Phi();

        if (analysis->varyJES) {
          TLorentzVector vpfUZUp = vpfUWUp+vObj2; gt->pfUZmagUp = vpfUZUp.Pt();
          TLorentzVector vpfUZDown = vpfUWDown+vObj2; gt->pfUZmagDown = vpfUZDown.Pt();
        }

        vpuppiU = vpuppiUZ; vpfU = vpfUZ;
        gt->whichRecoil = 2;

      } else if ( gt->nLooseLep>1 && ( (looseLep1PdgId==11 && looseLep2PdgId==-13)
                                       || ( looseLep1PdgId==-11 && looseLep2PdgId==13 )
                                       || ( looseLep1PdgId==13 && looseLep2PdgId==-11 )
                                       || ( looseLep1PdgId==-13 && looseLep2PdgId==11 )
                                       ) ){
	
	std::cout<<"two W with opposite flavor"<<std::endl;                                                                                                                                                                       
	panda::Lepton *lep2 = looseLeps.at(1);
        vObj2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());

	vpuppiUWW=vpuppiUW+vObj2; gt->puppiUWWmag=vpuppiUWW.Pt(); gt->puppiUWWphi=vpuppiUWW.Phi();
        vpfUWW=vpfUW+vObj2; gt->pfUWWmag=vpfUWW.Pt(); gt->pfUWWphi=vpfUWW.Phi();

        if (analysis->varyJES) {
	  TLorentzVector vpfUWWUp = vpfUWUp+vObj2; gt->pfUWWmagUp = vpfUWWUp.Pt();
	  TLorentzVector vpfUWWDown = vpfUWDown+vObj2; gt->pfUWWmagDown = vpfUWWDown.Pt();
	}

        vpuppiU = vpuppiUWW; vpfU = vpfUWW;
        gt->whichRecoil = 3;
      } else {
        vpuppiU = vpuppiUW; vpfU = vpfUW;
        gt->whichRecoil = 1;
      }
    }
    if (gt->nLoosePhoton>0) {
      panda::Photon *pho = loosePhos.at(0);
      vObj1.SetPtEtaPhiM(pho->pt(),pho->eta(),pho->phi(),0.);

      vpuppiUA=vPuppiMET+vObj1; gt->puppiUAmag=vpuppiUA.Pt(); gt->puppiUAphi=vpuppiUA.Phi();
      vpfUA=vPFMET+vObj1; gt->pfUAmag=vpfUA.Pt(); gt->pfUAphi=vpfUA.Phi();

      if (analysis->varyJES) {
        TLorentzVector vpfUAUp = vpfUp+vObj1; gt->pfUAmagUp = vpfUAUp.Pt();
        TLorentzVector vpfUADown = vpfDown+vObj1; gt->pfUAmagDown = vpfUADown.Pt();
      }

      if (gt->nLooseLep==0) {
        vpuppiU = vpuppiUA; vpfU = vpfUA;
        gt->whichRecoil = -1;
      }
    }
    if (gt->nLooseLep==0 && gt->nLoosePhoton==0) {
      vpuppiU = vPuppiMET;
      vpfU = vPFMET;
      gt->whichRecoil = 0;
    }
    gt->puppiUmag = vpuppiU.Pt();
    gt->puppiUphi = vpuppiU.Phi();
    gt->pfUmag = vpfU.Pt();
    gt->pfUphi = vpfU.Phi();

    tr->TriggerEvent("recoils");
}

void PandaAnalyzer::HeavyFlavorCounting() 
{
  // Simple B and C counting stored in nB, nHF
  for (auto& gen : event.genParticles) {
    float pt = gen.pt();
    int pdgid = gen.pdgid;
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
  }
  // Gen B jet counting stored in nBGenJets
  for (auto &gen : event.ak4GenJets) {
    if (gen.pt() > 20 && std::abs(gen.eta()) < 2.4 && (gen.numB != 0 || abs(gen.pdgid)==5))
      gt->nBGenJets++;
  }
}

void PandaAnalyzer::GetMETSignificance()
{
  gt->pfmetsig = event.pfMet.significance;
  gt->puppimetsig = event.puppiMet.significance;

  tr->TriggerEvent("MET significance");
}

