from PandaCore.Tools.Misc import *
from re import sub


metTrigger='(trigger&1)!=0'
eleTrigger='(trigger&2)!=0'
phoTrigger='(trigger&4)!=0'

cuts = {}
weights = {}
triggers = {}

baseline = 'metFilter==1 && nTau==0 && bosonpt>150 && Sum$(jetPt>30)>1 && !(fj1Pt>200) && jet1IsTight==1 '

#cuts for specific regions
cuts['signal'] = tAND(baseline,'nLooseLep==0 && nLooseElectron==0 && nLoosePhoton==0 && pfmet>200 && dphipfmet>0.8')
cuts['tmn'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && nTightMuon==1 && pfUWmag>200 && dphipfUW>0.8')
cuts['ten'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && nTightElectron==1 && pfUWmag>200 && dphipfUW>0.8 && pfmet>50')
cuts['wmn'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && nTightMuon==1 && pfUWmag>200 && dphipfUW>0.8')
cuts['wen'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && nTightElectron==1 && pfUWmag>200 && dphipfUW>0.8 && pfmet>50')
cuts['zmm'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==2 && nTightMuon==1 && pfUZmag>200 && dphipfUZ>0.8 && diLepMass>80 && diLepMass<100')
cuts['zee'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==2 && nTightElectron==1 && pfUZmag>200 && dphipfUZ>0.8 && diLepMass>80 && diLepMass<100')
cuts['pho'] = tAND(baseline,'(nLooseMuon+nLooseElectron+nTau)==0 && nLoosePhoton==1 && loosePho1IsTight==1 && pfUAmag>200 && dphipfUA>0.8')


for r in ['signal','tmn','ten','zmm','zee','wen','wmn','pho']:
#	cuts[r] = tAND(cuts[r],'(tanh(atanh((2*jetCSV[bosonjtidx[0]])-1)+atanh((2*jetCSV[bosonjtidx[1]])-1))+1)/2>0.48') 
        cuts[r+'_fail'] = tAND(cuts[r],'!(jetCSV[bosonjtidx[0]]>0.8484 && jetEta[bosonjtidx[0]]<2.4)||(jetCSV[bosonjtidx[1]]>0.8484 && jetEta[bosonjtidx[1]]<2.4)')
        cuts[r] = tAND(cuts[r],'(jetCSV[bosonjtidx[0]]>0.8484 && jetEta[bosonjtidx[0]]<2.4)||(jetCSV[bosonjtidx[1]]>0.8484 && jetEta[bosonjtidx[1]]<2.4)')

for r in ['signal','zmm','zee','wen','wmn','pho','signal_fail','wmn_fail','wen_fail','zmm_fail','zee_fail','pho','pho_fail']:
	if 'fail' in r:
          cuts[r] = tAND(cuts[r],'Sum$(jetCSV>0.8484 && jetEta<2.4)==0')
	else:
          cuts[r] = tAND(cuts[r],'Sum$(jetCSV>0.8484 && jetEta<2.4)<3')
#                         '(((jetCSV[bosonjtidx[0]]>0.8484 && jetEta[bosonjtidx[0]]<2.4)&&!(jetCSV[bosonjtidx[1]]>0.8484 && jetEta[bosonjtidx[1]]<2.4) && Sum$(jetCSV>0.8484 && jetEta<2.4)==1) ||
#                         (!(jetCSV[bosonjtidx[0]]>0.8484 && jetEta[bosonjtidx[0]]<2.4)&&(jetCSV[bosonjtidx[1]]>0.8484 && jetEta[bosonjtidx[1]]<2.4) && Sum$(jetCSV>0.8484 && jetEta<2.4)==1) ||
#                         ((jetCSV[bosonjtidx[0]]>0.8484 && jetEta[bosonjtidx[0]]<2.4)&&(jetCSV[bosonjtidx[1]]>0.8484 && jetEta[bosonjtidx[1]]<2.4) && Sum$(jetCSV>0.8484 && jetEta<2.4)==2))')

for r in ['tmn','ten','tmn_fail','ten_fail']:
	if 'fail' in r:
          cuts[r] = tAND(cuts[r],'Sum$(jetCSV>0.8484 && jetEta<2.4)>0')
	else:
          cuts[r] = tAND(cuts[r],'Sum$(jetCSV>0.8484 && jetEta<2.4)>2')
#                         '(((jetCSV[bosonjtidx[0]]>0.8484 && jetEta[bosonjtidx[0]]<2.4)&&!(jetCSV[bosonjtidx[1]]>0.8484 && jetEta[bosonjtidx[1]]<2.4) && Sum$(jetCSV>0.8484 && jetEta<2.4)>1) ||
#                         (!(jetCSV[bosonjtidx[0]]>0.8484 && jetEta[bosonjtidx[0]]<2.4)&&(jetCSV[bosonjtidx[1]]>0.8484 && jetEta[bosonjtidx[1]]<2.4) && Sum$(jetCSV>0.8484 && jetEta<2.4)>1) ||
#                         ((jetCSV[bosonjtidx[0]]>0.8484 && jetEta[bosonjtidx[0]]<2.4)&&(jetCSV[bosonjtidx[1]]>0.8484 && jetEta[bosonjtidx[1]]<2.4) && Sum$(jetCSV>0.8484 && jetEta<2.4)>2))')

#weights for specific regions
weights = {
  'signal'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_metTrig',
  'control'        : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'pho'            : '%f*sf_pu*normalizedWeight*sf_ewkV*sf_qcdV*sf_pho*sf_phoTrig *sf_qcdV2j', # add the additional 2-jet kfactor
}
weights['signal_fail'] = weights['signal']
weights['pho_fail'] = weights['pho']

for x in ['tmn','ten','tmn_fail','ten_fail']:
	if 'e' in x:
	  weights[x] = tTIMES(weights['control'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['control'],'sf_metTrig')
for x in ['wmn','wen','wmn_fail','wen_fail']:
	if 'e' in x:
	  weights[x] = tTIMES(weights['control'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['control'],'sf_metTrig')
for x in ['zmm','zee','zmm_fail','zee_fail']:
	if 'e' in x:
	  weights[x] = tTIMES(weights['control'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['control'],'sf_metTrig')

for x in ['signal','signal_fail','wmn','wen','zmm','zee','wmn_fail','wen_fail','zmm_fail','zee_fail','pho','pho_fail']:
	if 'fail' in x:
	  weights[x] = tTIMES(weights[x],'sf_btag0')
        else:
	  weights[x] = tTIMES(weights[x],'sf_btag2')
   
for x in ['tmn','ten','tmn_fail','ten_fail']:
	if 'fail' in x:
	  weights[x] = tTIMES(weights[x],'sf_btag1')
        else:
	  weights[x] = tTIMES(weights[x],'sf_btag3')


for r in ['signal','tmn','ten','wmn','wen','zmm','zee','signal_fail','tmn_fail','ten_fail','wmn_fail','wen_fail','zmm_fail','zee_fail','pho','pho_fail']:
    for shift in ['BUp','BDown','MUp','MDown']:
        for cent in ['sf_btag']:
            if 'btag0' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,weights[r])
            if 'btag1' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'1',cent+'1'+shift,weights[r])
            if 'btag2' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'2',cent+'2'+shift,weights[r])
            if 'btag3' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'3',cent+'3'+shift,weights[r])
