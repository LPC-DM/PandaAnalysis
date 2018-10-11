from PandaCore.Tools.Misc import *
from re import sub

metTrigger='(trigger&1)!=0'
eleTrigger='(trigger&2)!=0'
phoTrigger='(trigger&4)!=0'

metFilter='metFilter==1'# && egmFilter==1'
presel = '!(nFatjet==1 && fj1Pt>200) && jet1Pt>100 && abs(jet1Eta)<2.4 && jet1IsTight==1 && nTau==0'

cuts = {
    'signal' : tAND(metFilter,tAND(presel,'nLooseMuon==0 && nLooseElectron==0 && nLoosePhoton==0 && pfmet>200 && dphipfmet>0.5')),
    'mn'    : tAND(metFilter,tAND(presel,'nLoosePhoton==0 && nLooseElectron==0 && nLooseMuon==1 && nTightMuon==1 && pfUWmag>200 && dphipfUW>0.5')),
    'en'    : tAND(metFilter,tAND(presel,'nLoosePhoton==0 && nLooseMuon==0 && nLooseElectron==1 && nTightElectron==1 && pfmet>50 && pfUWmag>200 && dphipfUW>0.5')),
    'zmm'    : tAND(metFilter,tAND(presel,'pfUZmag>200 && dphipfUZ>0.5 && nLooseElectron==0 && nLoosePhoton==0 && nLooseMuon==2 && nTightLep>0 && 60<diLepMass && diLepMass<120')),
    'zee'    : tAND(metFilter,tAND(presel,'pfUZmag>200 && dphipfUZ>0.5 && nLoosePhoton==0 && nLooseMuon==0 && nLooseElectron==2 && nTightLep>0 && 60<diLepMass && diLepMass<120')),
    'pho'    : tAND(metFilter,tAND(presel,'pfUAmag>200 && dphipfUA>0.5 && nLooseLep==0 && nLoosePhoton==1 && loosePho1IsTight==1 && fabs(loosePho1Eta)<1.4442')),
    }

for r in ['signal','zmm','zee','mn','en','pho']:
    cuts[r+'_0tag'] = tAND(cuts[r],'Sum$(jetPt>20 && jetCSV>0.8 && abs(jetEta)<2.4)==0')
    cuts[r+'_1tag'] = tAND(cuts[r],'Sum$(jetPt>20 && jetCSV>0.8 && abs(jetEta)<2.4)==1')
    cuts[r+'_2tag'] = tAND(cuts[r],'Sum$(jetPt>20 && jetCSV>0.8 && abs(jetEta)<2.4)==2')

weights = {
  'signal'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_metTrig',
  'control'        : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'pho'            : '%f*sf_pu*normalizedWeight*sf_ewkV*sf_qcdV*sf_pho*sf_phoTrig *sf_qcdV2j', # add the additional 2-jet kfactor
}

for x in ['mn','en','zee','zmm']:
	if 'en' in x or 'ee' in x:
	  weights[x] = tTIMES(weights['control'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['control'],'sf_metTrig')

for x in ['signal','mn','en','zee','zmm','pho']:
    weights[x+'_0tag'] = tTIMES(weights[x],'sf_btag0')
    weights[x+'_1tag'] = tTIMES(weights[x],'sf_btag1')
    weights[x+'_2tag'] = tTIMES(weights[x],'sf_btag2')

for x in ['signal','mn','en','zee','zmm','pho']:
    for y in ['0tag','1tag','2tag']:
        r = x+'_'+y
        for shift in ['BUp','BDown','MUp','MDown']:
            for cent in ['sf_btag']:
               weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,sub(cent+'1',cent+'1'+shift,sub(cent+'2',cent+'2'+shift,weights[r])))
