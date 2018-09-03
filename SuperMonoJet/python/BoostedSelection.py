from PandaCore.Tools.Misc import *
from re import sub

metTrigger='(trigger&1)!=0'
eleTrigger='(trigger&2)!=0'
phoTrigger='(trigger&4)!=0'

#metFilter='metFilter==1 && egmFilter==1'
metFilter='metFilter==1'
#presel = 'nFatjet==1 && fj1Pt>250 && nTau==0 && Sum$(jetPt>30 && jetIso)<2'
presel = 'nFatjet>0 && fj1Pt>200 && nTau==0  && jet1IsTight==1'

cuts = {
 'signal' : tAND(metFilter,tAND(presel,'nLooseLep==0 && nLooseElectron==0 && nLoosePhoton==0 && pfmet>250 && dphipfmet>0.8')), 
 'mn'      : tAND(metFilter,tAND(presel,'nLoosePhoton==0 && nTau==0 && nLooseLep==1 && nTightMuon==1 && pfUWmag>250 && dphipfUW>0.8 && mT<160')),
 'en'      : tAND(metFilter,tAND(presel,'nLoosePhoton==0 && nTau==0 && nLooseLep==1 && nTightElectron==1 && pfmet>50 && pfUWmag>250 && dphipfUW>0.4 && mT<160')),
 'zmm'    : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.8 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && nLooseMuon==2 && nTightLep>0 && 60<diLepMass && diLepMass<120')),
 'zee'    : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.8 && nLooseMuon==0 && nLoosePhoton==0 && nTau==0 && nLooseElectron==2 && nTightLep>0 && 60<diLepMass && diLepMass<120')),
 'tme'    : tAND(metFilter,tAND(presel,'pfUWWmag>250 && dphipfUWW>0.8 && nLoosePhoton==0 && nLooseElectron==1 && nTightMuon==1')),# && Sum$(jetPt>30 && jetCSV>0.54 && abs(jetEta)<2.4 && jetIso)==2')),
 'tem'    : tAND(metFilter,tAND(presel,'pfUWWmag>250 && dphipfUWW>0.8 && nLoosePhoton==0 && nLooseMuon==1 && nTightElectron==1 && nTightMuon==0')),
 'pho'    : tAND(metFilter,tAND(presel,'(nLooseMuon+nLooseElectron+nTau)==0 && nLoosePhoton==1 && loosePho1IsTight==1 && pfUAmag>250 && dphipfUA>0.8'))

}

for r in ['mn','en']:
	cuts['w'+r] = tAND(cuts[r],'Sum$(jetPt>30 && jetCSV>0.54 && abs(jetEta)<2.4 && jetIso)==0')
	cuts['t'+r]     = tAND(cuts[r],'Sum$(jetPt>30 && jetCSV>0.54 && abs(jetEta)<2.4 && jetIso)>0')
#        cuts['t'+r]     = tAND(cuts[r],'isojetNBtags>-1')
#	cuts['t'+r]     =  cuts[r]
for r in ['signal','zmm','zee','pho']:
        cuts[r] = tAND(cuts[r],'Sum$(jetPt>30 && jetCSV>0.54 && abs(jetEta)<2.4 && jetIso)==0')

for r in ['signal','wmn','tmn','wen','ten','zmm','zee','pho']:#,'tem','tme']:
	cuts[r+'_fail'] = tAND(cuts[r],'fj1DoubleCSV<=0.9')
	cuts[r] = tAND(cuts[r],'fj1DoubleCSV>0.9')

weights = {
  'signal'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_metTrig*sf_btag0',
  'top'            : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_btag1',
  'opposite'       : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'w'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_btag0',
  'z'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_btag0',
  'pho'         : '%f*sf_pu*normalizedWeight*sf_ewkV*sf_qcdV*sf_pho*sf_phoTrig *sf_qcdV2j*sf_btag0', # add the additional 2-jet kfactor
}
weights['qcd'] = weights['signal']
weights['signal_fail'] = weights['signal']
weights['pho_fail'] = weights['pho']

for x in ['tmn','ten','tmn_fail','ten_fail']:
        if 'e' in x:
          weights[x] = tTIMES(weights['top'],'sf_eleTrig')
        else:
          weights[x] = tTIMES(weights['top'],'sf_metTrig')
for x in ['tme','tem','tme_fail','tem_fail']:
	if 'em' in x:
	  weights[x] = tTIMES(weights['opposite'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['opposite'],'sf_metTrig')
for x in ['wmn','wen','wmn_fail','wen_fail']:
	if 'e' in x:
	  weights[x] = tTIMES(weights['w'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['w'],'sf_metTrig')
for x in ['zmm','zee','zmm_fail','zee_fail']:
	if 'e' in x:
	  weights[x] = tTIMES(weights['z'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['z'],'sf_metTrig')

for r in ['signal','top','w','tmn','ten','tme','tem','wmn','wen','zmm','zee','z','signal_fail','tmn_fail','ten_fail','tme_fail','tem_fail','wmn_fail','wen_fail','zmm_fail','zee_fail','pho','pho_fail']:
  for shift in ['BUp','BDown','MUp','MDown']:
    for cent in ['sf_btag']:
      weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,sub(cent+'1',cent+'1'+shift,weights[r]))
