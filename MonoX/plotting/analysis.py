#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse,sys

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir#.replace('0_4','0_4_egfix')
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default='.')
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--tt',metavar='tt',type=str,default='')
parser.add_argument('--bdtcut',type=float,default=None)
parser.add_argument('--masscut',type=float,default=None)
args = parser.parse_args()
lumi = 35800.
blind=True
region = args.region
sname = argv[0]

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
#import PandaAnalysis.Monotop.MonojetSelection as sel
#import PandaAnalysis.Monotop.LooseSelection as sel
#import PandaAnalysis.Monotop.TightSelection as sel
#import PandaAnalysis.Monotop.OneFatJetSelection as sel
#import PandaAnalysis.Monotop.CombinedBVetoSelection as sel
import PandaAnalysis.MonoX.Selection_bb as sel
#import PandaAnalysis.MonoX.MonoXSelection as sel                                                                                                                                                         
from PandaCore.Drawers.plot_utility import *

#reload(sys.modules['plot_utility'])

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)
if args.bdtcut:
    cut = tAND(cut,'top_ecf_bdt>%f'%args.bdtcut)
if args.masscut:
    cut = tAND(cut,'fj1MSD>%f'%args.masscut)

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(0.4)
if 'qcd' in region:
    plot.FixRatio(1)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.cut = cut
plot.SetEvtNum("eventNumber")
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights[region]%lumi
plot.mc_weight = weight

if args.bdtcut:
    plot.AddPlotLabel('BDT > %.2f'%args.bdtcut,.18,.7,False,42,.04)
if args.masscut:
    plot.AddPlotLabel('%i < m_{SD} < 210 GeV'%(int(args.masscut)),.18,.7,False,42,.04)

PInfo('cut',plot.cut)
PInfo('weight',plot.mc_weight)

#plot.add_systematic('QCD scale','scaleUp','scaleDown',root.kRed+2)
#plot.add_systematic('PDF','pdfUp','pdfDown',root.kBlue+2)

### DEFINE PROCESSES ###
zjets         = Process('Z+jets',root.kZjets)
wjets         = Process('W+jets',root.kWjets)
diboson       = Process('Diboson',root.kDiboson)
ttbar         = Process('t#bar{t}',root.kTTbar)
ttg           = Process('t#bar{t}#gamma',root.kTTbar)
singletop     = Process('Single t',root.kST)
singletopg    = Process('t#gamma',root.kST)
qcd           = Process("QCD",root.kQCD)
gjets         = Process('#gamma+jets',root.kGjets)
data          = Process("Data",root.kData)
signal        = Process('m_{V}=1.75 TeV, m_{#chi}=1 GeV',root.kSignal)
#processes = [qcd,diboson,singletop,ttbar,wewk,zewk,wjets,zjets]
#processes = [qcd,diboson,singletop,wjets,ttbar,zjets]
processes = [qcd,diboson,singletop,wjets,ttbar,zjets]

if 'qcd' in region:
    processes = [diboson,singletop,wjets,ttbar,zjets,qcd]

### ASSIGN FILES TO PROCESSES ###
if 'signal' in region or 'qcd' in region:
    zjets.add_file(baseDir+'ZtoNuNu.root')
    signal.add_file(baseDir+'Vector_MonoTop_NLO_Mphi-1750_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph.root')
else:
    zjets.add_file(baseDir+'ZJets.root')
    #zjets.add_file(baseDir+'ZJets_nlo.root')
wjets.add_file(baseDir+'WJets.root')
diboson.add_file(baseDir+'Diboson.root')
ttbar.add_file(baseDir+'TTbar%s.root'%(args.tt));
singletop.add_file(baseDir+'SingleTop.root')
#ttg.add_file(baseDir+'TTbar_Photon.root');
#singletopg.add_file(baseDir+'SingleTop_tG.root')
if 'pho' in region:
    #processes = [qcd,singletopg,ttg,gjets]
    processes = [qcd,gjets]
    gjets.add_file(baseDir+'GJets.root')
    qcd.add_file(baseDir+'SinglePhoton.root')
    qcd.additional_cut = sel.triggers['pho']
    qcd.use_common_weight = False
    qcd.additional_weight = 'sf_phoPurity'
else:
    qcd.add_file(baseDir+'QCD.root')

if any([x in region for x in ['singlemuonw','singleelectronw']]):
    processes = [qcd,diboson,singletop,zjets,ttbar,wjets,]
if any([x in region for x in ['singlemuontop','singleelectrontop']]):
    processes = [qcd,diboson,singletop,zjets,wjets,ttbar]
if any([x in region for x in ['signal','muon','qcd']]):
    data.additional_cut = sel.triggers['met']
    data.add_file(dataDir+'MET.root')
    lep='#mu'
elif 'electron' in region:
    if 'di' in region:
        data.additional_cut = tOR(sel.triggers['ele'],sel.triggers['pho'])
    else:
        data.additional_cut = sel.triggers['ele']
    data.add_file(dataDir+'SingleElectron.root')
    lep='e'
elif region=='photon':
    data.additional_cut = sel.triggers['pho']
    data.add_file(dataDir+'SinglePhoton.root')
elif region=='ZmmCR':
    data.additional_cut = sel.triggers['met']
    print dataDir+'MET.root'
    data.add_file(dataDir+'MET.root')
    lep='#mu'


processes.append(data)

for p in processes:
    plot.add_process(p)

recoilBins = [250,280,310,350,400,450,600,1000]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###
if 'signal' in region or 'qcd' in region:
    recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","Events/GeV")
elif any([x in region for x in ['singlemuonw','singleelectronw','singlemuontop','singleelectrontop','singlemuon','singleelectron']]):
    recoil=VDistribution("pfUWmag",recoilBins,"PF U(%s) [GeV]"%(lep),"Events/GeV")
    plot.add_distribution(FDistribution('looseLep1Pt',0,1000,20,'Leading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('looseLep1Eta',-2.5,2.5,20,'Leading %s #eta'%lep,'Events/bin'))
elif any([x in region for x in ['dielectron','dimuon']]):
    recoil=VDistribution("pfUZmag",recoilBins,"PF U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
    plot.add_distribution(FDistribution('diLepMass',60,120,20,'m_{ll} [GeV]','Events/3 GeV'))
    plot.add_distribution(FDistribution('looseLep1Pt',0,1000,20,'Leading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('looseLep1Eta',-2.5,2.5,20,'Leading %s #eta'%lep,'Events/bin'))
    plot.add_distribution(FDistribution('looseLep2Pt',0,1000,20,'Subleading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('looseLep2Eta',-2.5,2.5,20,'Subleading %s #eta'%lep,'Events/bin'))
elif region=='photon':
    recoil=VDistribution("pfUAmag",recoilBins,"PF U(#gamma) [GeV]","Events/GeV")
    plot.add_distribution(FDistribution('loosePho1Pt',0,1000,20,'Leading #gamma p_{T} [GeV]','Events/40 GeV'))
    plot.add_distribution(FDistribution('loosePho1Eta',-2.5,2.5,20,'Leading #gamma #eta','Events/bin'))
elif region=='ZmmCR':
    recoil=VDistribution("pfUZmag",recoilBins,"PF U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
    plot.add_distribution(FDistribution('diLepMass',60,120,20,'m_{ll} [GeV]','Events/3 GeV'))
    #plot.add_distribution(FDistribution('looseLep1Pt',0,1000,20,'Leading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    #plot.add_distribution(FDistribution('looseLep1Eta',-2.5,2.5,20,'Leading %s #eta'%lep,'Events/bin'))
    #plot.add_distribution(FDistribution('looseLep2Pt',0,1000,20,'Subleading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    #plot.add_distribution(FDistribution('looseLep2Eta',-2.5,2.5,20,'Subleading %s #eta'%lep,'Events/bin'))


#recoil.calc_chi2 = True
#plot.add_distribution(recoil)

#plot.add_distribution(FDistribution('nJet',0.5,8.5,8,'N_{jet}','Events'))
#plot.add_distribution(FDistribution('npv',0,45,45,'N_{PV}','Events'))
#plot.add_distribution(FDistribution('fj1MSD',50,250,10,'fatjet m_{SD} [GeV]','Events'))
#plot.add_distribution(FDistribution('fj1Pt',200,1000,20,'fatjet p_{T} [GeV]','Events'))
#plot.add_distribution(FDistribution('top_ecf_bdt',-1,1,20,'Top BDT','Events'))
#plot.add_distribution(FDistribution('fj1MaxCSV',0,1,20,'fatjet max CSV','Events'))
#plot.add_distribution(FDistribution('fj1Tau32',0,1,20,'fatjet #tau_{32}','Events'))
#plot.add_distribution(FDistribution('fj1Tau32SD',0,1,20,'fatjet #tau_{32}^{SD}','Events'))
#plot.add_distribution(FDistribution('jet1CSV',0,1,20,'jet 1 CSV','Events',filename='jet1CSV'))
#plot.add_distribution(FDistribution('dphipfmet',0,3.14,20,'min#Delta#phi(jet,E_{T}^{miss})','Events'))
#plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
if args.bdtcut:
    region += ('_bdt%.2f'%(args.bdtcut)).replace('.','p').replace('-','m')
if args.masscut:
    region += ('_mass%i'%(int(args.masscut)))
plot.draw_all(args.outdir+'/'+region+'_')
