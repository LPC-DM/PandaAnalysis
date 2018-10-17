#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse,sys

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default='.')
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--analysis',metavar='analysis',type=str,default=None)
parser.add_argument('--tt',metavar='tt',type=str,default='')
parser.add_argument('--bdtcut',type=float,default=None)
parser.add_argument('--masscut1',type=float,default=None)
parser.add_argument('--masscut2',type=float,default=None)
#parser.add_argument('--fromlimit',type=bool,default=False)
parser.add_argument('--fromlimit',action='store_true',help='enable direct plotting from fitting ntuple')
args = parser.parse_args()
lumi = 35800.
blind=False
SIGNAL=False
region = args.region
processes = []
sname = argv[0]

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
from PandaCore.Drawers.plot_utility import *
import PandaAnalysis.SuperMonoJet.BoostedSelection as sel

if args.analysis == "boosted":
    import PandaAnalysis.SuperMonoJet.BoostedSelection as sel                            
if args.analysis == "resolved":
    import PandaAnalysis.SuperMonoJet.ResolvedSelection as sel
if args.analysis == "monojet":
    import PandaAnalysis.SuperMonoJet.MonoJetSelection as sel

### SET GLOBAL VARIABLES ###
if not args.fromlimit:
    baseDir = getenv('PANDA_FLATDIR')+'/'
    dataDir = baseDir
    ### DEFINE REGIONS ###
    cut = tAND(sel.cuts[args.region],args.cut)
else:
    baseDir = getenv('PANDA_FLATDIR')+'/limits/'
#    baseDir = getenv('PANDA_FLATDIR')
    dataDir = baseDir
    cut = '1==1'

print "PLOTTER: input directory is:", baseDir

if args.bdtcut:
    cut = tAND(cut,'top_ecf_bdt>%f'%args.bdtcut)
if args.masscut1 and args.masscut2:
    cut = tAND(cut,'fj1MSD>%f'%args.masscut1)
    cut = tAND(cut,'fj1MSD<%f'%args.masscut2)
if args.masscut1:
    cut = tAND(cut,'fj1MSD>%f'%args.masscut1)

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(0.4)
if 'qcd' in region:
    plot.FixRatio(1)
plot.SetTDRStyle()
plot.InitLegend()
plot.cut = cut
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.SetEvtNum("eventNumber")
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

if args.bdtcut:
    plot.AddPlotLabel('BDT > %.2f'%args.bdtcut,.18,.7,False,42,.04)

if args.masscut1 and args.masscut2:
    plot.AddPlotLabel('%i GeV < m_{SD} < %i GeV'%(int(args.masscut1),int(args.masscut2)),.18,.7,False,42,.04)

def f(x):
    return dataDir + 'fittingForest_' + x + '.root'

def fillProcess(region, plotFromLimits):
    region=region
    global processes
    temp_processes = []

    limitLabel = None
    znunu         = Process('Z(#nu#nu)+jets',root.kZjets,'Zvv_'+region,root.kCyan-9, plotFromLimits)
    zjets         = Process('Z+jets',root.kZjets,'Zll_'+region,root.kCyan-9, plotFromLimits)
    wjets         = Process('W+jets',root.kWjets,'Wlv_'+region,root.kGreen-10, plotFromLimits)
    diboson       = Process('Diboson',root.kDiboson,'Diboson_'+region,root.kYellow-9, plotFromLimits)
    ttbar         = Process('t#bar{t}',root.kTTbar,'ttbar_'+region,root.kOrange-4, plotFromLimits)
    singletop     = Process('Single t',root.kST,'ST_'+region,root.kRed-9, plotFromLimits)
    qcd           = Process('QCD',root.kQCD,'QCD_'+region,root.kMagenta-10, plotFromLimits)
    gjets         = Process('#gamma+jets',root.kGjets,None,root.kBlue, plotFromLimits)
    data          = Process("Data",root.kData,'Data_'+region, root.kBlack, plotFromLimits)
    signal        = Process('m_{Zp}=1.0 TeV,m_{h}=90 GeV, m_{#chi}=300 GeV',root.kSignal,None, root.kBlack, plotFromLimits)
    
    if 'zee' in region or 'zmm' in region:
        temp_processes = [diboson,ttbar,zjets]
        print 'adding diboson, ttbar, and zjets to temp. Temp type =',type(temp_processes)
    if 'wen' in region or 'wmn' in region:
        temp_processes = [qcd,diboson,singletop,zjets,ttbar,wjets]
    if 'signal' in region or 'qcd' in region:
        temp_processes = [qcd,zjets,singletop,ttbar,diboson,wjets,znunu]
    if not blind:
        temp_processes.append(data)

    print 'Type of temp = ',type(temp_processes)

    if not plotFromLimits:
        zjets.add_file(baseDir+'ZJets.root')
        diboson.add_file(baseDir+'Diboson.root')
        ttbar.add_file(baseDir+'TTbar%s.root'%(args.tt)); print 'TTbar%s.root'%(args.tt)
        singletop.add_file(baseDir+'SingleTop.root')
        wjets.add_file(baseDir+'WJets.root')
        
        if 'signal' in region:
            SIGNAL=True
            signal.add_file(baseDir+'BBbarDM_MZprime-1000_Mhs-90_Mchi-300.root')
            temp_processes.append(signal)
        if 'signal' in region or 'qcd' in region:
            znunu.add_file(baseDir+'ZtoNuNu.root')
        if 'signal' in region or 'wen' in region or 'wmn' in region or 'ten' in region or 'tmn' in region:
            qcd.add_file(baseDir+'QCD.root')
    
        if any([x in region for x in ['signal','wmn','zmm','tmn','qcd','tme']]):
            if not blind:
                data.add_file(dataDir+'MET.root')
                data.additional_cut = sel.metTrigger
        elif any([x in region for x in ['wen','zee','ten','tem']]):
            data.additional_cut = sel.eleTrigger
            data.add_file(dataDir+'SingleElectron.root')
        elif 'pho' in region:
            temp_processes = [qcd,gjets]
            gjets.add_file(baseDir+'GJets.root')
            qcd.add_file(baseDir+'SinglePhoton.root')
            qcd.additional_cut = sel.phoTrigger
            qcd.use_common_weight = False
            qcd.additional_weight = 'sf_phoPurity'
            data.additional_cut = sel.phoTrigger
            data.add_file(dataDir+'SinglePhoton.root')
    print 'Type of temp = ',type(temp_processes)
    for p in temp_processes:
        print p.name
    processes = temp_processes

def normalPlotting(region):
    print 'Plotting from PandaAnalyzed  ntuple: ',region, ' region'
    global processes
    fillProcess(region, False)

    region=region
    weight = sel.weights[region]%lumi
    plot.mc_weight = weight

    PInfo('cut',plot.cut)
    PInfo('weight',plot.mc_weight)

    for p in processes:
        plot.add_process(p)

    if args.analysis == "monojet" or args.analysis == "resolved":
        recoilBins = [250,280,310,340,370,400,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250,1400]
    if args.analysis == "boosted":
        recoilBins = [250,270,350,475,1000]
    fatjetBins = [25,75,100,150,600]
    nRecoilBins = len(recoilBins)-1

    if any([x in region for x in ['signal','wmn','zmm','tmn','qcd','tme']]):
        lep='#mu'
    elif any([x in region for x in ['wen','zee','ten','tem']]):
        lep='e'

    ### CHOOSE DISTRIBUTIONS, LABELS ###
    if 'signal' in region or 'qcd' in region:
        recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","Events/GeV")
        plot.add_distribution(FDistribution('dphipfmet',0,3.14,10,'min#Delta#phi(AK4 jet,E_{T}^{miss})','Events'))
    elif any([x in region for x in ['wen','wmn','ten','tmn']]):
        recoil=VDistribution("pfUWmag",recoilBins,"PF U(%s) [GeV]"%(lep),"Events/GeV")
        #plot.add_distribution(FDistribution('mT',0,500,25,'Transverse Mass of W [GeV]','Events'))
        if not lep=="e":
            plot.add_distribution(FDistribution('muonPt[0]',0,400,15,'Leading %s p_{T} [GeV]'%lep,'Events/25 GeV'))
            plot.add_distribution(FDistribution('muonEta[0]',-2.5,2.5,10,'%s #eta'%lep,'Events'))
        else:
            plot.add_distribution(FDistribution('electronPt[0]',0,400,15,'Leading %s p_{T} [GeV]'%lep,'Events/25 GeV'))
            plot.add_distribution(FDistribution('electronEta[0]',-2.5,2.5,10,'%s #eta'%lep,'Events'))
        plot.add_distribution(FDistribution('dphipfUW',0,3.14,10,'min#Delta#phi(AK4 jet,E_{T}^{miss})','Events'))
        plot.add_distribution(FDistribution('jetNMBtags',0,5,5,'nmbtag jets','Events'))
        plot.add_distribution(FDistribution('jetNBtags',0,5,5,'nbtag jets','Events'))
        
    elif any([x in region for x in ['zee','zmm']]):
        recoil=VDistribution("pfUZmag",recoilBins,"PF U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
        #plot.add_distribution(FDistribution('diLepMass',60,120,20,'m_{ll} [GeV]','Events/3 GeV'))
        if not lep=="e":
            plot.add_distribution(FDistribution('muonPt[0]',0,400,15,'Leading %s p_{T} [GeV]'%lep,'Events/25 GeV'))
            plot.add_distribution(FDistribution('muonEta[0]',-2.5,2.5,10,'%s #eta'%lep,'Events'))
        else:
            plot.add_distribution(FDistribution('electronPt[0]',0,400,15,'Leading %s p_{T} [GeV]'%lep,'Events/25 GeV'))
            #plot.add_distribution(FDistribution('electronEta[0]',-2.5,2.5,10,'%s #eta'%lep,'Events'))
        plot.add_distribution(FDistribution('dphipfUZ',0,3.14,10,'min#Delta#phi(AK4 jet,E_{T}^{miss})','Events'))

    elif any([x in region for x in ['tem','tme']]):
        recoil=VDistribution("pfUWWmag",recoilBins,"PF U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
        #plot.add_distribution(FDistribution('diLepMass',60,120,20,'m_{ll} [GeV]','Events/3 GeV'))
        if not lep=="e":
            plot.add_distribution(FDistribution('muonPt[0]',0,400,15,'Leading %s p_{T} [GeV]'%lep,'Events/25 GeV'))
            plot.add_distribution(FDistribution('muonEta[0]',-2.5,2.5,10,'%s #eta'%lep,'Events'))
        else:
            plot.add_distribution(FDistribution('electronPt[0]',0,400,15,'Leading %s p_{T} [GeV]'%lep,'Events/25 GeV'))
            plot.add_distribution(FDistribution('electronEta[0]',-2.5,2.5,10,'%s #eta'%lep,'Events'))
        plot.add_distribution(FDistribution('dphipfUWW',0,3.14,10,'min#Delta#phi(AK4 jet,E_{T}^{miss})','Events'))

    elif 'pho' in region:
        recoil=VDistribution("pfUAmag",recoilBins,"PF U(#gamma) [GeV]","Events/GeV")
        plot.add_distribution(FDistribution('loosePho1Pt',0,1000,20,'Leading #gamma p_{T} [GeV]','Events/50 GeV'))
        plot.add_distribution(FDistribution('loosePho1Eta',-2.5,2.5,10,'Leading #gamma #eta','Events/bin'))
        plot.add_distribution(FDistribution('dphipfUA',0,3.14,10,'min#Delta#phi(jet,E_{T}^{miss})','Events'))
  
    plot.add_distribution(recoil)

    if args.analysis == "boosted":
        plot.add_distribution(FDistribution('fj1DoubleCSV',0,1,20,'fatjet 1 DoubleCSV','Events'))
   
    plot.add_distribution(FDistribution('fj1Pt',200,700,15,'fatjet p_{T} [GeV]','Events/25 GeV'))
    plot.add_distribution(FDistribution('fj1Eta',-2.5,2.5,10,'fatjet #eta [GeV]','Events'))

    plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

    system('mkdir -p %s/%s/%s' %(args.outdir,args.analysis,region))
    plot.draw_all(args.outdir+'/'+args.analysis+'/'+args.region+'/')

def fromLimit(region):
    print 'Plotting from fitting ntuple: ',region, ' region'
    fillProcess(region, True)
    global processes

    region=region
    
    for p in processes:
        p.add_file(f(region))
        p.additional_weight='weight'
        plot.add_process(p)
    
    if args.analysis == "monojet" or args.analysis == "resolved":
        recoilBins = [250,280,310,340,370,400,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250,1400]
    if args.analysis == "boosted":
        recoilBins = [250,280,310,350,400,450,600,1000]
    nRecoilBins = len(recoilBins)-1
    plot.add_distribution(FDistribution('fjpt',200,700,15,'fatjet p_{T} [GeV]','Events/25 GeV'))
    plot.add_distribution(FDistribution('fjmass',0,600,20,'fatjet m_{SD} [GeV]','Events'))
    if args.analysis == "boosted":
        plot.add_distribution(FDistribution('doubleb',0,1,20,'fatjet 1 DoubleCSV','Events'))
    plot.add_distribution(FDistribution('n2',0,0.5,10,'n2','Events'))
    plot.add_distribution(FDistribution('n2ddt56',0,0.5,10,'n2ddt56','Events'))
    plot.add_distribution(FDistribution('n2ddt53',0,0.5,10,'n2ddt53','Events'))
    recoil=VDistribution("met",recoilBins,"PF MET [GeV]","Events/GeV")
    plot.add_distribution(recoil)
    
    system('mkdir -p %s/%s/fromLimits/%s' %(args.outdir,args.analysis,region))
    plot.draw_all(args.outdir+'/'+args.analysis+'/fromLimits/'+args.region+'/')

if not args.fromlimit:
    normalPlotting(region)
else:
    fromLimit(region)
