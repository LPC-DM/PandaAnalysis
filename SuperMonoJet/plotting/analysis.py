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
    import PandaAnalysis.SuperMonoJet.MonojetSelection as sel

### SET GLOBAL VARIABLES ###
if not args.fromlimit:
    baseDir = getenv('PANDA_FLATDIR')+'/'
    dataDir = baseDir
    ### DEFINE REGIONS ###
    cut = tAND(sel.cuts[args.region],args.cut)
else:
    baseDir = getenv('PANDA_FLATDIR')+'/limits/'
    dataDir = baseDir
    cut = '1==1'

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

def normalPlotting(region):

    print 'Plotting from PandaAnalyzed  ntuple: ',region, ' region'
    region=region
    weight = sel.weights[region]%lumi
    plot.mc_weight = weight

    PInfo('cut',plot.cut)
    PInfo('weight',plot.mc_weight)

    #plot.add_systematic('QCD scale','scaleUp','scaleDown',root.kRed+2)
    #plot.add_systematic('PDF','pdfUp','pdfDown',root.kBlue+2)

    ### DEFINE PROCESSES ###
    znunu         = Process('Z(#nu#nu)+jets',root.kZjets,None,root.kCyan-9)
    zjets         = Process('Z+jets',root.kZjets,None,root.kCyan-9)
    wjets         = Process('W+jets',root.kWjets,None,root.kGreen-10)
    diboson       = Process('Diboson',root.kDiboson,None,root.kYellow-9)
    ttbar         = Process('t#bar{t}',root.kTTbar,None,root.kOrange-4)
    #ttg           = Process('t#bar{t}#gamma',root.kTTbar,root.kOrange-4)
    singletop     = Process('Single t',root.kST,None,root.kRed-9)
    #singletopg    = Process('t#gamma',root.kST,root.kRed-9)
    qcd           = Process("QCD",root.kQCD,None,root.kMagenta-10)
    #gjets         = Process('#gamma+jets',root.kGjets)
    data          = Process("Data",root.kData)
    #signal        = Process('m_{V}=1.75 TeV, m_{#chi}=1 GeV',root.kSignal)
    processes = [qcd,diboson,singletop,wjets,ttbar,zjets]

    if 'qcd' in region:
        processes = [diboson,singletop,wjets,ttbar,zjets,qcd]

    ### ASSIGN FILES TO PROCESSES ###
    if 'qcd' in region:
        processes = [diboson,singletop,wjets,ttbar,zjets,qcd]
    if 'zee' or 'zmm' in region:
        processes = [qcd,diboson,singletop,ttbar,wjets,zjets]
    if 'wen' or 'wmn' in region:
        processes = [qcd,diboson,singletop,zjets,ttbar,wjets]
    ### ASSIGN FILES TO PROCESSES ###
    if 'signal' in region or 'qcd' in region:
        processes = [qcd,zjets,singletop,ttbar,diboson,wjets,znunu]
        znunu.add_file(baseDir+'ZtoNuNu.root')

    zjets.add_file(baseDir+'ZJets.root')
    wjets.add_file(baseDir+'WJets.root')
    diboson.add_file(baseDir+'Diboson.root')
    ttbar.add_file(baseDir+'TTbar%s.root'%(args.tt)); print 'TTbar%s.root'%(args.tt)
    singletop.add_file(baseDir+'SingleTop.root')
    if 'pho' in region:
       #processes = [qcd,singletopg,ttg,gjets]
        processes = [qcd,gjets]
        gjets.add_file(baseDir+'GJets.root')
        qcd.add_file(baseDir+'SinglePhoton.root')
        qcd.additional_cut = sel.phoTrigger
        qcd.use_common_weight = False
        qcd.additional_weight = 'sf_phoPurity'
    else:
        qcd.add_file(baseDir+'QCD.root')

    if any([x in region for x in ['signal','wmn','zmm','tmn','qcd']]):
        if not blind:
            data.add_file(dataDir+'MET.root')

        data.additional_cut = sel.metTrigger
        lep='#mu'
    elif any([x in region for x in ['wen','zee','ten']]):
        data.additional_cut = sel.eleTrigger
        data.add_file(dataDir+'SingleElectron.root')
        lep='e'
    elif region=='pho':
        data.additional_cut = sel.phoTrigger
        data.add_file(dataDir+'SinglePhoton.root')

    if not blind:
        processes.append(data)
    if SIGNAL:
        processes.append(signal)
    #    #processes.append(signal1)

    for p in processes:
       #print "processess considered -> ", p
        plot.add_process(p)

    recoilBins = [250,270,350,475,1000]
    fatjetBins = [25,75,100,150,600]
    nRecoilBins = len(recoilBins)-1

    ### CHOOSE DISTRIBUTIONS, LABELS ###
    if 'signal' in region or 'qcd' in region:
        recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","Events/GeV")

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
        
    elif any([x in region for x in ['zee','zmm']]):
        recoil=VDistribution("pfUZmag",recoilBins,"PF U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
        #plot.add_distribution(FDistribution('diLepMass',60,120,20,'m_{ll} [GeV]','Events/3 GeV'))
        if not lep=="e":
            plot.add_distribution(FDistribution('muonPt[0]',0,400,15,'Leading %s p_{T} [GeV]'%lep,'Events/25 GeV'))
            plot.add_distribution(FDistribution('muonEta[0]',-2.5,2.5,10,'%s #eta'%lep,'Events'))
        else:
            plot.add_distribution(FDistribution('electronPt[0]',0,400,15,'Leading %s p_{T} [GeV]'%lep,'Events/25 GeV'))
            plot.add_distribution(FDistribution('electronEta[0]',-2.5,2.5,10,'%s #eta'%lep,'Events'))
        plot.add_distribution(FDistribution('dphipfUZ',0,3.14,10,'min#Delta#phi(AK4 jet,E_{T}^{miss})','Events'))

    elif region=='pho':
        recoil=VDistribution("pfUAmag",recoilBins,"PF U(#gamma) [GeV]","Events/GeV")
        plot.add_distribution(FDistribution('loosePho1Pt',0,1000,20,'Leading #gamma p_{T} [GeV]','Events/50 GeV'))
        plot.add_distribution(FDistribution('loosePho1Eta',-2.5,2.5,10,'Leading #gamma #eta','Events/bin'))
        plot.add_distribution(FDistribution('dphipfUA',0,3.14,10,'min#Delta#phi(jet,E_{T}^{miss})','Events'))
  
    #recoil.calc_chi2 = True
    plot.add_distribution(recoil)

    #global variable
    #plot.add_distribution(FDistribution('nJet',-0.5,9.5,9,'N_{jet}','Events'))
    #plot.add_distribution(FDistribution('npv',0,45,45,'N_{PV}','Events'))
    #plot.add_distribution(FDistribution('dphipfmet',0,3.14,20,'min#Delta#phi(jet,E_{T}^{miss})','Events'))

    #global hadronic
    #plot.add_distribution(FDistribution('jetPt[0]',0,1000,20,'Leading Jet p_{T} [GeV]','Events/50 GeV'))
    #plot.add_distribution(FDistribution('jetPt[1]',0,1000,20,'Subleading Jet p_{T} [GeV]','Events/50 GeV'))
    #plot.add_distribution(FDistribution('jet1Eta',-2.5,2.5,20,'Leading Jet #eta','Events/bin'))
    #plot.add_distribution(FDistribution('jet2Eta',-2.5,2.5,20,'Sub-Leading Jet #eta','Events/bin'))
    #plot.add_distribution(FDistribution('jetCSV[0]',0,1,20,'jet 1 CSV','Events'))
    #plot.add_distribution(FDistribution('jetCSV[1]',0,1,20,'jet 2 CSV','Events'))
    plot.add_distribution(FDistribution('fj1DoubleCSV',0,1,20,'fatjet 1 DoubleCSV','Events'))
    #plot.add_distribution(FDistribution('isojetNBtags',-0.5,9.5,9,'N_{isoBtagjet}','Events'))
    #global lepton
    #plot.add_distribution(FDistribution('nTightLep',-0.5,4.5,5,'Number of tight lepton','Events/bin'))
    #plot.add_distribution(FDistribution('nLooseLep',-0.5,4.5,5,'Number of loose lepton','Events/bin'))

    #plot.add_distribution(FDistribution('nLooseElectron',-0.5,4.5,5,'Number of loose Electron','Events/bin'))
    #plot.add_distribution(FDistribution('nLooseMuon',-0.5,4.5,5,'Number of loose Muon','Events/bin'))

    #fatjet
    plot.add_distribution(FDistribution('fj1MSD',0,600,20,'fatjet m_{SD} [GeV]','Events'))

    #fjmass=VDistribution("fj1MSD",fatjetBins,"fatjet m_{SD} [GeV]","Events")
    #plot.add_distribution(fjmass)

    ####
    plot.add_distribution(FDistribution('fj1Pt',200,700,15,'fatjet p_{T} [GeV]','Events/25 GeV'))
    plot.add_distribution(FDistribution('fj1Eta',-2.5,2.5,10,'fatjet #eta [GeV]','Events'))
    #plot.add_distribution(FDistribution('top_ecf_bdt',-1,1,20,'Top BDT','Events'))
    #plot.add_distribution(FDistribution('fj1MaxCSV',0,1,20,'fatjet max CSV','Events'))
    #plot.add_distribution(FDistribution('fj1Tau32',0,1,20,'fatjet #tau_{32}','Events'))
    #plot.add_distribution(FDistribution('fj1Tau32SD',0,1,20,'fatjet #tau_{32}^{SD}','Events'))
    #Cutflow
    plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))
    system('mkdir -p %s/%s/' %(args.outdir,region))
    plot.draw_all(args.outdir+'/'+region+'/')

def fromLimit(region):
    region=region
    print 'Plotting from fitting ntuple: ',region, ' region'
    znunu = Process('Z(#nu#nu)+jets',root.kZjets,'Zvv_'+region,root.kCyan-9)
    zjets         = Process('Z+jets',root.kZjets,'Zll_'+region,root.kCyan-9)
    wjets         = Process('W+jets',root.kWjets,'Wlv_'+region,root.kGreen-10)
    diboson       = Process('Diboson',root.kDiboson,'Diboson_'+region,root.kYellow-9)
    ttbar         = Process('t#bar{t}',root.kTTbar,'ttbar_'+region,root.kOrange-4)
    #ttg           = Process('t#bar{t}#gamma',root.kTTbar,'ttbarg_'+args.region)
    singletop     = Process('Single t',root.kST,'ST_'+region,root.kRed-9)
    #singletopg    = Process('t#gamma',root.kST,args.region)
    qcd           = Process('QCD',root.kQCD,'QCD_'+region,root.kMagenta-10)
    #gjets         = Process('#gamma+jets',root.kGjets,args.region)
    data          = Process("Data",root.kData,'Data_'+region)
    #signal        = Process('m_{V}=1.75 TeV, m_{#chi}=1 GeV',root.kSignal,args.region)
    processes = [qcd,diboson,singletop,wjets,ttbar,zjets]

    
    if 'qcd' in region:
        processes = [diboson,singletop,wjets,ttbar,zjets,qcd]
    if 'zee' or 'zmm' in region:
        processes = [qcd,diboson,singletop,ttbar,wjets,zjets]
    if 'wen' or 'wmn' in region:
        processes = [qcd,diboson,singletop,zjets,ttbar,wjets]
    ### ASSIGN FILES TO PROCESSES ###
    if 'signal' in region or 'qcd' in region:
        processes = [qcd,zjets,singletop,ttbar,diboson,wjets,znunu]

    for p in processes:
        p.add_file(f(region))
        p.additional_weight='weight'
        plot.add_process(p)

    #Inclusion of data
    data.add_file(f(region))
    data.additional_weight='weight'
    processes.append(data)
    plot.add_process(data)

    recoilBins = [250,280,310,350,400,450,600,1000]
    nRecoilBins = len(recoilBins)-1
    plot.add_distribution(FDistribution('fjpt',200,700,15,'fatjet p_{T} [GeV]','Events/25 GeV'))
    plot.add_distribution(FDistribution('fjmass',0,600,20,'fatjet m_{SD} [GeV]','Events'))
    plot.add_distribution(FDistribution('doubleb',0,1,20,'fatjet 1 DoubleCSV','Events'))
    plot.add_distribution(FDistribution('n2',0,0.5,10,'n2','Events'))
    plot.add_distribution(FDistribution('n2ddt56',0,0.5,10,'n2ddt56','Events'))
    plot.add_distribution(FDistribution('n2ddt53',0,0.5,10,'n2ddt53','Events'))
    recoil=VDistribution("met",recoilBins,"PF MET [GeV]","Events/GeV")
    plot.add_distribution(recoil)

    system('mkdir -p %s/%s/' %(args.outdir,region))
    plot.draw_all(args.outdir+'/'+region+'/')

if not args.fromlimit:
    normalPlotting(region)
else:
    fromLimit(region)
