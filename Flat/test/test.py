#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv
import json

debug_level = 0
torun = argv[1]
#output = 'testskim_resolved.root'
output = 'testskim_boosted.root'
if len(argv)>2:
    debug_level = int(argv[2])
    if len(argv)>3:
        output = argv[3]

argv = []

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
from PandaAnalysis.Flat.analysis import *

Load('PandaAnalyzer')

skimmer = root.PandaAnalyzer(debug_level)

analysis = resolved(True)
#analysis = resolved(True)
#analysis = monojet(True)
#analysis = monojet(True)
skimmer.SetAnalysis(analysis)
skimmer.firstEvent=0
skimmer.lastEvent=1000
skimmer.isData=False
#skimmer.SetPreselectionBit(root.PandaAnalyzer.kBoosted)
#skimmer.SetPreselectionBit(root.PandaAnalyzer.kResolved)
#skimmer.SetPreselectionBit(root.PandaAnalyzer.kMonojet)
skimmer.SetPreselectionBit(root.PandaAnalyzer.kMET)
if skimmer.isData:
    with open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt') as jsonFile:
        payload = json.load(jsonFile)
        for run,lumis in payload.iteritems():
            for l in lumis:
                skimmer.AddGoodLumiRange(int(run),l[0],l[1])
fin = root.TFile.Open(torun)
tree = fin.FindObjectAny("events")
hweights = fin.FindObjectAny("hSumW")
weights = fin.FindObjectAny('weights')
if not weights:
    weights = None
skimmer.SetDataDir(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/')
skimmer.Init(tree,hweights,weights)
skimmer.SetOutputFile(output)
skimmer.Run()
skimmer.Terminate()
