#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv,path
from time import clock,time
import json

which = int(argv[1])
submit_id = int(argv[2])
sname = argv[0]
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
import PandaCore.Tools.job_config as cb
import PandaAnalysis.LPC_T3.job_utilities as utils
from PandaAnalysis.Flat.analysis import monojet

Load('PandaAnalyzer')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

def fn(input_name, isData, full_path):
    
    PInfo(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    skimmer = root.PandaAnalyzer()
    analysis = monojet(True)
    analysis.processType = utils.classify_sample(full_path, isData)
    analysis.genOnly = False
    skimmer.SetAnalysis(analysis)
    skimmer.isData = isData
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kPho)

    return utils.run_PandaAnalyzer(skimmer, isData, input_name)

if __name__ == "__main__":
    sample_list = cb.read_sample_config('local.cfg',as_dict=False)
    to_run = None 
    for s in sample_list:
        if which==s.get_id():
            to_run = s
            break
    if not to_run:
        PError(sname,'Could not find a job for PROCID=%i'%(which))
        exit(3)

    outdir = getenv('SUBMIT_OUTDIR')
    lockdir = getenv('SUBMIT_LOCKDIR')
    outfilename = to_run.name+'_%i.root'%(submit_id)
    processed = {}
    
    utils.main(to_run, processed, fn)

    utils.hadd(processed.keys())
    utils.print_time('hadd')

    ret = utils.stageout(outdir,outfilename)
    utils.cleanup('*.root')
    utils.print_time('stageout and cleanup')
    if not ret:
        utils.write_lock(lockdir,outfilename,processed)
        utils.cleanup('*.lock')
        utils.print_time('create lock')
    else:
        exit(-1*ret)

    exit(0)

