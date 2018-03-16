#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv,path
from time import clock,time
import json
from glob import glob

which = int(argv[1])
submit_id = int(argv[2])
sname = argv[0]
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
import PandaCore.Tools.job_config as cb
import PandaAnalysis.LPC_T3.job_utilities as utils
import PandaAnalysis.Deep.job_deepgen_utilities as deep_utils
from PandaAnalysis.Flat.analysis import deepgen

deep_utils.STORE = True
deep_utils.SAVE = True
deep_utils.INFER = False
deep_utils.NORM = False # temporary, need to recalculate normalizations

Load('PandaAnalyzer')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

def fn(input_name, isData, full_path):
    
    PInfo(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    skimmer = root.PandaAnalyzer()

    processType = utils.classify_sample(full_path, isData)
    if processType == root.kSignal:
        processType = root.kTop
    analysis = deepgen() 
    analysis.processType=processType 
#    analysis.deepAntiKtSort = True
    analysis.dump()
    skimmer.SetAnalysis(analysis)
    skimmer.isData=isData
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kGenFatJet)

    outpath = utils.run_PandaAnalyzer(skimmer, isData, input_name)
    if not outpath:
        return False 
    deep_utils.run_model(outpath.replace('.root','_gen_%i.root'), outpath)
    return True


if __name__ == "__main__":
    sample_list = cb.read_sample_config('local.cfg',as_dict=False)
    to_run = None #sample_list[which]
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
    if deep_utils.STORE and False:
        utils.stageout(outdir,outfilename.replace('.root','_arrays.root'),'arrays.root')
    utils.cleanup('*.root')
    if deep_utils.SAVE:
        if not ret:
            data = {}
            for f in glob('*npz'):
                f_data = deep_utils.np.load(f)
                for k,v in f_data.iteritems():
                    if k not in data:
                        data[k] = []
                    if v.shape[0] > 0:
                        data[k].append(v)
            if len(data['pt']) > 0:
                merged_data = {k : deep_utils.np.concatenate(v) for k,v in data.iteritems() if (k != 'singleton_branches')}
                merged_data['singleton_branches'] = data['singleton_branches'][0] 
                deep_utils.np.savez('merged_arrays.npz', **merged_data)
                utils.print_time('merging npz')
                ret = max(ret, utils.stageout(outdir, outfilename.replace('.root', '.npz'), 'merged_arrays.npz'))
        utils.cleanup('*.npz')
    utils.print_time('stageout and cleanup')
    if not ret:
        utils.write_lock(lockdir,outfilename,processed)
        utils.cleanup('*.lock')
        utils.print_time('create lock')
    else:
        exit(-1*ret)

    exit(0)
