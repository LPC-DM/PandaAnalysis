#!/usr/bin/env python

# imports and load libraries
from array import array
from glob import glob
from re import sub
from sys import argv,exit
import sys
from os import environ,system,path,remove
from argparse import ArgumentParser

sname = argv[0]
parser = ArgumentParser()
parser.add_argument('--silent', action='store_true')
parser.add_argument('--cfg', type=str, default='common')
parser.add_argument('--skip_missing', action='store_true')
parser.add_argument('arguments', type=str, nargs='+')
args = parser.parse_args()
arguments = args.arguments
VERBOSE = not args.silent
skip_missing = args.skip_missing
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load

if args.cfg == 'leptonic':
    from PandaCore.Tools.process_leptonic import *
    xsecscale = 1000
else:
    from PandaCore.Tools.process import *
    xsecscale = 1

sys.path.append(environ['CMSSW_BASE'] + '/src/PandaAnalysis/LPC_T3/merging/configs/')
cfg = __import__(args.cfg)

Load('Normalizer')

# global variables
pds = {}
for k,v in processes.iteritems():
    if v[1]=='MC':
        pds[v[0]] = (k,v[2])  
    else:
        pds[v[0]] = (k,-1)

user = environ['USER']
system('mkdir -p /uscmst1b_scratch/lpc1/3DayLifetime/%s/split'%user) # tmp dir
system('mkdir -p /uscmst1b_scratch/lpc1/3DayLifetime/%s/merged'%user) # tmp dir

#inbase = environ['SUBMIT_OUTDIR']
inbase = "/eos/uscms/store/group/lpcmetx/panda/80X-v1dot1/lepmonotop/batch/"
#outbase = environ['PANDA_FLATDIR']
#outbase = "/uscms_data/d3/lpcmetx/panda/80X-v1dot1/lepmonotop/flat/"
outbase = "/uscms/home/dsilveri/nobackup/PandaAnalysis_2018/CMSSW_8_0_29/src/PandaAnalysis/LPC_T3/merging/"

hadd_cmd = 'hadd -k -ff -n 100 -f '

if VERBOSE:
    suffix = ''
else:
    suffix = ' > /dev/null '

# helper functions
def hadd(inpath,outpath):
    if type(inpath)==type('str'):
        infiles = glob(inpath)
        #if len(infiles) > 1: # if 1 file, use mv
        PInfo(sname,'hadding %s into %s'%(inpath,outpath))
        cmd = '%s %s %s %s'%(hadd_cmd, outpath,inpath,suffix)
        system(cmd)
        return True
    else:
        infiles = inpath
    if len(infiles)==0:
        PWarning(sname,'nothing hadded into '+outpath)
        return False
    elif len(infiles)==1:
        PInfo(sname,'moving %s to %s'%(inpath[0],outpath))
        cmd = 'mv -v %s %s'%(infiles[0],outpath)
    else:
        cmd = '%s %s '%(hadd_cmd, outpath)
        for f in infiles:
            if path.isfile(f):
                cmd += '%s '%f
        PInfo(sname,'hadding into %s'%(outpath))
    if VERBOSE: PInfo(sname,cmd)
    system(cmd+suffix)
    return True

def normalizeFast(fpath,opt):
    xsec=-1
    if type(opt)==type(1.) or type(opt)==type(1):
        xsec = opt
    else:
        try:
            xsec = processes[proc][2]
        except KeyError:
            for k,v in processes.iteritems():
                if proc in k:
                    xsec = v[2]
    if xsec<0:
        PError(sname,'could not find xsec, skipping %s!'%opt)
        return
    xsec *= xsecscale
    PInfo(sname,'normalizing %s (%s) ...'%(fpath,opt))
    n = root.Normalizer();
    n.NormalizeTree(fpath,xsec)

def merge(shortnames,mergedname):
    to_skip = []
    for shortname in shortnames:
        if 'monotop' in shortname:
            pd = shortname
            xsec = 1
        elif 'vector' in shortname:
            tmp_ = shortname
            replacements = {
                'Vector_MonoTop_Leptonic_Mphi_':'',
                '_Mchi_':'_',
                }
            for k,v in replacements.iteritems():
                tmp_ = tmp_.replace(k,v)
            m_V,m_DM = [int(x) for x in tmp_.split('_')]
            params = read_nr_model(m_V,m_DM)
            if params:
                xsec = params.sigma
            else:
                xsec = 1
        elif 'Scalar' in shortname:
            tmp_ = shortname
            replacements = {
                'Scalar_MonoTop_LO_Mphi-':'',
                '_13TeV-madgraph':'',
                '_Mchi-':'_',
                }
            for k,v in replacements.iteritems():
                tmp_ = tmp_.replace(k,v)
            m_V,m_DM = [int(x) for x in tmp_.split('_')]
            params = read_r_model(m_V,m_DM)
            if params:
                xsec = params.sigma
            else:
                xsec = 1
        elif shortname in pds:
            pd = pds[shortname][0]
            xsec = pds[shortname][1]
        else:
            for shortname_ in [shortname.split('_')[0],shortname.split('_')[-1]]:
                if shortname_ in pds:
                    pd = pds[shortname_][0]
                    xsec = pds[shortname_][1]
                    break
        inpath = '`xrdfs root://cmseos.fnal.gov ls -u ' + inbase + ' | grep \'' + shortname + '_\'`'
        success = hadd(inpath,'/uscmst1b_scratch/lpc1/3DayLifetime/%s/split/%s.root'%(user,shortname))
        if xsec>0 and success:
            normalizeFast('/uscmst1b_scratch/lpc1/3DayLifetime/%s/split/%s.root'%(user,shortname),xsec)
        if not success:
            if not skip_missing:
                PError(sname, 'Could not merge %s, exiting!'%shortname)
                exit(1)
            else:
                to_skip.append(shortname)
    to_hadd = ['/uscmst1b_scratch/lpc1/3DayLifetime/%s/split/%s.root'%(user,x) for x in shortnames if x not in to_skip]
    hadd(to_hadd, '/uscmst1b_scratch/lpc1/3DayLifetime/%s/merged/%s.root'%(user,mergedname))
    for f in to_hadd:
        system('rm -f %s'%f)


args = {}

for pd in arguments:
    if pd in cfg.d:
        args[pd] = cfg.d[pd]
    else:
        args[pd] = [pd]

for pd in args:
    merge(args[pd],pd)
    merged_file = '/uscmst1b_scratch/lpc1/3DayLifetime/%s/merged/%s.root'%(user,pd)
    hadd(merged_file ,outbase+pd+'.root') # really an mv
    system('rm -f %s'%merged_file)
    PInfo(sname,'finished with '+pd)

