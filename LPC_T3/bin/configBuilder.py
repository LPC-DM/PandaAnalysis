#!/usr/bin/env python

import sys
import argparse
import subprocess
from re import sub
from os import getenv
from PandaCore.Tools.Misc import PInfo
from PandaCore.Tools.job_management import DataSample,convert_catalog

workdir = getenv('SUBMIT_WORKDIR')
parser = argparse.ArgumentParser(description='convert configuration')
parser.add_argument('--infile',type=str,default=None)
parser.add_argument('--outfile',type=str,default=None)
parser.add_argument('--nfiles',type=int,default=None)
args = parser.parse_args()

fin = open(args.infile)
samples = convert_catalog(list(fin),as_dict=True)

fout = open(args.outfile,'w')
keys = sorted(samples)
#print keys
counter=0
for k in keys:
	sample = samples[k]
	configs = sample.get_config(args.nfiles,suffix='_%i')
	if 'MET' in k or 'Znunu' in k or 'ZJets' in k or 'WJets' in k: configs = sample.get_config(args.nfiles*2,suffix='_%i')
	if 'SingleElectron' in k: configs = sample.get_config(args.nfiles*3,suffix='_%i')
	for c in configs:
		fout.write(c%(counter,counter))
		counter += 1

PInfo('configBuilder.py','Submission will have %i jobs'%(counter))

fout.close()
