#!/usr/bin/env python

from glob import glob
from os import stat,getenv,system,path
from multiprocessing import Pool
from PandaCore.Tools.Misc import *
from re import sub, match
from sys import argv
import argparse

parser = argparse.ArgumentParser(description='make config file')
parser.add_argument('--catalog',type=str,default='/home/cmsprod/catalog/t2mit/pandaf/009')
parser.add_argument('--user_catalog', action='store_true')
parser.add_argument('--mc_catalog',type=str,default=None)
parser.add_argument('--data_catalog',type=str,default=None)
parser.add_argument('--cfg',type=str,default='common')
parser.add_argument('--outfile',type=str,default=None)
parser.add_argument('--include',nargs='+',type=str,default=None)
parser.add_argument('--exclude',nargs='+',type=str,default=None)
parser.add_argument('--smartcache',action='store_true')
parser.add_argument('--force',action='store_true')
args = parser.parse_args()

if not args.mc_catalog:
    args.mc_catalog = args.catalog
if not args.data_catalog:
    args.data_catalog = args.catalog

if args.cfg == 'leptonic':
    from PandaCore.Tools.process_leptonic import *
else:
    from PandaCore.Tools.process import *

user = getenv('USER')

class CatalogSample:
    def __init__(self,name,dtype,xsec):
        self.name = name
        self.dtype = dtype
        self.xsec = xsec
        self.files = []
    def add_file(self,f):
        self.files.append(f)
    def get_lines(self,smartcache_args=None):
        lines = []
        nickname = self.name+'_%i'
        for f in self.files:
            ds_ = f.split('/')[-2]
            f_ = f.split('/')[-1]
            book_ = '/'.join(args.catalog.split('/')[-2:])
            lines.append('{0:<25} {2:<10} {3:<15} {1}\n'.format(nickname,f,self.dtype,self.xsec)) 
            if smartcache_args is not None:
                smartcache_args.append('/cms/store/user/paus/%s/%s'%(book_,ds_))
        return lines

def smartcache(arguments):
    #arguments = ' '.join(arguments)
    for a in arguments:
        cmd = ('dynamoCache request --datasets %s'%(a))
        #print cmd
        system(cmd)

def checkDS(nickname,include,exclude):
  included=False
  if include:
    for i in include:
      if i in nickname:
        included=True
        break
  else:
    included=True
  excluded=False
  if exclude:
    for e in exclude:
      if e in nickname:
        excluded=True
        break
  else:
    excluded=False
  return (included and not(excluded))

samples = {}

could_not_find = []

def cat(catalog, condition=lambda x : True): 
    global samples, could_not_find
    for d in sorted(glob(catalog+'/*')):
        dirname = d.split('/')[-1]
        if not condition(dirname):
            continue
        shortname = dirname.split('+')[0]
        try:
            properties = processes[shortname]
            found = True
        except KeyError:
            if args.force:
              found = False
              properties = (shortname,'MC',1)
            else:
              continue
        dtype = 'MC' if 'MINIAODSIM' in dirname else 'Data'
        if not checkDS(properties[0],args.include,args.exclude):
            continue
        if not found:
          could_not_find.append(shortname)
        if properties[0] not in samples:
            samples[properties[0]] = CatalogSample(*properties)
        PInfo(argv[0], 'Selecting %s'%properties[0])
        sample = samples[properties[0]]
        for rfpath in glob(d+'/RawFiles.*'):
            rawfile = open(rfpath)
            for line in rawfile:
                sample.add_file(line.split()[0])

cat(args.mc_catalog, lambda x : bool(match('.*SIM$', x)))
cat(args.data_catalog, lambda x : bool(match('.*AOD$', x)))
if args.user_catalog:
    cat(args.catalog.replace('cmsprod',user))
    cat(args.catalog.replace('cmsprod',user).replace('t2','t3'))

if len(could_not_find)>0:
    PWarning(argv[0],"Could not properly catalog following datasets (force=%s)"%('True' if args.force else 'False'))
    for c in could_not_find:
        PWarning(argv[0],'\t'+c)

cfg_file = open(args.outfile,'w')
lines = []
smartcache_args = [] if args.smartcache else None
for k in sorted(samples):
    sample = samples[k]
    lines += sample.get_lines(smartcache_args)
for iL in xrange(len(lines)):
    cfg_file.write(lines[iL]%(iL))

cfg_file.close()
PInfo(argv[0],'Cataloged %i files for %i datasets'%(len(lines),len(samples)))
PInfo(argv[0],'Output written to '+args.outfile)

if args.smartcache:
    smartcache_datasets = list(set(smartcache_args))
    PInfo(argv[0],'Making smartcache requests for files')
    smartcache(smartcache_datasets)
