#!/usr/bin/env python

# imports and load libraries
from array import array
from glob import glob
from re import sub
from sys import argv,exit
import os
import sys
import subprocess
from os import environ,system,path,remove
from argparse import ArgumentParser

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load

parser = ArgumentParser(description='Parsing dataset folder')
parser.add_argument('--folder',type=str,default=None)
args = parser.parse_args()

inbase="/store/user/paus/pandaf/009/"
inpath = 'xrdfs root://xrootd.cmsaf.mit.edu ls ' + inbase
files = subprocess.check_output(inpath,shell=True)
file = files.split()

def runAll():
    for f in file:
        split = f.split('/')
        direc =  split[-1]
        print direc
        path = "/uscms_data/d3/lpcmetx/catalog/80x-v1/"
        #path = "/uscms_data/d3/lpcmetx/test-shoh/"
        os.makedirs(path + direc)

        subfiles=subprocess.check_output("xrdfs root://xrootd.cmsaf.mit.edu ls "+f,shell=True)
        subfile = subfiles.split()
        with open(path + direc +"/RawFiles.00","w") as f:
            for s in subfile:
                rootfile = root.TFile.Open("root://xrootd.cmsaf.mit.edu/" + s)
                tree = rootfile.Get("events")
                nevents = tree.GetEntriesFast()
                rootfile.Close()
                print "root://xrootd.cmsaf.mit.edu/" + s + " " + str(nevents) + " " + str(nevents) + " 1 1 1 1"
                f.write("root://xrootd.cmsaf.mit.edu/" + s + " " + str(nevents) + " " + str(nevents) + " 1 1 1 1" + '\n')
            f.close()

def runOneFolder(folderName):
    f=inbase+'/'+folderName
    if f in file:
        split = f.split('/')
        direc =  split[-1]
        print direc
        #path = "/uscms_data/d3/lpcmetx/test-shoh/"
        path = "/uscms_data/d3/lpcmetx/catalog/80x-v1/"
        os.makedirs(path + direc)
        subfiles=subprocess.check_output("xrdfs root://xrootd.cmsaf.mit.edu ls "+f,shell=True)
        subfile = subfiles.split()
        with open(path + direc +"/RawFiles.00","w") as f:
            for s in subfile:
                rootfile = root.TFile.Open("root://xrootd.cmsaf.mit.edu/" + s)
                tree = rootfile.Get("events")
                nevents = tree.GetEntriesFast()
                rootfile.Close()
                print "root://xrootd.cmsaf.mit.edu/" + s + " " + str(nevents) + " " + str(nevents) + " 1 1 1 1"
                f.write("root://xrootd.cmsaf.mit.edu/" + s + " " + str(nevents) + " " + str(nevents) + " 1 1 1 1" + '\n')
            f.close()


if not args.folder:
    print 'folder not specify, running over %s' %inbase
    runAll()
else:
    print 'folder is specified: %s' %args.folder
    runOneFolder(args.folder)
