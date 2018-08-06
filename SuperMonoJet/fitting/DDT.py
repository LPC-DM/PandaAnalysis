#
import os
import ROOT
from ROOT import *
from array import array
import math
from math import *
import sys
import pdb
import os
import argparse

parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--output',metavar='output',type=str,default='')
parser.add_argument('--input',metavar='input',type=str,default='/data/t3home000/mcremone/lpc/jorgem/skim/v_8026_0_4/monohiggs_boosted/')
parser.add_argument('--process',metavar='process',type=str,default='')
basedir = parser.parse_args().input
foredir = parser.parse_args().output
process = parser.parse_args().process


def ComputeDDT(name, point, nPtBins, nRhoBins, H):
	DDT = TH2F(name, "", nRhoBins, -6, -1.5, nPtBins, 200, 1600)
	DDT.SetStats(0)
	nXb = H.GetXaxis().GetNbins()
	nYb = H.GetYaxis().GetNbins()
	for x in range(nXb):
		for y in range(nYb):
			proj = H.ProjectionZ("H3"+str(x)+str(y),x+1,x+1,y+1,y+1)
			print str(x+1) + "," + str(y+1) + ":    "+ str(proj.Integral())
			p = array('d', [point])
			q = array('d', [0.0]*len(p))
			proj.GetQuantiles(len(p), q, p)
			DDT.SetBinContent( x+1, y+1, q[0] );
	return DDT

def DisplayDDT(DDT, Title, SaveName):
	C = TCanvas("TempCanvas", "Title", 575, 500)
	C.cd()
	DDT.SetStats(0)
	DDT.GetXaxis().SetTitle("jet #rho")
	DDT.GetYaxis().SetTitle("jet p_{T}")
	DDT.Draw("COLZ")
	C.Print("MAP_"+SaveName+".gif")
	C.SaveAs("DDT_"+SaveName+".pdf")

def setlist(path,process):
	f = []
	if os.path.isfile(path):
		if not "root" in path:
			print "ERROR PLEASE INSERT A THE PATH OF A FOLDER WITH ROOT FILES OR THE PATH TO A ROOT FILE"
			sys.exit("Error message")
		if "/" in path:	
			f_ = path.split("/")[-1]
			if "." in f_:
				f_ = f_.split(".")[0]
		f.append((f_,path))	
	if os.path.isdir(path):
		#print "hello"
		lroot = os.listdir(path)
		#print lroot
		for lr in lroot:
			if 'DDT' in lr: continue
			if ".root" in lr:
				f_ = lr.split(".")[0]
				if f_==process: f.append((f_,path+'/'+lr))
	#print f
	return f
def outname(path):
	if os.path.isdir(path):
		f = path.split("/")[-1]
		return path,f
	if os.path.isfile(path):
		f = path.split("/")[0:-1]
		fs = path.split("/")[-2]
		for g in f:
			fg = g+'/'
		return fg,fs
	else: return sys.env('SKIM_MONOHIGGS_BOOSTED_FLATDIR'),"random"

path_ = basedir
#print path_
if path_ :
	Bkgs_tags = setlist(path_,process)
#print Bkgs_tags
H3={}


for bks,B in Bkgs_tags:

	H3[bks] = TH3F("H3_%s"%(bks), "H3_%s"%(bks), 9, -6, -1.5, 12, 200, 1600, 500, 0, 0.5)
if  foredir == '': 
	print 'No output: Taking default one'
	foredir = path_
	
for bks,B in Bkgs_tags:

	print 'Starting with '+bks+'-------------------------- :)'

	H3[bks].SetStats(0)
	F = TFile(B)
	#print B
	tree = "events"
	if "test" in B:
		tree = "Diboson_test"
	T = F.Get("%s"%tree)
	#print T
	n = T.GetEntries()
	print n
	for j in range(0, n): # Here is where we loop over all events.
		'''
		print  j % (1 * n/100)
                if(j % (1 * n/100) == 0):
                        sys.stdout.write("\r[" + "="*int(20*j/n) + " " + str(round(100.*j/n,0)) + "% done")
                        sys.stdout.flush()
		'''
		T.GetEntry(j)
		weight = T.normalizedWeight*T.sf_pu*T.sf_ewkV*T.sf_qcdV
		#weight = 1
		PT = T.fj1Pt
		preRHO = T.fj1MSD_corr*T.fj1MSD_corr/T.fj1Pt/T.fj1Pt
		if preRHO > 0. and T.fj1ECFN_1_2_10 != 0.:
			RHO = math.log(preRHO)
			jtN2b1sd_8 = T.fj1ECFN_2_3_10/math.pow(T.fj1ECFN_1_2_10,2.00)
			if PT > 200 and RHO < -1.5 and RHO > -6.0 and jtN2b1sd_8 > 0.:
				H3[bks].Fill(RHO, PT, jtN2b1sd_8, weight)
DDT_5by6={}
DDT_5by3={}
print "%s/DDT.root"%(foredir)
Fout = TFile("%s/DDT.root"%foredir, "recreate")
Fout.cd()
for key in H3:
	DDT_5by6[key]=ComputeDDT('DDT_5by6', 0.6, 12, 9, H3[key])
	#DisplayDDT(DDT_5by6[key], "DDT vals at 5%", "DDT_5by6")
	DDT_5by6[key].Write()
	DDT_5by3[key]=ComputeDDT('DDT_5by3', 0.6, 3, 9, H3[key])
	#DisplayDDT(DDT_5by3[key], "DDT vals at 5%", "DDT_5by3")
	DDT_5by3[key].Write()
Fout.Close()
#DDT_5by6 = ComputeDDT("DDT_5by6", 0.05, 12, 9, H3)
#DisplayDDT(DDT_5by[]6, "DDT vals at 5%", "DDT_5by6")
#DDT_5by3 = ComputeDDT("DDT_5by3", 0.05, 3, 9, H3)
