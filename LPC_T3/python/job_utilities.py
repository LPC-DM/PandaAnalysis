#!/usr/bin/env python

import json
import socket
from re import sub
from sys import exit
from random import choice
from time import clock,time,sleep
from os import system,getenv,path,environ

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load
import PandaCore.Tools.job_config as cb

_sname = 'T3.job_utilities'                                    # name of this module
_data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'  # data directory
_host = socket.gethostname()                                   # where we're running
_IS_T3 = (_host[:2] == 't3')                                   # are we on the T3?
REMOTE_READ = True                                             # should we read from hadoop or copy locally?
local_copy = bool(smart_getenv('SUBMIT_LOCALACCESS', True))    # should we always xrdcp from T2?

stageout_protocol = None                                       # what stageout should we use?
if _IS_T3:
    stageout_protocol = 'cp' 
elif system('which gfal-copy') == 0:
    stageout_protocol = 'gfal'
elif system('which lcg-cp') == 0:
    stageout_protocol = 'lcg'
else:
    try:
        ret = system('wget http://t3serv001.mit.edu/~snarayan/misc/lcg-cp.tar.gz')
        ret = max(ret, system('tar -xvf lcg-cp.tar.gz'))
        if ret:
            raise RuntimeError
        environ['PATH'] = '$PWD/lcg-cp:'+environ['PATH']
        environ['LD_LIBRARY_PATH'] = '$PWD/lcg-cp:'+environ['LD_LIBRARY_PATH']
        stageout_protocol = 'lcg'
    except Exception as e:
        PError(_sname,
               'Could not install lcg-cp in absence of other protocols!')
        raise e


# derived from t3serv006.mit.edu:/etc/bestman2/conf/bestman2.rc
_gsiftp_doors = [
        't3btch000.mit.edu',
        't3btch001.mit.edu',
        't3btch003.mit.edu',
        't3btch004.mit.edu',
        't3btch005.mit.edu',
        't3btch006.mit.edu',
        't3btch010.mit.edu',
        't3btch013.mit.edu',
        't3btch014.mit.edu',
        't3btch018.mit.edu',
        't3btch021.mit.edu',
        't3btch025.mit.edu',
        't3btch026.mit.edu',
        't3btch027.mit.edu',
        't3btch028.mit.edu',
        't3btch029.mit.edu',
        't3btch030.mit.edu',
        ]



# global to keep track of how long things take
_stopwatch = time() 
def print_time(label):
    global _stopwatch
    now_ = time()
    PDebug(_sname+'.print_time:',
           '%.1f s elapsed performing "%s"'%((now_-_stopwatch),label))
    _stopwatch = now_


# convert an input name to an output name
def input_to_output(name):
    if 'input' in name:
        return name.replace('input','output')
    else:
        return 'output_' + name.split('/')[-1]


# find data and bring it into the job somehow 
#  - if local_copy and it exists locally, then:
#    - if REMOTE_READ, read from hadoop
#    - else copy it locally
#  - else xrdcp it locally
def copy_local(long_name):
    full_path = long_name
    PInfo(_sname,full_path)

    panda_id = long_name.split('/')[-1].split('_')[-1].replace('.root','')
    input_name = 'input_%s.root'%panda_id
    copied = False

    # if the file is cached locally, why not use it?
    if 'scratch' in full_path:
        local_path = full_path.replace('root://t3serv006.mit.edu/','/mnt/hadoop')
    else:
        local_path = full_path.replace('root://xrootd.cmsaf.mit.edu/','/mnt/hadoop/cms')
    PInfo(_sname+'.copy_local','Local access is configured to be %s'%('on' if local_copy else 'off'))
    if local_copy and path.isfile(local_path): 
        # apparently SmartCached files can be corrupted...
        ftest = root.TFile(local_path)
        if ftest and not(ftest.IsZombie()):
            PInfo(_sname+'.copy_local','Opting to read locally')
            if REMOTE_READ:
                return local_path
            else:
                cmd = 'cp %s %s'%(local_path, input_name)
                PInfo(_sname+'.copy_local',cmd)
                system(cmd)
                copied = True

    if not copied:
        cmd = "xrdcp %s %s"%(full_path,input_name)
        PInfo(_sname+'.copy_local',cmd)
        system(cmd)
        copied = True
            
    if path.isfile(input_name):
        PInfo(_sname+'.copy_local','Successfully copied to %s'%(input_name))
        return input_name
    else:
        PError(_sname+'.copy_local','Failed to copy %s'%input_name)
        return None


# wrapper around rm -f. be careful!
def cleanup(fname):
    ret = system('rm -f %s'%(fname))
    if ret:
        PError(_sname+'.cleanup','Removal of %s exited with code %i'%(fname,ret))
    else:
        PInfo(_sname+'.cleanup','Removed '+fname)
    return ret


# wrapper around hadd
def hadd(good_inputs, output='output.root'):
    good_outputs = ' '.join([input_to_output(x) for x in good_inputs])
    cmd = 'hadd -f ' + output + ' ' + good_outputs
    print cmd
    ret = system(cmd)    
    if not ret:
        PInfo(_sname+'.hadd','Merging exited with code %i'%ret)
    else:
        PError(_sname+'.hadd','Merging exited with code %i'%ret)


# remove any irrelevant branches from the final tree.
# this MUST be the last step before stageout or you 
# run the risk of breaking something
def drop_branches(to_drop=None, to_keep=None):
    if not to_drop and not to_keep:
        return 0

    if to_drop and to_keep:
        PError(_sname+'.drop_branches','Can only provide to_drop OR to_keep')
        return 0

    f = root.TFile('output.root','UPDATE')
    t = f.FindObjectAny('events')
    n_entries = t.GetEntriesFast() # to check the file wasn't corrupted
    if to_drop:
        if type(to_drop)==str:
            t.SetBranchStatus(to_drop,False)
        else:
            for b in to_drop:
                t.SetBranchStatus(b,False)
    elif to_keep:
        t.SetBranchStatus('*',False)
        if type(to_keep)==str:
            t.SetBranchStatus(to_keep,True)
        else:
            for b in to_keep:
                t.SetBranchStatus(b,True)
    t_clone = t.CloneTree()
    f.WriteTObject(t_clone,'events','overwrite')
    f.Close()

    # check that the write went okay
    f = root.TFile('output.root')
    if f.IsZombie():
        PError(_sname+'.drop_branches','Corrupted file trying to drop '+str(to_drop))
        return 1 
    t_clone = f.FindObjectAny('events')
    if (n_entries==t_clone.GetEntriesFast()):
        return 0
    else:
        PError(_sname+'.drop_branches','Corrupted tree trying to drop '+str(to_drop))
        return 2



# stageout a file (e.g. output or lock)
#  - if _IS_T3, execute a simple cp
#  - else, use lcg-cp
# then, check if the file exists:
#  - if _IS_T3, use os.path.isfile
#  - else, use lcg-ls
def stageout(outdir,outfilename,infilename='output.root',n_attempts=10,ls=None):
    gsiftp_doors = _gsiftp_doors[:]
    if ls is None:
        # if it's not on hadoop, copy the file to test it exists, to force nfs refresh
        ls = not('hadoop' in outdir) 
    ls = False # override, I don't trust network-mounted filesystems anymore 
    if stageout_protocol is None:
        PError(_sname+'.stageout',
               'Stageout protocol has not been satisfactorily determined! Cannot proceed.')
        return -2
    timeout = 300
    ret = -1
    for i_attempt in xrange(n_attempts):
        door = choice(gsiftp_doors); gsiftp_doors.remove(door)
        failed = False
        if stageout_protocol == 'cp':
            cpargs =     ' '.join(['cp',
                                   '-v', 
                                   '$PWD/%s'%infilename,
                                   '%s/%s'%(outdir,outfilename)])
            if ls:
                lsargs = ' '.join(['ls',
                                   '%s/%s'%(outdir,outfilename)])
            else:
                lsargs = ' '.join(['cp',
                                   '-v',
                                   '%s/%s'%(outdir,outfilename),
                                   '$PWD/testfile'])
        elif stageout_protocol == 'gfal':
            cpargs =     ' '.join(['gfal-copy',
                                   '-f', 
                                   '--transfer-timeout %i'%timeout,
                                   '$PWD/%s'%infilename,
                                   'gsiftp://%s:2811//%s/%s'%(door,outdir,outfilename)])
            if ls:
                lsargs = ' '.join(['gfal-ls',
                                   'gsiftp://%s:2811//%s/%s'%(door,outdir,outfilename)])
            else:
                lsargs = ' '.join(['gfal-copy',
                                   '-f', 
                                   '--transfer-timeout %i'%timeout,
                                   'gsiftp://%s:2811//%s/%s'%(door,outdir,outfilename),
                                   '$PWD/testfile'])
        elif stageout_protocol == 'lcg':
            cpargs =     ' '.join(['lcg-cp',
                                   '-v -D srmv2 -b', 
                                   'file://$PWD/%s'%infilename,
                                   'gsiftp://%s:2811//%s/%s'%(door,outdir,outfilename)])
            if ls:       
                lsargs = ' '.join(['lcg-ls',
                                   '-v -D srmv2 -b', 
                                   'gsiftp://%s:2811//%s/%s'%(door,outdir,outfilename)])
            else:
                lsargs = ' '.join(['lcg-cp',
                                   '-v -D srmv2 -b', 
                                   'gsiftp://%s:2811//%s/%s'%(door,outdir,outfilename),
                                   'file://$PWD/testfile'])
        PInfo(_sname+'.stageout',cpargs)
        ret = system(cpargs)
        if not ret:
            PInfo(_sname+'.stageout','Move exited with code %i'%ret)
            sleep(10) # give the filesystem a chance to respond
        else:
            PError(_sname+'.stageout','Move exited with code %i'%ret)
            failed = True
        if not failed:
            PInfo(_sname+'.stageout',lsargs)
            ret = system(lsargs)
            if ret:
                PError(_sname+'.stageout','Output file is missing!')
                failed = True
        if not failed:
            PInfo(_sname+'.stageout', 'Copy succeeded after %i attempts'%(i_attempt+1))
            return ret
        else:
            timeout = int(timeout * 1.5)
        system('rm -f testfile')
    PError(_sname+'.stagoeut', 'Copy failed after %i attempts'%(n_attempts))
    return ret


# write a lock file, based on what succeeded,
# and then stage it out to a lock directory
def write_lock(outdir,outfilename,processed):
    outfilename = outfilename.replace('.root','.lock')
    flock = open(outfilename,'w')
    for k,v in processed.iteritems():
        flock.write(v+'\n')
    flock.close()
    stageout(outdir,outfilename,outfilename,ls=False)


# make a record in the primary output of what
# inputs went into it
def record_inputs(outfilename,processed):
    fout = root.TFile.Open(outfilename,'UPDATE')
    names = root.TNamed('record',
                        ','.join(processed.values()))
    fout.WriteTObject(names)
    fout.Close()


# classify a sample based on its name
def classify_sample(full_path, isData):
    if not isData:
        if any([x in full_path for x in ['Vector_','Scalar_']]):
            return root.kSignal
        elif any([x in full_path for x in ['ST_','ZprimeToTT']]):
            return root.kTop
        elif 'EWKZ2Jets' in full_path:
            return root.kZEWK
        elif 'EWKW' in full_path:
            return root.kWEWK
        elif 'ZJets' in full_path or 'DY' in full_path:
            return root.kZ
        elif 'WJets' in full_path:
            return root.kW
        elif 'GJets' in full_path:
            return root.kA
        elif 'TTJets' in full_path or 'TT_' in full_path:
            return root.kTT
        elif 'HTo' in full_path:
            return root.kH
    return root.kNoProcess


# read a CERT json and add it to the skimmer
def add_json(skimmer, json_path):
    with open(json_path) as jsonFile:
        payload = json.load(jsonFile)
        for run_str,lumis in payload.iteritems():
            run = int(run_str)
            for l in lumis:
                skimmer.AddGoodLumiRange(run,l[0],l[1])


# some common stuff that doesn't need to be configured
def run_PandaAnalyzer(skimmer, isData, input_name):
    # read the inputs
    try:
        fin = root.TFile.Open(input_name)
        tree = fin.FindObjectAny("events")
        weight_table = fin.FindObjectAny('weights')
        hweights = fin.FindObjectAny("hSumW")
    except:
        PError(_sname+'.run_PandaAnalyzer','Could not read %s'%input_name)
        return False # file open error => xrootd?
    if not tree:
        PError(_sname+'.run_PandaAnalyzer','Could not recover tree in %s'%input_name)
        return False
    if not hweights:
        PError(_sname+'.run_PandaAnalyzer','Could not recover hweights in %s'%input_name)
        return False
    if not weight_table:
        weight_table = None

    output_name = input_to_output(input_name)
    skimmer.SetDataDir(_data_dir)
    if isData:
        add_json(skimmer, _data_dir+'/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt')

    rinit = skimmer.Init(tree,hweights,weight_table)
    if rinit:
        PError(_sname+'.run_PandaAnalyzer','Failed to initialize %s!'%(input_name))
        return False 
    skimmer.SetOutputFile(output_name)

    # run and save output
    skimmer.Run()
    skimmer.Terminate()

    ret = path.isfile(output_name)
    if ret:
        PInfo(_sname+'.run_PandaAnalyzer','Successfully created %s'%(output_name))
        return output_name 
    else:
        PError(_sname+'.run_PandaAnalyzer','Failed in creating %s!'%(output_name))
        return False


# main function to run a skimmer, customizable info 
# can be put in fn
def main(to_run, processed, fn):
    print_time('loading')
    for f in to_run.files:
        input_name = copy_local(f)
        print_time('copy %s'%input_name)
        if input_name:
            success = fn(input_name,(to_run.dtype!='MC'),f)
            print_time('analyze %s'%input_name)
            if success:
                processed[input_name] = f
            if input_name[:5] == 'input': # if this is a local copy
                cleanup(input_name)
            print_time('remove %s'%input_name)
    
    if len(processed)==0:
        PWarning(_sname+'.main', 'No successful outputs!')
        exit(1)


