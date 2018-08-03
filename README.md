# Panda Analysis

Throughout this readme, I will refer to several user-defined environment variables. 
These are typically defined in `LPC_T3/lpc_setup.sh`, but the user can define things elswhere.

## Installation

```bash
cmsrel CMSSW_8_0_29
cd CMSSW_8_0_29/src
cmsenv
git clone https://github.com/LPC-DM/PandaTree
git clone https://github.com/LPC-DM/PandaCore
git clone https://github.com/LPC-DM/PandaAnalysis
PandaCore/bin/genDict -f -j8                                          # I typically add PandaCore/bin to my $PATH
scram b -j8
```

Explanation: `PandaTree` is the data format, `PandaCore` is some core utilities for dealing with ROOT and python, and `PandaAnalysis` implements the analyses.
Most of your interaction should only be with `PandaAnalysis`, unless you find a bug, which you will.

numpy_root installation: (you may need this dependency for makefittingforest.py)
```bash
cd ~
```
install from the link: http://scikit-hep.org/root_numpy/install.html


## Producing a flat tree for analysis

There is a core analysis tool (`PandaAnalyzer`) that outputs a data format (`GeneralTree`).
`PandaAnalyzer` is configured using an `Analysis` class, which can turn on and off analysis-specific calculations and corresponding branches in the output tree.
The functions are implemented in `PandaAnalysis/Flat/src/Modules*cc`.

For what follows, I'll assume you're inside `$CMSSW_BASE/src/PandaAnalysis`.

To see what is produced in the output tree, you can look in `Flat/config/GeneralTree.cfg`. 
If you want to add a variable, just put it in the config and then run:
```bash
./config/generateTreeClass.py --config config/GeneralTree.cfg
```
Any custom code in the class definition will be preserved, and the new variables will be added.

### Defining your analysis

For what follows, I assume you're in `$CMSSW_BASE/src/PandaAnalysis/LPC_T3`.

First, open up `lpc_setup.sh`. 
This defines all the environment variables we'll need.
Make sure anything that is a path (`PANDA_FLATDIR`,`SUBMIT_LOGDIR`,etc) is writable by you.
The `SUBMIT*` environment variables are used for running the jobs, and the `PANDA*` variables are for running things on the outputs of the jobs.
`SUBMIT_TMPL` points to a script inside `inputs/`, that will define your analysis.
For a straightforward example, do `export SUBMIT_TMPL=skim_monojet_tmpl.py`.

Now, let's open up `inputs/skim_monojet_tmpl.py` as a concrete example and look at it.
The main function you have to worry about is `fn`.
The rest is all fluff.
For convenience, here's the (annotated) content of `fn`:
```python
def fn(input_name, isData, full_path):

    PInfo(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    skimmer = root.PandaAnalyzer()
    analysis = monojet(True)                                          # this is imported from PandaAnalysis.Flat.Analysis, where the defaults are set
    analysis.processType = utils.classify_sample(full_path, isData)  # set the type of the process
    skimmer.SetAnalysis(analysis)
    skimmer.isData=isData
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kRecoil)             # set the preselection
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kPassTrig)         # only save data events that trip a trigger

    return utils.run_PandaAnalyzer(skimmer, isData, input_name)      # run the analysis 
```

### Testing your analyzer

Inside `Flat/test`, there is a testing script that runs as:
```bash
./test.py /path/to/input/panda.root [DEBUG_LEVEL]
```
You have to open up the script and modify the number of events you want to run, the type of file it is (data, W+jets MC, etc), and what flags are on.

## Running on the grid

### Cataloging inputs

```bash
./catalogT2Prod.py --outfile /path/to/config.cfg [ --include datasets to include ] [ --exclude datasets to skip ] [ --force ] [--smartcache]
```

I recommend you put the config in a web-accessible place for use in later steps. For example:
```bash
./catalogT2Prod.py --force --outfile /publicweb/<initial-username>/<username>/$(date +%Y%m%d).cfg --include TT --exclude TTbarDM --smartcache
```
(To create a webpage at Fermilab server, you can follow this link -
https://fermi.service-now.com/kb_view_customer.do?sysparm_article=KB0011889)

The above command will do the following things:

- It will only check datasets that contain `TT` and do not contain `TTbarDM` in the dataset's nickname

- If a dataset is not found in `PandaCore.Tools.process`, it will guess the nickname of the dataset, give it a xsec of 1, and write it to the catalog (`--force`)

- If the file does not exist locally on the T3, a smartcache request will be made

- The output will be a timestamped web-facing config file


### Building the work environment

First, make sure the following environment variables are defined (some examples shown):
```bash
export PANDA_CFG="http://t3serv001.mit.edu/~mcremone/eoscatalog/test2_009.cfg"  # location of config file from previous section
export SUBMIT_TMPL="skim_monojet_tmpl.py"  # name of template script in LPC_T3/inputs
export SUBMIT_NAME="v_80X_v1_May15"  # name for this job
export SUBMIT_WORKDIR="${scratch_area}/${USER}/condor/"${SUBMIT_NAME}"/work/"  # staging area for submission
export SUBMIT_LOGDIR="${scratch_area}/${USER}/condor/"${SUBMIT_NAME}"/logs/"  # log directory
export SUBMIT_LOGDIR="${scratch_area}/${USER}/condor/"${SUBMIT_NAME}"/locks/"  # lock directory
export SUBMIT_OUTDIR="/store/user/${USER}/panda/"${SUBMIT_NAME}"/batch/"  # location of unmerged files
export PANDA_FLATDIR="${scratch_area}/${USER}/panda/"${SUBMIT_NAME}"/flat/"   # merged output
eosmkdir -p $SUBMIT_OUTDIR

export SUBMIT_CONFIG=T2  # allow running on T3 or T2. if $SUBMIT_CONFIG==T3, then only run on T3
```

`LPC_T3/inputs/$SUBMIT_TMPL` should be the skimming configuration you wish to run your files through. 

### Configure condor submission
```bash
cd ~
wget http://shoh.web.cern.ch/shoh/public/Panda/condor-8.6.3-x86_64_RedHat6-stripped.tar.gz .
tar zxvf condor-8.6.3-x86_64_RedHat6-stripped.tar.gz
```
Please add the line below in you your bashrc
```bash
export PYTHONPATH=/uscms/home/USERNAME/condor-8.6.3-x86_64_RedHat6-stripped/lib/python:$PYTHONPATH
```
after,
```bash
source ~/.bashrc
```

### Submitting the Condor Jobs
The Analyser make use of lpc computing infrastructure such as uscms_data space for storing processed ntuple; eos space for condor job output and scratch_3day space for tmp.
the environment variables are configured in PandaAnalyzer/LPC_T3/lpc_setup.sh

To submit jobs, simply do
```bash 
source lpc_setup.sh $ANALYSIS $REGION
```
```bash 
$ANALYSIS can be one of [boosted resolved monojet]
```
and
```
$REGION can be one of [met singleele singlemu diele dimu pho] 
```
e.g. in order to run the job for monojet analysis in met region do the setup as follows:
```
source lpc_setup.sh monojet met
```
In order to prepare the job with grid authentication, go to PandaAnalyzer/LPC_T3/bin/
```bash
sh buildMergedInputs.sh -t -n 40
```
where n = filesetSize and t = doTar

To submit the job:
```bash
python submit.py
```
To check the job progress:
```bash
condor_q USERNAME | tail
```

Output files reside in LPC EOS area. 

## Merging

Make sure `$PANDA_FLATDIR` exists. Then, go into `LPC_T3/merging` and do:
```bash
./merge.sh
```
or
```bash
./merge.py [--cfg CONFIG] [--silent] TTbar_Powheg
```
to merge the Powheg TT sample, for example. 
If provided, `CONFIG` is the module that is imported from `configs/`. 
The default is `common`, but there are others, like `leptonic`.
To merge en-masse (e.g. many many signal outputs), you can do something like:
```bash
submit --exec merge.py --arglist list_of_signals.txt
```
This assumes that you have `PandaCore/bin` in your `$PATH`
The `submit` command will print a cache directory it created to keep track of the jobs.
You can check the status of the jobs by doing
```bash
check --cache <cache_directory> [--resubmit_failed]
```
The last flag is optional and will resubmit anything that failed (exited without code 0).
NB: the `submit` and `check` executables are very generic and don't know anything about the code they are running.
For them, "success" simply means "exited with code 0".
So it is important to check that the output looks sane to you.
