#!/usr/bin/env python

from PandaCore.Tools.Load import Load 
from PandaCore.Tools.Misc import PInfo, PError
import ROOT as root

Load('PandaAnalyzer')

def _dump(a):
    PInfo('PandaAnalysis.Flat.analysis','Summary of analysis %s:'%(a.name))
    for k in dir(a):
        if k[0] == '_':
            continue
        if type(getattr(a, k)) != int:
            continue
        PInfo('PandaAnalysis.Flat.analysis','    %20s = %s'%(k, 'True' if getattr(a, k) else 'False'))



def _analysis(name, verbose, **kwargs):
    a = root.Analysis(name)
    for k,v in kwargs.iteritems():
        if not hasattr(a, k):
            PError('PandaAnalysis.Flat.analysis','Could not set property %s'%k)
            return None 
        setattr(a, k, bool(v))
    setattr(a, 'dump', lambda : _dump(a))
    if verbose:
        a.dump()
    return a

def analysis(name, **kwargs):
    return _analysis(name, verbose=True, **kwargs)


# predefined!
boosted = lambda v=False : _analysis(
        name = 'boosted',
        verbose = v,
        boosted = True,
    )

resolved = lambda v=False : _analysis(
        name = 'resolved',
        verbose = v,
    )

monojet = lambda v=False : _analysis(
        name = 'monojet',
        verbose = v,
    )

lepmonotop = lambda v=False : _analysis(
        name = 'lepmonotop',
        verbose = v,
    )

monotop = lambda v=False : _analysis(
        name = 'monotop',
        verbose = v,
    )

vbf = lambda v=False : _analysis(
        name = 'vbf',
        verbose = v,
        vbf = True,
        fatjet = False,
        btagSFs = False,
        puppi_jets = False
    )

gghbb = lambda v=False : _analysis(
        name = 'gghbb',
        verbose = v,
        recoil = False,
        ak8 = True,
    )

deep = lambda v=False : _analysis(
        name = 'deep',
        verbose = v,
        ak8 = True,
        deep = True,
        deepTracks = True,
        deepSVs = True,
        deepAntiKtSort = True,
        btagSFs = False,
        jetFlavorPartons = False,
    )

deepgen = lambda v=False : _analysis(
        name = 'deepgen',
        verbose = v,
        ak8 = True,
        btagSFs = False,
        deepGen = True,
        genOnly = True,
    )

wlnhbb = lambda v=False : _analysis(
        name = 'wlnhbb',
        verbose = v,
        ak8 = True,
        hbb = True,
        monoh = False,
        recoil = True,
        fatjet = True,
        btagSFs = True,
        btagWeights = True,
        useCMVA = True,
        complicatedLeptons = True,
        hfCounting = True,
        reclusterGen = True,
        bjetRegression = True,
        varyJES = True,
        rerunJES = True,
        jetFlavorPartons = False,
        jetFlavorJets = True,
    )
wlnhbb_ca15 = lambda v=False : _analysis(
        name = 'wlnhbb',
        verbose = v,
        ak8 = False,
        hbb = True,
        monoh = False,
        recoil = True,
        fatjet = True,
        btagSFs = True,
        btagWeights = True,
        useCMVA = True,
        complicatedLeptons = True,
        hfCounting = True,
        reclusterGen = True,
        bjetRegression = True,
        varyJES = True,
        rerunJES = True,
        jetFlavorPartons = False,
        jetFlavorJets = True,
    )

vv = lambda v=False : _analysis(
        name = 'vv',
        verbose = v,
        recoil = False,
        fatjet = False,
        btagSFs = True,
        complicatedLeptons = True,
        complicatedPhotons = True,
        applyMCTriggers = True,
        varyJES = True,
        jetFlavorPartons = False,
        jetFlavorJets = True,
    )
