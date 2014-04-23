#!/usr/bin/env python

# generates datacards for combine

import operator, sys, itertools

# using the RooDataHists (e.g. for comparison), can't scale these with lambda
# rootFile = "workspace.root"; mode = "RooDataHist"

# objects with lambda scaling enabled
rootFile = "workspace.root"; mode = "RooHistPdf"

#----------------------------------------------------------------------

def getBins():
    """ finds all bins to loop over (excludes those with no background) """

    retval = []

    import ROOT

    # open the workspace file
    fin = ROOT.TFile.Open(rootFile)
    assert fin != None and fin.IsOpen(), "failed to open workspace file " + rootFile

    ws = fin.Get("CMS_hh_workspace")
    assert ws != None

    for finalState in [ "mm", "em", "ee" ]:

        numNonEmptyCategories = 0

        for cat in itertools.count():

            binName = "{finalState}_cat{cat}".format(**locals())

            # get the RooDataHist objects
            sigDS  = ws.obj("sig_lambdaXX_" + binName)
            bgDS   = ws.obj("bkg_lambdaXX_{finalState}_cat{cat}".format(**locals()))
            dataDS = ws.obj("data_lambdaXX_{finalState}_cat{cat}".format(**locals()))
            
            if sigDS == None or bgDS == None or dataDS == None:
                # this category seems not to exist
                print >> sys.stderr,"found %d categories (%d non-empty) for final state %s" % (cat, numNonEmptyCategories, finalState)
                break

            # check if the background RooDataHist is not empty
            if bgDS.sumEntries() == 0:
                # note that numEntries(..) does not work for this purpose (probably returns the number of bins...)
                continue

            numNonEmptyCategories += 1

            retval.append(dict(cat = cat,
                               catName = "cat%d" % cat,
                               finalState = finalState,
                               binName = finalState + "_cat%d" % cat,
                               ))

    return retval

#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------

allBins = getBins()

if mode == 'RooDataHist':
    sigPdfName = "sig_lambdaXX_$CHANNEL"
    bkgPdfName = "bkg_lambdaXX_$CHANNEL"
elif mode == "RooHistPdf":
    sigPdfName = "pdf_sig_lambdaXX_$CHANNEL"
    bkgPdfName = "pdf_bkg_lambdaXX_$CHANNEL"
else:
    raise Exception("unknown mode " + mode)

print "# datacards"
print "#"
print "---------------------------------------------"
print "imax *"
print "jmax *"
print "kmax *"
print "---------------------------------------------"



print "shapes data_obs * {rootFile}     CMS_hh_workspace:data_lambdaXX_$CHANNEL".format(**globals())
print "shapes bkg *      {rootFile}     CMS_hh_workspace:{bkgPdfName}".format(**globals())
print "shapes sig *      {rootFile}     CMS_hh_workspace:{sigPdfName}".format(**globals())

# observed data
print "---------------------------------------------"
print "bin"," ".join(bin['binName'] for bin in allBins)
print "observation"," ".join([ "-1" ] * len(allBins))

print "---------------------------------------------"

processes = [
    dict(name = "sig", number = 0), # signals must be <= 0
    dict(name = "bkg", number = 1), # backgrounds must be >= 1
    ]

# cartesian product of categories and processes
columns = list(itertools.product(allBins, processes))

print "bin"," ".join(x[0]['finalState'] + "_" + x[0]['catName'] for x in columns)

print "process"," ".join(x[1]['name'] for x in columns)
print "process"," ".join(str(x[1]['number']) for x in columns)

if mode == 'RooDataHist':
    # auto-determination of the number of events by combine
    print "rate"," ".join([ "-1" if x[1]['number'] <= 0 else "-1" for x in columns ])
else:
    print "rate"," ".join([ "1" if x[1]['number'] <= 0 else "1" for x in columns ]) 

#----------------------------------------
# nuisances
#----------------------------------------

# for the moment affects signal AND background
# DEBUG commented out

print "lumi_8TeV lnN"," ".join([ "1.044000" ] * len(columns) )

