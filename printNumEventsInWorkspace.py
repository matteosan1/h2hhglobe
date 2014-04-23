#!/usr/bin/env python

import sys, itertools


#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------
ARGV = sys.argv[1:]

assert len(ARGV) > 0

import ROOT; gcs = []

for fname in ARGV:
    fin = ROOT.TFile.Open(fname)
    assert fin != None
    assert fin.IsOpen()

    ws = fin.Get("CMS_hh_workspace")

    for finalState in [ "mm", "em", "ee" ]:
    
        for cat in itertools.count():
            binName = "{finalState}_cat{cat}".format(**locals())

            sigDS = ws.obj("sig_lambdaXX_" + binName)
            bgDS = ws.obj("bkg_lambdaXX_{finalState}_cat{cat}".format(**locals()))
            dataDS = ws.obj("data_lambdaXX_{finalState}_cat{cat}".format(**locals()))

            if sigDS == None or bgDS == None or dataDS == None:
                # this category seems not to exist
                break

            numSig = sigDS.sumEntries()

            # DEBUG
            # numSig = ws.obj("pdf_sig_lambdaXX_" + binName + "_norm").getVal()

            numBkg = bgDS.sumEntries()
            
            print "%-20s: sig: %7.1f bkg: %7.1f obs: %7f" % ("{finalState}_cat{cat}".format(**locals()),
                                                             numSig,
                                                             numBkg,
                                                             dataDS.sumEntries()
                                                             ),

            if numBkg > 0:
                print "s/b: %7.2f" % (numSig/numBkg), 
            print

            
        # end of loop over categories    
    # end of loop over final states    
