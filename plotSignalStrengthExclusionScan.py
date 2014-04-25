#!/usr/bin/env python

import sys, os


#----------------------------------------------------------------------
def makeGraph(xvalues, yvalues):

    import ROOT, array

    numPoints = len(xvalues)
    assert len(yvalues) == numPoints

    return ROOT.TGraph(numPoints,

                       array.array('f',xvalues),
                       array.array('f',yvalues),
                       )

#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------

ARGV = sys.argv[1:]

assert len(ARGV) == 1

execfile(ARGV[0])

import ROOT; gcs = []

mg = ROOT.TMultiGraph(); gcs.append(mg)

# sort the values by lambda
limits.sort(key = lambda x: x['lambdaScale'])

lambdaScalings = [ x['lambdaScale'] for x in limits ]


# bands
for quantiles, color in (
    (['2.5%',  '97.5%'], ROOT.kYellow),
    (['16.0%', '84.0%'], ROOT.kGreen),
    ):

    yvalues = [ x['limits']['Expected %5s' % quantiles[0]] for x in limits]
    yvalues += [ x['limits']['Expected %5s' % quantiles[1]] for x in limits[::-1]]

    gr = makeGraph(lambdaScalings + lambdaScalings[::-1],
                   yvalues); gcs.append(gr)

    gr.SetFillColor(color)
    mg.Add(gr, "F")


# median
grMedian = makeGraph(lambdaScalings,
                     [ x['limits']['Expected 50.0%'] for x in limits]); gcs.append(grMedian)
grMedian.SetLineWidth(2)

mg.Add(grMedian,"L")


# draw the multigraph
mg.Draw("A")

histo = mg.GetHistogram()

# draw the line at mu = 1
line = ROOT.TLine(histo.GetXaxis().GetXmin(), 1,
                  histo.GetXaxis().GetXmax(), 1); gcs.append(line)
line.SetLineWidth(2)
line.SetLineStyle(ROOT.kDashed)
line.Draw()

# axis titles
histo.SetXTitle("#lambda / #lambda_{SM}")
histo.SetYTitle("signal strength modifier #mu")

ROOT.gPad.SetGrid()
ROOT.gPad.SaveAs("signalStrengthExclusionScan.png")
