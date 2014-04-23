#!/usr/bin/env python

import sys, os, tempfile, commands, re
from pprint import pprint

# run combine to determine the limit on the signal strength
# modifier at various points of lambda/lambda_SM

# lambdaScales = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] + \
#                [1,2,3,4,5,6,7,8,9,10] + \
#                [20,30,40,50,60,70]

lambdaScales = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] + \
               [1,1.2,1.4,1.6,1.8,
                2,2.2, 2.4, 2.6, 2.8,
                3,3.2, 3.4, 3.6, 3.8,
                4,4.2, 4.4, 4.6, 4.8,
                5, 5.2, 5.4, 5.6, 5.8]

workspaceFname = "workspace.root"

datacardFname = "datacards.txt"

#----------------------------------------------------------------------
def runCommand(cmd):
    print >> sys.stderr,"running",cmd
    res =  os.system(cmd)
    assert res == 0, "failed to run command " + cmd


#----------------------------------------------------------------------

def parseCombineOutput(lines):
    
    lines = lines[:]

    while lines:
        line = lines.pop(0).strip()
        if line == "-- Asymptotic --":
            break

    if not lines:
        return None

    retval = {}

    for line in lines:

        mo = re.match("\s*([^:]+): r < (\S+)\s*$",line)

        if not mo:
            break
        
        retval[mo.group(1)] = float(mo.group(2))

    return retval    


#----------------------------------------------------------------------
# old-fashioned way: just copy the workspace file, modify the lambda scaling by hand
# (but at least we know what is going on...)

def runCombineStandard(lambdaScale):
    # create a working directory
    workdir = tempfile.mkdtemp()

    # the output of text2workspace.py
    compactWorkspace = os.path.join(workdir, "ws.root")

    # run text2workspace first

    cmdParts = [
        "text2workspace.py",
        datacardFname,
        '-m 125',
        "-o " + compactWorkspace
        ]

    runCommand(" ".join(cmdParts))

    #----------
    # set the lambda scaling value
    #----------
    import ROOT
    ROOT.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
    fin = ROOT.TFile(compactWorkspace,"UPDATE")
    ws = fin.Get("w")
    assert ws != None

    lambdaScalingVar = ws.obj("lambdaScaling")
    assert lambdaScalingVar != None

    lambdaScalingVar.setVal(lambdaScale)

    # also set the range of r
    # rVar = ws.obj("r")
    # assert rVar != None
    # 
    # rVar.setMax(0.1)

    fin.cd()
    ws.Write()
    fin.Close()
    ROOT.gROOT.cd()
    #----------    

    # only run the expected limit for the moment
    cmdParts = [
        'cd ' + workdir,
        " && ",
        'combine',
        '-M Asymptotic',
        
        compactWorkspace,

        "-m 125",
        "--expectSignal 0",  # assuming background only (zero signal injected)
        '-t -1',              # create an Asimov dataset

        '> log 2>&1',
        ]
    
    runCommand(" ".join(cmdParts))

    result = parseCombineOutput(open(os.path.join(workdir, "log")).read().splitlines())

    return result
                 
                 

#----------------------------------------------------------------------
# main
#----------------------------------------------------------------------

# this does not give errors as the CLs is already < 0.05 ?! 
# runCombineStandard(70)
# runCombineStandard(50)

# produces the error values
# runCombineStandard(20)

# median expected r excluded at r = 0.46
# runCombineStandard(5)

# median expected r excluded at r = 2.68
# runCombineStandard(2)

allResults = []

for lambdaScale in lambdaScales:
    thisResult = runCombineStandard(lambdaScale)

    allResults.append(dict(lambdaScale = lambdaScale,
                           limits = thisResult))



# write the results to an output file
import time
outputFname = time.strftime("signal-strength-scan-%Y-%m-%d:%H:%M.py")
fout = open(outputFname,"w")
print >> fout,"limits=",allResults
fout.close()
print >> sys.stderr,"wrote results to",outputFname 
