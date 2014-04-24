import ROOT
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--input", default="", help="Input ROOT file")
parser.add_option("-o", "--output", default="", help="Output ROOT file")
parser.add_option("-d", "--dirname", default="", help="Dirname of input tree")
parser.add_option("-t", "--treename", default="opttree", help="Output tree name")
parser.add_option("-r", "--remove", default="", help="List of tags to remove")
(options, arg) = parser.parse_args()

def mergeOther(inputFile, listOfTrees):
    global options
    outFileName = options.output
    
    print "Adding ", listOfTrees[0]
    firstTree = inputFile.Get(listOfTrees[0])
    firstTree.SetBranchStatus('*',1)
    
    outFile = ROOT.TFile(outFileName, 'recreate')
    outTree = firstTree.CloneTree()
    outTree.SetName(options.treename)
    
    for elem in listOfTrees[1:]:
        print "Adding ", elem
        ttree = inputFile.Get(elem)
        ttree.SetBranchStatus('*', 1)
        outTree.CopyEntries(ttree, -1, "FAST")
    
    outFile.cd()
    outTree.Write()
    outFile.Close()

tags_to_remove = []
if (options.remove):
    tags_to_remove = options.remove.split(",")
dirname = options.dirname

list = []
f = ROOT.TFile(options.input, "read")
ROOT.gROOT.cd()
mydir = 0

if(dirname != ""):
    f.cd(dirname)
    mydir = f.Get(dirname)
  
keys = f.GetListOfKeys()
if (dirname != ""):
    keys = mydir.GetListOfKeys()

for k in keys:
    name = k.GetName()
    className = k.GetClassName()

    if (len(tags_to_remove) > 0):
        toSkip = False
        for t in tags_to_remove:
            if (t in name):
                toSkip = True
                break
        if (toSkip):
            continue

    if (className == "TTree"):
        if ((name == "lumi") or (name == "plotvariables") or (name == "inputfiles")):
            continue
        if (dirname == ""):
            list.append(name)
        else:
            list.Add(name)

mergeOther(f, list)
f.cd()
f.Close()
