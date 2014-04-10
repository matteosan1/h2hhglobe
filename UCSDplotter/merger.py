import ROOT, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--input", default="", help="Input ROOT file")
parser.add_option("-o", "--output", default="", help="Output ROOT file")
parser.add_option("-d", "--dirname", default="", help="Dirname of input tree")
parser.add_option("-t", "--treename", default="opttree", help="Output tree name")
parser.add_option("-r", "--remove", default="", help="List of tags to remove")
(options, arg) = parser.parse_args()

tags_to_remove = options.remove.split(",")
dirname = options.dirname

#tree = []
list = ROOT.TList()
firstName = ""

f = ROOT.TFile(options.input)
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
            if (toSkip and len(tags_to_remove) > 0):
                continue

    if (className == "TTree"):
        if ((name == "lumi") or (name == "plotvariables") or (name == "inputfiles")):
            continue
        print "Adding ", name

        if (dirname == ""):
            list.Add(f.Get(name))
        else:
            list.Add(mydir.Get(name))

if (list.GetEntries() > 0):
    out = ROOT.TFile(options.output, "recreate")
    out.cd()
    opttree = ROOT.TTree.MergeTrees(list)
    opttree.SetName(options.treename)
    f.Close()
    opttree.Write()
    out.Close()
#  
#
#  out = new TFile(outfilename, "update");
#  TIter nextKey(out->GetListOfKeys());
#  TKey* key;
#  
#  while (key = (TKey*)nextKey()) {
#    TString name(key->GetName());
#    TString className(key->GetClassName());
#    if (className.CompareTo("TTree") == 0) {
#      if (name.CompareTo(treename) != 0) {
#	std::cout << "Cleaning " << key->GetName() << std::endl;
#	key->Delete();
#      } else {
#	key->SetTitle("Opttree");
#      }
#    }
#  }
#
#  out->Close();
#  
#  return 0;
#}
