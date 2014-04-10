import ROOT
import array, sys
from  makeWS import WSProducer
from optparse import OptionParser

hMET  = []
hMass = []
hBtag = []
colors = (ROOT.kBlack, ROOT.kRed, ROOT.kBlue-10, ROOT.kBlue, ROOT.kViolet, ROOT.kGreen, ROOT.kYellow, ROOT.kOrange)

def processCats(itype):
    Ztypes     = [1, 2]
    ZZtypes    = [20, 21]
    TTtypes    = [10, 11, 12, 13, 14, 15]
    WHtypes    = [40, 41, 42]
    Wtypes     = range(30, 40)
    WJetstypes = [50, 51, 52, 53]
    
    if (itype == 0):
        return 0

    if (itype < 0):
        return 1

    if (itype in Ztypes):
        return 2

    if (itype in ZZtypes):
        return 3
    
    if (itype in TTtypes):
        return 4

    if (itype in WHtypes):
        return 5

    if (itype in Wtypes):
        return 6

    if (itype in WJetstypes):
        return 7
    
def histograms(cats):
    global hMET, hMass, hBtag
    for i in xrange(cats):
        for j in xrange(3):
            hMET.append(ROOT.TH1F("hmet"+str(i)+"_typ"+str(j), "hmet"+str(i)+"_typ"+str(j), 100, 0, 200))
            hMass.append(ROOT.TH1F("hmass"+str(i)+"_typ"+str(j), "hmass"+str(i)+"_typ"+str(j), 50, 0, 200))
            hBtag.append(ROOT.TH1F("hbtag"+str(i)+"_typ"+str(j), "hbtag"+str(i)+"_typ"+str(j), 100, -1, 1))

    
def makePlots(processCats, pairCats, isBlind):
    global hMET, hMass, hBtag
    input = ROOT.TFile("sslept_output.root")
    for i in xrange(processCats):
        for j in xrange(pairCats):
            hMET.append(input.Get("hmet"+str(i)+"_typ"+str(j)))
            hMET[-1].SetFillColor(colors[i])
            hMass.append(input.Get("hmass"+str(i)+"_typ"+str(j)))
            hMass[-1].SetFillColor(colors[i])
            hBtag.append(input.Get("hbtag"+str(i)+"_typ"+str(j)))
            hBtag[-1].SetFillColor(colors[i])

    type_range = range(2, processCats)

    canvases = []
    stacks = []
    for j in xrange(pairCats):
        canvases.append(ROOT.TCanvas("c"+str(j), "c"+str(j)))
        hMass[j+3].Scale(10)
        hMass[j+3].SetFillStyle(0)
        hMass[j+3].SetLineWidth(2)
        hMass[j+3].SetLineColor(colors[1])        
        hMass[j+3].Draw("HIST")

        stacks.append(ROOT.THStack("s"+str(j), "s"+str(j)))
        for i in type_range:
            stacks[-1].Add(hMass[j+i*3])

        stacks[-1].Draw("SAME")
        hMass[j+3].Draw("HISTSAME")
        canvases[-1].SaveAs("mass"+str(j)+".png")
        
        canvases.append(ROOT.TCanvas("cmet"+str(j), "cmet"+str(j)))
        stacks.append(ROOT.THStack("smet"+str(j), "smet"+str(j)))
        for i in type_range:
            stacks[-1].Add(hMET[j+i*3])

        stacks[-1].Draw()
        hMET[j+3].Scale(10)
        hMET[j+3].SetFillStyle(0)
        hMET[j+3].SetLineWidth(2)
        hMET[j+3].SetLineColor(colors[1])        
        hMET[j+3].Draw("HISTSAME")
        canvases[-1].SaveAs("met"+str(j)+".png")

        canvases.append(ROOT.TCanvas("cbtag"+str(j), "cbtag"+str(j)))
        stacks.append(ROOT.THStack("sbtag"+str(j), "sbtag"+str(j)))
        for i in type_range:
            stacks[-1].Add(hBtag[j+i*3])

        stacks[-1].Draw()
        hBtag[j+3].Scale(10)
        hBtag[j+3].SetFillStyle(0)
        hBtag[j+3].SetLineWidth(2)
        hBtag[j+3].SetLineColor(colors[1])        
        hBtag[j+3].Draw("HISTSAME")
        #canvases[-1].SaveAs("btag"+str(j)+".png")
    input.Close()
    
parser = OptionParser()
parser.add_option("-p", "--plot", default=False, action="store_true", help="Plot the histograms")
parser.add_option("-b", "--blind", default=False, action="store_true", help="Do not plot data")
(options, arg) = parser.parse_args()

if (options.plot):
    makePlots(8, 3, options.blind)
    sys.exit()

wsProducer = WSProducer()
wsProducer.prepareDataSets(3)
histograms(8)

file = ROOT.TFile("UCSDplotter/sslep_v4.root")
tree = file.Get("opttree")
tree.SetBranchStatus("*",0)
tree.SetBranchStatus("itype",1)
tree.SetBranchStatus("weight",1)
tree.SetBranchStatus("pairs", 1)
tree.SetBranchStatus("mass",1)
tree.SetBranchStatus("type",1)
tree.SetBranchStatus("met",1)
tree.SetBranchStatus("njets",1)
tree.SetBranchStatus("btag", 1)

entries = tree.GetEntries()

for z in xrange(entries):
    tree.GetEntry(z)
    processCat = processCats(tree.itype)
    if (tree.njets < 2):
        continue
    pairs = tree.pairs
    for pair in xrange(pairs):
        hMET[tree.type[pair]+processCat*3].Fill(tree.met, tree.weight)
    
        if (tree.type[pair]==0 and tree.met < 50.):
            continue
        if (tree.type[pair]>0 and tree.met < 50.):
            continue 

        if (tree.btag[0] > 0.6 or tree.btag[1] > 0.6):
            continue

        if (tree.mass[pair] > 8.):
            hMass[tree.type[pair]+processCat*3].Fill(tree.mass[pair], tree.weight)
            hBtag[tree.type[pair]+processCat*3].Fill(tree.btag[0], tree.weight)
            hBtag[tree.type[pair]+processCat*3].Fill(tree.btag[1], tree.weight)
            wsProducer.fillDataset(tree.itype, tree.type[pair], tree.cat[pair], tree.mass[pair], tree.weight)

wsProducer.saveWS()

output = ROOT.TFile("sslept_output.root", "recreate")
for h in hMET:
    h.Write()

for h in hMass:
    h.Write()

for h in hBtag:
    h.Write()

output.Close()
