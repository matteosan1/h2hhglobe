import ROOT
#import sys

class WSProducer:
    mass = ROOT.RooRealVar("CMS_hh_mass", "CMS_hh_mass", 800, 0, 400)
    workspace = ROOT.RooWorkspace("CMS_hh_workspace") 
    datasets  = {}
    datahists = {}
    set = ROOT.RooArgSet("set")
    set.add(mass)
    process_type = ["mm", "ee", "em"]
    
    def prepareDataSets(self, cats):
        for i in xrange(cats):
            for p in self.process_type:
                name = "sig_lambdaXX_"+p+"_cat"+str(i)
                self.addDataHist(name)
                name = "bkg_lambdaXX_"+p+"_cat"+str(i)
                self.addDataHist(name)
                name = "data_lambdaXX_"+p+"_cat"+str(i)
                self.addDataSet(name)
        
    def addDataSet(self, name):
        self.datasets[name] = ROOT.RooDataHist(name, name, self.set, "weight")

    def addDataHist(self, name):
        self.datahists[name] = ROOT.RooDataHist(name, name, self.set)
        
    def saveWS(self):
        getattr(self.workspace, 'import')(self.mass)
        for d in self.datasets.values():
            getattr(self.workspace, 'import')(d)

        for h in self.datahists.values():
            getattr(self.workspace, 'import')(h)

        self.workspace.writeToFile("workspace.root")

    def fillDataset(self, itype, process, cat, mass, weight):
        self.mass.setVal(mass)

        if (itype == 0):
            name = "data_lambdaXX_"+self.process_type[process]+"_cat"+str(cat)
            self.datasets[name].add(self.set, weight)
        elif (itype > 0):
            name = "bkg_lambdaXX_"+self.process_type[process]+"_cat"+str(cat)
            self.datahists[name].add(self.set, weight)
        else:
            name = "sig_lambdaXX_"+self.process_type[process]+"_cat"+str(cat)
            self.datahists[name].add(self.set, weight)

            
