import ROOT

# the garbage collection saver (to avoid python deleting ROOT objects which ROOT still uses internally)
gcs = []
#import sys

#----------------------------------------------------------------------
# parameters
#----------------------------------------------------------------------

# coefficients for cross section polynomial

# 0.0025155441561416935 * lambda_ratio^2 -0.011835060325668232 * lambda_ratio + 0.017223395762987696
lambdaPolyCoeffs = [ 0.0025155441561416935, # lambda_ratio^2
                     -0.011835060325668232, # lambda_ratio
                     0.017223395762987696,  # 1
                     ]

# branching ratio at 125.0 GeV, taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3
brHtoWW =  2.15E-01 

# at lambda/lambda_SM = 70, HPAIR gives: 11.5149189 pb
# the signal cross section used to create this workspace 
# (depends on lambda)
usedXsectBR = 11.5149189 * brHtoWW**2

#----------------------------------------------------------------------
def lambdaPoly(lambdaScaling):
    """ returns the cross section in pb at a given lambda scaling"""

    return sum([
        lambdaScaling ** 2 * lambdaPolyCoeffs[0],
        lambdaScaling      * lambdaPolyCoeffs[1],
                             lambdaPolyCoeffs[2],
        ])
#----------------------------------------------------------------------

class WSProducer:
    mass = ROOT.RooRealVar("CMS_hh_mass", "CMS_hh_mass", 800, 0, 400)
    workspace = ROOT.RooWorkspace("CMS_hh_workspace") 
    datasets  = {}
    datahists = {}

    # other objects to import into the workspace
    otherObjectsToImport = []

    set = ROOT.RooArgSet("set")
    set.add(mass)
    process_type = ["mm", "ee", "em"]

    #----------------------------------------
    
    def prepareDataSets(self, cats):
        # keep the number of categories
        self.numCategories = cats

        for i in xrange(cats):
            for p in self.process_type:
                name = "sig_lambdaXX_"+p+"_cat"+str(i)
                self.addDataHist(name)
                name = "bkg_lambdaXX_"+p+"_cat"+str(i)
                self.addDataHist(name)
                name = "data_lambdaXX_"+p+"_cat"+str(i)
                self.addDataSet(name)

    #----------------------------------------
        
    def addDataSet(self, name):
        self.datasets[name] = ROOT.RooDataHist(name, name, self.set, "weight")

    #----------------------------------------

    def addDataHist(self, name):
        self.datahists[name] = ROOT.RooDataHist(name, name, self.set)

    #----------------------------------------
        
    def saveWS(self):
        getattr(self.workspace, 'import')(self.mass)
        for d in self.datasets.values():
            getattr(self.workspace, 'import')(d)

        for h in self.datahists.values():
            getattr(self.workspace, 'import')(h)

        for obj in self.otherObjectsToImport:
            getattr(self.workspace, 'import')(obj, ROOT.RooFit.RecycleConflictNodes())

        self.workspace.writeToFile("workspace.root")

    #----------------------------------------

    def finalize(self):
        """ adds more logic (e.g. scaling of cross section with lambda)
        to the workspace. Should be called after filling is complete
        and before saving the workspace """

        #----------
        # signal normalization
        #----------

        # create the scaled lambda variable
        # we should also let it go to negative values in the future
        self.lambdaScalingVar = ROOT.RooRealVar("lambdaScaling", "lambda / lambda_SM", 70, 0, 200) 

        # avoid that the fit in combine moves this around
        self.lambdaScalingVar.setConstant(True)

        # BR(H->WW) constant
        self.brHWWvar = ROOT.RooConstVar("brHWW", "BR(H->WW)", brHtoWW)

        # expression for the scaled cross section
        tmp = ROOT.RooArgList(self.lambdaScalingVar)
        self.xsectVar = ROOT.RooFormulaVar("xsect", "lambda dependent cross section",
                              "%f * lambdaScaling * lambdaScaling %+f * lambdaScaling %+f" % tuple(lambdaPolyCoeffs),
                              tmp)

        # cross section times BR(H->WW)^2
        tmp = ROOT.RooArgList(self.xsectVar, self.brHWWvar)
        self.xsbrVar = ROOT.RooFormulaVar("xsbr","cross section times BR(H->WW)^2",
                                          "@0 * @1 * @1",
                                          tmp)

        # variable for the cross section which was used to create this workspace
        self.usedXsectBRVar = ROOT.RooConstVar("usedXsectBR",
                                             "cross section  * BR(H->WW)**2 to produce the initial workspace",
                                             usedXsectBR)

        # multiply with the (lambda dependent) cross section, divide by the used cross section
        # (note that this does NOT depend on the channel)
        tmp = ROOT.RooArgList(self.xsbrVar, self.usedXsectBRVar)
        self.scaleFactorVar = ROOT.RooFormulaVar("scaleFactor",
                                                 "scale factor for cross section",
                                                 "xsbr / usedXsectBR",
                                                 tmp
                                                 )

        for catNum in range(self.numCategories):
            cat = "cat%d" % catNum 

            for finalState in self.process_type:

                binName = "{finalState}_{cat}".format(**locals())

                print "binName=",binName

                #--------------------
                # create scaled RooHistPdfs for signal
                #--------------------
                name = "sig_lambdaXX_" + binName
                sigHisto = self.datahists[name]

                # create a RooHistPdf from that

                tmp = ROOT.RooArgSet(self.mass)
                sigHistPdf = ROOT.RooHistPdf("pdf_" + name,
                                             "pdf_" + name,
                                             tmp,
                                             sigHisto,
                                             )

                # create a variable with the actual normalization
                # (which is probably lost when converting the RooDataHist)

                numEventsVar = ROOT.RooConstVar("origNumEvents_" + binName,
                                                "original num signal events " + binName,
                                                sigHisto.sumEntries()); gcs.append(numEventsVar)

                tmp = ROOT.RooArgList(numEventsVar, self.scaleFactorVar)                                        
                scaledNumEventsVar = ROOT.RooFormulaVar(sigHistPdf.GetName() + "_norm", # combine naming convention
                                                        "scaled number of expected events " + binName,
                                                        "@0 * @1",
                                                        tmp); gcs.append(scaledNumEventsVar)

                self.otherObjectsToImport.append(sigHistPdf)
                self.otherObjectsToImport.append(scaledNumEventsVar)

                #--------------------
                # create RooHistPdfs for background
                #--------------------
                name = "bkg_lambdaXX_" + binName
                bkgHisto = self.datahists[name]

                assert bkgHisto != None

                tmp = ROOT.RooArgSet(self.mass)
                bkgHistPdf = ROOT.RooHistPdf("pdf_" + name,
                                             "pdf_" + name,
                                             tmp,
                                             bkgHisto,
                                             )

                self.otherObjectsToImport.append(bkgHistPdf)

                # create a norm variable
                normVar = ROOT.RooConstVar(bkgHistPdf.GetName() + "_norm",
                                           "background normalization " + binName,
                                           bkgHisto.sumEntries()); gcs.append(normVar)

                self.otherObjectsToImport.append(normVar)

            # end of loop over final states
        # end of loop over categories

    #----------------------------------------

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

    #----------------------------------------            
