#ifndef __SSLeptonAnalysis__
#define __SSLeptonAnalysis__

#include "PhotonAnalysis/interface/StatAnalysis.h"
#include "SSLeptonAnalysis/interface/MuScleFitCorrector.h"
//#include "JetAnalysis/interface/JetHandler.h"
#include "TRegexp.h"

#include <string>

class SSLeptonAnalysis : public StatAnalysis {
 public:
  SSLeptonAnalysis();
  virtual ~SSLeptonAnalysis();
  
  virtual const std::string & name() const { return name_; };
  
  // LoopAll analysis interface implementation
  virtual void Init(LoopAll&);
  virtual void Term(LoopAll&);
  
  virtual void ReducedOutputTree(LoopAll &l, TTree *);
  virtual void GetBranches(TTree *, std::set<TBranch *>& );
  
  virtual void ResetAnalysis();
  
  virtual void FillReductionVariables(LoopAll& l, int jentry);   
  virtual bool SelectEventsReduction(LoopAll&, int);
  
  virtual bool SkimEvents(LoopAll&, int);
  virtual bool SelectEvents(LoopAll&, int);
  virtual bool Analysis(LoopAll&, Int_t);
  
  bool checkEventHLT(LoopAll& l, std::vector<std::string> paths);
  void Tree(LoopAll& l, Int_t lept1, Int_t lept2, const TLorentzVector & Higgs, Int_t cat, Int_t vbfcat, Float_t weight, Float_t pu_weight, bool isSyst, std::string name1, bool* jetid_flags);
  void FillRooContainer(LoopAll& l, int cur_type, float mass, int category, float weight);
  void FillSignalLabelMap(LoopAll & l);
  void buildBkgModel(LoopAll& l, const std::string & postfix);
  std::string GetSignalLabel(int id);
  
  int nCategories_;
  float massMin,massMax;
  int nDataBins;
  std::vector<int> bkgPolOrderByCat;
  std::map<int,std::string> signalLabels;
  std::string muFitParametersFile;
  bool doBlinding;
  bool dataIs2011;
  
  std::string jetHandlerCfg;
  //JetHandler* jethandler_;
  
 protected:
  std::string name_;
  std::vector<std::string> hltSelection;
  MuScleFitCorrector* muCorrector_;

};

#endif
