#include "SSLeptonAnalysis/interface/SSLeptonAnalysis.h"
#include <iostream>

#define SSLeptonAnalysisDEBUG 0

using namespace std;

// ----------------------------------------------------------------------------------------------------
SSLeptonAnalysis::SSLeptonAnalysis()  : 
  name_("SSLeptonAnalysis")
{}

// ----------------------------------------------------------------------------------------------------
SSLeptonAnalysis::~SSLeptonAnalysis() 
{}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::Term(LoopAll& l) 
{}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::Init(LoopAll& l) {
  
  if (l.typerun != l.kReduce) {
    if (puHist != "" && puHist != "auto" ) {
      if(DEBUG)
        cout << "Opening PU file"<<endl;
      TFile* puFile = TFile::Open( puHist );
      if (puFile) {
        TH1 * target = 0;
	
        if( puTarget != "" ) {
          TFile * puTargetFile = TFile::Open( puTarget );
          assert( puTargetFile != 0 );
          target = (TH1*)puTargetFile->Get("pileup");
          if( target == 0 ) { target = (TH1*)puTargetFile->Get("target_pu"); }
          target->Scale( 1. / target->Integral() );
        }
	
        if( puMap != "" ) {
          loadPuMap(puMap, puFile, target);
        } else {
	  loadPuWeights(0, puFile, target);
        }
        puFile->Close();
      }
      else {
        cout<<"Error opening " <<puHist<<" pileup reweighting histogram, using 1.0"<<endl;
        weights[0].resize(50);
        for (unsigned int i=0; i<weights[0].size(); i++) weights[0][i] = 1.0;
      }
      if(DEBUG)
        cout << "Opening PU file END"<<endl;
    } else if ( puHist == "auto" ) {
      TFile * puTargetFile = TFile::Open( puTarget );
      assert( puTargetFile != 0 );
      puTargetHist = (TH1*)puTargetFile->Get("pileup");
      if( puTargetHist == 0 ) {
        puTargetHist = (TH1*)puTargetFile->Get("target_pu");
      }
      puTargetHist = (TH1*)puTargetHist->Clone();
      puTargetHist->SetDirectory(0);
      puTargetHist->Scale( 1. / puTargetHist->Integral() );
      puTargetFile->Close();
    }
  }
  
  
  if (l.typerun == l.kReduce)
    muCorrector_ = new MuScleFitCorrector(muFitParametersFile);
  
  if (l.typerun != l.kReduce) {
    // setup roocontainer
    FillSignalLabelMap(l);
    
    l.rooContainer->BlindData(doBlinding);
    l.rooContainer->AddGlobalSystematic("lumi",1.044,1.00);
    l.rooContainer->SetNCategories(nCategories_);
    l.rooContainer->AddObservable("CMS_hll_mass" ,massMin,massMax);
    l.rooContainer->AddConstant("IntLumi",l.intlumi_);
    
    //// SM Model
    //l.rooContainer->AddConstant("XSBR_ggh_125",0.0350599);
    //l.rooContainer->AddConstant("XSBR_vbf_125",0.00277319);
    //l.rooContainer->AddConstant("XSBR_wzh_125",0.002035123);
    //l.rooContainer->AddConstant("XSBR_tth_125",0.000197718);
    
    l.rooContainer->CreateDataSet("CMS_hll_mass","data_mass"    ,nDataBins);
    l.rooContainer->CreateDataSet("CMS_hll_mass","bkg_mass"     ,nDataBins);
    
    for(size_t isig=0; isig<sigPointsToBook.size(); ++isig) {
      int sig = sigPointsToBook[isig];
      l.rooContainer->CreateDataSet("CMS_hll_mass",Form("sig_ggh_mass_m%d",sig),nDataBins);
      l.rooContainer->CreateDataSet("CMS_hll_mass",Form("sig_vbf_mass_m%d",sig),nDataBins);
      l.rooContainer->CreateDataSet("CMS_hll_mass",Form("sig_rsg_mass_m%d",sig),nDataBins);
    }
    
    std::string postfix=(dataIs2011?"":"_8TeV");
    // build the model
    buildBkgModel(l, postfix);
  }
  
  ////
  //// Jet Handler for sorting out jet energies
  ////
  //
  //if( recomputeBetas || recorrectJets || rerunJetMva || recomputeJetWp || applyJer || applyJecUnc || l.typerun != l.kFill ) {
  //std::cout << "JetHandler: \n"
  //	    << "recomputeBetas " << recomputeBetas << "\n"
  //	    << "recorrectJets " << recorrectJets << "\n"
  //	    << "rerunJetMva " << rerunJetMva << "\n"
  //	    << "recomputeJetWp " << recomputeJetWp
  //	    << std::endl;
  //jetHandler_ = new JetHandler(jetHandlerCfg, l);
  //}
  
  hltSelection.push_back("HLT_Mu17_Mu8_v*");
  hltSelection.push_back("HLT_IsoMu24_v*");
  hltSelection.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
  hltSelection.push_back("HLT_DoubleMu8_Mass8_v*");
  hltSelection.push_back("HLT_DoubleMu14_Mass8_v*");
  hltSelection.push_back("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass_v*");
  hltSelection.push_back("HLT_DoubleEle14_Mass8_CaloIdT_TrkIdVL_pfMHT40_v*");
  hltSelection.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
  hltSelection.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoV_v*");
  hltSelection.push_back("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass*_v*");
  hltSelection.push_back("HLT_Mu14_Ele14_Mass8_CaloIdT_CaloIsoVL_TrkIdVL_pfMHT40_v*");
}
	      
void SSLeptonAnalysis::buildBkgModel(LoopAll& l, const std::string& postfix) {
  // sanity check
  if( bkgPolOrderByCat.size() != nCategories_ ) {
    std::cout << "Number of categories not consistent with specified background model " << nCategories_ << " " << bkgPolOrderByCat.size() << std::endl;
    assert( 0 );
  }
  
  
  l.rooContainer->AddRealVar("CMS_hll_pol6_0"+postfix,-0.1,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_pol6_1"+postfix,-0.1,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_pol6_2"+postfix,-0.1,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_pol6_3"+postfix,-0.01,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_pol6_4"+postfix,-0.01,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_pol6_5"+postfix,-0.01,-1.0,1.0);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol6_0"+postfix,"@0*@0","CMS_hll_pol6_0"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol6_1"+postfix,"@0*@0","CMS_hll_pol6_1"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol6_2"+postfix,"@0*@0","CMS_hll_pol6_2"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol6_3"+postfix,"@0*@0","CMS_hll_pol6_3"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol6_4"+postfix,"@0*@0","CMS_hll_pol6_4"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol6_5"+postfix,"@0*@0","CMS_hll_pol6_4"+postfix);
  l.rooContainer->AddRealVar("CMS_hll_pol5_0"+postfix,-0.1,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_pol5_1"+postfix,-0.1,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_pol5_2"+postfix,-0.1,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_pol5_3"+postfix,-0.01,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_pol5_4"+postfix,-0.01,-1.0,1.0);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol5_0"+postfix,"@0*@0","CMS_hll_pol5_0"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol5_1"+postfix,"@0*@0","CMS_hll_pol5_1"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol5_2"+postfix,"@0*@0","CMS_hll_pol5_2"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol5_3"+postfix,"@0*@0","CMS_hll_pol5_3"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modpol5_4"+postfix,"@0*@0","CMS_hll_pol5_4"+postfix);

  l.rooContainer->AddRealVar("CMS_hll_quartic0"+postfix,-0.1,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_quartic1"+postfix,-0.1,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_quartic2"+postfix,-0.1,-1.0,1.0);
  l.rooContainer->AddRealVar("CMS_hll_quartic3"+postfix,-0.01,-1.0,1.0);
  l.rooContainer->AddFormulaVar("CMS_hll_modquartic0"+postfix,"@0*@0","CMS_hll_quartic0"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modquartic1"+postfix,"@0*@0","CMS_hll_quartic1"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modquartic2"+postfix,"@0*@0","CMS_hll_quartic2"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modquartic3"+postfix,"@0*@0","CMS_hll_quartic3"+postfix);

  l.rooContainer->AddRealVar("CMS_hll_quad0"+postfix,-0.1,-1.5,1.5);
  l.rooContainer->AddRealVar("CMS_hll_quad1"+postfix,-0.01,-1.5,1.5);
  l.rooContainer->AddFormulaVar("CMS_hll_modquad0"+postfix,"@0*@0","CMS_hll_quad0"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modquad1"+postfix,"@0*@0","CMS_hll_quad1"+postfix);

  l.rooContainer->AddRealVar("CMS_hll_cubic0"+postfix,-0.1,-1.5,1.5);
  l.rooContainer->AddRealVar("CMS_hll_cubic1"+postfix,-0.1,-1.5,1.5);
  l.rooContainer->AddRealVar("CMS_hll_cubic2"+postfix,-0.01,-1.5,1.5);
  l.rooContainer->AddFormulaVar("CMS_hll_modcubic0"+postfix,"@0*@0","CMS_hll_cubic0"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modcubic1"+postfix,"@0*@0","CMS_hll_cubic1"+postfix);
  l.rooContainer->AddFormulaVar("CMS_hll_modcubic2"+postfix,"@0*@0","CMS_hll_cubic2"+postfix);

  l.rooContainer->AddRealVar("CMS_hll_lin0"+postfix,-0.01,-1.5,1.5);
  l.rooContainer->AddFormulaVar("CMS_hll_modlin0"+postfix,"@0*@0","CMS_hll_lin0"+postfix);

  l.rooContainer->AddRealVar("CMS_hll_plaw0"+postfix,0.01,-10,10);

  l.rooContainer->AddRealVar("CMS_hll_voigtexp_0"+postfix, 91.186, 80.0, 100.0);
  l.rooContainer->AddRealVar("CMS_hll_voigtexp_1"+postfix, 2.125, 0.0, 5.0);
  l.rooContainer->AddRealVar("CMS_hll_voigtexp_2"+postfix, 5, 0, 10.0);
  l.rooContainer->AddRealVar("CMS_hll_voigtexp_3"+postfix,-1.,-10.0,0.0);
  l.rooContainer->AddRealVar("CMS_hll_voigtexp_4"+postfix, .5, .0 , 1.0);

  //l.rooContainer->AddRealVar("CMS_hll_voigt_0"+postfix, 91.186, 88.0, 94.0);
  //l.rooContainer->AddRealVar("CMS_hll_voigt_1"+postfix, 2.495, 0.0, 0.0);
  //l.rooContainer->AddRealVar("CMS_hll_voigt_2"+postfix, 5., 0., 10.0);

  // prefix for models parameters
  std::map<int,std::string> parnames;
  parnames[1] = "modlin";
  parnames[2] = "modquad";
  parnames[3] = "modcubic";
  parnames[4] = "modquartic";
  parnames[5] = "modpol5_";
  //parnames[5] = "voigt_";
  parnames[6] = "modpol6_";
  parnames[7] = "voigtexp_";
  parnames[-1] = "plaw";


  // map order to categories flags + parameters names
  std::map<int, std::pair<std::vector<int>, std::vector<std::string> > > catmodels;
  // fill the map
  for(int icat=0; icat<nCategories_; ++icat) {
    // get the poly order for this category
    int catmodel = bkgPolOrderByCat[icat];
    if (catmodel == 7)
      catmodel = 5;
    int parnameIndex = bkgPolOrderByCat[icat];
    std::vector<int> & catflags = catmodels[catmodel].first;
    std::vector<std::string> & catpars = catmodels[catmodel].second;
    // if this is the first time we find this order, build the parameters
    if( catflags.empty() ) {
      assert( catpars.empty() );
      // by default no category has the new model
      catflags.resize(nCategories_, 0);
      std::string & parname = parnames[parnameIndex];
      if( catmodel > 0 ) {
        for(int iorder = 0; iorder<catmodel; ++iorder) {
          std::cout << "HtoLL " <<   Form( "CMS_hll_%s%d%s", parname.c_str(), iorder, +postfix.c_str() ) << std::endl;
          catpars.push_back( Form( "CMS_hll_%s%d%s", parname.c_str(), iorder, +postfix.c_str() ) );
        }
      } else {
        if( catmodel != -1 ) {
          std::cout << "The only supported negative bkg poly order is -1, ie 1-parmeter power law" << std::endl;
          assert( 0 );
        }
        catpars.push_back( Form( "CMS_hll_%s%d%s", parname.c_str(), 0, +postfix.c_str() ) );
      }
    } else if ( catmodel != -1 ) {
      
      assert( catflags.size() == nCategories_ && catpars.size() == catmodel );
    }
    // chose category order
    catflags[icat] = 1;
  }
  
  // now loop over the models and allocate the pdfs
  for(std::map<int, std::pair<std::vector<int>, std::vector<std::string> > >::iterator modit = catmodels.begin();
      modit!=catmodels.end(); ++modit ) {
    std::vector<int> & catflags = modit->second.first;
    std::vector<std::string> & catpars = modit->second.second;
    
    if( modit->first > 0 ) {
      l.rooContainer->AddSpecificCategoryPdf(&catflags[0],"data_voigtexp_model"+postfix,
                                             "0","CMS_hll_mass", catpars, 7);
      // >= 71 means RooBernstein of order >= 1
    } else {
      l.rooContainer->AddSpecificCategoryPdf(&catflags[0],"data_pol_model"+postfix,
                                             "0","CMS_hll_mass",catpars,6);
      // 6 is power law
    }
  }
}

int SSLeptonAnalysis::categories(TLorentzVector* p1, TLorentzVector* p2, bool mixed) {
  
  if (fabs(p1->Eta())<1.479 and fabs(p2->Eta()) < 1.479)
    return 0;
  else {
    if (!mixed)
      return 1;
    else
      if (fabs(p1->Eta()) > 1.479)
	return 2;
      else
	return 1;
  }
}

// ----------------------------------------------------------------------------------------------------
bool SSLeptonAnalysis::Analysis(LoopAll& l, Int_t jentry) {

  //apply pileup reweighting
  float weight = 1.;
  float pu_weight = 1.;
  int cur_type = l.itype[l.current];

  if (cur_type != 0) {
    unsigned int n_pu = l.pu_n;
    pu_weight = getPuWeight(l.pu_n, cur_type, &(l.sampleContainer[l.current_sample_index]), jentry == 1);
    weight = pu_weight * l.sampleContainer[l.current_sample_index].weight();
  }

  std::vector<unsigned int> goodEl, goodMu;
  
  for (int i=0; i<l.el_std_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)(l.el_std_p4_corr->At(i));
    if (l.el_std_mva_trig[i] < eleIDCut and p4->Et() > 20.)
      goodEl.push_back(i);
  }

  for (int i=0; i<l.mu_glo_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)(l.mu_glo_p4_corr->At(i));
    if ((muIDTight and l.mu_glo_id_tight[i] == 1 and p4->Pt() > 20.) or
	(!muIDTight and l.mu_glo_id_loose[i] == 1 and p4->Pt() > 20.))
      goodMu.push_back(i);
  }

  Float_t mass[100];
  Int_t type[100], cat[100];

  Int_t pairs = 0;
  if (goodMu.size() != 0) {
    for (unsigned int i=0; i<goodMu.size()-1; i++) {
      TLorentzVector* p1 = (TLorentzVector*)(l.mu_glo_p4_corr->At(goodMu[i]));
      for (unsigned int j=i+1; j<goodMu.size(); j++) {
	if (l.mu_glo_charge[i] == l.mu_glo_charge[j]) {
	  TLorentzVector* p2 = (TLorentzVector*)(l.mu_glo_p4_corr->At(goodMu[j]));
	  mass[pairs] = (*p1+*p2).M();
	  type[pairs] = 0;	
	  cat[pairs] = categories(p1, p2);
	  pairs++;
	}
      }
    }
  }

  if (goodEl.size() != 0) {	
  for (unsigned int i=0; i<goodEl.size()-1; i++) {
    TLorentzVector* p1 = (TLorentzVector*)(l.el_std_p4_corr->At(goodEl[i]));
    for (unsigned int j=i+1; j<goodEl.size(); j++) {
      if (l.el_std_charge[i] == l.el_std_charge[j]) {
	TLorentzVector* p2 = (TLorentzVector*)(l.el_std_p4_corr->At(goodEl[j]));
	mass[pairs] = (*p1+*p2).M();
	type[pairs] = 1;
	cat[pairs] = categories(p1, p2);
	pairs++;
      }
    }
  }
  }

  if (goodMu.size() != 0 and goodEl.size() != 0) {
  for (unsigned int i=0; i<goodEl.size(); i++) {
    TLorentzVector* p1 = (TLorentzVector*)(l.el_std_p4_corr->At(goodEl[i]));
    for (unsigned int j=0; j<goodMu.size(); j++) {
      if (l.el_std_charge[i] == l.mu_glo_charge[j]) {
	TLorentzVector* p2 = (TLorentzVector*)(l.mu_glo_p4_corr->At(goodMu[j]));
	mass[pairs] = (*p1+*p2).M();
	type[pairs] = 2;
	cat[pairs] = categories(p1, p2, true);
	pairs++;
      }
    }
  }
  }

  if (pairs > 0) {
    Tree(l, pairs, type, mass, cat, weight, pu_weight);
    return true;
  }

  return false;
}

void SSLeptonAnalysis::FillRooContainer(LoopAll& l, int cur_type, float mass, int category, float weight) {

  if (cur_type == 0 ) {
    l.rooContainer->InputDataPoint("data_mass",category,mass);
  } else if (cur_type > 0 ) {
    if( doMcOptimization ) {
      l.rooContainer->InputDataPoint("data_mass",category,mass,weight);
    } else if ( cur_type != 3 && cur_type != 4 ) {
      l.rooContainer->InputDataPoint("bkg_mass",category,mass,weight);
    }
  } else if (cur_type < 0) {
    l.rooContainer->InputDataPoint("sig_"+GetSignalLabel(cur_type),category,mass,weight);
  }
}

void SSLeptonAnalysis::FillSignalLabelMap(LoopAll & l) {

  //signalLabels[-113]="ggh_mass_m125";
  //signalLabels[-213]="vbf_mass_m125";
  //signalLabels[-313]="wzh_mass_m125";
  //signalLabels[-413]="tth_mass_m125";
  //signalLabels[-513]="rsg_mass_m125";
}

std::string SSLeptonAnalysis::GetSignalLabel(int id) {
  // For the lazy man, can return a memeber of the map rather than doing it yourself
  std::map<int,std::string>::iterator it = signalLabels.find(id);
  
  if (it!=signalLabels.end()){
    return it->second;
    
  } else {
    
    std::cerr << "No Signal Type defined in map with id - " << id << std::endl;
    return "NULL";
  }
}


// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::FillReductionVariables(LoopAll& l, int jentry) {

  l.mu_glo_p4_corr->Clear();
  
  for (int i=0; i<l.mu_glo_n; i++) {
    TLorentzVector* p = (TLorentzVector*)(l.mu_glo_p4->At(i)->Clone());
    if (l.itype[l.current] == 0)
      muCorrector_->applyPtCorrection(*p, l.mu_glo_charge[i]);
    else
      muCorrector_->applyPtSmearing(*p, l.mu_glo_charge[i]);
    
    new((*(l.mu_glo_p4_corr))[i]) TLorentzVector();
    ((TLorentzVector*)l.mu_glo_p4_corr->At(i))->SetXYZM(p->X(), p->Y(), p->Z(), p->M());

    if (l.MuonLooseID2012(i) && l.MuonIsolation2012(i, p->Pt(), false))
      l.mu_glo_id_loose[i] = 1;
    else
      l.mu_glo_id_loose[i] = 0;

    if (l.MuonTightID2012(i) && l.MuonIsolation2012(i, p->Pt(), true))
      l.mu_glo_id_tight[i] = 1;
    else
      l.mu_glo_id_tight[i] = 0;
  }

//  if (l.itype[l.current] < 0) {
//    for (int i=0; i<l.gp_n; i++) {
//      if (l.gp_pdgid[i] == 25 && l.gp_status[i] == 3) {
//        TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(i);
//        *(l.higgs) = *p4;
//        break;
//      }
//    }
//    for (int i=0; i<l.mu_glo_n; i++) {
//      Int_t mc1=-1;
//      Int_t mc2=-1;
//      Int_t pho=-1;
//      
//      l.FindMCLeptons(i, mc1, mc2, pho, 13);
//      // et, eta, phi
//      
//      if (mc1 != -1) {
//	TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(mc1);
//	l.mc_et[i] = p4->Et();
//	l.mc_phi[i] = p4->Phi();
//	l.mc_eta[i] = p4->Eta();
//      } else {
//	l.mc_et[i] = -999.;
//	l.mc_phi[i] = -999.;
//	l.mc_eta[i] = -999.;
//      }
//      
//      if (pho != -1) {
//	TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(pho);
//	l.fsr_et[i] = p4->Et();
//	l.fsr_phi[i] = p4->Phi();
//	l.fsr_eta[i] = p4->Eta();
//      } else {
//	l.fsr_et[i] = -999.;
//	l.fsr_phi[i] = -999.;
//	l.fsr_eta[i] = -999.;
//      }
//    }
//    for (int i=0; i<l.el_std_n; i++) {
//      Int_t mc1=-1;
//      Int_t mc2=-1;
//      Int_t pho=-1;
//      
//      l.FindMCLeptons(i, mc1, mc2, pho, 11);
//      // et, eta, phi
//      
//      if (mc1 != -1) {
//	TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(mc1);
//	l.mc_et[i] = p4->Et();
//	l.mc_phi[i] = p4->Phi();
//	l.mc_eta[i] = p4->Eta();
//      } else {
//	l.mc_et[i] = -999.;
//	l.mc_phi[i] = -999.;
//	l.mc_eta[i] = -999.;
//      }
//      
//      if (pho != -1) {
//	TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(pho);
//	l.fsr_et[i] = p4->Et();
//	l.fsr_phi[i] = p4->Phi();
//	l.fsr_eta[i] = p4->Eta();
//      } else {
//	l.fsr_et[i] = -999.;
//	l.fsr_phi[i] = -999.;
//	l.fsr_eta[i] = -999.;
//      }
//    }
//  }    
}

// ----------------------------------------------------------------------------------------------------
bool SSLeptonAnalysis::SelectEventsReduction(LoopAll& l, int jentry) {
  
  if (l.itype[l.current] == 0 && !checkEventHLT(l, hltSelection))
    if (!checkEventHLT(l, hltSelection))
      return false;
  
  // Two muons/electrons with pT > 20
  int goodMu = 0;
  for (int i=0; i<l.mu_glo_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)l.mu_glo_p4->At(i);
    if (p4->Pt() > 20.)
      goodMu++;
  }
  
  int goodEl = 0;
  for(int i =0; i<l.el_std_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)l.el_std_p4->At(i);
    if (p4->Et() > 20.)
      goodEl++;
  }
  
  return ((goodMu+goodEl) > 1);
}

// ----------------------------------------------------------------------------------------------------
bool SSLeptonAnalysis::SkimEvents(LoopAll& l, int jentry) {
  return true;
}

// ----------------------------------------------------------------------------------------------------
bool SSLeptonAnalysis::SelectEvents(LoopAll& l, int jentry) {
  return true;
}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) {

  //l.higgs = new TLorentzVector(0,0,0,0);
  l.mu_glo_p4_corr = new TClonesArray("TLorentzVector", 100);
  
  l.Branch_mu_glo_p4_corr(outputTree);
  l.Branch_mu_glo_id_tight(outputTree);
  l.Branch_mu_glo_id_loose(outputTree);
  l.Branch_mc_et(outputTree);
  l.Branch_mc_eta(outputTree);
  l.Branch_mc_phi(outputTree);
  //l.Branch_fsr_et(outputTree);
  //l.Branch_fsr_eta(outputTree);
  //l.Branch_fsr_phi(outputTree);
  //l.Branch_higgs(outputTree);
}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::ResetAnalysis()
{}

void SSLeptonAnalysis::Tree(LoopAll& l, Int_t pairs, Int_t* type, Float_t* mass, Int_t* cat,
			    Float_t weight, Float_t pu_weight) {

  l.FillTree("run", l.run);
  l.FillTree("lumis", l.lumis);
  l.FillTree("event", (double)l.event);
  l.FillTree("itype", (float)l.itype[l.current]);
  l.FillTree("nvtx", (float)l.vtx_std_n);
  l.FillTree("rho", (float)l.rho_algo1);
  
  l.FillTree("pairs", pairs);
  l.FillTree("type",  type, pairs);
  l.FillTree("mass",  mass, pairs);
  l.FillTree("cat",    cat, pairs);
  l.FillTree("weight", (float)weight);
  l.FillTree("pu_weight", (float)pu_weight);
    
}

bool SSLeptonAnalysis::checkEventHLT(LoopAll& l, std::vector<std::string> paths) {
  
  bool result = false;

  std::vector<unsigned short> hltNumbers;
  
  for (unsigned int i=0; i<paths.size(); i++) {
    //std::cout << i << " " << paths[i] << std::endl;
    //std::cout << "_______________________" << std::endl;
    TRegexp e(TString(paths[i].c_str()), true);
    for (unsigned int j=0; j<l.hlt_path_names_HLT->size(); j++) {
      TString str1((*l.hlt_path_names_HLT)[j].c_str());
      //std::cout << (*l.hlt_path_names_HLT)[j] << std::endl;
      if (str1.Contains(e)) {
        //std::cout << (*l.hlt_path_names_HLT)[j] << std::endl;
        hltNumbers.push_back(j);
      }
    }
  }
  
  //system ("sleep 100000000");
  for (int j=0; j< hltNumbers.size(); j++) {
    for (int i=0; i<(*l.hlt_bit).size(); i++) {
      if (hltNumbers[j] == (*l.hlt_bit)[i]) {
        result = true;
        break;
      }
    }
  }

  return result;
}
