/***************************************************************************
 * Project: CMS detector at the LHC, CERN
 * Package: RooFit
 *    File: $Id$
 * Authors:
 *   Giordano Cerizza
 *   
 * History:
 *   08-15-2008 TS started
 *

 *****************************************************************************/
#include <map>
#include <string>
#include "TLatex.h"
#include "RooFit.h"
#include "Riostream.h"
#include "RooNLLVar.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "Riostream.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooFormula.h"
#include "TNtupleD.h"
#include "RooStringVar.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TMath.h"
#include "RooFit.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooRandom.h"
#include "Riostream.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "RooAddPdf.h"
#include "RooAddModel.h"
#include "RooAbsReal.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooVoigtian.h"
#include "RooBreitWigner.h"
#include "RooBifurGauss.h"
#include "RooArgusBG.h"
#include "RooCBShape.h"
#include "RooGaussModel.h"
#include "RooNovosibirsk.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooJpsiJpsiFit/RooJpsiJpsiFit.hh"
#include "RooGaussModel.h"
#include "RooGExpModel.h"
#include "RooLandau.h"
#include "RooDecay.h"
#include "RooMCStudy.h"
#include "RooGlobalFunc.h"
#include "RooChebychev.h"
#include "TApplication.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooParametricStepFunction.h"
#include "RooGenericPdf.h"
using namespace RooFit;

RooJpsiJpsiFit::RooJpsiJpsiFit()  { 

  cout << " //////////////////////////////////////////////////// " << endl;
  cout << " //                                                // " << endl;
  cout << " //              J/Psi J/Psi Fitter                // " << endl;
  cout << " //                                                // " << endl;
  cout << " // Author:    Giordano Cerizza                    // " << endl; 
  cout << " //                                                // " << endl;
  cout << " // Institute: University of Tennessee             // " << endl;
  cout << " //            401 Nielsen Physics Building        // " << endl;
  cout << " //            1408 Circle Drive                   // " << endl;
  cout << " //            Knoxville, TN 37996-1200            // " << endl;
  cout << " //                                                // " << endl;
  cout << " // email:     gcerizza@utk.edu                    // " << endl;   
  cout << " //                                                // " << endl;
  cout << " //////////////////////////////////////////////////// " << endl;

  initDataVars();

  gROOT->SetStyle("Plain");
}

RooJpsiJpsiFit::~RooJpsiJpsiFit()  { }

void 
RooJpsiJpsiFit::initDataVars() {

  cout << "RooJpsiJpsiFit: Initialization Data Variables" << endl ;

  FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",0.,999.);
  Psi1To2Significance = new RooRealVar("Psi1To2Significance","Psi1To2Significance",0.,8);
  Psi1_Mass = new RooRealVar("Psi1_Mass","Psi1_Mass",2.85,3.35);
  Psi2_Mass = new RooRealVar("Psi2_Mass","Psi2_Mass",2.85,3.35);
  Psi1_CTxy = new RooRealVar("Psi1_CTxy","Psi1_CTxy",-0.03,0.1);
  Psi2_CTxy = new RooRealVar("Psi2_CTxy","Psi2_CTxy",-0.05,0.1);

  // bins
  //  FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",6.,8.);
  //  FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",8.,13.);
  //  FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",13.,22.);
  //  FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",22.,35.);
  //  FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",35.,999.);

  // bins
  //  FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",6.,8.82);
  //  FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",8.82,18.42);
  //  FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",18.42,80.);


  //  dataVars.add(RooArgSet(*Psi1_Mass,*Psi2_Mass,*Psi1_CTxy,*Psi2_CTxy,*Psi1To2Significance));
  dataVars.add(RooArgSet(*FourMu_Mass,*Psi1_Mass,*Psi2_Mass,*Psi1_CTxy,*Psi2_CTxy,*Psi1To2Significance));
  _initDataVars = kTRUE ;
  
}

RooDataSet*
RooJpsiJpsiFit::readData(std::map<std::string, double> & dataFile) {
  
  std::map<std::string, double>::iterator iter = dataFile.begin();
 
  RooRealVar weight("weight","weight",0.1,-10.,10.);
  RooArgSet vars(dataVars);
  vars.add(weight);
  
  RooDataSet *data_nw(0);
  
  while (iter!=dataFile.end()) {
    
    cout << "RooJpsiJpsiFit: Reading data " << iter->first << endl;

    RooDataSet *subdata;
    
    // check for root files as input
    
    if (iter->first.rfind(".root")==iter->first.size()-5) {

      TFile *f = new TFile(iter->first.c_str());
      TTree* tree = (TTree*)f->Get("PATEventTree");
      subdata = new RooDataSet((iter->first+"_subdata").c_str(),"subdata",tree,dataVars);

    } else {

      // txt mode
      subdata = RooDataSet::read(iter->first.c_str(),dataVars);
      subdata->SetNameTitle((iter->first+"_subdata").c_str(),"subdata");
      
    }
    
    weight.setVal(iter->second);
    subdata->addColumn(weight);
      
    if (data_nw==0) {
      data_nw = new RooDataSet("data_nw","data_nw",vars,Import(*subdata));
    }
    else {
      data_nw->append(*subdata);
    }

    delete subdata;
    ++iter;
  }

  if (!data_nw) cout << "Error reading data file " << endl;
  RooDataSet*  data = new RooDataSet("data","data",vars,Import(*data_nw),WeightVar("weight"));
  delete data_nw;
  cout << "RooJpsiJpsiFit: sumEntries = " << data->sumEntries() << endl;
  cout << "RooJpsiJpsiFit: Number of Entries --> " << data->numEntries() << endl;
  
  return data;
}

void
RooJpsiJpsiFit::PDFmaker(std::map<std::string, double>& filename, TString pdf) {


  //////////////////////////////
  // read the dataset
  //////////////////////////////
  RooDataSet* data;
  data = readData(filename);

  cout << "RooJpsiJpsiFit: Number of Entries --> " << data->numEntries() << endl;

  /////////////////////
  // PDF parameters
  /////////////////////

  frac_1 = new RooRealVar("frac_1","",0.65,0.,1.);
  frac_2 = new RooRealVar("frac_2","",0.65,0.,1.);

  // Jpsi masses (signal)

  jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.1,2.8,3.4);
  jpsi1_mass_2 = new RooRealVar("jpsi1_mass_2","",3.1,2.8,3.4);
  jpsi1_mass_3 = new RooRealVar("jpsi1_mass_3","",3.1,2.8,3.4);
  jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",0.05,0.,0.9);
  jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",0.03,0.,0.9);
  jpsi1_width_3 = new RooRealVar("jpsi1_width_3","",0.03,0.,0.9);
  jpsi1_width_a = new RooFormulaVar("jpsi1_width_a","","@0*@1",RooArgList(*jpsi1_width_1,*jpsi1_width_2));
  jpsi1_width_b = new RooFormulaVar("jpsi1_width_b","","@0*@1",RooArgList(*jpsi1_width_1,*jpsi1_width_3));
  RooRealVar* jpsi1_CBalpha = new RooRealVar("jpsi1_CBalpha","",0.5,0.1,5.);
  RooRealVar* jpsi1_CBenne = new RooRealVar("jpsi1_CBenne","",10.,1.,60.);

  jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.1,2.8,3.4);
  jpsi2_mass_2 = new RooRealVar("jpsi2_mass_2","",3.1,2.8,3.4);
  jpsi2_mass_3 = new RooRealVar("jpsi2_mass_3","",3.1,2.8,3.4);
  jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",0.05,0.,0.9);
  jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",0.03,0.,0.9);
  jpsi2_width_3 = new RooRealVar("jpsi2_width_3","",0.03,0.,0.9);
  jpsi2_width_a = new RooFormulaVar("jpsi2_width_a","","@0*@1",RooArgList(*jpsi2_width_1,*jpsi2_width_2));
  jpsi2_width_b = new RooFormulaVar("jpsi2_width_b","","@0*@1",RooArgList(*jpsi2_width_1,*jpsi2_width_3));
  RooRealVar* jpsi2_CBalpha = new RooRealVar("jpsi2_CBalpha","",0.5,0.1,5.);
  RooRealVar* jpsi2_CBenne = new RooRealVar("jpsi2_CBenne","",10.,1.,60.);

  // Jpsi mass (bkg)

  bkg_p0_jpsi1 = new RooRealVar("bkg_p0_jpsi1","",0.,-5.,5.);
  bkg_p1_jpsi1 = new RooRealVar("bkg_p1_jpsi1","",0.,-5.,5.);
  bkg_p2_jpsi1 = new RooRealVar("bkg_p2_jpsi1","",0.,-5.,5.);


  bkg_p0_jpsi2 = new RooRealVar("bkg_p0_jpsi2","",0.,-5.,5.);
  bkg_p1_jpsi2 = new RooRealVar("bkg_p1_jpsi2","",0.,-5.,5.);

  //  ct (Bbkg)

  Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","",0.04,0,1);
  Bbkg_lambda_CT2 = new RooRealVar("Bbkg_lambda_CT2","",0.04,0,1);
  etab_lambda = new RooRealVar("etab_lambda","",9,0,100);

  // Eta_b ct (bkg)

  bkg_lambda1 = new RooRealVar("bkg_lambda1","",0.04,0.,0.1);
  bkg_lambda2 = new RooRealVar("bkg_lambda2","",0.04,0.,0.1);

  // Resolution function

  R_mean_core = new RooRealVar("R_mean_core","",0.,-0.04,0.04);
  R_sigma_core = new RooRealVar("R_sigma_core","",0.4,0.,5);
  R_mean_tail = new RooRealVar("R_mean_tail","",0.,-0.04,0.04);
  R_sigma_tail = new RooRealVar("R_sigma_tail","",0.8,0.,5);
  R_sigma_tot = new RooFormulaVar("R_sigma_tot","","@0*@1",RooArgList(*R_sigma_tail,*R_sigma_core));

  // For significance

  bkg_co0 =  new RooRealVar("bkg_co0","", 0.5,0.,1.);
  bkg_co1 =  new RooRealVar("bkg_co1","", 0.2,0.,0.95);
  bkg_meanlandau = new RooRealVar("bkg_meanlandau","",1.5, 0., 5.);
  bkg_sigmalandau = new RooRealVar("bkg_sigmalandau","",0.8, 0., 2.);
  bkg_flau = new RooRealVar("bkg_flau","",0.6,0.,1.);

  RooRealVar* Sig_mean_core = new RooRealVar("Sig_mean_core","",0.5,0.,2.);
  RooRealVar* Sig_sigma_core = new RooRealVar("Sig_sigma_core","",0.5,0.,2.);
  RooRealVar* Sig_lambda = new RooRealVar("Sig_lambda","",1.,0.,8.);

  RooRealVar* SigB_mean_core = new RooRealVar("SigB_mean_core","",2.,0.,8.);
  RooRealVar* SigB_sigma_core = new RooRealVar("SigB_sigma_core","",1.,0.,4.);
  RooRealVar* SigB_lambda = new RooRealVar("SigB_lambda","",1.,0.,100.);


  bkg_p0_distT = new RooRealVar("bkg_p0_distT","",-3.85548e-01,-5.,5.);
  bkg_p1_distT = new RooRealVar("bkg_p1_distT","",-7.58832e-01,-5.,5.);
  bkg_p2_distT = new RooRealVar("bkg_p2_distT","",6.53086e-01,-5.,5.);
  bkg_p3_distT = new RooRealVar("bkg_p3_distT","",1.90211e-01,-1.,1.); 
  bkg_p4_distT = new RooRealVar("bkg_p4_distT","",2.17327e-01,0.,1.);          

  /////////////////////  
  // PDF definition
  /////////////////////

  if ( pdf == "plotM_4mu"){
    
    RooPlot* frame = FourMu_Mass->frame(Title("#mu^{+}#mu^{-}#mu^{+}#mu^{-} invariant mass"));
    frame = FourMu_Mass->frame(30);
    data->plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-}#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();
  } 

  if ( pdf == "plotM_jpsi1"){
    
    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    data->plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();
  } 

  if ( pdf == "plotM_jpsi2"){
    
    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(20);
    data->plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();
  } 
  
  if ( pdf == "plotCt"){
    
    RooPlot* frame = Psi1To2Significance->frame(Title("ct (cm)"));
    frame = Psi1To2Significance->frame(40);
    data->plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("ct (cm)");
    frame->Draw();
  } 

  if ( pdf == "G_jpsi1"){
    
    RooGaussian jpsimass_1("jpsimass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    
    jpsimass_1.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    data->plotOn(frame);
    jpsimass_1.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 5;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi1_Mass_G.pdf");
    c1->Close();
  } 

  if ( pdf == "G_jpsi2"){
    
    RooGaussian jpsimass_2("jpsimass_2","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_1);
    
    jpsimass_2.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(20);
    data->plotOn(frame);
    jpsimass_2.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 5;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi2_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi2_Mass_G.pdf");
    c1->Close();
  } 

  if ( pdf == "Pol2G1"){
    
    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooGaussian jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_2,*jpsi1_width_a);
    RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));

    RooChebychev bkg_distT_Pol("bkg_distT_Pol","",*Psi1_Mass,RooArgList(*bkg_p0_distT,*bkg_p1_distT));
    RooAddPdf bkg_distT("bkg_distT","",RooArgList(jpsi1mass,bkg_distT_Pol),RooArgList(*frac_2));    

    bkg_distT.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    data->plotOn(frame);
    bkg_distT.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();
  } 

  if ( pdf == "2G_jpsi1"){
    
    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooGaussian jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_a);
    RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));

    jpsi1mass.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    data->plotOn(frame);
    jpsi1mass.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 5;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi1_Mass_2G.pdf");
    c1->Close();
  }

  if ( pdf == "CBGaus_jpsi1"){

    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooCBShape jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_2,*jpsi1_width_a,*jpsi1_CBalpha,*jpsi1_CBenne);
    RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));

    jpsi1mass.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    //data->plotOn(frame,DataError(RooAbsData::SumW2));
    data->plotOn(frame);
    jpsi1mass.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 7;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi1_Mass_CBGaus.pdf");
    c1->Close();
  }

  if ( pdf == "2G_jpsi1_a"){
    
    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooGaussian jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_2);
    RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));

    jpsi1mass.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    data->plotOn(frame);
    jpsi1mass.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();
  } 

  if ( pdf == "2G_jpsi1_b"){

    // From SPS signal parameterization
    //jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.09373e+00);
    //jpsi1_mass_2 = new RooRealVar("jpsi1_mass_2","",3.09553e+00);
    //jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",5.80659e-02);
    //jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",4.56585e-01);
    //frac_1 = new RooRealVar("frac_1","frac_1",1.78448e-01);

    // From DPS signal parameterization
    //jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.09373e+00);
    //jpsi1_mass_2 = new RooRealVar("jpsi1_mass_2","",3.09553e+00);
    //jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",5.80659e-02);
    //jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",4.56585e-01);
    //frac_1 = new RooRealVar("frac_1","frac_1",1.78448e-01);

    // From letting float on data
    //jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.09121e+00);
    //jpsi1_mass_2 = new RooRealVar("jpsi1_mass_2","",3.09359e+00);
    //jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",4.98331e-02);
    //jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",3.70876e-01);
    //frac_1 = new RooRealVar("frac_1","frac_1",6.25067e-01);

    RooFormulaVar* jpsi1_width_calc = new RooFormulaVar("jpsi1_width_calc","","@0*@1",RooArgList(*jpsi1_width_1,*jpsi1_width_2));
    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooGaussian jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_2,*jpsi1_width_calc);
    RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));
    RooChebychev Bkglinear("Bkglinear","",*Psi1_Mass,RooArgList(*bkg_p0_jpsi1,*bkg_p1_jpsi1,*bkg_p2_jpsi1));
    RooAddPdf jpsi1massPol("jpsi1massPol","jpsi1 mass distribution",RooArgList(jpsi1mass,Bkglinear),RooArgList(*frac_2));

    jpsi1massPol.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    data->plotOn(frame);
    jpsi1massPol.plotOn(frame);
    jpsi1massPol.plotOn(frame,RooFit::Components("Bkg,Bkg*"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 9;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi1_Mass_2GPol.pdf");
    c1->Close();

  } 

  if ( pdf == "2G_jpsi2_b"){

    // From SPS signal parameterization
    //jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.09637e+00);
    //jpsi2_mass_2 = new RooRealVar("jpsi2_mass_2","",3.09285e+00);
    //jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",4.95406e-02);
    //jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",5.71976e-01);
    //frac_1 = new RooRealVar("frac_1","frac_1",4.97356e-01);

    // From DPS signal parameterization
    //jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.09637e+00);
    //jpsi2_mass_2 = new RooRealVar("jpsi2_mass_2","",3.09285e+00);
    //jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",4.95406e-02);
    //jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",5.71976e-01);
    //frac_1 = new RooRealVar("frac_1","frac_1",4.97356e-01);

    // From letting float on data
    //jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.08892e+00);
    //jpsi2_mass_2 = new RooRealVar("jpsi2_mass_2","",3.09571e+00);
    //jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",3.33911e-02);
    //jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",8.31105e-02e-01);
    //frac_1 = new RooRealVar("frac_1","frac_1",9.88042e-01);

    RooFormulaVar* jpsi2_width_calc = new RooFormulaVar("jpsi2_width_calc","","@0*@1",RooArgList(*jpsi2_width_1,*jpsi2_width_2));
    RooGaussian jpsi2mass_1("jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_1);
    RooGaussian jpsi2mass_2("jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_2,*jpsi2_width_calc);
    RooAddPdf jpsi2mass("jpsi2mass","jpsi2 mass distribution",RooArgList(jpsi2mass_1,jpsi2mass_2),RooArgList(*frac_1));
    RooChebychev Bkglinear("Bkglinear","",*Psi2_Mass,RooArgList(*bkg_p0_jpsi1,*bkg_p1_jpsi1,*bkg_p2_jpsi1));
    RooAddPdf jpsi2massPol("jpsi2massPol","jpsi2 mass distribution",RooArgList(jpsi2mass,Bkglinear),RooArgList(*frac_2));

    jpsi2massPol.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(20);
    data->plotOn(frame);
    jpsi2massPol.plotOn(frame);
    jpsi2massPol.plotOn(frame,RooFit::Components("Bkg,Bkg*"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 9;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi2_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi2_Mass_2GPol.pdf");
    c1->Close();

  } 

  if ( pdf == "2G_jpsi2_a"){
    
    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi2_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooGaussian jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi2_Mass,*jpsi1_mass_2,*jpsi1_width_2);
    RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));

    jpsi1mass.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(20);
    data->plotOn(frame);
    jpsi1mass.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();
  } 

  if ( pdf == "2G_jpsi2"){
    
    RooGaussian jpsi2mass_1("jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_1);
    RooGaussian jpsi2mass_2("jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_a);
    RooAddPdf jpsi2mass("jpsi2mass","jpsi2 mass distribution",RooArgList(jpsi2mass_1,jpsi2mass_2),RooArgList(*frac_1));
    
    jpsi2mass.fitTo(*data,SumW2Error(kTRUE));
    Psi2_Mass->setBins(22); 
    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(20);
    data->plotOn(frame);
    jpsi2mass.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();
    
    int nFitPar = 5;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi2_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl; 

    c1->SaveAs("pic/JPsi2_Mass_2G.pdf");
    c1->Close();
  } 

  if ( pdf == "CBGaus_jpsi2"){

    RooGaussian jpsi2mass_1("jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_1);
    RooCBShape jpsi2mass_2("jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_2,*jpsi2_width_a,*jpsi2_CBalpha,*jpsi2_CBenne);
    RooAddPdf jpsi2mass("jpsi2mass","jpsi2 mass distribution",RooArgList(jpsi2mass_1,jpsi2mass_2),RooArgList(*frac_1));

    jpsi2mass.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(20);
    data->plotOn(frame);
    jpsi2mass.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 7;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi2_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi2_Mass_CBGaus.pdf");
    c1->Close();
  }

  if ( pdf == "3G_jpsi1"){
    
    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooGaussian jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_2,*jpsi1_width_a);
    RooGaussian jpsi1mass_3("jpsi1mass_3","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_3,*jpsi1_width_b);
    RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2,jpsi1mass_3),RooArgList(*frac_1,*frac_2));
    
    jpsi1mass.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    data->plotOn(frame);
    frame->SetMinimum(0.002);
    jpsi1mass.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();
  } 

  if ( pdf == "3G_jpsi2"){
    
    RooGaussian jpsi2mass_1("jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_1);
    RooGaussian jpsi2mass_2("jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_2,*jpsi2_width_a);
    RooGaussian jpsi2mass_3("jpsi2mass_3","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_3,*jpsi2_width_b);
    RooAddPdf jpsi2mass("jpsi2mass","jpsi2 mass distribution",RooArgList(jpsi2mass_1,jpsi2mass_2,jpsi2mass_3),RooArgList(*frac_1,*frac_2));
    
    jpsi2mass.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(20);
    data->plotOn(frame);
    frame->SetMinimum(0.002);
    jpsi2mass.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();
  } 

  if ( pdf == "2G_Bjpsi1"){

    jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.09374e+00);
    jpsi1_mass_2 = new RooRealVar("jpsi1_mass_2","",3.09653e+00);
    jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",4.97088e-02);
    jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",4.76911e-01);
    jpsi1_width_a = new RooFormulaVar("jpsi1_width_a","","@0*@1",RooArgList(*jpsi1_width_1,*jpsi1_width_2));
    frac_1 = new RooRealVar("frac_1","",2.98315e-01);

    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooGaussian jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_2,*jpsi1_width_a);
    RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));

    //jpsi1mass.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    data->plotOn(frame);
    jpsi1mass.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 5;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi1_Mass_2G.pdf");
    c1->Close();
  }

  if ( pdf == "2G_Bjpsi2"){
    jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.09907e+00);
    jpsi2_mass_2 = new RooRealVar("jpsi2_mass_2","",3.09559e+00);
    jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",7.08343e-02);
    jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",4.97124e-01);
    jpsi2_width_a = new RooFormulaVar("jpsi2_width_a","","@0*@1",RooArgList(*jpsi2_width_1,*jpsi2_width_2));
    frac_2 = new RooRealVar("frac_2","",1.32151e-01);

    RooGaussian jpsi2mass_1("jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_1);
    RooGaussian jpsi2mass_2("jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_2,*jpsi2_width_a);
    RooAddPdf jpsi2mass("jpsi2mass","jpsi2 mass distribution",RooArgList(jpsi2mass_1,jpsi2mass_2),RooArgList(*frac_1));

    //  jpsi2mass.fitTo(*data,SumW2Error(kTRUE));
    Psi2_Mass->setBins(22);
    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));

    frame = Psi2_Mass->frame(20);
    data->plotOn(frame);
    jpsi2mass.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 5;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi2_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi1_Mass_CBGaus.pdf");
    c1->Close();
  }

  if ( pdf == "Pol1_jpsi1") {
    
    RooChebychev Bkglinear("Bkglinear","",*Psi1_Mass,*bkg_p0_jpsi1);
    
    Bkglinear.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(16);
    data->plotOn(frame);
    Bkglinear.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

  } 

  if ( pdf == "Pol1_jpsi2") {
    
    RooChebychev Bkglinear("Bkglinear","",*Psi2_Mass,*bkg_p0_jpsi2);
    
    Bkglinear.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(16);
    data->plotOn(frame);
    Bkglinear.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

  } 

  if ( pdf == "Pol3_distTot") {
    
    RooChebychev Bkglinear("Bkglinear","",*Psi1To2Significance,RooArgList(*bkg_p0_jpsi1,*bkg_p1_jpsi1,*bkg_p2_jpsi1));
    
    Bkglinear.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1To2Significance->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1To2Significance->frame(25);
    data->plotOn(frame);
    Bkglinear.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("Psi1To2DistTot (cm)");
    frame->Draw();

  } 

  if ( pdf == "Pol2_jpsi1") {
    
    RooChebychev Bkglinear("Bkglinear","",*Psi1_Mass,RooArgList(*bkg_p0_jpsi1,*bkg_p1_jpsi1));
    
    Bkglinear.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(16);
    data->plotOn(frame);
    Bkglinear.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

  } 

  if ( pdf == "Pol2_jpsi2") {
    
    RooChebychev Bkglinear("Bkglinear","",*Psi2_Mass,RooArgList(*bkg_p0_jpsi2,*bkg_p1_jpsi2));
    
    Bkglinear.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(16);
    data->plotOn(frame);
    Bkglinear.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

  } 

  if ( pdf == "GExpCTxy1") {

    RooGaussModel resolution_core("resolution_core","",*Psi1_CTxy,*R_mean_core,*R_sigma_core);
    RooDecay modelC("modelC","",*Psi1_CTxy,*Bbkg_lambda_CT1,resolution_core,RooDecay::SingleSided);

    modelC.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi1_CTxy->frame(Title("#eta_{b} ct"));
    frame = Psi1_CTxy->frame(20);
    data->plotOn(frame,DataError(RooAbsData::SumW2));
    modelC.plotOn(frame);
    frame->SetMinimum(0.1);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    c1->SetLogy();
    frame->SetMaximum(500);
    frame->SetTitle("");
    frame->SetXTitle("Psi1_CTxy (cm)");
    frame->Draw();

    int nFitPar = 3;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1_CTxy->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi1_CTxy_GExp.pdf");
    c1->Close();
  }

  if ( pdf == "GExpCTxy2") {

    RooGaussModel resolution_core("resolution_core","",*Psi2_CTxy,*R_mean_core,*R_sigma_core);
    RooDecay modelC("modelC","",*Psi2_CTxy,*Bbkg_lambda_CT2,resolution_core,RooDecay::SingleSided);

    modelC.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi2_CTxy->frame(Title("#eta_{b} ct"));
    frame = Psi2_CTxy->frame(20);
    data->plotOn(frame,DataError(RooAbsData::SumW2));
    //data->plotOn(frame);
    modelC.plotOn(frame);
    frame->SetMinimum(0.1);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    c1->SetLogy();
    frame->SetMaximum(500);
    frame->SetTitle("");
    frame->SetXTitle("Psi2_CTxy (cm)");
    frame->Draw();

    int nFitPar = 3;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi2_CTxy->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi2_CTxy_GExp.pdf");
    c1->Close();
  }

  if ( pdf == "2GCTxy1"){

    RooGaussian ct_1("ct_1","ct distribution",*Psi1_CTxy,*R_mean_core,*R_sigma_core);
    RooGaussian ct_2("ct_2","ct distribution",*Psi1_CTxy,*R_mean_tail,*R_sigma_tot);
    RooAddPdf ct("ct","ct distribution",RooArgList(ct_1,ct_2),RooArgList(*frac_1));

    ct.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi1_CTxy->frame(Title("#mu^{+}#mu^{-}#mu^{+}#mu^{-} ct"));
    frame = Psi1_CTxy->frame(25);
    data->plotOn(frame, DataError(RooAbsData::SumW2));
    ct.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    c1->SetLogy();
    frame->SetMaximum(3000);
    //    frame->SetMaximum(500);
    frame->SetTitle("");
    frame->SetXTitle("Psi1_CTxy (cm)");
    frame->SetMinimum(0.1);
    frame->Draw();

    int nFitPar = 5;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1_CTxy->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi1_CTxy_2G.pdf");
    c1->Close();
  }

  if ( pdf == "GCTxy1"){

    RooGaussian ct("ct","ct distribution",*Psi1_CTxy,*R_mean_core,*R_sigma_core);

    ct.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi1_CTxy->frame(Title("#mu^{+}#mu^{-}#mu^{+}#mu^{-} ct"));
    frame = Psi1_CTxy->frame(25);
    data->plotOn(frame);
    ct.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    frame->SetMaximum(500);
    frame->SetTitle("");
    frame->SetXTitle("Psi1_CTxy (cm)");
    frame->SetMinimum(0.1);
    frame->Draw();
  }

  if ( pdf == "2GCTxy2"){

    RooGaussian ct_1("ct_1","ct distribution",*Psi2_CTxy,*R_mean_core,*R_sigma_core);
    RooGaussian ct_2("ct_2","ct distribution",*Psi2_CTxy,*R_mean_tail,*R_sigma_tot);
    RooAddPdf ct("ct","ct distribution",RooArgList(ct_1,ct_2),RooArgList(*frac_1));

    ct.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi2_CTxy->frame(Title("#mu^{+}#mu^{-}#mu^{+}#mu^{-} ct"));
    frame = Psi2_CTxy->frame(25);
    
    data->plotOn(frame, DataError(RooAbsData::SumW2) );
    ct.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    c1->SetLogy();
    frame->SetMaximum(3000);
    //    frame->SetMaximum(500);
    frame->SetMinimum(0.1);
    frame->SetTitle("");
    frame->SetXTitle("Psi2_CTxy (cm)");
    frame->Draw();

    int nFitPar = 5;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1_CTxy->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;
 
    c1->SaveAs("pic/JPsi2_CTxy_2G.pdf");
    c1->Close();
  }

  if ( pdf == "GCTxy2"){

    RooGaussian ct("ct","ct distribution",*Psi2_CTxy,*R_mean_core,*R_sigma_core);

    ct.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi2_CTxy->frame(Title("#mu^{+}#mu^{-}#mu^{+}#mu^{-} ct"));
    frame = Psi2_CTxy->frame(25);
    data->plotOn(frame);
    ct.plotOn(frame);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    frame->SetMaximum(500);
    frame->SetTitle("");
    frame->SetXTitle("Psi2_CTxy (cm)");
    frame->SetMinimum(0.1);
    frame->Draw();
  }

  if ( pdf == "SB_GPol3_1"){

    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooChebychev Bkglinear("Bkglinear","",*Psi1_Mass,RooArgList(*bkg_p0_jpsi1,*bkg_p1_jpsi1,*bkg_p2_jpsi1));
    RooAddPdf jpsi1massPol("jpsi1massPol","jpsi1 mass distribution",RooArgList(jpsi1mass_1,Bkglinear),RooArgList(*frac_2));

    jpsi1massPol.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi1_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi1_Mass->frame(20);
    data->plotOn(frame,DataError(RooAbsData::SumW2));
    jpsi1massPol.plotOn(frame);
    jpsi1massPol.plotOn(frame,RooFit::Components("Bkg,Bkg*"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 6;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi1_Mass_SB1_2G_GPol3.pdf");
    c1->Close(); 
  }

  if ( pdf == "SB_GPol3_2"){

    RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi2_Mass,*jpsi1_mass_1,*jpsi1_width_1);
    RooChebychev Bkglinear("Bkglinear","",*Psi2_Mass,RooArgList(*bkg_p0_jpsi1,*bkg_p1_jpsi1,*bkg_p2_jpsi1));
    RooAddPdf jpsi1massPol("jpsi1massPol","jpsi1 mass distribution",RooArgList(jpsi1mass_1,Bkglinear),RooArgList(*frac_2));

    jpsi1massPol.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi2_Mass->frame(Title("#mu^{+}#mu^{-} invariant mass"));
    frame = Psi2_Mass->frame(20);
    data->plotOn(frame,DataError(RooAbsData::SumW2));
    jpsi1massPol.plotOn(frame);
    jpsi1massPol.plotOn(frame,RooFit::Components("Bkg,Bkg*"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd();
    c1->SetFillColor(kWhite);
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
    frame->Draw();

    int nFitPar = 6;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi2_Mass->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/JPsi2_Mass_SB2_2G_GPol3.pdf");
    c1->Close();
  }


  if ( pdf == "2GExp") {

    RooGaussModel resolution_core("resolution_core","",*Psi1To2Significance,*R_mean_core,*R_sigma_core);
    RooGaussModel resolution_tail("resolution_tail","",*Psi1To2Significance,*R_mean_tail,*R_sigma_tot);

    RooAddModel modelR("modelR","double exp distribution",RooArgList(resolution_core,resolution_tail),RooArgList(*frac_1));

    RooDecay modelC("modelC","",*Psi1To2Significance,*etab_lambda,modelR,RooDecay::SingleSided);
    
    modelC.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1To2Significance->frame(Title("#eta_{b} ct"));
    frame = Psi1To2Significance->frame(10);
    data->plotOn(frame);
    modelC.plotOn(frame);
    frame->SetMinimum(0.1);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    //    frame->SetMaximum(1000);
    frame->SetTitle("");
    frame->SetLabelOffset(-0.003,"Y");
    frame->SetXTitle("Transverse #mu^{+}#mu^{-}#mu^{+}#mu^{-} ct (cm)");
    frame->Draw();

  } 

  if ( pdf == "Bbkgnew"){
    RooGaussian ct_2("ct_2","ct distribution",*Psi1To2Significance,*R_mean_tail,*R_sigma_tail);
    RooGaussModel resolution_core("resolution_core","",*Psi1To2Significance,*R_mean_core,*R_sigma_core);
    RooDecay modelC("modelC","",*Psi1To2Significance,*etab_lambda,resolution_core,RooDecay::SingleSided);

    RooAddPdf bkg_distT_G("bkg_distT_G","ct distribution",RooArgList(modelC,ct_2),RooArgList(*frac_1)); 

    bkg_distT_G.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1To2Significance->frame(Title("#eta_{b} ct"));
    frame = Psi1To2Significance->frame(16);
    data->plotOn(frame);
    bkg_distT_G.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("Psi1To2Significance");
    frame->Draw();

  }

  if ( pdf == "Significance_GExp") {

    RooGaussModel resolution_core("resolution_core","",*Psi1To2Significance,*Sig_mean_core,*Sig_sigma_core);
    //RooDecay modelC("modelC","",*Psi1To2Significance,*Sig_lambda,resolution_core,RooDecay::SingleSided);
    
    //modelC.fitTo(*data,SumW2Error(kTRUE));
    // constrain pdf to zero at zero
    Int_t nbins(2);
    TArrayD limits(nbins+1);
    limits[0] = -0.5; limits[1] = 0.0; limits[2] = 8.0; 
    //RooArgList* list = new RooArgList("list");
    RooRealVar* binHeight0 = new RooRealVar("binHeight0","bin 0 Value",10.);
    //list->add(binHeight0);
    RooRealVar* binHeight1 = new RooRealVar("binHeight1","bin 1 Value",1.0);
    //list->add(binHeight1);
    RooParametricStepFunction  force_val("force_val","PSF",*Psi1To2Significance,RooArgList(*binHeight0,*binHeight1),limits,nbins);    
    RooDecay modelC("modelC","",*Psi1To2Significance,*Sig_lambda,resolution_core,RooDecay::SingleSided);
    //modelC.fitTo(*data,ExternalConstraints(force_val),SumW2Error(kTRUE));
    
    // alternate way
    RooDecay model_temp("model_temp","",*Psi1To2Significance,*Sig_lambda,resolution_core,RooDecay::SingleSided);
    //RooProdPdf modelC("modelC","model with constraint",RooArgSet(model_temp,force_val)) ;
    //modelC.fitTo(*data,Constrain(force_val),SumW2Error(kTRUE));
    modelC.fitTo(*data,SumW2Error(kTRUE));

    RooPlot* frame = Psi1To2Significance->frame(Title("Psi1To2Significance"));
    frame = Psi1To2Significance->frame(16);
    data->plotOn(frame);
    modelC.plotOn(frame);
    frame->SetMinimum(0.0);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("Psi1To2Significance");
    frame->Draw();

    int nFitPar = 3;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1To2Significance->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/Significance_GExp.pdf");
    c1->Close();

  } 

  if ( pdf == "SignificanceBbkg_GExp") {

    RooGaussModel resolution_core("resolution_core","",*Psi1To2Significance,*SigB_mean_core,*SigB_sigma_core);
    RooDecay modelC("modelC","",*Psi1To2Significance,*SigB_lambda,resolution_core,RooDecay::SingleSided);
    
    modelC.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1To2Significance->frame(Title("Psi1To2Significance"));
    frame = Psi1To2Significance->frame(16);
    data->plotOn(frame);
    modelC.plotOn(frame);
    frame->SetMinimum(0.1);

    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("Psi1To2Significance");
    frame->Draw();

    int nFitPar = 3;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1To2Significance->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/Significance_GExp.pdf");
    c1->Close();

  } 

  if ( pdf == "Significance_LandauPoly"){

    // Landau + Poly
    RooLandau bkg_landau("bkg_landau", "bkg_landau", *Psi1To2Significance, *bkg_meanlandau, *bkg_sigmalandau);
    RooChebychev bkg_polyshape("bkg_polyshape","",*Psi1To2Significance,RooArgList(*bkg_co0,*bkg_co1));  
    RooAddPdf modelC("modelC","", RooArgList(bkg_landau,bkg_polyshape),*bkg_flau);

    modelC.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1To2Significance->frame(Title("Psi1To2Significance"));
    //frame = Psi1To2Significance->frame(35);
    frame = Psi1To2Significance->frame(16);
    data->plotOn(frame, DataError(RooAbsData::SumW2));
    modelC.plotOn(frame);
    //bkg_polyshape->plotOn(frame);
    //bkg_landau.plotOn(frame);
    //    modelD.plotOn(frame,RooFit::LineColor(kRed));
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    //    frame->SetMaximum(3000);
    frame->SetXTitle("Psi1To2Significance");
    frame->Draw();

    int nFitPar = 4;
    double chi2 = frame->chiSquare(nFitPar);
    int nDOF = Psi1To2Significance->getBinning().numBins() - nFitPar;
    double prob = TMath::Prob(chi2,nDOF);
    cout<<"Reduced #chi^{2}="<<chi2<<", #chi^{2} probability="<<prob<<endl;

    c1->SaveAs("pic/Significance_LandauPoly.pdf");
    c1->Close();

  }

  if ( pdf == "Pol3G"){

    //    RooGaussian ct_1("ct_1","ct distribution",*Psi1To2DistTot,*R_mean_core,*R_sigma_core);
    //    RooGaussian ct_2("ct_2","ct distribution",*Psi1To2DistTot,*R_mean_tail,*R_sigma_tot);
    //    RooAddPdf bkg_distT_G("bkg_distT_G","ct distribution",RooArgList(ct_1,ct_2),RooArgList(*frac_1));

    RooGaussian bkg_distT_G("bkg_distT_G","ct distribution",*Psi1To2Significance,*bkg_p3_distT,*bkg_p4_distT);
    RooChebychev bkg_distT_Pol("bkg_distT_Pol","",*Psi1To2Significance,RooArgList(*bkg_p0_distT,*bkg_p1_distT,*bkg_p2_distT));
    RooAddPdf bkg_distT("bkg_distT","",RooArgList(bkg_distT_G,bkg_distT_Pol),RooArgList(*frac_2));

    bkg_distT.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1To2Significance->frame(Title("#eta_{b} ct"));
    frame = Psi1To2Significance->frame(25);
    data->plotOn(frame);
    bkg_distT.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("Psi1To2DistTot (cm)");
    frame->Draw();

  }

  if ( pdf == "Gct"){
    
    RooGaussian ct("ct","ct distribution",*Psi1To2Significance,*R_mean_core,*R_sigma_core);
    
    ct.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1To2Significance->frame(Title("#eta_{b} ct"));
    frame = Psi1To2Significance->frame(35);
    data->plotOn(frame);
    ct.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetTitle("");
    frame->SetXTitle("#eta_{b} ct (cm)");
    frame->Draw();
  } 

  if ( pdf == "2Gct"){
    
    RooGaussian ct_1("ct_1","ct distribution",*Psi1To2Significance,*R_mean_core,*R_sigma_core);
    RooGaussian ct_2("ct_2","ct distribution",*Psi1To2Significance,*R_mean_tail,*R_sigma_tot);
    RooAddPdf ct("ct","ct distribution",RooArgList(ct_1,ct_2),RooArgList(*frac_1));
    
    ct.fitTo(*data,SumW2Error(kTRUE));
    
    RooPlot* frame = Psi1To2Significance->frame(Title("#mu^{+}#mu^{-}#mu^{+}#mu^{-} ct"));
    frame = Psi1To2Significance->frame(25);
    data->plotOn(frame);
    ct.plotOn(frame);
    
    TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
    c1->cd(); 
    c1->SetFillColor(kWhite); 
    frame->SetMaximum(1000);
    frame->SetTitle("");
    frame->SetXTitle("Psi1To2DistTot (cm)");
    frame->Draw();
  } 

}


void
RooJpsiJpsiFit::pureMCtoysGeneration() {

  /////////////////////////
  // Signal PDFs parameters
  /////////////////////////

  jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.09313e+00);
  //jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.09245e+00);
  jpsi1_mass_2 = new RooRealVar("jpsi1_mass_2","",3.09553e+00);
  jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",5.79266e-02);
  jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",4.56886e-01);
  jpsi1_width_a = new RooFormulaVar("jpsi1_width_a","","@0*@1",RooArgList(*jpsi1_width_1,*jpsi1_width_2));
  frac_1 = new RooRealVar("frac_1","",1.80473e-01);

  jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.08933e+00);
  //jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.08904e+00);
  jpsi2_mass_2 = new RooRealVar("jpsi2_mass_2","",3.09285e+00);
  jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",5.06183e-02);
  jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",5.82996e-01);
  jpsi2_width_a = new RooFormulaVar("jpsi2_width_a","","@0*@1",RooArgList(*jpsi2_width_1,*jpsi2_width_2));
  frac_2 = new RooRealVar("frac_2","",4.47328e-01);

  R_mean_coreCT1 = new RooRealVar("R_mean_coreCT1","",2.63860e-05);
  R_mean_tailCT1 = new RooRealVar("R_mean_tailCT1","",1.19874e-03);
  R_sigma_coreCT1 = new RooRealVar("R_sigma_coreCT1","",2.72201e-03);
  R_sigma_tailCT1 = new RooRealVar("R_sigma_tailCT1","",2.72258e+00);
  R_sigma_totCT1 = new RooFormulaVar("R_sigma_totCT1","","@0*@1",RooArgList(*R_sigma_coreCT1,*R_sigma_tailCT1));
  R_fracCT1 = new RooRealVar("R_fracCT1","",8.47797e-01);

  R_mean_core = new RooRealVar("R_mean_core","",4.08278e-01);
  R_sigma_core = new RooRealVar("R_sigma_core","",2.17748e-01);
  etab_lambda = new RooRealVar("etab_lambda","",6.76288e-01);

  RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
  RooGaussian jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_a);
  RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));

  RooGaussian jpsi2mass_1("jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_1);
  RooGaussian jpsi2mass_2("jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_a);
  RooAddPdf jpsi2mass("jpsi2mass","jpsi2 mass distribution",RooArgList(jpsi2mass_1,jpsi2mass_2),RooArgList(*frac_2));

  RooGaussian ct_1a("ct_1a","ct distribution",*Psi1_CTxy,*R_mean_coreCT1,*R_sigma_coreCT1);
  RooGaussian ct_1b("ct_1b","ct distribution",*Psi1_CTxy,*R_mean_tailCT1,*R_sigma_totCT1);
  RooAddPdf ct_1("ct_1","ct distribution",RooArgList(ct_1a,ct_1b),RooArgList(*R_fracCT1));

  RooGaussModel resolution_core("resolution_core","",*Psi1To2Significance,*R_mean_core,*R_sigma_core);
  RooDecay sig_distT("sig_distT","",*Psi1To2Significance,*etab_lambda,resolution_core,RooDecay::SingleSided);

  RooProdPdf sig_model("sig_model","",RooArgList(jpsi1mass,jpsi2mass,ct_1,sig_distT));

  // B bkg PDF 

  Bbkg_jpsi1_mass_1 = new RooRealVar("Bbkg_jpsi1_mass_1","",3.09313e+00);
  Bbkg_jpsi1_mass_2 = new RooRealVar("Bbkg_jpsi1_mass_2","",3.09553e+00);
  Bbkg_jpsi1_width_1 = new RooRealVar("Bbkg_jpsi1_width_1","",5.79266e-02);
  Bbkg_jpsi1_width_2 = new RooRealVar("Bbkg_jpsi1_width_2","",4.56886e-01);
  Bbkg_jpsi1_width_a = new RooFormulaVar("Bbkg_jpsi1_width_a","","@0*@1",RooArgList(*Bbkg_jpsi1_width_1,*Bbkg_jpsi1_width_2));
  Bbkg_frac_1 = new RooRealVar("Bbkg_frac_1","",1.80473e-01);

  Bbkg_jpsi2_mass_1 = new RooRealVar("Bbkg_jpsi2_mass_1","",3.08933e+00);
  Bbkg_jpsi2_mass_2 = new RooRealVar("Bbkg_jpsi2_mass_2","",3.09285e+00);
  Bbkg_jpsi2_width_1 = new RooRealVar("Bbkg_jpsi2_width_1","",5.06183e-02);
  Bbkg_jpsi2_width_2 = new RooRealVar("Bbkg_jpsi2_width_2","",5.82996e-01);
  Bbkg_jpsi2_width_a = new RooFormulaVar("Bbkg_jpsi2_width_a","","@0*@1",RooArgList(*Bbkg_jpsi2_width_1,*Bbkg_jpsi2_width_2));
  Bbkg_frac_2 = new RooRealVar("Bbkg_frac_2","",4.47328e-01);

  Bbkg_mean_CT1 = new RooRealVar("Bbkg_mean_CT1","",6.51060e-04);
  Bbkg_width_CT1 = new RooRealVar("Bbkg_width_CT1","",3.77146e-03);
  Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","Bbkg_lambda_CT1",1.59424e-02); // for eff_cut dataset
  //Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","Bbkg_lambda_CT1",1.60933e-02); // for eff_cut dataset
  //Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","Bbkg_lambda_CT1",1.59728e-02); // for all dataset

  Bbkg_p3_distT = new RooRealVar("Bbkg_p3_distT","",1.18073e+00);
  Bbkg_p4_distT = new RooRealVar("Bbkg_p4_distT","",5.15922e-01);
  Bbkg_lambda1 = new RooRealVar("Bbkg_lambda1","",9.99999e+01);

  RooGaussian Bbkg_jpsi1mass_1("Bbkg_jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*Bbkg_jpsi1_mass_1,*Bbkg_jpsi1_width_1);
  RooGaussian Bbkg_jpsi1mass_2("Bbkg_jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*Bbkg_jpsi1_mass_1,*Bbkg_jpsi1_width_a);
  RooAddPdf Bbkg_jpsi1mass("Bbkg_jpsi1mass","jpsi1 mass distribution",RooArgList(Bbkg_jpsi1mass_1,Bbkg_jpsi1mass_2),RooArgList(*Bbkg_frac_1));

  RooGaussian Bbkg_jpsi2mass_1("Bbkg_jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*Bbkg_jpsi2_mass_1,*Bbkg_jpsi2_width_1);
  RooGaussian Bbkg_jpsi2mass_2("Bbkg_jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*Bbkg_jpsi2_mass_1,*Bbkg_jpsi2_width_a);
  RooAddPdf Bbkg_jpsi2mass("Bbkg_jpsi2mass","jpsi1 mass distribution",RooArgList(Bbkg_jpsi2mass_1,Bbkg_jpsi2mass_2),RooArgList(*Bbkg_frac_2));

  RooGaussModel resolution_R1("resolution_R1","",*Psi1_CTxy,*Bbkg_mean_CT1,*Bbkg_width_CT1);
  RooDecay Bbkg_CT1("Bbkg_CT1","",*Psi1_CTxy,*Bbkg_lambda_CT1,resolution_R1,RooDecay::SingleSided);

  RooGaussModel resolution_core2("resolution_core2","",*Psi1To2Significance,*Bbkg_p3_distT,*Bbkg_p4_distT);
  RooDecay Bbkg_distT("Bbkg_distT","",*Psi1To2Significance,*Bbkg_lambda1,resolution_core2,RooDecay::SingleSided);

  RooProdPdf Bbkg_model("Bbkg_model","",RooArgList(Bbkg_jpsi1mass,Bbkg_jpsi2mass,Bbkg_CT1,Bbkg_distT));

  // bkg PDF 

  // jpsi-sideband case
  bkg_jpsi1_mass_1 = new RooRealVar("bkg_jpsi1_mass_1","",3.09313e+00);
  //bkg_jpsi1_mass_1 = new RooRealVar("bkg_jpsi1_mass_1","",3.09250e+00);
  //RooRealVar*    bkg_jpsi1_mass_2 = new RooRealVar("bkg_jpsi1_mass_2","",3.09313e+00);
  bkg_jpsi1_width_1 = new RooRealVar("bkg_jpsi1_width_1","",5.79266e-02);
  bkg_jpsi1_width_2 = new RooRealVar("bkg_jpsi1_width_2","",4.56886e-01);
  bkg_jpsi1_width_a = new RooFormulaVar("bkg_jpsi1_width_a","","@0*@1",RooArgList(*bkg_jpsi1_width_1,*bkg_jpsi1_width_2));
  bkg_frac_6 = new RooRealVar("bkg_frac_6","",1.80473e-01);

  bkg_p3 = new RooRealVar("bkg_p3","",-2.10268e-01);
  bkg_p4 = new RooRealVar("bkg_p4","",-1.95504e-01);
  bkg_p5 = new RooRealVar("bkg_p5","",5.00755e-02);

  bkg1_R_mean_core = new RooRealVar("bkg1_R_mean_core","",1.04093e-02);
  bkg1_R_mean_tail = new RooRealVar("bkg1_R_mean_tail","",-2.51869e-04);
  bkg1_R_sigma_core = new RooRealVar("bkg1_R_sigma_core","",1.24609e-02);
  bkg1_R_sigma_tail = new RooRealVar("bkg1_R_sigma_tail","",2.63519e-01);
  bkg_frac_1 = new RooRealVar("bkg_frac_1","",7.07283e-01);
  //bkg1_R_mean_core = new RooRealVar("bkg1_R_mean_core","",2.63860e-05);
  //bkg1_R_mean_tail = new RooRealVar("bkg1_R_mean_tail","",1.19874e-03);
  //bkg1_R_sigma_core = new RooRealVar("bkg1_R_sigma_core","",2.72201e-03);
  //bkg1_R_sigma_tail = new RooRealVar("bkg1_R_sigma_tail","",2.72258e+00);
  //bkg_frac_1 = new RooRealVar("bkg_frac_1","",8.47797e-01);
  bkg1_R_sigma_tot = new RooFormulaVar("bkg1_R_sigma_tot","","@0*@1",RooArgList(*bkg1_R_sigma_tail,*bkg1_R_sigma_core));

  bkg_p3_distT = new RooRealVar("bkg_p3_distT","",2.28997e-01);
  bkg_p4_distT = new RooRealVar("bkg_p4_distT","",7.22911e-02);
  bkg_lambda1 = new RooRealVar("bkg_lambda1","",4.94307e+00);
  RooGaussModel resolution_core3("resolution_core3","",*Psi1To2Significance,*bkg_p3_distT,*bkg_p4_distT);
  RooDecay bkg_distT("bkg_distT","",*Psi1To2Significance,*bkg_lambda1,resolution_core3,RooDecay::SingleSided);

  bkg_co0 =  new RooRealVar("bkg_co0","",9.99840e-01);
  bkg_co1 =  new RooRealVar("bkg_co1","",2.00471e-06);
  bkg_flau = new RooRealVar("bkg_flau","",6.48517e-01);
  bkg_meanlandau = new RooRealVar("bkg_meanlandau","",1.00181e+00);
  bkg_sigmalandau = new RooRealVar("bkg_sigmalandau","",4.35740e-01);
  RooLandau bkg_landau("bkg_landau", "bkg_landau", *Psi1To2Significance, *bkg_meanlandau, *bkg_sigmalandau);
  RooChebychev bkg_polyshape("bkg_polyshape","",*Psi1To2Significance,RooArgList(*bkg_co0,*bkg_co1));
  //RooAddPdf bkg_distT("bkg_distT","", RooArgList(bkg_landau,bkg_polyshape),RooArgList(*bkg_flau));

  //RooGaussian bkg_jpsi1mass("bkg_jpsi1mass","jpsi mass distribution",*Psi1_Mass,*bkg_jpsi1_mass_1,*bkg_jpsi1_width_1);
  RooGaussian bkg_jpsi1mass_1("bkg_jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*bkg_jpsi1_mass_1,*bkg_jpsi1_width_1);
  RooGaussian bkg_jpsi1mass_2("bkg_jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*bkg_jpsi1_mass_1,*bkg_jpsi1_width_a);
  RooAddPdf bkg_jpsi1mass("bkg_jpsi1mass","jpsi1 mass distribution",RooArgList(bkg_jpsi1mass_1,bkg_jpsi1mass_2),RooArgList(*frac_1));
  RooChebychev bkg_jpsi1mass_Pol("bkg_jpsi1mass_Pol","",*Psi2_Mass,RooArgList(*bkg_p3,*bkg_p4,*bkg_p5));

  RooGaussian bkg_ctSB1_1a("bkg_ctSB1_1a","ct distribution",*Psi1_CTxy,*bkg1_R_mean_core,*bkg1_R_sigma_core);
  RooGaussian bkg_ctSB1_1b("bkg_ctSB1_1b","ct distribution",*Psi1_CTxy,*bkg1_R_mean_tail,*bkg1_R_sigma_tot);
  RooAddPdf bkg_CT1_SB1("bkg_CT1_SB1","ct distribution",RooArgList(bkg_ctSB1_1a,bkg_ctSB1_1b),RooArgList(*bkg_frac_1));

  //RooProdPdf bkg_mass1("bkg_mass1","",RooArgList(bkg_jpsi1mass,bkg_jpsi1mass_Pol,bkg_CT1_SB1,bkg_distT));
  RooProdPdf bkg_mass1("bkg_mass1","",RooArgList(bkg_jpsi1mass,bkg_jpsi1mass_Pol,bkg_CT1_SB1,bkg_distT));

  // sideband-jpsi case
  bkg_jpsi2_mass_1 = new RooRealVar("bkg_jpsi2_mass_1","",3.08933e+00);
  //bkg_jpsi2_mass_2 = new RooRealVar("bkg_jpsi2_mass_1","",3.08904e+00);
  bkg_jpsi2_width_1 = new RooRealVar("bkg_jpsi2_width_1","",5.06183e-02);
  bkg_jpsi2_width_2 = new RooRealVar("bkg_jpsi2_width_2","",5.82996e-01);
  bkg_jpsi2_width_a = new RooFormulaVar("bkg_jpsi1_width_a","","@0*@1",RooArgList(*bkg_jpsi2_width_1,*bkg_jpsi2_width_2));
  bkg_frac_7 = new RooRealVar("bkg_frac_7","",4.47328e-01);

  bkg_p0 = new RooRealVar("bkg_p0","",-2.93132e-01);
  bkg_p1 = new RooRealVar("bkg_p1","",-3.89092e-01);
  bkg_p2 = new RooRealVar("bkg_p2","",1.94808e-01);
  
  bkg3_R_mean_core = new RooRealVar("bkg3_R_mean_core","",3.60489e-02);
  bkg3_R_mean_tail = new RooRealVar("bkg3_R_mean_tail","",4.32342e-03);
  bkg3_R_sigma_core = new RooRealVar("bkg3_R_sigma_core","",2.89854e-02);
  bkg3_R_sigma_tail = new RooRealVar("bkg3_R_sigma_tail","",3.60637e-01);
  bkg3_R_sigma_tot = new RooFormulaVar("bkg3_R_sigma_tot","","@0*@1",RooArgList(*bkg3_R_sigma_tail,*bkg3_R_sigma_core));
  bkg_frac_3 = new RooRealVar("bkg_frac_3","",6.27677e-02);

  //RooGaussian bkg_jpsi2mass("bkg_jpsi2mass","jpsi mass distribution",*Psi2_Mass,*bkg_jpsi2_mass_1,*bkg_jpsi2_width_1);
  RooGaussian bkg_jpsi2mass_1("bkg_jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*bkg_jpsi2_mass_1,*bkg_jpsi2_width_1);
  RooGaussian bkg_jpsi2mass_2("bkg_jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*bkg_jpsi2_mass_1,*bkg_jpsi2_width_a);
  RooAddPdf bkg_jpsi2mass("bkg_jpsi2mass","jpsi2 mass distribution",RooArgList(bkg_jpsi2mass_1,bkg_jpsi2mass_2),RooArgList(*bkg_frac_7));
  RooChebychev bkg_jpsi2mass_Pol("bkg_jpsi2mass_Pol","",*Psi1_Mass,RooArgList(*bkg_p0,*bkg_p1,*bkg_p2));

  RooGaussian bkg_ctSB2_1a("bkg_ctSB2_1a","ct distribution",*Psi1_CTxy,*bkg3_R_mean_core,*bkg3_R_sigma_core);
  RooGaussian bkg_ctSB2_1b("bkg_ctSB2_1b","ct distribution",*Psi1_CTxy,*bkg3_R_mean_tail,*bkg3_R_sigma_tot);
  RooAddPdf bkg_CT1_SB2("bkg_CT1_SB2","ct distribution",RooArgList(bkg_ctSB2_1a,bkg_ctSB2_1b),RooArgList(*bkg_frac_3));

  bkg_p5_distT = new RooRealVar("bkg_p5_distT","",2.49205e-01);
  bkg_p6_distT = new RooRealVar("bkg_p6_distT","",1.04013e-01);
  bkg_lambda2 = new RooRealVar("bkg_lambda2","",6.54912e+00);
  RooGaussModel resolution_core4("resolution_core4","",*Psi1To2Significance,*bkg_p5_distT,*bkg_p6_distT);
  //RooDecay bkg_distT2("bkg_distT2","",*Psi1To2Significance,*bkg_lambda2,resolution_core4,RooDecay::SingleSided);

  bkg_co02 =  new RooRealVar("bkg_co02","", 6.58858e-01);
  bkg_co12 =  new RooRealVar("bkg_co12","", 2.48596e-05);
  bkg_flau2 = new RooRealVar("bkg_flau2","",5.33780e-01);
  bkg_meanlandau2 = new RooRealVar("bkg_meanlandau2","",1.09999e+00);
  bkg_sigmalandau2 = new RooRealVar("bkg_sigmalandau2","",4.64046e-01);
  RooChebychev bkg_polyshape2("bkg_polyshape2","",*Psi1To2Significance,RooArgList(*bkg_co02,*bkg_co12));
  RooLandau bkg_landau2("bkg_landau2", "bkg_landau2", *Psi1To2Significance, *bkg_meanlandau2, *bkg_sigmalandau2);
  RooAddPdf bkg_distT2("bkg_distT2","", RooArgList(bkg_landau2,bkg_polyshape2),RooArgList(*bkg_flau2));

  RooProdPdf bkg_mass2("bkg_mass2","",RooArgList(bkg_jpsi2mass,bkg_jpsi2mass_Pol,bkg_CT1_SB2,bkg_distT2));

  // fraction for the J/psi-flat flat-J/psi ratio (Andrew: I think this is defined as (# J/psi-flat)/(# J/psi-flat + # flat-J/psi))
  bkg_frac_5 = new RooRealVar("bkg_frac_5","",5.95167e-01); // for eff_cut dataset
  //bkg_frac_5 = new RooRealVar("bkg_frac_5","",4.80678e-01); // for eff_cut dataset
  //bkg_frac_5 = new RooRealVar("bkg_frac_5","",6.47958e-01); // for all dataset
  RooAddPdf bkg_model("bkg_model","",RooArgList(bkg_mass1,bkg_mass2),RooArgList(*bkg_frac_5));

  // bkg flat flat

  bkg7_R_mean_core = new RooRealVar("bkg7_R_mean_core","",3.45913e-03);
  bkg7_R_sigma_core = new RooRealVar("bkg7_R_sigma_core","",1.10565e-02);

  RooGaussian bkg2_CT1("bkg2_CT1","ct distribution",*Psi1_CTxy,*bkg7_R_mean_core,*bkg7_R_sigma_core);

  bkg_p7_distT = new RooRealVar("bkg_p7_distT","",3.89967e-01);
  bkg_p8_distT = new RooRealVar("bkg_p8_distT","",3.31232e-01);
  bkg_lambda3 = new RooRealVar("bkg_lambda3","",4.67109e+00);
  RooGaussModel resolution_core5("resolution_core5","",*Psi1To2Significance,*bkg_p7_distT,*bkg_p8_distT);
  //RooDecay bkg2_distT("bkg2_distT","",*Psi1To2Significance,*bkg_lambda3,resolution_core5,RooDecay::SingleSided);

  bkg2_co0 =  new RooRealVar("bkg2_co0","", 3.86180e-01);
  bkg2_co1 =  new RooRealVar("bkg2_co1","", 9.49975e-01);
  bkg2_flau = new RooRealVar("bkg2_flau","",7.28332e-01);
  bkg2_meanlandau = new RooRealVar("bkg2_meanlandau","",1.56581e+00);
  bkg2_sigmalandau = new RooRealVar("bkg2_sigmalandau","",5.51089e-01);
  RooLandau bkg2_landau("bkg2_landau", "bkg2_landau", *Psi1To2Significance, *bkg2_meanlandau, *bkg2_sigmalandau);
  RooChebychev bkg2_polyshape("bkg2_polyshape","",*Psi1To2Significance,RooArgList(*bkg2_co0,*bkg2_co1));
  RooAddPdf bkg2_distT("bkg2_distT","", RooArgList(bkg2_landau,bkg2_polyshape),RooArgList(*bkg2_flau));
  
  RooProdPdf bkg2_model("bkg2_model","",RooArgList(bkg_jpsi2mass_Pol,bkg_jpsi1mass_Pol,bkg2_CT1,bkg2_distT));
  
  /////////////////////////////////////////////

  RooRealVar nsig("nsig","number of signal events",4.45625e+02,1,1000); // for eff_cut dataset
  RooRealVar nBbkg("nBbkg","number of B background events",1.81904e+02,1,1000); // for eff_cut dataset
  RooRealVar nbkg("nbkg","number of background events",3.21429e+02,1,2000); // for eff_cut dataset
  RooRealVar nbkg2("nbkg2","number of background events",9.40499e+01,1,1000); // for eff_cut dataset
  //RooRealVar nsig("nsig","number of signal events",4.51246e+02,1,1000); // for all dataset
  //RooRealVar nBbkg("nBbkg","number of B background events",1.82217e+02,1,1000); // for all dataset
  //RooRealVar nbkg("nbkg","number of background events",3.40374e+02,1,2000); // for all dataset
  //RooRealVar nbkg2("nbkg2","number of background events",1.22149e+02,1,1000); // for all dataset

  RooAddPdf model("model","model",RooArgList(sig_model,Bbkg_model,bkg_model,bkg2_model),RooArgList(nsig,nBbkg,nbkg,nbkg2));

  RooMCStudy *mgr = new RooMCStudy(model,RooArgSet(*Psi1_Mass,*Psi2_Mass,*Psi1_CTxy,*Psi1To2Significance),Silence(),Extended(),
                                   FitOptions(Save(kTRUE),PrintEvalErrors(0)));
  
  mgr->generateAndFit(30000);

  RooPlot* frame0 = mgr->plotNLL(-12300,-9500,40);

  RooPlot* frame1a = mgr->plotParam(nsig,Bins(40));                                                                                   
  RooPlot* frame1b = mgr->plotError(nsig,Bins(40));                                                                                       
  RooPlot* frame1c = mgr->plotPull(nsig,FrameRange(-3,3),FrameBins(20),FitGauss(kTRUE));
  //RooPlot* frame0 = mgr->plotNLL(FrameRange(-12300,-9500),Bins(40),FitGauss(kTRUE));
  //RooPlot* frame1a = mgr->plotParam(nsig,FrameRange(325,575),Bins(40),FitGauss(kTRUE));
  //RooPlot* frame1b = mgr->plotError(nsig,FrameRange(20.5,26.3),Bins(40),FitGauss(kTRUE));
  RooPlot* frame2a = mgr->plotParam(nBbkg,Bins(40));                                                                                 
  RooPlot* frame2b = mgr->plotError(nBbkg,Bins(40));
  RooPlot* frame2c = mgr->plotPull(nBbkg,FrameRange(-3,3),FrameBins(20),FitGauss(kTRUE));   

  RooPlot* frame3a = mgr->plotParam(nbkg,Bins(40));                                                                                 
  RooPlot* frame3b = mgr->plotError(nbkg,Bins(40));
  RooPlot* frame3c = mgr->plotPull(nbkg,FrameRange(-3,3),FrameBins(20),FitGauss(kTRUE));   

  RooPlot* frame4a = mgr->plotParam(nbkg2,Bins(40));                                                                                 
  RooPlot* frame4b = mgr->plotError(nbkg2,Bins(40));
  RooPlot* frame4c = mgr->plotPull(nbkg2,FrameRange(-3,3),FrameBins(20),FitGauss(kTRUE));   

  TCanvas* c1 = new TCanvas("c1","",0,0,700,500);
  c1->cd(); c1->SetFillColor(kWhite); 
  frame0->SetTitle("");
  frame0->SetStats(0);
  frame0->Draw();                                                                                      
  c1->SaveAs("pic/plotNLL.pdf");
  //c1->Close();

  TCanvas* c1a = new TCanvas("c1a","",10,10,700,500);
  c1a->cd(); c1a->SetFillColor(kWhite); 
  frame1a->SetTitle("");
  frame1a->SetStats(0);
  frame1a->Draw();                                                                                      
  c1a->SaveAs("pic/sigYields.pdf");
  //c1a->Close();

  TCanvas* c1b = new TCanvas("c1b","",20,20,700,500);
  c1b->cd(); c1b->SetFillColor(kWhite); 
  frame1b->SetTitle("");
  frame1b->SetStats(0);
  frame1b->Draw();                                                                                       
  c1b->SaveAs("pic/sigErrors.pdf");
  //c1b->Close();
  
  TCanvas* c1c = new TCanvas("c1c","",30,30,700,500);
  c1c->cd(); c1c->SetFillColor(kWhite); 
  frame1c->SetTitle("");
  frame1c->GetXaxis()->SetTitle("Signal Pull");
  frame1c->SetStats(0);
  frame1c->Draw();
  c1c->SaveAs("pic/sigPulls.pdf");
  //c1c->Close();

  
  /*TCanvas* c2a = new TCanvas("c2a","",40,40,700,500);
  c2a->cd(); c2a->SetFillColor(kWhite); 
  frame2a->SetTitle("");
  frame2a->Draw();                                                                                      

  TCanvas* c2b = new TCanvas("c2b","",50,50,700,500);
  c2b->cd(); c2b->SetFillColor(kWhite); 
  frame2b->SetTitle("");
  frame2b->Draw();                                                                                       
  
  TCanvas* c2c = new TCanvas("c2c","",60,60,700,500);
  c2c->cd(); c2c->SetFillColor(kWhite); 
  frame2c->SetTitle("");
  frame2c->GetXaxis()->SetTitle("B Background Pull");
  frame2c->SetStats(0);
  frame2c->Draw();

  TCanvas* c3a = new TCanvas("c3a","",70,70,700,500);
  c3a->cd(); c3a->SetFillColor(kWhite); 
  frame3a->SetTitle("");
  frame3a->Draw();                                                                                      

  TCanvas* c3b = new TCanvas("c3b","",80,80,700,500);
  c3b->cd(); c3b->SetFillColor(kWhite); 
  frame3b->SetTitle("");
  frame3b->Draw();                                                                                       
  
  TCanvas* c3c = new TCanvas("c3c","",90,90,700,500);
  c3c->cd(); c3c->SetFillColor(kWhite); 
  frame3c->SetTitle("");
  frame3c->GetXaxis()->SetTitle("Background Pull");
  frame3c->SetStats(0);
  frame3c->Draw();

  TCanvas* c4a = new TCanvas("c4a","",70,70,700,500);
  c4a->cd(); c4a->SetFillColor(kWhite); 
  frame4a->SetTitle("");
  frame4a->Draw();                                                                                      

  TCanvas* c4b = new TCanvas("c4b","",80,80,700,500);
  c4b->cd(); c4b->SetFillColor(kWhite); 
  frame4b->SetTitle("");
  frame4b->Draw();                                                                                       
  
  TCanvas* c4c = new TCanvas("c4c","",90,90,700,500);
  c4c->cd(); c4c->SetFillColor(kWhite); 
  frame4c->SetTitle("");
  frame4c->GetXaxis()->SetTitle("Background 2 Pull");
  frame4c->SetStats(0);
  frame4c->Draw();*/


}


void
RooJpsiJpsiFit::fitData(std::map<std::string, double> & dataFile) {

  std::map<std::string, double>::iterator iter = dataFile.begin();

  RooRealVar weight("weight","weight",0.1,-10.,10.);
  RooArgSet vars(dataVars);
  vars.add(weight);

  RooDataSet *data_nw(0);

  while (iter!=dataFile.end()) {

    cout << "RooJpsiJpsiFit: Reading data " << iter->first << endl;

    RooDataSet *subdata;

    // check for root files as input

    if (iter->first.rfind(".root")==iter->first.size()-5) {

      TFile *f = new TFile(iter->first.c_str());
      TTree* tree = (TTree*)f->Get("PATEventTree");
      subdata = new RooDataSet((iter->first+"_subdata").c_str(),"subdata",tree,dataVars);

    } else {

      // txt mode
      subdata = RooDataSet::read(iter->first.c_str(),dataVars);
      subdata->SetNameTitle((iter->first+"_subdata").c_str(),"subdata");

    }

    weight.setVal(iter->second);
    subdata->addColumn(weight);

    if (data_nw==0) {
      data_nw = new RooDataSet("data_nw","data_nw",vars,Import(*subdata));
    }
    else {
      data_nw->append(*subdata);
    }

    delete subdata;
    ++iter;
  }

  if (!data_nw) cout << "Error reading data file " << endl;
  RooDataSet*  data = new RooDataSet("data","data",vars,Import(*data_nw),WeightVar("weight"));
  delete data_nw;
  cout << "RooJpsiJpsiFit: sumEntries = " << data->sumEntries() << endl;
  cout << "RooJpsiJpsiFit: Number of Entries --> " << data->numEntries() << endl;

  // modify PDFs near origin with expontential pdf
  // 1/(1+exp(-1000*(x-0.005)))
  //RooGenericPdf pdf("pdf", "my pdf","1.0/(1.0+exp(-1000.0*(x-0.005)))", RooArgSet(x, y, a, b);
  //RooExponential bkg2("bkg2","Background 2",*Psi1To2Significance,alpha);

  /////////////////////////
  // Signal PDFs parameters
  /////////////////////////
  //jpsi1_mass_1 = new RooRealVar("jpsi1_mass_2","",3.09313e+00);
  //jpsi2_mass_1 = new RooRealVar("jpsi2_mass_2","",3.08933e+00);

  //jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.09300);
  jpsi1_mass_1 = new RooRealVar("jpsi1_mass_2","",3.09, 3.07, 3.11);
  // from SPS
  //jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.09315e+00);
  //jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",4.90541e-02);
  //jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",5.25564e-01);
  //frac_1 = new RooRealVar("frac_1","",2.63474e-01);
  // from DPS
  //jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.09373e+00);
  //jpsi1_mass_2 = new RooRealVar("jpsi1_mass_2","",3.09553e+00);
  //jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",5.79266e-02);
  //jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",4.56886e-01);
  //frac_1 = new RooRealVar("frac_1","",1.80473e-01);
  //jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",0.046952);
  jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",0.05,0.,0.9);
  jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",0.05,0.,0.9);
  frac_1 = new RooRealVar("frac_1","",0.547353);
  jpsi1_width_a = new RooFormulaVar("jpsi1_width_a","","@0*@1",RooArgList(*jpsi1_width_1,*jpsi1_width_2));

  //jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.08938);
  jpsi2_mass_1 = new RooRealVar("jpsi2_mass_2","",3.09, 3.07, 3.11);
  // from SPS
  //jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.08926e+00);
  //jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",4.92804e-02);
  //jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",5.59273e-01);
  //frac_2 = new RooRealVar("frac_2","",3.71463e-01);
  // from DPS
  //jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.09637e+00);
  //jpsi2_mass_2 = new RooRealVar("jpsi2_mass_2","",3.09285e+00);
  //jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",5.06183e-02);
  //jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",5.82996e-01);
  //frac_2 = new RooRealVar("frac_2","",4.47328e-01);
  //jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",0.056935);
  jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",0.05,0.,2.9);
  jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",0.05,0.,0.9);
  frac_2 = new RooRealVar("frac_2","",0.362924);
  jpsi2_width_a = new RooFormulaVar("jpsi2_width_a","","@0*@1",RooArgList(*jpsi2_width_1,*jpsi2_width_2));

  // from SPS
  //R_mean_coreCT1 = new RooRealVar("R_mean_coreCT1","",1.05389e-04);
  //R_mean_tailCT1 = new RooRealVar("R_mean_tailCT1","",4.56787e-04);
  //R_sigma_coreCT1 = new RooRealVar("R_sigma_coreCT1","",3.05011e-03);
  //R_sigma_tailCT1 = new RooRealVar("R_sigma_tailCT1","",3.08677e+00);
  //R_sigma_totCT1 = new RooFormulaVar("R_sigma_totCT1","","@0*@1",RooArgList(*R_sigma_coreCT1,*R_sigma_tailCT1));
  //R_fracCT1 = new RooRealVar("R_fracCT1","",8.93546e-01);
  // from DPS
  R_mean_coreCT1 = new RooRealVar("R_mean_coreCT1","",2.63860e-05);
  R_mean_tailCT1 = new RooRealVar("R_mean_tailCT1","",1.19874e-03);
  R_sigma_coreCT1 = new RooRealVar("R_sigma_coreCT1","",2.72201e-03);
  R_sigma_tailCT1 = new RooRealVar("R_sigma_tailCT1","",2.72258e+00);
  R_sigma_totCT1 = new RooFormulaVar("R_sigma_totCT1","","@0*@1",RooArgList(*R_sigma_coreCT1,*R_sigma_tailCT1));
  R_fracCT1 = new RooRealVar("R_fracCT1","",8.47797e-01);

  // from SPS
  //R_mean_core = new RooRealVar("R_mean_core","",4.20320e-01);
  //R_sigma_core = new RooRealVar("R_sigma_core","",2.11550e-01);
  //etab_lambda = new RooRealVar("etab_lambda","",7.27597e-01);
  // from DPS
  R_mean_core = new RooRealVar("R_mean_core","",4.08278e-01);
  R_sigma_core = new RooRealVar("R_sigma_core","",2.17748e-01);
  etab_lambda = new RooRealVar("etab_lambda","",6.76288e-01);
  //R_sigma_core = new RooRealVar("R_sigma_core","",1.7e-01);

  RooGaussian jpsi1mass_1("jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_1);
  RooGaussian jpsi1mass_2("jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*jpsi1_mass_1,*jpsi1_width_a);
  RooAddPdf jpsi1mass("jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));

  RooGaussian jpsi2mass_1("jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_1);
  RooGaussian jpsi2mass_2("jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*jpsi2_mass_1,*jpsi2_width_a);
  RooAddPdf jpsi2mass("jpsi2mass","jpsi2 mass distribution",RooArgList(jpsi2mass_1,jpsi2mass_2),RooArgList(*frac_2));

  RooGaussian ct_1a("ct_1a","ct distribution",*Psi1_CTxy,*R_mean_coreCT1,*R_sigma_coreCT1);
  RooGaussian ct_1b("ct_1b","ct distribution",*Psi1_CTxy,*R_mean_tailCT1,*R_sigma_totCT1);
  RooAddPdf ct_1("ct_1","ct distribution",RooArgList(ct_1a,ct_1b),RooArgList(*R_fracCT1));

  RooGaussModel resolution_core("resolution_core","",*Psi1To2Significance,*R_mean_core,*R_sigma_core);
  RooDecay sig_distT("sig_distT","",*Psi1To2Significance,*etab_lambda,resolution_core,RooDecay::SingleSided);

  RooProdPdf sig_model("sig_model","",RooArgList(jpsi1mass,jpsi2mass,ct_1,sig_distT));

  // B bkg PDF 

  Bbkg_jpsi1_mass_1 = new RooRealVar("Bbkg_jpsi1_mass_1","",3.09313e+00);
  Bbkg_jpsi1_mass_2 = new RooRealVar("Bbkg_jpsi1_mass_2","",3.09553e+00);
  Bbkg_jpsi1_width_1 = new RooRealVar("Bbkg_jpsi1_width_1","",5.79266e-02);
  Bbkg_jpsi1_width_2 = new RooRealVar("Bbkg_jpsi1_width_2","",4.56886e-01);
  Bbkg_jpsi1_width_a = new RooFormulaVar("Bbkg_jpsi1_width_a","","@0*@1",RooArgList(*Bbkg_jpsi1_width_1,*Bbkg_jpsi1_width_2));
  Bbkg_frac_1 = new RooRealVar("Bbkg_frac_1","",1.80473e-01);

  Bbkg_jpsi2_mass_1 = new RooRealVar("Bbkg_jpsi2_mass_1","",3.08933e+00);
  Bbkg_jpsi2_mass_2 = new RooRealVar("Bbkg_jpsi2_mass_2","",3.09285e+00);
  Bbkg_jpsi2_width_1 = new RooRealVar("Bbkg_jpsi2_width_1","",5.06183e-02);
  Bbkg_jpsi2_width_2 = new RooRealVar("Bbkg_jpsi2_width_2","",5.82996e-01);
  Bbkg_jpsi2_width_a = new RooFormulaVar("Bbkg_jpsi2_width_a","","@0*@1",RooArgList(*Bbkg_jpsi2_width_1,*Bbkg_jpsi2_width_2));
  Bbkg_frac_2 = new RooRealVar("Bbkg_frac_2","",4.47328e-01);

  Bbkg_mean_CT1 = new RooRealVar("Bbkg_mean_CT1","",6.51060e-04);
  Bbkg_width_CT1 = new RooRealVar("Bbkg_width_CT1","",3.77146e-03);
  Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","",1.59424e-02);
  //Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","Bbkg_lambda_CT1",1.59424e-02, 0.010, 0.020);  // for initial fit
  //Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","Bbkg_lambda_CT1",1.59424e-02); // for eff_cut dataset
  //Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","Bbkg_lambda_CT1",1.59728e-02); // for all dataset

  Bbkg_p3_distT = new RooRealVar("Bbkg_p3_distT","",1.18073e+00);
  Bbkg_p4_distT = new RooRealVar("Bbkg_p4_distT","",5.15922e-01);
  Bbkg_lambda1 = new RooRealVar("Bbkg_lambda1","",9.99999e+01);

  RooGaussian Bbkg_jpsi1mass_1("Bbkg_jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*Bbkg_jpsi1_mass_1,*Bbkg_jpsi1_width_1);
  RooGaussian Bbkg_jpsi1mass_2("Bbkg_jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*Bbkg_jpsi1_mass_1,*Bbkg_jpsi1_width_a);
  RooAddPdf Bbkg_jpsi1mass("Bbkg_jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));

  RooGaussian Bbkg_jpsi2mass_1("Bbkg_jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*Bbkg_jpsi2_mass_1,*Bbkg_jpsi2_width_1);
  RooGaussian Bbkg_jpsi2mass_2("Bbkg_jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*Bbkg_jpsi2_mass_1,*Bbkg_jpsi2_width_a);
  RooAddPdf Bbkg_jpsi2mass("Bbkg_jpsi2mass","jpsi1 mass distribution",RooArgList(jpsi2mass_1,jpsi2mass_2),RooArgList(*frac_2));

  RooGaussModel resolution_R1("resolution_R1","",*Psi1_CTxy,*Bbkg_mean_CT1,*Bbkg_width_CT1);
  RooDecay Bbkg_CT1("Bbkg_CT1","",*Psi1_CTxy,*Bbkg_lambda_CT1,resolution_R1,RooDecay::SingleSided);

  RooGaussModel resolution_core2("resolution_core2","",*Psi1To2Significance,*Bbkg_p3_distT,*Bbkg_p4_distT);
  RooDecay Bbkg_distT("Bbkg_distT","",*Psi1To2Significance,*Bbkg_lambda1,resolution_core2,RooDecay::SingleSided);

  RooProdPdf Bbkg_model("Bbkg_model","",RooArgList(Bbkg_jpsi1mass,Bbkg_jpsi2mass,Bbkg_CT1,Bbkg_distT));

  // bkg PDF 

  // jpsi-sideband case
  bkg_jpsi1_mass_1 = new RooRealVar("bkg_jpsi1_mass_1","",3.09313e+00);
  bkg_jpsi1_width_1 = new RooRealVar("bkg_jpsi1_width_1","",5.79266e-02);
  bkg_jpsi1_width_2 = new RooRealVar("bkg_jpsi1_width_2","",4.56886e-01);
  bkg_jpsi1_width_a = new RooFormulaVar("bkg_jpsi1_width_a","","@0*@1",RooArgList(*bkg_jpsi1_width_1,*bkg_jpsi1_width_2));
  RooRealVar*    bkg_frac_6 = new RooRealVar("bkg_frac_6","",1.80473e-01);

  bkg_p3 = new RooRealVar("bkg_p3","",-2.10268e-01);
  bkg_p4 = new RooRealVar("bkg_p4","",-1.95504e-01);
  bkg_p5 = new RooRealVar("bkg_p5","",5.00755e-02);

  bkg1_R_mean_core = new RooRealVar("bkg1_R_mean_core","",1.04093e-02);
  bkg1_R_mean_tail = new RooRealVar("bkg1_R_mean_tail","",-2.51869e-04);
  bkg1_R_sigma_core = new RooRealVar("bkg1_R_sigma_core","",1.24609e-02);
  bkg1_R_sigma_tail = new RooRealVar("bkg1_R_sigma_tail","",2.63519e-01);
  bkg_frac_1 = new RooRealVar("bkg_frac_1","",7.07283e-01);
  bkg1_R_sigma_tot = new RooFormulaVar("bkg1_R_sigma_tot","","@0*@1",RooArgList(*bkg1_R_sigma_tail,*bkg1_R_sigma_core));

  //RooGaussian bkg_jpsi1mass("bkg_jpsi1mass","jpsi mass distribution",*Psi1_Mass,*bkg_jpsi1_mass_1,*bkg_jpsi1_width_1);
  RooGaussian bkg_jpsi1mass_1("bkg_jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*bkg_jpsi1_mass_1,*bkg_jpsi1_width_1);
  RooGaussian bkg_jpsi1mass_2("bkg_jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*bkg_jpsi1_mass_1,*bkg_jpsi1_width_a);
  RooAddPdf bkg_jpsi1mass("bkg_jpsi1mass","jpsi1 mass distribution",RooArgList(jpsi1mass_1,jpsi1mass_2),RooArgList(*frac_1));
  RooChebychev bkg_jpsi1mass_Pol("bkg_jpsi1mass_Pol","",*Psi2_Mass,RooArgList(*bkg_p3,*bkg_p4,*bkg_p5));

  RooGaussian bkg_ctSB1_1a("bkg_ctSB1_1a","ct distribution",*Psi1_CTxy,*bkg1_R_mean_core,*bkg1_R_sigma_core);
  RooGaussian bkg_ctSB1_1b("bkg_ctSB1_1b","ct distribution",*Psi1_CTxy,*bkg1_R_mean_tail,*bkg1_R_sigma_tot);
  //RooGaussian bkg_ctSB1_1a("bkg_ctSB1_1a","ct distribution",*Psi1_CTxy,*bkg1_R_mean_coreCT1,*bkg1_R_sigma_coreCT1);
  //RooGaussian bkg_ctSB1_1b("bkg_ctSB1_1b","ct distribution",*Psi1_CTxy,*bkg1_R_mean_tailCT1,*bkg1_R_sigma_totCT1);
  RooAddPdf bkg_CT1_SB1("bkg_CT1_SB1","ct distribution",RooArgList(bkg_ctSB1_1a,bkg_ctSB1_1b),RooArgList(*bkg_frac_1));

  bkg_p3_distT = new RooRealVar("bkg_p3_distT","",2.28997e-01);
  bkg_p4_distT = new RooRealVar("bkg_p4_distT","",7.22911e-02);
  bkg_lambda1 = new RooRealVar("bkg_lambda1","",4.94307e+00);
  RooGaussModel resolution_core3("resolution_core3","",*Psi1To2Significance,*bkg_p3_distT,*bkg_p4_distT);
  //RooDecay bkg_distT("bkg_distT","",*Psi1To2Significance,*bkg_lambda1,resolution_core3,RooDecay::SingleSided);

  bkg_co0 =  new RooRealVar("bkg_co0","",9.99840e-01);
  bkg_co1 =  new RooRealVar("bkg_co1","",2.00471e-06);
  bkg_flau = new RooRealVar("bkg_flau","",6.48517e-01);
  bkg_meanlandau = new RooRealVar("bkg_meanlandau","",1.00181e+00);
  bkg_sigmalandau = new RooRealVar("bkg_sigmalandau","",4.35740e-01);
  RooLandau bkg_landau("bkg_landau", "bkg_landau", *Psi1To2Significance, *bkg_meanlandau, *bkg_sigmalandau);
  RooChebychev bkg_polyshape("bkg_polyshape","",*Psi1To2Significance,RooArgList(*bkg_co0,*bkg_co1));
  RooAddPdf bkg_distT("bkg_distT","", RooArgList(bkg_landau,bkg_polyshape),RooArgList(*bkg_flau));

  RooProdPdf bkg_mass1("bkg_mass1","",RooArgList(bkg_jpsi1mass,bkg_jpsi1mass_Pol,bkg_CT1_SB1,bkg_distT));


  // sideband-jpsi case
  bkg_jpsi2_mass_1 = new RooRealVar("bkg_jpsi2_mass_1","",3.08933e+00);
  bkg_jpsi2_width_1 = new RooRealVar("bkg_jpsi2_width_1","",5.06183e-02);
  bkg_jpsi2_width_2 = new RooRealVar("bkg_jpsi2_width_2","",5.82996e-01);
  bkg_jpsi2_width_a = new RooFormulaVar("bkg_jpsi1_width_a","","@0*@1",RooArgList(*bkg_jpsi2_width_1,*bkg_jpsi2_width_2));
  bkg_frac_7 = new RooRealVar("bkg_frac_7","",4.47328e-01);

  bkg_p0 = new RooRealVar("bkg_p0","",-2.93132e-01);
  bkg_p1 = new RooRealVar("bkg_p1","",-3.89092e-01);
  bkg_p2 = new RooRealVar("bkg_p2","",1.94808e-01);
  
  bkg3_R_mean_core = new RooRealVar("bkg3_R_mean_core","",3.60489e-02);
  bkg3_R_mean_tail = new RooRealVar("bkg3_R_mean_tail","",4.32342e-03);
  bkg3_R_sigma_core = new RooRealVar("bkg3_R_sigma_core","",2.89854e-02);
  bkg3_R_sigma_tail = new RooRealVar("bkg3_R_sigma_tail","",3.60637e-01);
  bkg3_R_sigma_tot = new RooFormulaVar("bkg3_R_sigma_tot","","@0*@1",RooArgList(*bkg3_R_sigma_tail,*bkg3_R_sigma_core));
  bkg_frac_3 = new RooRealVar("bkg_frac_3","",6.27677e-02);

  //RooGaussian bkg_jpsi2mass("bkg_jpsi2mass","jpsi mass distribution",*Psi2_Mass,*bkg_jpsi2_mass_1,*bkg_jpsi2_width_1);
  RooGaussian bkg_jpsi2mass_1("bkg_jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*bkg_jpsi2_mass_1,*bkg_jpsi2_width_1);
  RooGaussian bkg_jpsi2mass_2("bkg_jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*bkg_jpsi2_mass_1,*bkg_jpsi2_width_a);
  RooAddPdf bkg_jpsi2mass("bkg_jpsi2mass","jpsi2 mass distribution",RooArgList(jpsi2mass_1,jpsi2mass_2),RooArgList(*frac_2));
  RooChebychev bkg_jpsi2mass_Pol("bkg_jpsi2mass_Pol","",*Psi1_Mass,RooArgList(*bkg_p0,*bkg_p1,*bkg_p2));

  RooGaussian bkg_ctSB2_1a("bkg_ctSB2_1a","ct distribution",*Psi1_CTxy,*bkg3_R_mean_core,*bkg3_R_sigma_core);
  RooGaussian bkg_ctSB2_1b("bkg_ctSB2_1b","ct distribution",*Psi1_CTxy,*bkg3_R_mean_tail,*bkg3_R_sigma_tot);
  RooAddPdf bkg_CT1_SB2("bkg_CT1_SB2","ct distribution",RooArgList(bkg_ctSB2_1a,bkg_ctSB2_1b),RooArgList(*bkg_frac_3));

  bkg_p5_distT = new RooRealVar("bkg_p5_distT","",2.49205e-01);
  bkg_p6_distT = new RooRealVar("bkg_p6_distT","",1.04013e-01);
  bkg_lambda2 = new RooRealVar("bkg_lambda2","",6.54912e+00);
  RooGaussModel resolution_core4("resolution_core4","",*Psi1To2Significance,*bkg_p5_distT,*bkg_p6_distT);
  //RooDecay bkg_distT2("bkg_distT2","",*Psi1To2Significance,*bkg_lambda2,resolution_core4,RooDecay::SingleSided);

  bkg_co02 =  new RooRealVar("bkg_co02","", 6.58858e-01);
  bkg_co12 =  new RooRealVar("bkg_co12","", 2.48596e-05);
  bkg_flau2 = new RooRealVar("bkg_flau2","",5.33780e-01);
  bkg_meanlandau2 = new RooRealVar("bkg_meanlandau2","",1.09999e+00);
  bkg_sigmalandau2 = new RooRealVar("bkg_sigmalandau2","",4.64046e-01);
  RooChebychev bkg_polyshape2("bkg_polyshape2","",*Psi1To2Significance,RooArgList(*bkg_co02,*bkg_co12));
  RooLandau bkg_landau2("bkg_landau2", "bkg_landau2", *Psi1To2Significance, *bkg_meanlandau2, *bkg_sigmalandau2);
  RooAddPdf bkg_distT2("bkg_distT2","", RooArgList(bkg_landau2,bkg_polyshape2),RooArgList(*bkg_flau2));

  RooProdPdf bkg_mass2("bkg_mass2","",RooArgList(bkg_jpsi2mass,bkg_jpsi2mass_Pol,bkg_CT1_SB2,bkg_distT2));

  // fraction for the J/psi-flat flat-J/psi ratio (Andrew: I think this is defined as (# J/psi-flat)/(# J/psi-flat + # flat-J/psi))
  //bkg_frac_5 = new RooRealVar("bkg_frac_5","",7.72600e-01,0.,1.); // for initial fit
  bkg_frac_5 = new RooRealVar("bkg_frac_5","",5.95167e-01); // for eff_cut dataset
  //bkg_frac_5 = new RooRealVar("bkg_frac_5","",6.47958e-01); // for all dataset
  RooAddPdf bkg_model("bkg_model","",RooArgList(bkg_mass1,bkg_mass2),RooArgList(*bkg_frac_5));

  // bkg flat flat

  bkg7_R_mean_core = new RooRealVar("bkg7_R_mean_core","",3.45913e-03);
  bkg7_R_sigma_core = new RooRealVar("bkg7_R_sigma_core","",1.10565e-02);

  RooGaussian bkg2_CT1("bkg2_CT1","ct distribution",*Psi1_CTxy,*bkg7_R_mean_core,*bkg7_R_sigma_core);

  bkg_p7_distT = new RooRealVar("bkg_p7_distT","",3.89967e-01);
  bkg_p8_distT = new RooRealVar("bkg_p8_distT","",3.31232e-01);
  bkg_lambda3 = new RooRealVar("bkg_lambda3","",4.67109e+00);
  RooGaussModel resolution_core5("resolution_core5","",*Psi1To2Significance,*bkg_p7_distT,*bkg_p8_distT);
  //RooDecay bkg2_distT("bkg2_distT","",*Psi1To2Significance,*bkg_lambda3,resolution_core5,RooDecay::SingleSided);

  bkg2_co0 =  new RooRealVar("bkg2_co0","", 3.86180e-01);
  bkg2_co1 =  new RooRealVar("bkg2_co1","", 9.49975e-01);
  bkg2_flau = new RooRealVar("bkg2_flau","",7.28332e-01);
  bkg2_meanlandau = new RooRealVar("bkg2_meanlandau","",1.56581e+00);
  bkg2_sigmalandau = new RooRealVar("bkg2_sigmalandau","",5.51089e-01);
  RooLandau bkg2_landau("bkg2_landau", "bkg2_landau", *Psi1To2Significance, *bkg2_meanlandau, *bkg2_sigmalandau);
  RooChebychev bkg2_polyshape("bkg2_polyshape","",*Psi1To2Significance,RooArgList(*bkg2_co0,*bkg2_co1));
  RooAddPdf bkg2_distT("bkg2_distT","", RooArgList(bkg2_landau,bkg2_polyshape),RooArgList(*bkg2_flau));
  
  RooProdPdf bkg2_model("bkg2_model","",RooArgList(bkg_jpsi2mass_Pol,bkg_jpsi1mass_Pol,bkg2_CT1,bkg2_distT));

  /////////////////////////////////////////////

  RooRealVar nsig("nsig","number of signal events",700,1,5000);
  RooRealVar nBbkg("nBbkg","number of B background events",1000,1,5000);
  RooRealVar nbkg("nbkg","number of background events",200,1,5000);
  RooRealVar nbkg2("nbkg2","number of background events",200,1,5000);

  RooAddPdf model("model","model",RooArgList(sig_model,Bbkg_model,bkg_model,bkg2_model),RooArgList(nsig,nBbkg,nbkg,nbkg2));
  model.fitTo(*data,Extended(kTRUE),Hesse(kTRUE),Minos(kFALSE));


  RooPlot* frame1 = Psi1_Mass->frame(Title("#mu^{+}#mu^{-}"));
  frame1 = Psi1_Mass->frame(20);
  data->plotOn(frame1);
  model.plotOn(frame1);
  model.plotOn(frame1,RooFit::Components("sig_*"),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
  model.plotOn(frame1,RooFit::Components("bkg_model"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
  model.plotOn(frame1,RooFit::Components("bkg2_model"),RooFit::LineColor(kBlack), RooFit::LineStyle(kDashed));
  model.plotOn(frame1,RooFit::Components("Bbkg_*"),RooFit::LineColor(6), RooFit::LineStyle(kDashed));
  TCanvas* c1 = new TCanvas("c1","",10,10,700,500);
  c1->cd();
  c1->SetFillColor(kWhite);
  frame1->SetTitle("");
  frame1->SetXTitle("#mu^{+}#mu^{-} 1 Invariant Mass [GeV/c^{2}]");
  frame1->Draw();
	frame1->SetStats(1);
  c1->SaveAs("pic/Psi1_mass.pdf");
  //c1->Close();

  RooPlot* frame2 = Psi2_Mass->frame(Title("#mu^{+}#mu^{-}"));
  frame2 = Psi2_Mass->frame(20);
  data->plotOn(frame2);
  model.plotOn(frame2);
  model.plotOn(frame2,RooFit::Components("Sig,sig_*"),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
  model.plotOn(frame2,RooFit::Components("bkg_model"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
  model.plotOn(frame2,RooFit::Components("bkg2_model"),RooFit::LineColor(kBlack), RooFit::LineStyle(kDashed));
  model.plotOn(frame2,RooFit::Components("BBkg,Bbkg_*"),RooFit::LineColor(6), RooFit::LineStyle(kDashed));
  TCanvas* c2 = new TCanvas("c2","",10,10,700,500);
  c2->cd();
  c2->SetFillColor(kWhite);
  frame2->SetTitle("");
  frame2->SetXTitle("#mu^{+}#mu^{-} 2 Invariant Mass [GeV/c^{2}]");
  frame2->Draw();
  c2->SaveAs("pic/Psi2_mass.pdf");
  //c2->Close();

  RooPlot* frame5 = Psi1_CTxy->frame(Title("#mu^{+}#mu^{-}#mu^{+}#mu^{-} ct"));
  frame5 = Psi1_CTxy->frame(30);
  data->plotOn(frame5);
  model.plotOn(frame5);
  model.plotOn(frame5,RooFit::Components("sig_*"),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
  model.plotOn(frame5,RooFit::Components("bkg_model"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
  model.plotOn(frame5,RooFit::Components("bkg2_model"),RooFit::LineColor(kBlack), RooFit::LineStyle(kDashed));
  model.plotOn(frame5,RooFit::Components("Bbkg_*"),RooFit::LineColor(6), RooFit::LineStyle(kDashed));
  TCanvas* c5 = new TCanvas("c5","",10,10,700,500);
  c5->cd();
  c5->SetFillColor(kWhite);
  frame5->SetTitle("");
  frame5->SetMaximum(3000);
  frame5->SetXTitle("J/#psi^{1} ct_{xy} [cm]");
  frame5->Draw();
  c5->SetLogy();
  c5->SaveAs("pic/Psi1_CTxy.pdf");
  //c5->Close();

  RooPlot* frame3 = Psi1To2Significance->frame(Title("#mu^{+}#mu^{-}#mu^{+}#mu^{-} ct"));
  frame3 = Psi1To2Significance->frame(20);
  data->plotOn(frame3);
  model.plotOn(frame3);
  model.plotOn(frame3,RooFit::Components("sig_*"),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
  model.plotOn(frame3,RooFit::Components("bkg_model"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
  model.plotOn(frame3,RooFit::Components("bkg2_model"),RooFit::LineColor(kBlack), RooFit::LineStyle(kDashed));
  model.plotOn(frame3,RooFit::Components("Bbkg_*"),RooFit::LineColor(6), RooFit::LineStyle(kDashed));
  TCanvas* c3 = new TCanvas("c3","",10,10,700,500);
  c3->cd();
  c3->SetFillColor(kWhite);
  frame3->SetTitle("");
  frame3->SetXTitle("J/#psi Distance Significance");
  frame3->Draw();
  c3->SaveAs("pic/Psi1To2Significance.pdf");
  //c3->Close();

  //std::cout << endl << " B-lifetime (Bbkg_lambda_CT1) found from floating: ";
  //Bbkg_lambda_CT1->printValue(std::cout);
  //std::cout << endl << " bkg1 fraction (bkg_frac_5) found from floating: ";
  //bkg_frac_5->printValue(std::cout);
  //std::cout << endl << endl;
  

  ///////////////////////////
  // LIKELIHOOD PLOT
  ///////////////////////////

  /*
  ProfileLikelihoodCalculator pl(*data,model,nsig);
  //  pl.SetConfidenceLevel(0.95);
  LikelihoodInterval* interval = pl.GetInterval();

  LikelihoodIntervalPlot plot(interval);
  TCanvas* c = new TCanvas("c","",200,200,700,500);
  c->cd();
  plot.Draw();
  */

}

