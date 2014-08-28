#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include <TPaveText.h>

using namespace RooFit ;
using namespace RooStats ;

void AddModel(RooWorkspace*);
void DoSPlot(RooWorkspace*);
void MakePlots(RooWorkspace*);

void SplotJpsiJpsi_2012()
{

  gROOT->SetStyle("Plain");

  RooWorkspace* wspace = new RooWorkspace("myWS");

  AddModel(wspace);
  DoSPlot(wspace);
  MakePlots(wspace);

  delete wspace;
  
}

 
//____________________________________
void AddModel(RooWorkspace* ws){

  // Variables
  RooRealVar* run = new RooRealVar("run","run number",100000,300000);
  RooRealVar* event = new RooRealVar("event","event number",0.,1e10);
  //  RooRealVar* FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",6.,100.);
  RooRealVar* FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",0.,999.);
  //  RooRealVar* FourMu_Mass = new RooRealVar("FourMu_Mass","FourMu_Mass",6.,20.);
  RooRealVar* Psi1_Mass = new RooRealVar("Psi1_Mass","Psi1_Mass",2.85,3.35);
  RooRealVar* Psi2_Mass = new RooRealVar("Psi2_Mass","Psi2_Mass",2.85,3.35);
  RooRealVar* Psi1_CTxy = new RooRealVar("Psi1_CTxy","Psi1_CTxy",-0.05,0.1);
  RooRealVar* Psi2_CTxy = new RooRealVar("Psi2_CTxy","Psi2_CTxy",-0.05,0.1);
  RooRealVar* Psi1To2_dY = new RooRealVar("Psi1To2_dY","Psi1To2_dY",0.,4.8);
  RooRealVar* Psi1To2Significance = new RooRealVar("Psi1To2Significance","Psi1To2Significance",0.,8);

  RooArgSet dataVars;
  dataVars.add(RooArgSet(*run,*event,*FourMu_Mass,*Psi1_Mass,*Psi2_Mass,*Psi1_CTxy,*Psi2_CTxy,*Psi1To2_dY,*Psi1To2Significance));
  RooDataSet *data;
  //TFile* f= (TFile*)gROOT->FindObject("./data.root"); if (f) f->Close();
  //f = new TFile ("./data.root","UPDATE");
  //TFile* f= (TFile*)gROOT->FindObject("./Input_To_Fit_pT_Sort.root"); if (f) f->Close();
  //TFile* f = new TFile ("./Input_To_Fit_pT_Sort_iter_7_Data_Good_Eff.root","READ");
  TFile* f = new TFile ("./Input_To_Fit_pT_Sort_2012_Data_04_tight_mass.root","READ");
  //TFile* f = new TFile ("./Input_To_Fit_pT_Sort_iter_7_Data.root","READ");
  TTree *myTree = (TTree*)f->Get("PATEventTree");
  data = new RooDataSet("data","data",myTree,dataVars);

  /////////////////////////
  // Signal PDFs parameters
  /////////////////////////

  RooRealVar*  jpsi1_mass_1 = new RooRealVar("jpsi1_mass_1","",3.08988e+00);
  RooRealVar*  jpsi1_mass_2 = new RooRealVar("jpsi1_mass_2","",3.09553e+00);
  RooRealVar*  jpsi1_width_1 = new RooRealVar("jpsi1_width_1","",7.00000e-02);
  RooRealVar*  jpsi1_width_2 = new RooRealVar("jpsi1_width_2","",4.30143e-01);
  RooFormulaVar*  jpsi1_width_a = new RooFormulaVar("jpsi1_width_a","","@0*@1",RooArgList(*jpsi1_width_1,*jpsi1_width_2));
  RooRealVar*  frac_1 = new RooRealVar("frac_1","",3.76426e-01);

  RooRealVar*  jpsi2_mass_1 = new RooRealVar("jpsi2_mass_1","",3.08680e+00);
  RooRealVar*  jpsi2_mass_2 = new RooRealVar("jpsi2_mass_2","",3.09285e+00);
  RooRealVar*  jpsi2_width_1 = new RooRealVar("jpsi2_width_1","",7.00000e-02);
  RooRealVar*  jpsi2_width_2 = new RooRealVar("jpsi2_width_2","",4.79895e-01);
  RooFormulaVar*  jpsi2_width_a = new RooFormulaVar("jpsi2_width_a","","@0*@1",RooArgList(*jpsi2_width_1,*jpsi2_width_2));
  RooRealVar*  frac_2 = new RooRealVar("frac_2","",5.00000e-01);  

  RooRealVar*  R_mean_coreCT1 = new RooRealVar("R_mean_coreCT1","",2.63860e-05);
  RooRealVar*  R_mean_tailCT1 = new RooRealVar("R_mean_tailCT1","",1.19874e-03);
  RooRealVar*  R_sigma_coreCT1 = new RooRealVar("R_sigma_coreCT1","",2.72201e-03);
  RooRealVar*  R_sigma_tailCT1 = new RooRealVar("R_sigma_tailCT1","",2.72258e+00);
  RooFormulaVar*  R_sigma_totCT1 = new RooFormulaVar("R_sigma_totCT1","","@0*@1",RooArgList(*R_sigma_coreCT1,*R_sigma_tailCT1));
  RooRealVar*  R_fracCT1 = new RooRealVar("R_fracCT1","",8.47797e-01);

  RooRealVar*  R_mean_core = new RooRealVar("R_mean_core","",4.08278e-01);
  RooRealVar*  R_sigma_core = new RooRealVar("R_sigma_core","",2.17748e-01);
  RooRealVar*  etab_lambda = new RooRealVar("etab_lambda","",6.76288e-01);

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

  /////////////////////////
  // B bkg PDF 
  /////////////////////////

  RooRealVar*   Bbkg_jpsi1_mass_1 = new RooRealVar("Bbkg_jpsi1_mass_1","",3.08988e+00);
  RooRealVar*   Bbkg_jpsi1_mass_2 = new RooRealVar("Bbkg_jpsi1_mass_2","",3.09553e+00);
  RooRealVar*   Bbkg_jpsi1_width_1 = new RooRealVar("Bbkg_jpsi1_width_1","",7.00000e-02);
  RooRealVar*   Bbkg_jpsi1_width_2 = new RooRealVar("Bbkg_jpsi1_width_2","",4.30143e-01);
  RooFormulaVar*   Bbkg_jpsi1_width_a = new RooFormulaVar("Bbkg_jpsi1_width_a","","@0*@1",RooArgList(*Bbkg_jpsi1_width_1,*Bbkg_jpsi1_width_2));
  RooRealVar*   Bbkg_frac_1 = new RooRealVar("Bbkg_frac_1","",3.76426e-01);

  RooRealVar*   Bbkg_jpsi2_mass_1 = new RooRealVar("Bbkg_jpsi2_mass_1","",3.08680e+00);
  RooRealVar*   Bbkg_jpsi2_mass_2 = new RooRealVar("Bbkg_jpsi2_mass_2","",3.09285e+00);
  RooRealVar*   Bbkg_jpsi2_width_1 = new RooRealVar("Bbkg_jpsi2_width_1","",7.00000e-02);
  RooRealVar*   Bbkg_jpsi2_width_2 = new RooRealVar("Bbkg_jpsi2_width_2","",4.79895e-01);
  RooFormulaVar*   Bbkg_jpsi2_width_a = new RooFormulaVar("Bbkg_jpsi2_width_a","","@0*@1",RooArgList(*Bbkg_jpsi2_width_1,*Bbkg_jpsi2_width_2));
  RooRealVar*   Bbkg_frac_2 = new RooRealVar("Bbkg_frac_2","",5.00000e-01);

  RooRealVar*   Bbkg_mean_CT1 = new RooRealVar("Bbkg_mean_CT1","",6.51060e-04);
  RooRealVar*   Bbkg_width_CT1 = new RooRealVar("Bbkg_width_CT1","",3.77146e-03);
  //RooRealVar*   Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","Bbkg_lambda_CT1",1.59424e-02); // for eff_cut dataset
  //RooRealVar*   Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","Bbkg_lambda_CT1",1.60933e-02); // for eff_cut dataset
  RooRealVar*   Bbkg_lambda_CT1 = new RooRealVar("Bbkg_lambda_CT1","Bbkg_lambda_CT1",1.86124e-02); // for all dataset

  RooRealVar*   Bbkg_p3_distT = new RooRealVar("Bbkg_p3_distT","",1.18073e+00);
  RooRealVar*   Bbkg_p4_distT = new RooRealVar("Bbkg_p4_distT","",5.15922e-01);
  RooRealVar*   Bbkg_lambda1 = new RooRealVar("Bbkg_lambda1","",9.99999e+01);

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

  /////////////////////////
  // bkg PDF 
  /////////////////////////

  // jpsi-sideband case
  RooRealVar*    bkg_jpsi1_mass_1 = new RooRealVar("bkg_jpsi1_mass_1","",3.08988e+00);
  //RooRealVar*    bkg_jpsi1_mass_2 = new RooRealVar("bkg_jpsi1_mass_2","",3.08988e+00);
  RooRealVar*    bkg_jpsi1_width_1 = new RooRealVar("bkg_jpsi1_width_1","",7.00000e-02);
  RooRealVar*    bkg_jpsi1_width_2 = new RooRealVar("bkg_jpsi1_width_2","",4.30143e-01);
  RooFormulaVar* bkg_jpsi1_width_a = new RooFormulaVar("bkg_jpsi1_width_a","","@0*@1",RooArgList(*bkg_jpsi1_width_1,*bkg_jpsi1_width_2));
  RooRealVar*    bkg_frac_6 = new RooRealVar("bkg_frac_6","",3.76426e-01);

  RooRealVar*    bkg_p3 = new RooRealVar("bkg_p3","",-2.10268e-01);
  RooRealVar*    bkg_p4 = new RooRealVar("bkg_p4","",-1.95504e-01);
  RooRealVar*    bkg_p5 = new RooRealVar("bkg_p5","",5.00755e-02);

  RooRealVar*    bkg1_R_mean_core = new RooRealVar("bkg1_R_mean_core","",1.04093e-02);
  RooRealVar*    bkg1_R_mean_tail = new RooRealVar("bkg1_R_mean_tail","",-2.51869e-04);
  RooRealVar*    bkg1_R_sigma_core = new RooRealVar("bkg1_R_sigma_core","",1.24609e-02);
  RooRealVar*    bkg1_R_sigma_tail = new RooRealVar("bkg1_R_sigma_tail","",2.63519e-01);
  RooRealVar*    bkg_frac_1 = new RooRealVar("bkg_frac_1","",7.07283e-01);
  //RooRealVar*    bkg1_R_mean_core = new RooRealVar("bkg1_R_mean_core","",2.63860e-05);
  //RooRealVar*    bkg1_R_mean_tail = new RooRealVar("bkg1_R_mean_tail","",1.19874e-03);
  //RooRealVar*    bkg1_R_sigma_core = new RooRealVar("bkg1_R_sigma_core","",2.72201e-03);
  //RooRealVar*    bkg1_R_sigma_tail = new RooRealVar("bkg1_R_sigma_tail","",2.72258e+00);
  //RooRealVar*    bkg_frac_1 = new RooRealVar("bkg_frac_1","",8.47797e-01);

  RooFormulaVar* bkg1_R_sigma_tot = new RooFormulaVar("bkg1_R_sigma_tot","","@0*@1",RooArgList(*bkg1_R_sigma_tail,*bkg1_R_sigma_core));
  //RooGaussian bkg_jpsi1mass("bkg_jpsi1mass","jpsi mass distribution",*Psi1_Mass,*bkg_jpsi1_mass_1,*bkg_jpsi1_width_1);
  RooGaussian bkg_jpsi1mass_1("bkg_jpsi1mass_1","jpsi mass distribution",*Psi1_Mass,*bkg_jpsi1_mass_1,*bkg_jpsi1_width_1);
  RooGaussian bkg_jpsi1mass_2("bkg_jpsi1mass_2","jpsi mass distribution",*Psi1_Mass,*bkg_jpsi1_mass_1,*bkg_jpsi1_width_a);
  RooAddPdf bkg_jpsi1mass("bkg_jpsi1mass","jpsi1 mass distribution",RooArgList(bkg_jpsi1mass_1,bkg_jpsi1mass_2),RooArgList(*bkg_frac_6));
  RooChebychev bkg_jpsi1mass_Pol("bkg_jpsi1mass_Pol","",*Psi2_Mass,RooArgList(*bkg_p3,*bkg_p4,*bkg_p5));

  RooGaussian bkg_ctSB1_1a("bkg_ctSB1_1a","ct distribution",*Psi1_CTxy,*bkg1_R_mean_core,*bkg1_R_sigma_core);
  RooGaussian bkg_ctSB1_1b("bkg_ctSB1_1b","ct distribution",*Psi1_CTxy,*bkg1_R_mean_tail,*bkg1_R_sigma_tot);
  RooAddPdf bkg_CT1_SB1("bkg_CT1_SB1","ct distribution",RooArgList(bkg_ctSB1_1a,bkg_ctSB1_1b),RooArgList(*bkg_frac_1));

  //RooRealVar*    bkg_p3_distT = new RooRealVar("bkg_p3_distT","",2.87578e-01);
  //RooRealVar*    bkg_p4_distT = new RooRealVar("bkg_p4_distT","",1.26380e-01);
  //RooRealVar*    bkg_lambda1 = new RooRealVar("bkg_lambda1","",3.65276e+00);
  //RooGaussModel resolution_core3("resolution_core3","",*Psi1To2Significance,*bkg_p3_distT,*bkg_p4_distT);
  //RooDecay bkg_distT("bkg_distT","",*Psi1To2Significance,*bkg_lambda1,resolution_core3,RooDecay::SingleSided);

  RooRealVar* bkg_co0 =  new RooRealVar("bkg_co0","",9.99840e-01);
  RooRealVar* bkg_co1 =  new RooRealVar("bkg_co1","",2.00471e-06);
  RooRealVar* bkg_flau = new RooRealVar("bkg_flau","",6.48517e-01);
  RooRealVar* bkg_meanlandau = new RooRealVar("bkg_meanlandau","",1.00181e+00);
  RooRealVar* bkg_sigmalandau = new RooRealVar("bkg_sigmalandau","",4.35740e-01);
  RooLandau bkg_landau("bkg_landau", "bkg_landau", *Psi1To2Significance, *bkg_meanlandau, *bkg_sigmalandau);
  RooChebychev bkg_polyshape("bkg_polyshape","",*Psi1To2Significance,RooArgList(*bkg_co0,*bkg_co1));
  RooAddPdf bkg_distT("bkg_distT","", RooArgList(bkg_landau,bkg_polyshape),*bkg_flau);

  RooProdPdf bkg_mass1("bkg_mass1","",RooArgList(bkg_jpsi1mass,bkg_jpsi1mass_Pol,bkg_CT1_SB1,bkg_distT));

  // sideband-jpsi case
  RooRealVar*    bkg_jpsi2_mass_1 = new RooRealVar("bkg_jpsi2_mass_1","",3.08680e+00);
  //RooRealVar*    bkg_jpsi2_mass_2 = new RooRealVar("bkg_jpsi2_mass_2","",3.08680e+00);
  RooRealVar*    bkg_jpsi2_width_1 = new RooRealVar("bkg_jpsi2_width_1","",7.00000e-02);
  RooRealVar*    bkg_jpsi2_width_2 = new RooRealVar("bkg_jpsi2_width_2","",4.79895e-01);
  RooFormulaVar* bkg_jpsi2_width_a = new RooFormulaVar("bkg_jpsi1_width_a","","@0*@1",RooArgList(*bkg_jpsi2_width_1,*bkg_jpsi2_width_2));
  RooRealVar*    bkg_frac_7 = new RooRealVar("bkg_frac_7","",5.00000e-01);

  RooRealVar*  bkg_p0 = new RooRealVar("bkg_p0","",-2.93132e-01);
  RooRealVar*  bkg_p1 = new RooRealVar("bkg_p1","",-3.89092e-01);
  RooRealVar*  bkg_p2 = new RooRealVar("bkg_p2","",1.94808e-01);
  
  RooRealVar*  bkg3_R_mean_core = new RooRealVar("bkg3_R_mean_core","",3.60489e-02);
  RooRealVar*  bkg3_R_mean_tail = new RooRealVar("bkg3_R_mean_tail","",4.32342e-03);
  RooRealVar*  bkg3_R_sigma_core = new RooRealVar("bkg3_R_sigma_core","",2.89854e-02);
  RooRealVar*  bkg3_R_sigma_tail = new RooRealVar("bkg3_R_sigma_tail","",3.60637e-01);
  RooRealVar*  bkg_frac_3 = new RooRealVar("bkg_frac_3","",6.27677e-02);
  RooFormulaVar*  bkg3_R_sigma_tot = new RooFormulaVar("bkg3_R_sigma_tot","","@0*@1",RooArgList(*bkg3_R_sigma_tail,*bkg3_R_sigma_core));

  //RooGaussian bkg_jpsi2mass("bkg_jpsi2mass","jpsi mass distribution",*Psi2_Mass,*bkg_jpsi2_mass_1,*bkg_jpsi2_width_1);
  RooGaussian bkg_jpsi2mass_1("bkg_jpsi2mass_1","jpsi mass distribution",*Psi2_Mass,*bkg_jpsi2_mass_1,*bkg_jpsi2_width_1);
  RooGaussian bkg_jpsi2mass_2("bkg_jpsi2mass_2","jpsi mass distribution",*Psi2_Mass,*bkg_jpsi2_mass_1,*bkg_jpsi2_width_a);
  RooAddPdf bkg_jpsi2mass("bkg_jpsi2mass","jpsi2 mass distribution",RooArgList(bkg_jpsi2mass_1,bkg_jpsi2mass_2),RooArgList(*bkg_frac_7));
  RooChebychev bkg_jpsi2mass_Pol("bkg_jpsi2mass_Pol","",*Psi1_Mass,RooArgList(*bkg_p0,*bkg_p1,*bkg_p2));

  RooGaussian bkg_ctSB2_1a("bkg_ctSB2_1a","ct distribution",*Psi1_CTxy,*bkg3_R_mean_core,*bkg3_R_sigma_core);
  RooGaussian bkg_ctSB2_1b("bkg_ctSB2_1b","ct distribution",*Psi1_CTxy,*bkg3_R_mean_tail,*bkg3_R_sigma_tot);
  RooAddPdf bkg_CT1_SB2("bkg_CT1_SB2","ct distribution",RooArgList(bkg_ctSB2_1a,bkg_ctSB2_1b),RooArgList(*bkg_frac_3));

  //RooRealVar*  bkg_p5_distT = new RooRealVar("bkg_p5_distT","",2.41297e-01);
  //RooRealVar*  bkg_p6_distT = new RooRealVar("bkg_p6_distT","",7.29316e-02);
  //RooRealVar*  bkg_lambda2 = new RooRealVar("bkg_lambda2","",4.66802e+00);
  //RooGaussModel resolution_core4("resolution_core4","",*Psi1To2Significance,*bkg_p5_distT,*bkg_p6_distT);
  //RooDecay bkg_distT2("bkg_distT2","",*Psi1To2Significance,*bkg_lambda2,resolution_core4,RooDecay::SingleSided);

  RooRealVar* bkg_co02 =  new RooRealVar("bkg_co02","", 6.58858e-01);
  RooRealVar* bkg_co12 =  new RooRealVar("bkg_co12","", 2.48596e-05);
  RooRealVar* bkg_flau2 = new RooRealVar("bkg_flau2","",5.33780e-01);
  RooRealVar* bkg_meanlandau2 = new RooRealVar("bkg_meanlandau2","",1.09999e+00);
  RooRealVar* bkg_sigmalandau2 = new RooRealVar("bkg_sigmalandau2","",4.64046e-01);
  RooChebychev bkg_polyshape2("bkg_polyshape2","",*Psi1To2Significance,RooArgList(*bkg_co02,*bkg_co12));
  RooLandau bkg_landau2("bkg_landau2", "bkg_landau2", *Psi1To2Significance, *bkg_meanlandau2, *bkg_sigmalandau2);
  RooAddPdf bkg_distT2("bkg_distT2","", RooArgList(bkg_landau2,bkg_polyshape2),*bkg_flau2);

  RooProdPdf bkg_mass2("bkg_mass2","",RooArgList(bkg_jpsi2mass,bkg_jpsi2mass_Pol,bkg_CT1_SB2,bkg_distT2));

  // fraction for the J/psi-flat flat-J/psi ratio (Andrew: I think this is defined as (# J/psi-flat)/(# J/psi-flat + # flat-J/psi))
  //RooRealVar*  bkg_frac_5 = new RooRealVar("bkg_frac_5","",);
  //RooRealVar*  bkg_frac_5 = new RooRealVar("bkg_frac_5","",5.95167e-01); // for eff_cut dataset
  RooRealVar*  bkg_frac_5 = new RooRealVar("bkg_frac_5","",8.92796e-01); // for all dataset
  RooAddPdf bkg_model("bkg_model","",RooArgList(bkg_mass1,bkg_mass2),RooArgList(*bkg_frac_5));

  /////////////////////////
  // bkg flat flat
  /////////////////////////
  RooRealVar*  bkg7_R_mean_core = new RooRealVar("bkg7_R_mean_core","",3.45913e-03);
  RooRealVar*  bkg7_R_sigma_core = new RooRealVar("bkg7_R_sigma_core","",1.10565e-02);
  RooGaussian bkg2_CT1("bkg2_CT1","ct distribution",*Psi1_CTxy,*bkg7_R_mean_core,*bkg7_R_sigma_core);

  //RooRealVar*  bkg_p7_distT = new RooRealVar("bkg_p7_distT","",5.13453e-01);
  //RooRealVar*  bkg_p8_distT = new RooRealVar("bkg_p8_distT","",2.85147e-01);
  //RooRealVar*  bkg_lambda3 = new RooRealVar("bkg_lambda3","",3.54163e+00);
  //RooGaussModel resolution_core5("resolution_core5","",*Psi1To2Significance,*bkg_p7_distT,*bkg_p8_distT);
  //RooDecay bkg2_distT("bkg2_distT","",*Psi1To2Significance,*bkg_lambda3,resolution_core5,RooDecay::SingleSided);

  RooRealVar* bkg2_co0 =  new RooRealVar("bkg2_co0","", 3.86180e-01);
  RooRealVar* bkg2_co1 =  new RooRealVar("bkg2_co1","", 9.49975e-01);
  RooRealVar* bkg2_flau = new RooRealVar("bkg2_flau","",7.28332e-01);
  RooRealVar* bkg2_meanlandau = new RooRealVar("bkg2_meanlandau","",1.56581e+00);
  RooRealVar* bkg2_sigmalandau = new RooRealVar("bkg2_sigmalandau","",5.51089e-01);
  RooLandau bkg2_landau("bkg2_landau", "bkg2_landau", *Psi1To2Significance, *bkg2_meanlandau, *bkg2_sigmalandau);
  RooChebychev bkg2_polyshape("bkg2_polyshape","",*Psi1To2Significance,RooArgList(*bkg2_co0,*bkg2_co1));
  RooAddPdf bkg2_distT("bkg2_distT","", RooArgList(bkg2_landau,bkg2_polyshape),*bkg2_flau);

  RooProdPdf bkg2_model("bkg2_model","",RooArgList(bkg_jpsi2mass_Pol,bkg_jpsi1mass_Pol,bkg2_CT1,bkg2_distT));

  /////////////////////////////////////////////

  RooRealVar nsig("nsig","number of signal events",700,1,5000);
  RooRealVar nBbkg("nBbkg","number of B background events",1000,1,5000);
  RooRealVar nbkg("nbkg","number of background events",200,1,5000);
  RooRealVar nbkg2("nbkg2","number of background events",200,1,5000);

  RooAddPdf model("model","model",RooArgList(sig_model,Bbkg_model,bkg_model,bkg2_model),RooArgList(nsig,nBbkg,nbkg,nbkg2));

  std::cout << "import model" << std::endl;
  ws->import(model);
  std::cout << "import data" << std::endl;
  ws->import(*data);

  cout << "data entries: " << data->numEntries() << endl;

  // fitting dY with SPS and DPS

  RooRealVar*   m_dY = new RooRealVar("m_dY","",3.11781e+00);
  RooRealVar*   w_dY = new RooRealVar("w_dY","",4.99424e-01);
  /*
  RooRealVar*   p0_dY = new RooRealVar("p0_dY","",-8.06009e-01,-1.,1.);
  RooRealVar*   p1_dY = new RooRealVar("p1_dY","",3.19965e-01,-1.,1.);
  RooRealVar*   p2_dY = new RooRealVar("p2_dY","",-4.41787e-02,-1.,1.);
  */
  RooRealVar*   p0_dY = new RooRealVar("p0_dY","",-82.9767/3182.94);
  RooRealVar*   p1_dY = new RooRealVar("p1_dY","",-61.3831/3182.94);
  RooRealVar*   frac_dY = new RooRealVar("frac_dY","",1.61012e-01);
  RooGaussian g_dY("g_dY"," ",*Psi1To2_dY,*m_dY,*w_dY);
  //  RooPolynomial pol_dY("pol_dY"," ",*Psi1To2_dY,RooArgList(*p0_dY,*p1_dY,*p2_dY));
  //  RooAddPdf model_DPS("model_DPS"," ",RooArgList(g_dY,pol_dY),*frac_dY);
  //  RooPolynomial model_DPS("model_DPS"," ",*Psi1To2_dY,RooArgList(*p0_dY,*p1_dY,*p2_dY));
  RooPolynomial model_DPS("model_DPS"," ",*Psi1To2_dY,RooArgList(*p0_dY,*p1_dY));
  //  RooPolynomial model_DPS("model_DPS"," ",*Psi1To2_dY,RooArgList(*p0_dY));

  RooRealVar*  Rm_dY = new RooRealVar("Rm_dY"," " , -5.82342e-01);
  RooRealVar*  Rw_dY = new RooRealVar("Rw_dY"," " , 7.35732e-01);
  RooRealVar*  lambda_dY = new RooRealVar("lambda_dY"," " , 3.26722e-01);
  RooGaussModel res_dY("res_dY","",*Psi1To2_dY,*Rm_dY,*Rw_dY);
  RooDecay model_SPS("model_SPS","",*Psi1To2_dY,*lambda_dY,res_dY,RooDecay::SingleSided);

  RooRealVar nsigDPS("nsigDPS","number of signal events",300,1,500);
  RooRealVar nsigSPS("nsigSPS","number of signal events",200,1,500);

  RooAddPdf model_dY("model_dY","model_dY",RooArgList(model_DPS,model_SPS),RooArgList(nsigDPS,nsigSPS));
  std::cout << "import model for dY fitting" << std::endl;
  ws->import(model_dY);

}

//____________________________________
void DoSPlot(RooWorkspace* ws){
  std::cout << "Calculate sWeights" << std::endl;

  RooAbsPdf* model = ws->pdf("model");
  RooRealVar* nsig = ws->var("nsig");
  RooRealVar* nBbkg = ws->var("nBbkg");
  RooRealVar* nbkg = ws->var("nbkg");
  RooRealVar* nbkg2 = ws->var("nbkg2");
  RooDataSet* data = (RooDataSet*) ws->data("data");

  // fit the model to the data.
  model->fitTo(*data, Extended() );

  RooMsgService::instance().setSilentMode(true);

  // Now we use the SPlot class to add SWeights to our data set
  // based on our model and our yield variables
  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
					       *data, model, RooArgList(*nsig,*nBbkg,*nbkg,*nbkg2) );


  // Check that our weights have the desired properties

  std::cout << "Check SWeights:" << std::endl;

  std::cout << std::endl <<  "Yield of sig is " 
	    << nsig->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("nsig") << std::endl;

  std::cout << std::endl <<  "Yield of Bbkg is " 
	    << nBbkg->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("nBbkg") << std::endl;

  std::cout << std::endl <<  "Yield of bkg is " 
	    << nbkg->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("nbkg") << std::endl;

  std::cout << std::endl <<  "Yield of bkg2 is " 
	    << nbkg2->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("nbkg2") << std::endl;

  cout << endl;   cout << endl;   cout << endl;
  float sum20=0;
  float sum50=0;
  float sum100=0;
  float sum200=0;
  float sum300=0;
  float sum600=0;
  float sum900=0;
  float sum1200=0;
  float total=0;

  // saving weights into a file
  ofstream myfile;
  myfile.open ("weights.txt");
  // plot the weight event by event with the Sum of events values as cross-check
  for(Int_t i=0; i < data->numEntries(); i++) {
      //myfile << sData->GetSWeight(i,"nsig") << " " << sData->GetSWeight(i,"nBbkg") << " " << sData->GetSWeight(i,"nbkg") << " " << sData->GetSWeight(i,"nbkg2") << endl;  
      //myfile << sData->GetSWeight(i,"nsig") <<endl;
    myfile << (unsigned int) data->get(i)->getRealValue("run")
             << " " << (unsigned int) data->get(i)->getRealValue("event")
	   << " " << (float) data->get(i)->getRealValue("FourMu_Mass")
             << " " << sData->GetSWeight(i,"nsig")
             << endl;
     // std::cout << "nsig Weight   " << sData->GetSWeight(i,"nsig") 
     //		<< "   nBbkg Weight   " << sData->GetSWeight(i,"nBbkg")
     //		<< "   nbkg Weight   " << sData->GetSWeight(i,"nbkg")
     //		<< "   nbkg2 Weight  " << sData->GetSWeight(i,"nbkg2")
//		<< "   Total Weight   " << sData->GetSumOfEventSWeight(i) 
//		<< std::endl;
      total+=sData->GetSWeight(i,"nsig");         
      if(i<20) sum20+=sData->GetSWeight(i,"nsig");
      if(i<50) sum50+=sData->GetSWeight(i,"nsig");
      if(i<100) sum100+=sData->GetSWeight(i,"nsig");
      if(i<200) sum200+=sData->GetSWeight(i,"nsig");
      if(i<300) sum300+=sData->GetSWeight(i,"nsig");
      if(i<600) sum600+=sData->GetSWeight(i,"nsig");
      if(i<900) sum900+=sData->GetSWeight(i,"nsig");
      if(i<1200) sum1200+=sData->GetSWeight(i,"nsig");

    }
  myfile.close();

  std::cout << std::endl;

  std::cout<<"Sum of the sWeights is: "<<total<<std::endl;
  std::cout<<"Sum of the first 20 sWeights is: "<<sum20<<std::endl;
  std::cout<<"Sum of the first 50 sWeights is: "<<sum50<<std::endl;
  std::cout<<"Sum of the first 100 sWeights is: "<<sum100<<std::endl;
  std::cout<<"Sum of the first 200 sWeights is: "<<sum200<<std::endl;
  std::cout<<"Sum of the first 300 sWeights is: "<<sum300<<std::endl;
  std::cout<<"Sum of the first 600 sWeights is: "<<sum600<<std::endl;
  std::cout<<"Sum of the first 900 sWeights is: "<<sum900<<std::endl;
  std::cout<<"Sum of the first 1200 sWeights is: "<<sum1200<<std::endl;
  std::cout<<"Total # of events: "<<data->numEntries()<<std::endl;

  // import this new dataset with sWeights
  std::cout << "import new dataset with sWeights" << std::endl;
  ws->import(*data, Rename("dataWithSWeights"));

}

void MakePlots(RooWorkspace* ws){

  std::cout << "make plots" << std::endl;

  RooAbsPdf* model = ws->pdf("model");
  RooAbsPdf* model_dY = ws->pdf("model_dY");
  RooRealVar* nsig = ws->var("nsig");
  RooRealVar* nsigDPS = ws->var("nsigDPS");
  RooRealVar* nsigSPS = ws->var("nsigSPS");
  RooRealVar* nBbkg = ws->var("nBbkg");
  RooRealVar* nbkg = ws->var("nbkg");
  RooRealVar* nbkg2 = ws->var("nbkg2");
  RooRealVar* FourMu_Mass = ws->var("FourMu_Mass");
  RooRealVar* Psi1_Mass = ws->var("Psi1_Mass");
  RooRealVar* Psi2_Mass = ws->var("Psi2_Mass");
  RooRealVar* Psi1_CTxy = ws->var("Psi1_CTxy");
  RooRealVar* Psi2_CTxy = ws->var("Psi2_CTxy");
  RooRealVar* Psi1To2_dY = ws->var("Psi1To2_dY");
  RooRealVar* Psi1To2Significance = ws->var("Psi1To2Significance");

  // note, we get the dataset with sWeights
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
  model->fitTo(*data, Extended() );

  // make TTree with efficiency variation info
  //TFile *fFile = new TFile("Fit_Results_Eff_Cut.root","recreate");
  setTDRStyle();
  // make our canvas
  TCanvas* cmass1 = new TCanvas("sPlotMass1","sPlotMass1", 600, 600);
  cmass1->cd();
  cmass1->SetFillColor(kWhite);

  RooPlot* Mass1Plot = Psi1_Mass->frame(20);
  data->plotOn(Mass1Plot,Name("data"), DataError(RooAbsData::SumW2));
  model->plotOn(Mass1Plot,Name("all"));
  model->plotOn(Mass1Plot,Name("sig"),RooFit::Components("Sig,sig_*"),RooFit::LineColor(kRed), RooFit::LineStyle(2));
  model->plotOn(Mass1Plot,Name("bkg"),RooFit::Components("bkg_model"),RooFit::LineColor(kGreen), RooFit::LineStyle(3));
  model->plotOn(Mass1Plot,Name("bkg2"),RooFit::Components("bkg2_model"),RooFit::LineColor(kBlack), RooFit::LineStyle(4));
  model->plotOn(Mass1Plot,Name("Bbkg"),RooFit::Components("BBkg,Bbkg_*"),RooFit::LineColor(6), RooFit::LineStyle(7));
  Mass1Plot->SetTitle("");
  Mass1Plot->SetXTitle("#mu^{+}#mu^{-} 1 Invariant Mass (GeV/c^{2})");
  Mass1Plot->SetYTitle("Events / 0.025 GeV/c^{2}");
  Mass1Plot->SetLabelOffset(0.012);
  //Mass1Plot->SetTitleOffset(0.95);
  Mass1Plot->Draw();
  //cmass1->SaveAs("pic/Psi1_mass.pdf");
  //cmass1->Close();

  TCanvas* cmass2 = new TCanvas("sPlotMass2","sPlotMass2", 600, 600);
  cmass2->cd();
  cmass2->SetFillColor(kWhite);
  RooPlot* Mass2Plot = Psi2_Mass->frame(20); 
  data->plotOn(Mass2Plot,Name("data"), DataError(RooAbsData::SumW2)); 
  model->plotOn(Mass2Plot,Name("all"));   
  model->plotOn(Mass2Plot,Name("sig"),RooFit::Components("Sig,sig_*"),RooFit::LineColor(kRed), RooFit::LineStyle(2));
  model->plotOn(Mass2Plot,Name("bkg"),RooFit::Components("bkg_model"),RooFit::LineColor(kGreen), RooFit::LineStyle(3));
  model->plotOn(Mass2Plot,Name("bkg2"),RooFit::Components("bkg2_model"),RooFit::LineColor(kBlack), RooFit::LineStyle(4));
  model->plotOn(Mass2Plot,Name("Bbkg"),RooFit::Components("BBkg,Bbkg_*"),RooFit::LineColor(6), RooFit::LineStyle(7));
  Mass2Plot->SetTitle("");
  Mass2Plot->SetXTitle("#mu^{+}#mu^{-} 2 Invariant Mass (GeV/c^{2})");
  Mass2Plot->SetYTitle("Events / 0.025 GeV/c^{2}");
  Mass2Plot->SetLabelOffset(0.012);
  //Mass2Plot->SetTitleOffset(0.95);
  Mass2Plot->Draw();
  //cmass2->SaveAs("pic/Psi2_mass.pdf");
  //cmass2->Close();

  TCanvas* cctxy1 = new TCanvas("sPlotCTxy1","sPlotCTxy1", 600, 600);
  cctxy1->cd();
  cctxy1->SetFillColor(kWhite);
  RooPlot* CTxy1Plot = Psi1_CTxy->frame(30); 
  data->plotOn(CTxy1Plot,Name("data"), DataError(RooAbsData::SumW2)); 
  model->plotOn(CTxy1Plot,Name("all"));   
  model->plotOn(CTxy1Plot,Name("sig"),RooFit::Components("Sig,sig_*"),RooFit::LineColor(kRed), RooFit::LineStyle(2));
  model->plotOn(CTxy1Plot,Name("bkg"),RooFit::Components("bkg_model"),RooFit::LineColor(kGreen), RooFit::LineStyle(3));
  model->plotOn(CTxy1Plot,Name("bkg2"),RooFit::Components("bkg2_model"),RooFit::LineColor(kBlack), RooFit::LineStyle(4));
  model->plotOn(CTxy1Plot,Name("Bbkg"),RooFit::Components("BBkg,Bbkg_*"),RooFit::LineColor(6), RooFit::LineStyle(7));
  CTxy1Plot->SetTitle("");
  CTxy1Plot->SetXTitle("J/#psi^{1} ct_{xy} (cm)");
  CTxy1Plot->SetYTitle("Events / 0.005 cm");
  CTxy1Plot->SetMaximum(2000);
  CTxy1Plot->SetMinimum(0.1);
  CTxy1Plot->Draw();
  cctxy1->SetLogy();
  //cctxy1->SaveAs("pic/Psi1_CTxy.pdf");
  //cctxy1->Close();

  TCanvas* csig = new TCanvas("sPlotSig","sPlotSig", 600, 600);
  csig->cd();
  RooPlot* SigPlot = Psi1To2Significance->frame(20); 
  data->plotOn(SigPlot,Name("data"), DataError(RooAbsData::SumW2)); 
  model->plotOn(SigPlot,Name("all"));   
  model->plotOn(SigPlot,Name("sig"),RooFit::Components("Sig,sig_*"),RooFit::LineColor(kRed), RooFit::LineStyle(2));
  model->plotOn(SigPlot,Name("bkg"),RooFit::Components("bkg_model"),RooFit::LineColor(kGreen), RooFit::LineStyle(3));
  model->plotOn(SigPlot,Name("bkg2"),RooFit::Components("bkg2_model"),RooFit::LineColor(kBlack), RooFit::LineStyle(4));
  model->plotOn(SigPlot,Name("Bbkg"),RooFit::Components("BBkg,Bbkg_*"),RooFit::LineColor(6), RooFit::LineStyle(7));
  SigPlot->SetTitle("");
  SigPlot->SetYTitle("Events / 0.4");
  SigPlot->SetXTitle("J/#psi Distance Significance");
  SigPlot->Draw();
  //csig->SaveAs("pic/Psi1To2Significance.pdf");
  //csig->Close();

  // create weighted data set (signal-weighted)
  //RooDataSet * dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"nsig_sw");
  //  model_dY->fitTo(*data);
  //TCanvas* cdata2 = new TCanvas("sPlot2","sPlots2", 700, 500);
  //cdata2->cd();
  //RooPlot* frame2 = Psi1To2_dY->frame(20); 
  //dataw_sig->plotOn(frame2, DataError(RooAbsData::SumW2) ); 
  /*
  model_dY->fitTo(*dataw_sig, SumW2Error(kTRUE), Extended());
  model_dY->plotOn(frame2);
  model_dY->plotOn(frame2,RooFit::Components("*DPS"),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
  model_dY->plotOn(frame2,RooFit::Components("*SPS"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));
  */
  //frame2->SetTitle("DeltaY distribution for sig");
  //frame2->SetMinimum(1e-03);
  //frame2->Draw();

  //TCanvas* cdata3 = new TCanvas("sPlot3","sPlots3", 700, 500);
  //cdata3->cd();
  //RooPlot* frame3 = FourMu_Mass->frame(Bins(14),Range(6.,20.)); 
  //dataw_sig->plotOn(frame3, DataError(RooAbsData::SumW2));
  //frame3->SetTitle("FourMu_Mass distribution for sig");
  //frame3->SetMinimum(1e-03);
  //frame3->Draw();

  //fFile->Write();
  //fFile->Close();
  //delete fFile;

}
