//--*-C++-*--
/*****************************************************************************
 * Project: CMS detector at LHC, CERN
 * Package: RooFit
 *    File: $Id$
 * Authors:
 *   G.C. 
 * History:
 *   08-Aug-2008 TS started
 *
 * Copyright (C) 2008 UTK
 *****************************************************************************/
#ifndef ROOJPSIJPSIFIT
#define ROOJPSIJPSIFIT

#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategory.h"
#include "RooUnblindUniform.h"
#include "RooAbsPdf.h"
#include "TCanvas.h"

class RooRealVar;
class RooStringVar;
class RooAbsReal;
class RooAbsPdf;
class RooCategory;
class RooRealFrac;
class RooArgSet;


class RooJpsiJpsiFit {

public:

  RooJpsiJpsiFit();
  ~RooJpsiJpsiFit();

  RooDataSet* readData(TString dataFile);
  RooDataSet* readData(std::map<std::string, double> & dataFile);
  void plotData(std::map<std::string, double>& filename);
  void PDFmaker(std::map<std::string, double>& filename, TString pdf);
  void pureMCtoysGeneration();
  void fitData(std::map<std::string, double> & dataFile);
  void embedMCtoysGeneration();
  void PlotEmbedMCtoys();
  
protected:

  Bool_t _initDataVars;

  void initDataVars();

public:

  // variables

  RooArgSet dataVars;
  RooRealVar * Psi1_Mass;
  RooRealVar * Psi2_Mass;
  RooRealVar * FourMu_CTxy;
  RooRealVar * FourMu_pT;
  RooRealVar * Psi1To2_S;
  RooRealVar * Psi1To2_dY;
  RooRealVar * Psi1To2DistT;
  RooRealVar * Psi1To2DistTot;
  RooRealVar * Psi1To2Significance;
  RooRealVar * FourMu_Mass;
  RooRealVar * PvN;
  RooRealVar * Psi1_CTxy;
  RooRealVar * Psi2_CTxy;


  // parameters

  RooRealVar * frac_1;
  RooRealVar * frac_2;
  RooRealVar * Bbkg_frac_1;
  RooRealVar * Bbkg_frac_2;
  RooRealVar * bkg_frac_1;
  RooRealVar * bkg_frac_2;
  RooRealVar * bkg_frac_3;
  RooRealVar * bkg_frac_4;
  RooRealVar * bkg_frac_5;
  RooRealVar * bkg_frac_6;
  RooRealVar * bkg_frac_7;
  RooRealVar * bkg_frac_8;
  RooRealVar * jpsi1_f;
  RooRealVar * jpsi2_f;
  RooRealVar * bkg_jpsi1_f;
  RooRealVar * bkg_jpsi2_f;
  RooRealVar * R_f;

  RooRealVar *  jpsi1_mass_1;
  RooRealVar *  jpsi1_mass_2;
  RooRealVar *  jpsi1_mass_3;
  RooRealVar *  jpsi1_width_1;
  RooRealVar *  jpsi1_width_2;
  RooRealVar *  jpsi1_width_3;
  RooRealVar *  jpsi1_CBalpha;
  RooRealVar *  jpsi1_CBenne;
  RooFormulaVar *  jpsi1_width_a;
  RooFormulaVar *  jpsi1_width_b;

  RooRealVar *  jpsi2_mass_1;
  RooRealVar *  jpsi2_mass_2;
  RooRealVar *  jpsi2_mass_3;
  RooRealVar *  jpsi2_width_1;
  RooRealVar *  jpsi2_width_2;
  RooRealVar *  jpsi2_width_3;
  RooRealVar *  jpsi2_CBalpha;
  RooRealVar *  jpsi2_CBenne;
  RooFormulaVar *  jpsi2_width_a;
  RooFormulaVar *  jpsi2_width_b;

  RooRealVar *  bkg_jpsi1_mass_1;
  RooRealVar *  bkg_jpsi1_mass_2;
  RooRealVar *  bkg_jpsi1_mass_3;
  RooRealVar *  bkg_jpsi1_width_1;
  RooRealVar *  bkg_jpsi1_width_2;
  RooRealVar *  bkg_jpsi1_width_3;
  RooFormulaVar *  bkg_jpsi1_width_a;
  RooFormulaVar *  bkg_jpsi1_width_b;

  RooRealVar *  bkg_jpsi2_mass_1;
  RooRealVar *  bkg_jpsi2_mass_2;
  RooRealVar *  bkg_jpsi2_mass_3;
  RooRealVar *  bkg_jpsi2_width_1;
  RooRealVar *  bkg_jpsi2_width_2;
  RooRealVar *  bkg_jpsi2_width_3;
  RooFormulaVar *  bkg_jpsi2_width_a;
  RooFormulaVar *  bkg_jpsi2_width_b;

  RooRealVar *  Bbkg_jpsi1_mass_1;
  RooRealVar *  Bbkg_jpsi1_mass_2;
  RooRealVar *  Bbkg_jpsi1_mass_3;
  RooRealVar *  Bbkg_jpsi1_width_1;
  RooRealVar *  Bbkg_jpsi1_width_2;
  RooRealVar *  Bbkg_jpsi1_width_3;
  RooFormulaVar *  Bbkg_jpsi1_width_a;
  RooFormulaVar *  Bbkg_jpsi1_width_b;

  RooRealVar *  Bbkg_jpsi2_mass_1;
  RooRealVar *  Bbkg_jpsi2_mass_2;
  RooRealVar *  Bbkg_jpsi2_mass_3;
  RooRealVar *  Bbkg_jpsi2_width_1;
  RooRealVar *  Bbkg_jpsi2_width_2;
  RooRealVar *  Bbkg_jpsi2_width_3;
  RooFormulaVar *  Bbkg_jpsi2_width_a;
  RooFormulaVar *  Bbkg_jpsi2_width_b;

  RooRealVar *  bkg_p0;
  RooRealVar *  bkg_p1;
  RooRealVar *  bkg_p2;
  RooRealVar *  bkg_p3;
  RooRealVar *  bkg_p4;
  RooRealVar *  bkg_p5;

  RooRealVar *  bkg_p0_jpsi1;
  RooRealVar *  bkg_p1_jpsi1;
  RooRealVar *  bkg_p2_jpsi1;

  RooRealVar *  bkg_p0_jpsi2;
  RooRealVar *  bkg_p1_jpsi2;

  RooRealVar *  Bbkg_p0_jpsi1;
  RooRealVar *  Bbkg_p1_jpsi1;

  RooRealVar *  Bbkg_p0_jpsi2;
  RooRealVar *  Bbkg_p1_jpsi2;

  RooRealVar *  bkg_p0_distT;
  RooRealVar *  bkg_p1_distT;
  RooRealVar *  bkg_p2_distT;
  RooRealVar *  bkg_p3_distT;
  RooRealVar *  bkg_p4_distT;
  RooRealVar *  bkg_p5_distT;
  RooRealVar *  bkg_p6_distT;
  RooRealVar *  bkg_p7_distT;
  RooRealVar *  bkg_p8_distT;
  RooRealVar *  bkg_p9_distT;
  RooRealVar *  bkg_co0;
  RooRealVar *  bkg_co1;
  RooRealVar *  bkg_flau;
  RooRealVar *  bkg_meanlandau;
  RooRealVar *  bkg_sigmalandau;
  RooRealVar *  bkg_co02;
  RooRealVar *  bkg_co12;
  RooRealVar *  bkg_flau2;
  RooRealVar *  bkg_meanlandau2;
  RooRealVar *  bkg_sigmalandau2;
  RooRealVar *  bkg2_co0;
  RooRealVar *  bkg2_co1;
  RooRealVar *  bkg2_flau;
  RooRealVar *  bkg2_meanlandau;
  RooRealVar *  bkg2_sigmalandau;
  RooRealVar *  Bbkg_p0_distT;
  RooRealVar *  Bbkg_p1_distT;
  RooRealVar *  Bbkg_p2_distT;
  RooRealVar *  Bbkg_p3_distT;
  RooRealVar *  Bbkg_p4_distT;
  RooRealVar *  Bbkg_p5_distT;
  RooRealVar *  Bbkg_p6_distT;

  RooRealVar *  etab_lambda;
  RooRealVar *  bkg_lambda1;
  RooRealVar *  bkg_lambda2;
  RooRealVar *  bkg_lambda3;
  RooRealVar *  Bbkg_lambda1;
  RooRealVar *  Bbkg_lambda2;
  RooRealVar *  Bbkg_lambda;

  RooRealVar *  R_mean_core;
  RooRealVar *  R_sigma_core;
  RooRealVar *  R_mean_tail;
  RooRealVar *  R_sigma_tail;
  RooRealVar *  R_mean_coreCT1;
  RooRealVar *  R_sigma_coreCT1;
  RooRealVar *  R_mean_tailCT1;
  RooRealVar *  R_sigma_tailCT1;
  RooRealVar *  R_fracCT1;
  RooRealVar *  R_mean_coreCT2;
  RooRealVar *  R_sigma_coreCT2;
  RooRealVar *  R_mean_tailCT2;
  RooRealVar *  R_sigma_tailCT2;
  RooRealVar *  R_fracCT2;
  RooFormulaVar *  R_sigma_tot;
  RooFormulaVar * R_sigma_totCT1;
  RooFormulaVar * R_sigma_totCT2;

  RooRealVar *    bkg1_R_mean_core;
  RooRealVar *    bkg1_R_sigma_core;
  RooRealVar *    bkg1_R_mean_tail;
  RooRealVar *    bkg1_R_sigma_tail;
  RooFormulaVar *  bkg1_R_sigma_tot;

  RooRealVar *    bkg2_R_mean_core;
  RooRealVar *    bkg2_R_sigma_core;
  RooRealVar *    bkg2_R_mean_tail;
  RooRealVar *    bkg2_R_sigma_tail;
  RooFormulaVar *  bkg2_R_sigma_tot;

  RooRealVar *    bkg3_R_mean_core;
  RooRealVar *    bkg3_R_sigma_core;
  RooRealVar *    bkg3_R_mean_tail;
  RooRealVar *    bkg3_R_sigma_tail;
  RooFormulaVar *  bkg3_R_sigma_tot;

  RooRealVar *    bkg4_R_mean_core;
  RooRealVar *    bkg4_R_sigma_core;
  RooRealVar *    bkg4_R_mean_tail;
  RooRealVar *    bkg4_R_sigma_tail;
  RooFormulaVar *  bkg4_R_sigma_tot;

  RooRealVar *    bkg7_R_mean_core;
  RooRealVar *    bkg7_R_sigma_core;
  RooRealVar *    bkg7_R_mean_tail;
  RooRealVar *    bkg7_R_sigma_tail;
  RooFormulaVar *  bkg7_R_sigma_tot;

  RooRealVar *    bkg8_R_mean_core;
  RooRealVar *    bkg8_R_sigma_core;
  RooRealVar *    bkg8_R_mean_tail;
  RooRealVar *    bkg8_R_sigma_tail;
  RooFormulaVar *  bkg8_R_sigma_tot;

  RooRealVar *  Bbkg_mean_CT1;
  RooRealVar *  Bbkg_width_CT1;
  RooRealVar *  Bbkg_lambda_CT1;

  RooRealVar *  Bbkg_mean_CT2;
  RooRealVar *  Bbkg_width_CT2;
  RooRealVar *  Bbkg_lambda_CT2;

  RooRealVar *  bkg_mean_CT1_SB1;
  RooRealVar *  bkg_width_CT1_SB1;
  RooRealVar *  bkg_lambda_CT1_SB1;

  RooRealVar *  bkg_mean_CT2_SB1;
  RooRealVar *  bkg_width_CT2_SB1;
  RooRealVar *  bkg_lambda_CT2_SB1;

  RooRealVar *  bkg_mean_CT1_SB2;
  RooRealVar *  bkg_width_CT1_SB2;
  RooRealVar *  bkg_lambda_CT1_SB2;

  RooRealVar *  bkg_mean_CT2_SB2;
  RooRealVar *  bkg_width_CT2_SB2;
  RooRealVar *  bkg_lambda_CT2_SB2;

  RooRealVar *  bkg2_mean_CT1;
  RooRealVar *  bkg2_width_CT1;
  RooRealVar *  bkg2_lambda_CT1;

  RooRealVar *  bkg2_mean_CT2;
  RooRealVar *  bkg2_width_CT2;
  RooRealVar *  bkg2_lambda_CT2;

};

#endif

