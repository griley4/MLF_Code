#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <string>
#include <stdexcept>
#include <map>

using namespace std;

fit()
{

  
  gROOT->SetStyle("Plain");
  gSystem->Load("libRooFit");
  gSystem->Load("libGpad");
  gSystem->Load("tmp/libRooJpsiJpsiFit.so");
  
  using namespace RooFit;
  using namespace RooStats;

  RooJpsiJpsiFit fit;

  std::map<std::string, double> filename;

  //////////
  // sig
  //////////

  //filename["./Input_To_Fit_pT_Sort_iter_7_SPS_matched.root"] = 0.2008;
  //filename["./Input_To_Fit_pT_Sort_iter_7_DPS_matched.root"] = 1.;

  /////////////
  // Bkg
  /////////////

  //filename["./Input_To_Fit_pT_Sort_iter_7_Bbkg.root"] = 1.;

  //////////
  // Data
  //////////
	//filename["./Input_To_Fit_pT_Sort_iter_7_Data_2011_Good_Eff.root"] = 1.;
	//filename["./Input_To_Fit_pT_Sort_iter_7_Data_Good_Eff.root"] = 1.;
//filename["./Input_To_Fit_pT_Sort_2012_Data_04_tight_mass.root"] = 1.;
//filename["./Modified_Input_To_Fit_pT_Sort_2012_Data_04_tight_prob.root"] = 1.;
//filename["./Modified_Input_To_Fit_pT_Sort_2012_Data_pT_Select.root"] = 1.;
filename["./Input_To_Fit_pT_Sort_2012_Data_pT_Select.root"] = 1.;
//filename["./Input_To_Fit_NO_pT_Sort_2012_Data_04_tight_prob.root"] = 1.;



  //filename["./Input_To_Fit_pT_Sort_iter_7_2012_Data.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_2012_Data_02.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_2012_Data_03.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_iter_7_Data_Good_Eff.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_iter_7_Data.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_Good_Eff_dY_0_0.199.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_dY_0.645_4.4.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_Good_Eff_Mass_17.8_80.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_Mass_17.8_80.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_iter_7_jpsi_bkg.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_iter_7_bkg_jpsi.root"] = 1.;
  //filename["./Input_To_Fit_pT_Sort_iter_7_bkg_bkg.root"] = 1.;
  //filename[""] = 1.;

  ////////////////////////
  // PDFs
  ////////////////////////
  
  //  fit.PDFmaker(filename,"Pol3G_4muMass"); // Mass fit background
  //  fit.PDFmaker(filename,"GausPol3G_4muMass"); // Mass fit sig+bkgd
  //  fit.PDFmaker(filename,"Pol3G_dY"); // DPS
  //  fit.PDFmaker(filename,"GExpo_dY"); // SPS
  //  fit.PDFmaker(filename,"Expo_dY"); // B bkg / bkg-bkg
  //  fit.PDFmaker(filename,"Pol4_dY"); // J/psi-bkg  and bkg-J/psi

  //  fit.PDFmaker(filename,"plotM_4mu");
  //  fit.PDFmaker(filename,"plotM_jpsi1");

  //  J/psi Mass pdfs
  //  fit.PDFmaker(filename,"2G_jpsi1");
  //  fit.PDFmaker(filename,"2G_jpsi2");
  //  fit.PDFmaker(filename,"CBGaus_jpsi1");
  //  fit.PDFmaker(filename,"CBGaus_jpsi2");

  //  fit.PDFmaker(filename,"2G_Bjpsi1");
  //  fit.PDFmaker(filename,"2G_Bjpsi2");
  //  fit.PDFmaker(filename,"G_jpsi1");
  //  fit.PDFmaker(filename,"G_jpsi2");
  //  fit.PDFmaker(filename,"2G_jpsi1_b");
  //  fit.PDFmaker(filename,"2G_jpsi2_b");
  //  fit.PDFmaker(filename,"SB_GPol3_1");
  //  fit.PDFmaker(filename,"SB_GPol3_2");

  //  significance
  //  fit.PDFmaker(filename,"Significance_GExp");
  //  fit.PDFmaker(filename,"SignificanceBbkg_GExp");
  //  fit.PDFmaker(filename,"Significance_LandauPoly");

  //  proper decay length
  // signal
  //  fit.PDFmaker(filename,"2GCTxy1");
  //  fit.PDFmaker(filename,"2GCTxy2");

  // Bbkg
  //  fit.PDFmaker(filename,"GExpCTxy1");
  //  fit.PDFmaker(filename,"GExpCTxy2");
  // bkg
  //  fit.PDFmaker(filename,"2GCTxy1"); // for 4D
  //  fit.PDFmaker(filename,"GExpCTxy1"); 
  //  fit.PDFmaker(filename,"2GExpCTxy1");
  //  fit.PDFmaker(filename,"GCTxy1");

  //  fit.PDFmaker(filename,"2GCTxy2");  // for 4D
  //  fit.PDFmaker(filename,"GExpCTxy2"); 
  //  fit.PDFmaker(filename,"2GExpCTxy2");
  //  fit.PDFmaker(filename,"GCTxy2");

  //fit.pureMCtoysGeneration();
  fit.fitData(filename);

}
