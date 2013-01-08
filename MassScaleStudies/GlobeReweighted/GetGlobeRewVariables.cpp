// g++ -Wall -o GetGlobeRewVariables `root-config --cflags --glibs` GetGlobeRewVariables.cpp

//example
// ./GetGlobeRewVariables 2012

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TRandom.h"
#include "TVirtualFitter.h"
#include "TLatex.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

TH1F* templateHisto;
TF1* templateFunc;
std::vector<double>* mydata;


int main(int argc, char** argv)
{
//   // Set style options
//   setTDRStyle();
//   gStyle->SetPadTickX(1);
//   gStyle->SetPadTickY(1);
//   gStyle->SetOptTitle(0); 
//   gStyle->SetOptStat(1110); 
//   gStyle->SetOptFit(1); 
  
  // Set fitting options
  
  
  //-----------------
  // Input parameters
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
//   char* EBEE = argv[1];
//   char* LOWHIGH = argv[2];
  std::string string_year = argv[1];
  int year = atoi(argv[1]);

//   std::cout << "EBEE:         " << EBEE         << std::endl;
//   std::cout << "LOWHIGH:      " << LOWHIGH       << std::endl;
  std::cout << "year:      " << year       << std::endl;
  

  // Get trees
  std::cout << std::endl;

  std::string nameNtuplesMCZee = "DYJetsToLL";
  std::string nameNtuplesMCggh = "ggh_m125_8TeV";
  std::string nameNtuplesMCvbf = "vbf_m125_8TeV";
  std::string nameNtuplesMCwzh = "wzh_m125_8TeV";
  std::string nameNtuplesMCtth = "tth_m125_8TeV";

  if(year == 2011) {
    nameNtuplesMCggh = "gluglu_H_gg_125_pu2011";
    nameNtuplesMCvbf = "vbf_H_gg_125_pu2011";
    nameNtuplesMCwzh = "wz_H_gg_125_pu2011";
    nameNtuplesMCtth = "tt_H_gg_125_pu2011";
  }

  TChain* ntu_MCZee = new TChain(nameNtuplesMCZee.c_str());
  TChain* ntu_MCggh = new TChain(nameNtuplesMCggh.c_str());
  TChain* ntu_MCvbf = new TChain(nameNtuplesMCvbf.c_str());
  TChain* ntu_MCwzh = new TChain(nameNtuplesMCwzh.c_str());
  TChain* ntu_MCtth = new TChain(nameNtuplesMCtth.c_str());

  //2012                                                                                                                                    
  if(year == 2012){
    //with Smearing and Correction
    ntu_MCZee->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_zee.root");
    ntu_MCggh->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_hgg125.root");
    ntu_MCvbf->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_hgg125.root");
    ntu_MCwzh->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_hgg125.root");
    ntu_MCtth->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_hgg125.root");
  }
  //2011                                                                                                                                    
  if(year == 2011){
    //with Smearing and Correction
     ntu_MCZee->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_zee_7TeV.root");
     ntu_MCggh->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_hgg125_7TeV.root");
     ntu_MCvbf->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_hgg125_7TeV.root");
     ntu_MCwzh->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_hgg125_7TeV.root");
     ntu_MCtth->Add("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_hgg125_7TeV.root");
  }

  std::cout << "     REFERENCE ntu_MCZee: " << std::setw(8) << ntu_MCZee->GetEntries() << " entries" << std::endl;
  std::cout << "     REFERENCE ntu_MCggh: " << std::setw(8) << ntu_MCggh->GetEntries() << " entries" << std::endl;
  std::cout << "     REFERENCE ntu_MCvbf: " << std::setw(8) << ntu_MCvbf->GetEntries() << " entries" << std::endl;
  std::cout << "     REFERENCE ntu_MCwzh: " << std::setw(8) << ntu_MCwzh->GetEntries() << " entries" << std::endl;
  std::cout << "     REFERENCE ntu_MCtth: " << std::setw(8) << ntu_MCtth->GetEntries() << " entries" << std::endl;

  
  if(ntu_MCZee->GetEntries() == 0 || 
     ntu_MCggh->GetEntries() == 0 || ntu_MCvbf->GetEntries() == 0 || 
     ntu_MCwzh->GetEntries() == 0 || ntu_MCtth->GetEntries() == 0 )
  {
    std::cout << "Error: At least one file is empty" << std::endl; 
    return -1;
  }
  
  // Set branch addresses
  int nVtx;
  float npu;
  // Electron data
  float scEne1, scEneReg1, R9_pho1, scEta1, ES1, P1, scERaw1, e3x31, e5x51;
  float scEne2, scEneReg2, R9_pho2,  scEta2, ES2, P2, scERaw2, e3x32, e5x52;
  int isEB1,isEB2;
  double scPt1 = 0.;
  double scPt2 = 0.;
  double Energy1 = 0.;
  double Energy2 = 0.;

  
  ntu_MCZee->SetBranchStatus("*",0);
  ntu_MCZee->SetBranchStatus("weight",1);                    ntu_MCZee->SetBranchAddress("weight",&npu);

  ntu_MCZee->SetBranchStatus("pho1_sceta", 1);               ntu_MCZee->SetBranchAddress("pho1_sceta", &scEta1);    
  ntu_MCZee->SetBranchStatus("pho1_energy_regr", 1);         ntu_MCZee->SetBranchAddress("pho1_energy_regr", &scEneReg1);
  ntu_MCZee->SetBranchStatus("pho1_r9", 1);                  ntu_MCZee->SetBranchAddress("pho1_r9", &R9_pho1);
  ntu_MCZee->SetBranchStatus("pho1_pt", 1);                  ntu_MCZee->SetBranchAddress("pho1_pt", &scPt1);

  ntu_MCZee->SetBranchStatus("pho2_sceta", 1);               ntu_MCZee->SetBranchAddress("pho2_sceta", &scEta2);
  ntu_MCZee->SetBranchStatus("pho2_energy_regr", 1);         ntu_MCZee->SetBranchAddress("pho2_energy_regr", &scEneReg2);
  ntu_MCZee->SetBranchStatus("pho2_r9", 1);                  ntu_MCZee->SetBranchAddress("pho2_r9", &R9_pho2);
  ntu_MCZee->SetBranchStatus("pho2_pt", 1);                  ntu_MCZee->SetBranchAddress("pho2_pt", &scPt2);
  /////////// 
  
  ntu_MCggh->SetBranchStatus("*",0);
  ntu_MCggh->SetBranchStatus("weight",1);                    ntu_MCggh->SetBranchAddress("weight",&npu);

  ntu_MCggh->SetBranchStatus("pho1_sceta", 1);               ntu_MCggh->SetBranchAddress("pho1_sceta", &scEta1);
  ntu_MCggh->SetBranchStatus("pho1_energy_regr", 1);         ntu_MCggh->SetBranchAddress("pho1_energy_regr", &scEneReg1);
  ntu_MCggh->SetBranchStatus("pho1_r9", 1);                  ntu_MCggh->SetBranchAddress("pho1_r9", &R9_pho1);
  ntu_MCggh->SetBranchStatus("pho1_pt", 1);                  ntu_MCggh->SetBranchAddress("pho1_pt", &scPt1);
  ntu_MCggh->SetBranchStatus("pho2_sceta", 1);               ntu_MCggh->SetBranchAddress("pho2_sceta", &scEta2);
  ntu_MCggh->SetBranchStatus("pho2_energy_regr", 1);         ntu_MCggh->SetBranchAddress("pho2_energy_regr", &scEneReg2);
  ntu_MCggh->SetBranchStatus("pho2_r9", 1);                  ntu_MCggh->SetBranchAddress("pho2_r9", &R9_pho2);
  ntu_MCggh->SetBranchStatus("pho2_pt", 1);                  ntu_MCggh->SetBranchAddress("pho2_pt", &scPt2);
  /////////////
  
  ntu_MCvbf->SetBranchStatus("*",0);
  ntu_MCvbf->SetBranchStatus("weight",1);                    ntu_MCvbf->SetBranchAddress("weight",&npu);

  ntu_MCvbf->SetBranchStatus("pho1_sceta", 1);               ntu_MCvbf->SetBranchAddress("pho1_sceta", &scEta1);
  ntu_MCvbf->SetBranchStatus("pho1_energy_regr", 1);         ntu_MCvbf->SetBranchAddress("pho1_energy", &scEneReg1);
  ntu_MCvbf->SetBranchStatus("pho1_r9", 1);                  ntu_MCvbf->SetBranchAddress("pho1_r9", &R9_pho1);
  ntu_MCvbf->SetBranchStatus("pho1_pt", 1);                  ntu_MCvbf->SetBranchAddress("pho1_pt", &scPt1);
  ntu_MCvbf->SetBranchStatus("pho2_sceta", 1);               ntu_MCvbf->SetBranchAddress("pho2_sceta", &scEta2);
  ntu_MCvbf->SetBranchStatus("pho2_energy_regr", 1);         ntu_MCvbf->SetBranchAddress("pho2_energy_regr", &scEneReg2);
  ntu_MCvbf->SetBranchStatus("pho2_r9", 1);                  ntu_MCvbf->SetBranchAddress("pho2_r9", &R9_pho2);
  ntu_MCvbf->SetBranchStatus("pho2_pt", 1);                  ntu_MCvbf->SetBranchAddress("pho2_pt", &scPt2);
  /////////////

  ntu_MCwzh->SetBranchStatus("*",0);
  ntu_MCwzh->SetBranchStatus("weight",1);                    ntu_MCwzh->SetBranchAddress("weight",&npu);

  ntu_MCwzh->SetBranchStatus("pho1_sceta", 1);               ntu_MCwzh->SetBranchAddress("pho1_sceta", &scEta1);
  ntu_MCwzh->SetBranchStatus("pho1_energy_regr", 1);         ntu_MCwzh->SetBranchAddress("pho1_energy_regr", &scEneReg1);
  ntu_MCwzh->SetBranchStatus("pho1_r9", 1);                  ntu_MCwzh->SetBranchAddress("pho1_r9", &R9_pho1);
  ntu_MCwzh->SetBranchStatus("pho1_pt", 1);                  ntu_MCwzh->SetBranchAddress("pho1_pt", &scPt1);
  ntu_MCwzh->SetBranchStatus("pho2_sceta", 1);               ntu_MCwzh->SetBranchAddress("pho2_sceta", &scEta2);
  ntu_MCwzh->SetBranchStatus("pho2_energy_regr", 1);         ntu_MCwzh->SetBranchAddress("pho2_energy_regr", &scEneReg2);
  ntu_MCwzh->SetBranchStatus("pho2_r9", 1);                  ntu_MCwzh->SetBranchAddress("pho2_r9", &R9_pho2);
  ntu_MCwzh->SetBranchStatus("pho2_pt", 1);                  ntu_MCwzh->SetBranchAddress("pho2_pt", &scPt2);
  /////////////

  ntu_MCtth->SetBranchStatus("*",0);
  ntu_MCtth->SetBranchStatus("weight",1);                    ntu_MCtth->SetBranchAddress("weight",&npu);

  ntu_MCtth->SetBranchStatus("pho1_sceta", 1);               ntu_MCtth->SetBranchAddress("pho1_sceta", &scEta1);
  ntu_MCtth->SetBranchStatus("pho1_energy_regr", 1);         ntu_MCtth->SetBranchAddress("pho1_energy_regr", &scEneReg1);
  ntu_MCtth->SetBranchStatus("pho1_r9", 1);                  ntu_MCtth->SetBranchAddress("pho1_r9", &R9_pho1);
  ntu_MCtth->SetBranchStatus("pho1_pt", 1);                  ntu_MCtth->SetBranchAddress("pho1_pt", &scPt1);
  ntu_MCtth->SetBranchStatus("pho2_sceta", 1);               ntu_MCtth->SetBranchAddress("pho2_sceta", &scEta2);
  ntu_MCtth->SetBranchStatus("pho2_energy_regr", 1);         ntu_MCtth->SetBranchAddress("pho2_energy_regr", &scEneReg2);
  ntu_MCtth->SetBranchStatus("pho2_r9", 1);                  ntu_MCtth->SetBranchAddress("pho2_r9", &R9_pho2);
  ntu_MCtth->SetBranchStatus("pho2_pt", 1);                  ntu_MCtth->SetBranchAddress("pho2_pt", &scPt2);
  /////////////

  TH1F* h_Et_allZee_scE_reg = new TH1F("h_Et_allZee_scE_reg", "", 1000, 0., 1000.);
  TH1F* h_Et_allggh_scE_reg = new TH1F("h_Et_allggh_scE_reg", "", 1000, 0., 1000.);
  TH1F* h_Et_allvbf_scE_reg = new TH1F("h_Et_allvbf_scE_reg", "", 1000, 0., 1000.);
  TH1F* h_Et_allwzh_scE_reg = new TH1F("h_Et_allwzh_scE_reg", "", 1000, 0., 1000.);
  TH1F* h_Et_alltth_scE_reg = new TH1F("h_Et_alltth_scE_reg", "", 1000, 0., 1000.);

  TH1F* h_Eta_allZee_scE_reg = new TH1F("h_Eta_allZee_scE_reg", "", 1000, -3., 3.);
  TH1F* h_Eta_allggh_scE_reg = new TH1F("h_Eta_allggh_scE_reg", "", 1000, -3., 3.);
  TH1F* h_Eta_allvbf_scE_reg = new TH1F("h_Eta_allvbf_scE_reg", "", 1000, -3., 3.);
  TH1F* h_Eta_allwzh_scE_reg = new TH1F("h_Eta_allwzh_scE_reg", "", 1000, -3., 3.);
  TH1F* h_Eta_alltth_scE_reg = new TH1F("h_Eta_alltth_scE_reg", "", 1000, -3., 3.);

  TH1F* h_Eta_allZee_HR9 = new TH1F("h_Eta_allZee_HR9", "", 1000, -3., 3.);
  TH1F* h_Eta_allggh_HR9 = new TH1F("h_Eta_allggh_HR9", "", 1000, -3., 3.);
  TH1F* h_Eta_allvbf_HR9 = new TH1F("h_Eta_allvbf_HR9", "", 1000, -3., 3.);
  TH1F* h_Eta_allwzh_HR9 = new TH1F("h_Eta_allwzh_HR9", "", 1000, -3., 3.);
  TH1F* h_Eta_alltth_HR9 = new TH1F("h_Eta_alltth_HR9", "", 1000, -3., 3.);

  TH1F* h_Eta_allZee_LR9 = new TH1F("h_Eta_allZee_LR9", "", 1000, -3., 3.);
  TH1F* h_Eta_allggh_LR9 = new TH1F("h_Eta_allggh_LR9", "", 1000, -3., 3.);
  TH1F* h_Eta_allvbf_LR9 = new TH1F("h_Eta_allvbf_LR9", "", 1000, -3., 3.);
  TH1F* h_Eta_allwzh_LR9 = new TH1F("h_Eta_allwzh_LR9", "", 1000, -3., 3.);
  TH1F* h_Eta_alltth_LR9 = new TH1F("h_Eta_alltth_LR9", "", 1000, -3., 3.);


  TH1F* h_R9_allZee_scE_reg = new TH1F("h_R9_allZee_scE_reg", "", 1000, 0., 1.2);
  TH1F* h_R9_allggh_scE_reg = new TH1F("h_R9_allggh_scE_reg", "", 1000, 0., 1.2);
  TH1F* h_R9_allvbf_scE_reg = new TH1F("h_R9_allvbf_scE_reg", "", 1000, 0., 1.2);
  TH1F* h_R9_allwzh_scE_reg = new TH1F("h_R9_allwzh_scE_reg", "", 1000, 0., 1.2);
  TH1F* h_R9_alltth_scE_reg = new TH1F("h_R9_alltth_scE_reg", "", 1000, 0., 1.2);

  h_Et_allZee_scE_reg->Sumw2();
  h_Et_allggh_scE_reg->Sumw2(); 
  h_Et_allvbf_scE_reg->Sumw2(); 
  h_Et_allwzh_scE_reg->Sumw2(); 
  h_Et_alltth_scE_reg->Sumw2(); 

  h_Eta_allZee_scE_reg->Sumw2();
  h_Eta_allggh_scE_reg->Sumw2();
  h_Eta_allvbf_scE_reg->Sumw2();
  h_Eta_allwzh_scE_reg->Sumw2();
  h_Eta_alltth_scE_reg->Sumw2();

  h_Eta_allZee_HR9->Sumw2();
  h_Eta_allggh_HR9->Sumw2();
  h_Eta_allvbf_HR9->Sumw2();
  h_Eta_allwzh_HR9->Sumw2();
  h_Eta_alltth_HR9->Sumw2();

  h_Eta_allZee_LR9->Sumw2();
  h_Eta_allggh_LR9->Sumw2();
  h_Eta_allvbf_LR9->Sumw2();
  h_Eta_allwzh_LR9->Sumw2();
  h_Eta_alltth_LR9->Sumw2();

  h_R9_allZee_scE_reg->Sumw2(); 
  h_R9_allggh_scE_reg->Sumw2(); 
  h_R9_allvbf_scE_reg->Sumw2(); 
  h_R9_allwzh_scE_reg->Sumw2(); 
  h_R9_alltth_scE_reg->Sumw2(); 
  

  for(int ientry = 0; ientry<ntu_MCZee->GetEntries(); ++ientry)
  {
  	if( (ientry%100000 == 0) ) std::cout << "reading MC entry " << ientry << "\r" << std::flush;
        ntu_MCZee->GetEntry(ientry);  

	if( (fabs(scEta1) < 1.4442 || fabs(scEta1) > 1.566) && fabs(scEta1) < 2.7) { //no crack1;
	  float Rt1 = sin(2*atan(exp(-scEta1)) );
	  if(scEneReg1*Rt1 > 25.) {
	    h_Et_allZee_scE_reg->Fill(scEneReg1*Rt1, npu);
	    h_Eta_allZee_scE_reg->Fill(scEta1, npu);
	    h_R9_allZee_scE_reg->Fill(R9_pho1, npu);
	    if(R9_pho1 > 0.94) h_Eta_allZee_HR9->Fill(scEta1, npu);
	    if(R9_pho1 < 0.94) h_Eta_allZee_LR9->Fill(scEta1, npu);
	  }
	}

	///// second
	if( (fabs(scEta2) < 1.4442 || fabs(scEta2) > 1.566) && fabs(scEta2) < 2.7) { //no crack2
	  float Rt2 = sin(2*atan(exp(-scEta2)) );
	  if(scEneReg2*Rt2 > 25.)	{
	    h_Et_allZee_scE_reg->Fill(scEneReg2*Rt2, npu);
	    h_Eta_allZee_scE_reg->Fill(scEta2, npu);
	    h_R9_allZee_scE_reg->Fill(R9_pho2, npu);
	    if(R9_pho2 > 0.94) h_Eta_allZee_HR9->Fill(scEta2, npu);
	    if(R9_pho2 < 0.94) h_Eta_allZee_LR9->Fill(scEta2, npu);
	  }
	}
  }


  ///// Hgg
  for(int ientry = 0; ientry<ntu_MCggh->GetEntries(); ++ientry)
  {
  	if( (ientry%100000 == 0) ) std::cout << "reading MC entry " << ientry << "\r" << std::flush;
        ntu_MCggh->GetEntry(ientry);  

	if( (fabs(scEta1) < 1.4442 || fabs(scEta1) > 1.566) && fabs(scEta1) < 2.7) { //no crack1;
	  float Rt1 = sin(2*atan(exp(-scEta1)) );
	  if(scEneReg1*Rt1 > 25.){
	    h_Et_allggh_scE_reg->Fill(scEneReg1*Rt1, npu);
	    h_Eta_allggh_scE_reg->Fill(scEta1, npu);
	    h_R9_allggh_scE_reg->Fill(R9_pho1, npu);
	    if(R9_pho1 > 0.94) h_Eta_allggh_HR9->Fill(scEta1, npu);
	    if(R9_pho1 < 0.94) h_Eta_allggh_LR9->Fill(scEta1, npu);
	  }
	}

	///// second
	if( (fabs(scEta2) < 1.4442 || fabs(scEta2) > 1.566) && fabs(scEta2) < 2.7) { //no crack2
	  float Rt2 = sin(2*atan(exp(-scEta2)) );
	  if(scEneReg2*Rt2 > 25.){
	    h_Et_allggh_scE_reg->Fill(scEneReg2*Rt2, npu);
	    h_Eta_allggh_scE_reg->Fill(scEta2, npu);
	    h_R9_allggh_scE_reg->Fill(R9_pho2, npu);
	    if(R9_pho2 > 0.94) h_Eta_allggh_HR9->Fill(scEta2, npu);
	    if(R9_pho2 < 0.94) h_Eta_allggh_LR9->Fill(scEta2, npu);
	  }
	}
  }


  ///////////
  for(int ientry = 0; ientry<ntu_MCvbf->GetEntries(); ++ientry)
  {
  	if( (ientry%100000 == 0) ) std::cout << "reading MC entry " << ientry << "\r" << std::flush;
        ntu_MCvbf->GetEntry(ientry);  

	if( (fabs(scEta1) < 1.4442 || fabs(scEta1) > 1.566) && fabs(scEta1) < 2.7) { //no crack1
	  float Rt1 = sin(2*atan(exp(-scEta1)) );
	  if(scEneReg1*Rt1 > 25.)	{
	    h_Et_allvbf_scE_reg->Fill(scEneReg1*Rt1, npu);
	    h_Eta_allvbf_scE_reg->Fill(scEta1, npu);
	    h_R9_allvbf_scE_reg->Fill(R9_pho1, npu);
	    if(R9_pho1 > 0.94) h_Eta_allvbf_HR9->Fill(scEta1, npu);
	    if(R9_pho1 < 0.94) h_Eta_allvbf_LR9->Fill(scEta1, npu);
	  }
	}
	///// second
	if( (fabs(scEta2) < 1.4442 || fabs(scEta2) > 1.566) && fabs(scEta2) < 2.7) { //no crack2 
	  float Rt2 = sin(2*atan(exp(-scEta2)) );
	  if(scEneReg2*Rt2 > 25.)	{
	    h_Et_allvbf_scE_reg->Fill(scEneReg2*Rt2, npu);
	    h_Eta_allvbf_scE_reg->Fill(scEta2, npu);
	    h_R9_allvbf_scE_reg->Fill(R9_pho2, npu);
	    if(R9_pho2 > 0.94) h_Eta_allvbf_HR9->Fill(scEta2, npu);
	    if(R9_pho2 < 0.94) h_Eta_allvbf_LR9->Fill(scEta2, npu);
	  }
	}
  }

  /////////
  for(int ientry = 0; ientry<ntu_MCwzh->GetEntries(); ++ientry)
  {
  	if( (ientry%100000 == 0) ) std::cout << "reading MC entry " << ientry << "\r" << std::flush;
        ntu_MCwzh->GetEntry(ientry);  

	if( (fabs(scEta1) < 1.4442 || fabs(scEta1) > 1.566) && fabs(scEta1) < 2.7) { //no crack1
	  float Rt1 = sin(2*atan(exp(-scEta1)) );
	  if(scEneReg1*Rt1 > 25.)	{
	    h_Et_allwzh_scE_reg->Fill(scEneReg1*Rt1, npu);
	    h_Eta_allwzh_scE_reg->Fill(scEta1, npu);
	    h_R9_allwzh_scE_reg->Fill(R9_pho1, npu);
	    if(R9_pho1 > 0.94) h_Eta_allwzh_HR9->Fill(scEta1, npu);
	    if(R9_pho1 < 0.94) h_Eta_allwzh_LR9->Fill(scEta1, npu);
	  }
	}
	///// second
	if( (fabs(scEta2) < 1.4442 || fabs(scEta2) > 1.566) && fabs(scEta2) < 2.7) { //no crack2 
	  float Rt2 = sin(2*atan(exp(-scEta2)) );
	  if(scEneReg2*Rt2 > 25.)	{
	    h_Et_allwzh_scE_reg->Fill(scEneReg2*Rt2, npu);
	    h_Eta_allwzh_scE_reg->Fill(scEta2, npu);
	    h_R9_allwzh_scE_reg->Fill(R9_pho2, npu);
	    if(R9_pho2 > 0.94) h_Eta_allwzh_HR9->Fill(scEta2, npu);
	    if(R9_pho2 < 0.94) h_Eta_allwzh_LR9->Fill(scEta2, npu);
	  }
	}
  }

  /////////
  for(int ientry = 0; ientry<ntu_MCtth->GetEntries(); ++ientry)
  {
  	if( (ientry%100000 == 0) ) std::cout << "reading MC entry " << ientry << "\r" << std::flush;
        ntu_MCtth->GetEntry(ientry);  

	if( (fabs(scEta1) < 1.4442 || fabs(scEta1) > 1.566) && fabs(scEta1) < 2.7) { //no crack1
	  float Rt1 = sin(2*atan(exp(-scEta1)) );
	  if(scEneReg1*Rt1 > 25.)	{
	    h_Et_alltth_scE_reg->Fill(scEneReg1*Rt1, npu);
	    h_Eta_alltth_scE_reg->Fill(scEta1, npu);
	    h_R9_alltth_scE_reg->Fill(R9_pho1, npu);
	    if(R9_pho1 > 0.94) h_Eta_alltth_HR9->Fill(scEta1, npu);
	    if(R9_pho1 < 0.94) h_Eta_alltth_LR9->Fill(scEta1, npu);
	  }
	}
	///// second
	if( (fabs(scEta2) < 1.4442 || fabs(scEta2) > 1.566) && fabs(scEta2) < 2.7) { //no crack2
	  float Rt2 = sin(2*atan(exp(-scEta2)) );
	  if(scEneReg2*Rt2 > 25.)	{
	    h_Et_alltth_scE_reg->Fill(scEneReg2*Rt2, npu);
	    h_Eta_alltth_scE_reg->Fill(scEta2, npu);
	    h_R9_alltth_scE_reg->Fill(R9_pho2, npu);
	    if(R9_pho2 > 0.94) h_Eta_alltth_HR9->Fill(scEta2, npu);
	    if(R9_pho2 < 0.94) h_Eta_alltth_LR9->Fill(scEta2, npu);
	  }
	}
  }


  std::string plotFolderName = "results_Globe_"+string_year+".root";


  TFile pippo((plotFolderName).c_str(),"recreate");

  h_Et_allZee_scE_reg->Write();
  h_Et_allggh_scE_reg->Write();
  h_Et_allvbf_scE_reg->Write();
  h_Et_allwzh_scE_reg->Write();
  h_Et_alltth_scE_reg->Write();

  h_Eta_allZee_scE_reg->Write();
  h_Eta_allggh_scE_reg->Write();
  h_Eta_allvbf_scE_reg->Write();
  h_Eta_allwzh_scE_reg->Write();
  h_Eta_alltth_scE_reg->Write();

  h_Eta_allZee_HR9->Write();
  h_Eta_allggh_HR9->Write();
  h_Eta_allvbf_HR9->Write();
  h_Eta_allwzh_HR9->Write();
  h_Eta_alltth_HR9->Write();

  h_Eta_allZee_LR9->Write();
  h_Eta_allggh_LR9->Write();
  h_Eta_allvbf_LR9->Write();
  h_Eta_allwzh_LR9->Write();
  h_Eta_alltth_LR9->Write();

  h_R9_allZee_scE_reg->Write();
  h_R9_allggh_scE_reg->Write();
  h_R9_allvbf_scE_reg->Write();
  h_R9_allwzh_scE_reg->Write();
  h_R9_alltth_scE_reg->Write();

  pippo.Close();
 
  
  std::cout << " plottato tutto " << std::endl;
  //std::cout << "CREATI I FILES" << std::endl;
  return (0);
}
