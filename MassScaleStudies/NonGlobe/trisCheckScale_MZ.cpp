// g++ -Wall -o trisCheckScale_MZ `root-config --cflags --glibs` ../Utils/setTDRStyle.cc ../Utils/ntupleUtils.cc ../Utils/stabilityUtils.cc ../Utils/ConvoluteTemplate.cc ../Utils/histoFunc.h ../Utils/TPileupReweighting.h trisCheckScale_MZ.cpp

#include "../Utils/setTDRStyle.h"
#include "../Utils/histoFunc.h"
#include "../Utils/ConvoluteTemplate.h"
#include "../Utils/ntupleUtils.h"
#include "../Utils/stabilityUtils.h"
#include "../Utils/TPileupReweighting.h"
#include "../Utils/GetShervinCorrections.h"
#include "../Utils/GetSmearings.h"

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
  // Set style options
  setTDRStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110); 
  gStyle->SetOptFit(1); 
  
  // Set fitting options
  TVirtualFitter::SetDefaultFitter("Fumili2");
  

  // //   ////////////// vs Et                                                                                 
  TF1* Et_highR9_2011 = new TF1("Et_highR9_2011", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
  Et_highR9_2011->SetParameters(1.59984924630326465e-02, 4.14188316002253587e-02, -6.49126732859059939e-03);
  TF1* Et_lowR9_2011 = new TF1("Et_lowR9_2011", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
  Et_lowR9_2011->SetParameters(2.20638739628473586e-02, 6.98744642383235803e-02, -1.85601207959524978e-02);

  TF1* Et_highR9_2012 = new TF1("Et_highR9_2012", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
  Et_highR9_2012->SetParameters(1.76747992064786620e-02, 3.73408739026924591e-02, -7.82929065282905561e-03);
  TF1* Et_lowR9_2012 = new TF1("Et_lowR9_2012", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
  Et_lowR9_2012->SetParameters(1.97205016874162468e-02, 4.41133183909690751e-02, -1.58915655671104904e-02);



  //   bool UsePhotonRegression = false;
  bool UsePhotonRegression = true;


  //bool useShCorr = true;  
  bool useShCorr = false;


  //    bool correctEt = false; 
  bool correctEt = false;

  //  bool useMomentum = false;
  bool useMomentum = true;

  //-----------------
  // Input parameters
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  char* EBEE = argv[1];
  char* LOWHIGH = argv[2];
  char* ENE = argv[3];
  int PU = atoi(argv[4]);
  int evtsPerPoint = atoi(argv[5]);
  std::string string_year = argv[6];
  int year = atoi(argv[6]);
  std::string doVsEach = argv[7];


  std::cout << "EBEE:         " << EBEE         << std::endl;
  std::cout << "LOWHIGH:      " << LOWHIGH       << std::endl;
  std::cout << "ENE:          " << ENE           << std::endl;
  std::cout << "PU:           " << PU            << std::endl;
  std::cout << "evtsPerPoint: " << evtsPerPoint  << std::endl;
  std::cout << "year:      " << year       << std::endl;
  std::cout << "doVsEach:      " << doVsEach       << std::endl;
  
  TPileupReweighting* puReweighting;
  //2012 prompt           
  if(year == 2012) puReweighting =
 new TPileupReweighting("/afs/cern.ch/work/a/amartell/public/weights/PUweights_DYJetsToLL_Summer12_53X_ShSkim_ABC_TrueNumInteractions.root","hweights");
    //    new TPileupReweighting("/afs/cern.ch/work/a/amartell/public/weights/PUweights_DYJetsToLL_Summer12_ABC_TrueNumInteractions.root","pileup");
  //     new TPileupReweighting("/afs/cern.ch/work/a/amartell/public/weights/PUweights_DYJetsToLL_Summer12_Prompt_TrueNumInteractions.root","hweights");

  //2011                                                                                                                                              
  if(year == 2011) puReweighting =
    new TPileupReweighting("/afs/cern.ch/work/a/amartell/public/weights/PUweights_2011_DYJetsToLL_Fall2011_TrueNumInteractions.root", "hweights");

  
  std::string R9MOD = std::string(LOWHIGH);
  std::string ENERGY = std::string(ENE);
  
  //-------------------
  // Define in/outfiles
  std::string folderName;
  if(PU == 0)
  	folderName = std::string(EBEE)+"_"+std::string(LOWHIGH)+"_"+std::string(ENE)+"_noPU";
  if(PU == 1)
        folderName = std::string(EBEE)+"_"+std::string(LOWHIGH)+"_"+std::string(ENE);
  
  // Get trees
  std::cout << std::endl;

  std::string nameNtuples = "simpleNtupleEoverP/SimpleNtupleEoverP";
  std::string nameNtuplesMC = "simpleNtupleEoverP/SimpleNtupleEoverP";
  //  if(year == 2011) nameNtuples = "ntu";
  if(year == 2012) nameNtuplesMC = "simpleNtupleEoverPSh/SimpleNtupleEoverP";
  TChain* ntu_MC = new TChain(nameNtuplesMC.c_str());
  TChain* ntu_DA = new TChain(nameNtuples.c_str());

  if(year == 2012){
    ntu_MC->Add("/tmp/amartell/DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM.root");
    ntu_MC->Add("/tmp/amartell/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_2.root");
    ntu_DA->Add("/tmp/amartell/DoubleElectronAB_13Jul2012.root");
    ntu_DA->Add("/tmp/amartell/DoubleElectron_C_Prompt.root");
    //     ntu_MC->Add("/tmp/amartell/WJetsToLNu_START53_V7A.root");                                                                       
    //     ntu_DA->Add("/tmp/amartell/Single_AB_Prompt.root");                                                                             
    //     ntu_DA->Add("/tmp/amartell/Single_C_Prompt.root");                                                                              
  }
  if(year == 2011){
    ntu_DA->Add("/tmp/amartell/DoubleElectron-RUN2011AB.root");
    ntu_MC->Add("/tmp/amartell/DYJetsToLL_Fall11_START44_V9B.root");
  }

  std::cout << "     REFERENCE: " << std::setw(8) << ntu_MC->GetEntries() << " entries" << std::endl;
  std::cout << "     DATA: " << std::setw(8) << ntu_DA->GetEntries() << " entries" << std::endl;
  
  if(ntu_DA->GetEntries() == 0 || ntu_MC->GetEntries() == 0 )
  {
    std::cout << "Error: At least one file is empty" << std::endl; 
    return -1;
  }

  std::vector<int> run_DA, time_DA, Z_DA, PV_DA;
  std::vector<int> run_MC, time_MC, Z_MC, PV_MC;
  std::vector<float> scEt_reg1_DA, R91_DA, P1_DA, EoP1_DA,  scEta1_DA, isEB1_DA, e3x31_DA, scERaw1_DA;
  std::vector<float> scEt_reg2_DA, R92_DA, P2_DA, EoP2_DA,  scEta2_DA, isEB2_DA, e3x32_DA, scERaw2_DA;
  std::vector<float> scEt_reg1_MC, R91_MC, P1_MC, EoP1_MC,  scEta1_MC, isEB1_MC, e3x31_MC, scERaw1_MC;
  std::vector<float> scEt_reg2_MC, R92_MC, P2_MC, EoP2_MC,  scEta2_MC, isEB2_MC, e3x32_MC, scERaw2_MC;
  std::vector<float> Mass_DA, Mass_MC, Ht_DA, Ht_MC, puRe;

  // Set branch addresses
  float ele1ele2_scM,ele1ele2_scM_regression;
  int isZ,runId,timeStamp,nVtx;
  float npu;
  
  ntu_DA->SetBranchStatus("*",0);                         
  ntu_DA->SetBranchStatus("runId",1);                        ntu_DA->SetBranchAddress("runId", &runId);  
  ntu_DA->SetBranchStatus("timeStampHigh",1);                ntu_DA->SetBranchAddress("timeStampHigh", &timeStamp);  
  ntu_DA->SetBranchStatus("isZ",1);                          ntu_DA->SetBranchAddress("isZ", &isZ);
  ntu_DA->SetBranchStatus("PV_n",1);                         ntu_DA->SetBranchAddress("PV_n",&nVtx);

  ntu_MC->SetBranchStatus("*",0);
  ntu_MC->SetBranchStatus("PUit_TrueNumInteractions", 1);      ntu_MC->SetBranchAddress("PUit_TrueNumInteractions", &npu);
  ntu_MC->SetBranchStatus("runId",1);                          ntu_MC->SetBranchAddress("runId", &runId);
  ntu_MC->SetBranchStatus("timeStampHigh",1);                  ntu_MC->SetBranchAddress("timeStampHigh", &timeStamp);
  ntu_MC->SetBranchStatus("isZ",1);                            ntu_MC->SetBranchAddress("isZ", &isZ);
  ntu_MC->SetBranchStatus("PV_n",1);                           ntu_MC->SetBranchAddress("PV_n",&nVtx);


  ntu_DA->SetBranchStatus("ele1ele2_scM", 1);                ntu_DA->SetBranchAddress("ele1ele2_scM", &ele1ele2_scM);
  ntu_DA->SetBranchStatus("ele1ele2_scM_regression", 1);     ntu_DA->SetBranchAddress("ele1ele2_scM_regression", &ele1ele2_scM_regression);

  ntu_MC->SetBranchStatus("ele1ele2_scM", 1);                 ntu_MC->SetBranchAddress("ele1ele2_scM", &ele1ele2_scM);
  ntu_MC->SetBranchStatus("ele1ele2_scM_regression", 1);      ntu_MC->SetBranchAddress("ele1ele2_scM_regression", &ele1ele2_scM_regression);

  // Electron data
  float scEne1, scEneReg1, scEt1, scEta1, P1, scERaw1, e3x31;
  float scEne2, scEneReg2, scEt2, scEta2, P2, scERaw2, e3x32;
  int isEB1,isEB2;
  int ele1_charge, ele2_charge;
 
  ntu_DA->SetBranchStatus("ele1_scE", 1);       ntu_DA->SetBranchAddress("ele1_scE", &scEne1);
  ntu_DA->SetBranchStatus("ele1_scEt", 1);      ntu_DA->SetBranchAddress("ele1_scEt", &scEt1);
  ntu_DA->SetBranchStatus("ele1_scEta", 1);     ntu_DA->SetBranchAddress("ele1_scEta", &scEta1);
  ntu_DA->SetBranchStatus("ele1_scERaw",1);      ntu_DA->SetBranchAddress("ele1_scERaw",&scERaw1);

  ntu_DA->SetBranchStatus("ele2_scE", 1);        ntu_DA->SetBranchAddress("ele2_scE", &scEne2);
  ntu_DA->SetBranchStatus("ele2_scEt", 1);       ntu_DA->SetBranchAddress("ele2_scEt", &scEt2);
  ntu_DA->SetBranchStatus("ele2_scEta", 1);      ntu_DA->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_DA->SetBranchStatus("ele2_scERaw",1);      ntu_DA->SetBranchAddress("ele2_scERaw",&scERaw2);

  if(!UsePhotonRegression)  {
    ntu_DA->SetBranchStatus("ele1_scE_regression", 1);      ntu_DA->SetBranchAddress("ele1_scE_regression", &scEneReg1);
    ntu_DA->SetBranchStatus("ele2_scE_regression",1);       ntu_DA->SetBranchAddress("ele2_scE_regression", &scEneReg2);
  }
  else {
    ntu_DA->SetBranchStatus("ele1_scE_regression_PhotonTuned", 1);    ntu_DA->SetBranchAddress("ele1_scE_regression_PhotonTuned", &scEneReg1);
    ntu_DA->SetBranchStatus("ele2_scE_regression_PhotonTuned",1);     ntu_DA->SetBranchAddress("ele2_scE_regression_PhotonTuned", &scEneReg2);
  }

  ntu_DA->SetBranchStatus("ele1_e3x3",1);        ntu_DA->SetBranchAddress("ele1_e3x3", &e3x31);
  ntu_DA->SetBranchStatus("ele1_isEB",1);        ntu_DA->SetBranchAddress("ele1_isEB",&isEB1);
  ntu_DA->SetBranchStatus("ele1_tkP",1);         ntu_DA->SetBranchAddress("ele1_tkP", &P1);
  ntu_DA->SetBranchStatus("ele1_charge",1);      ntu_DA->SetBranchAddress("ele1_charge", &ele1_charge);

  ntu_DA->SetBranchStatus("ele2_e3x3",1);        ntu_DA->SetBranchAddress("ele2_e3x3", &e3x32);
  ntu_DA->SetBranchStatus("ele2_isEB",1);        ntu_DA->SetBranchAddress("ele2_isEB",&isEB2);
  ntu_DA->SetBranchStatus("ele2_tkP",1);         ntu_DA->SetBranchAddress("ele2_tkP", &P2);
  ntu_DA->SetBranchStatus("ele2_charge",1);      ntu_DA->SetBranchAddress("ele2_charge", &ele2_charge);
  ///////////////////////                                                                        

  ntu_MC->SetBranchStatus("ele1_scE", 1);       ntu_MC->SetBranchAddress("ele1_scE", &scEne1);
  ntu_MC->SetBranchStatus("ele1_scEt", 1);      ntu_MC->SetBranchAddress("ele1_scEt", &scEt1);
  ntu_MC->SetBranchStatus("ele1_scEta", 1);     ntu_MC->SetBranchAddress("ele1_scEta", &scEta1);
  ntu_MC->SetBranchStatus("ele1_scERaw",1);      ntu_MC->SetBranchAddress("ele1_scERaw",&scERaw1);

  ntu_MC->SetBranchStatus("ele2_scE", 1);       ntu_MC->SetBranchAddress("ele2_scE", &scEne2);
  ntu_MC->SetBranchStatus("ele2_scEt", 1);      ntu_MC->SetBranchAddress("ele2_scEt", &scEt2);
  ntu_MC->SetBranchStatus("ele2_scEta", 1);     ntu_MC->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_MC->SetBranchStatus("ele2_scERaw",1);      ntu_MC->SetBranchAddress("ele2_scERaw",&scERaw2);
  if(!UsePhotonRegression)  {
    ntu_MC->SetBranchStatus("ele1_scE_regression", 1);      ntu_MC->SetBranchAddress("ele1_scE_regression", &scEneReg1);
    ntu_MC->SetBranchStatus("ele2_scE_regression",1);       ntu_MC->SetBranchAddress("ele2_scE_regression", &scEneReg2);
  }
  else {
    ntu_MC->SetBranchStatus("ele1_scE_regression_PhotonTuned", 1);    ntu_MC->SetBranchAddress("ele1_scE_regression_PhotonTuned", &scEneReg1);
    ntu_MC->SetBranchStatus("ele2_scE_regression_PhotonTuned",1);     ntu_MC->SetBranchAddress("ele2_scE_regression_PhotonTuned", &scEneReg2);
  }

  ntu_MC->SetBranchStatus("ele1_e3x3",1);        ntu_MC->SetBranchAddress("ele1_e3x3", &e3x31);
  ntu_MC->SetBranchStatus("ele1_isEB",1);        ntu_MC->SetBranchAddress("ele1_isEB",&isEB1);
  ntu_MC->SetBranchStatus("ele1_tkP",1);         ntu_MC->SetBranchAddress("ele1_tkP", &P1);
  ntu_MC->SetBranchStatus("ele1_charge",1);      ntu_MC->SetBranchAddress("ele1_charge", &ele1_charge);

  ntu_MC->SetBranchStatus("ele2_e3x3",1);        ntu_MC->SetBranchAddress("ele2_e3x3", &e3x32);
  ntu_MC->SetBranchStatus("ele2_isEB",1);        ntu_MC->SetBranchAddress("ele2_isEB",&isEB2);
  ntu_MC->SetBranchStatus("ele2_tkP",1);         ntu_MC->SetBranchAddress("ele2_tkP", &P2);
  ntu_MC->SetBranchStatus("ele2_charge",1);      ntu_MC->SetBranchAddress("ele2_charge", &ele2_charge);
  //////////////////////

  //  std::cout << "  >>>> ntu_DA -> GetEntries() = " << ntu_DA -> GetEntries() << std::endl;  
  for(int ientry = 0; ientry < ntu_DA -> GetEntries(); ientry++)
  {
  	if( (ientry%100000 == 0) ) std::cout << "reading DATA entry " << ientry << "\r" << std::flush;
        ntu_DA->GetEntry(ientry);  

	if(isZ == 0) continue;

	if(useMomentum){
	  scEneReg1 = P1; scEneReg2 = P2;
	}

	if(scEneReg1/scEne1*scEt1 < 25. || scEneReg2/scEne2*scEt2 < 25.) continue;


	float corrEtR9_1 = 1.;
	float corrEtR9_2 = 1.;
	if(correctEt == true){
	  if(year == 2012){
	    if(e3x31/scERaw1 < 0.94 ) corrEtR9_1 = corrEtR9_1 / (1. + Et_lowR9_2012->Eval(scEneReg1));
	    if(e3x32/scERaw2 < 0.94 ) corrEtR9_2 = corrEtR9_2 / (1. + Et_lowR9_2012->Eval(scEneReg2));
	    if(e3x31/scERaw1 >= 0.94 ) corrEtR9_1 = corrEtR9_1 / (1. + Et_highR9_2012->Eval(scEneReg1));
	    if(e3x32/scERaw2 >= 0.94 ) corrEtR9_2 = corrEtR9_2 / (1. + Et_highR9_2012->Eval(scEneReg2));
	  }
	  if(year == 2011){
	    if(e3x31/scERaw1 < 0.94 ) corrEtR9_1 = corrEtR9_1 / (1. + Et_lowR9_2011->Eval(scEneReg1));
	    if(e3x32/scERaw2 < 0.94 ) corrEtR9_2 = corrEtR9_2 / (1. + Et_lowR9_2011->Eval(scEneReg2));
	    if(e3x31/scERaw1 >= 0.94 ) corrEtR9_1 = corrEtR9_1 / (1. + Et_highR9_2011->Eval(scEneReg1));
	    if(e3x32/scERaw2 >= 0.94 ) corrEtR9_2 = corrEtR9_2 / (1. + Et_highR9_2011->Eval(scEneReg2));
	  }
	}


	if(useShCorr == true){
	  corrEtR9_1 = corrEtR9_1 * GetShervingCorrections(scEta1, e3x31/scERaw1, runId);
	  corrEtR9_2 = corrEtR9_2 * GetShervingCorrections(scEta2, e3x32/scERaw2, runId);
	}

        run_DA.push_back(runId);

	Mass_DA.push_back(ele1ele2_scM * sqrt(scEneReg1/scEne1 * corrEtR9_1 * scEneReg2/scEne2 * corrEtR9_2));

        time_DA.push_back(timeStamp);
        Z_DA.push_back(isZ);
	PV_DA.push_back(nVtx);    

	scEta1_DA.push_back(scEta1);
	scEta2_DA.push_back(scEta2);

	isEB1_DA.push_back(isEB1);
	isEB2_DA.push_back(isEB2);

	R91_DA.push_back(e3x31/scERaw1);
	R92_DA.push_back(e3x32/scERaw2);

	Ht_DA.push_back(scEneReg1/scEne1*scEt1* corrEtR9_1 + scEneReg2/scEne2*scEt2*corrEtR9_2);

	e3x31_DA.push_back(e3x31);
	e3x32_DA.push_back(e3x32);
  }
  
  std::cout << std::endl;
  float ww = 1.;
  for(int ientry = 0; ientry < ntu_MC -> GetEntries(); ientry++)
  {
  	if( (ientry%100000 == 0) ) std::cout << "reading MC entry " << ientry << "\r" << std::flush;
        ntu_MC->GetEntry(ientry);  
  
	if(isZ == 0) continue;

	if(useMomentum){
	  scEneReg1 = P1; scEneReg2 = P2;
	}

	float R9_ele1 = e3x31/scERaw1;
 	if(year == 2012 && isEB1 == 1) R9_ele1 = 0.0010 + 1.0045 * e3x31/scERaw1;
 	if(year == 2012 && isEB1 == 0) R9_ele1 = -0.0007 + 1.0086 * e3x31/scERaw1;
 	if(year == 2011) R9_ele1 = 1.0035 * e3x31/scERaw1;
	float R9_ele2 = e3x32/scERaw2;
 	if(year == 2012 && isEB2 == 1) R9_ele2 = 0.0010 + 1.0045 * e3x32/scERaw2;
 	if(year == 2012 && isEB2 == 0) R9_ele2 = -0.0007 + 1.0086 * e3x32/scERaw2;
 	if(year == 2011) R9_ele2 = 1.0035 * e3x32/scERaw2;

	float energySmearing1 = 1.;
	float energySmearing2 = 1.;
	energySmearing1 = gRandom->Gaus(1., GetSmearings(scEta1, R9_ele1, year, isEB1));
	energySmearing2 = gRandom->Gaus(1., GetSmearings(scEta2, R9_ele2, year, isEB2));

	if(useMomentum){
	  energySmearing1 = 1.; energySmearing2 = 1.;	}


	if(scEneReg1/scEne1*scEt1*energySmearing1 < 25. || scEneReg2/scEne2*scEt2*energySmearing2 < 25.) continue;


	ww = puReweighting->GetWeight((int)npu);
        puRe.push_back(ww);
        run_MC.push_back(runId);
        time_MC.push_back(timeStamp);
        Z_MC.push_back(isZ);
	PV_MC.push_back(nVtx);    

	Mass_MC.push_back(ele1ele2_scM * sqrt(scEneReg1/scEne1*energySmearing1 * scEneReg2/scEne2*energySmearing2) );

	P1_MC.push_back(P1);
	P2_MC.push_back(P2);

	Ht_MC.push_back(scEneReg1/scEne1*scEt1*energySmearing1 + scEneReg2/scEne2*scEt2*energySmearing2);
	
	scEta1_MC.push_back(scEta1);
	scEta2_MC.push_back(scEta2);

	isEB1_MC.push_back(isEB1);
	isEB2_MC.push_back(isEB2);

	R91_MC.push_back(R9_ele1);
	R92_MC.push_back(R9_ele2);

 	e3x31_MC.push_back(R9_ele1*scERaw1);
 	e3x32_MC.push_back(R9_ele2*scERaw2);
  }

  
  // Loop and sort events
  std::cout << std::endl;
  std::cout << "***** Sort events and define bins *****" << std::endl;
  
  int nEntries = Mass_DA.size();
  std::cout << " >>>>>> nEntries = " << nEntries << std::endl;
  int nSavePts = 0;
  std::vector<bool> isSavedEntries(nEntries);
  std::vector<SorterLC> sortedEntries;

  TH1F* etaDistrib = new TH1F("etaDistrib", "", 500, -3., 3);
  TH1F* etaDistribEB = new TH1F("etaDistribEB", "", 500, -3., 3);
  TH1F* etaDistribEBcrack = new TH1F("etaDistribEBcrack", "", 500, -3., 3);
  
  for(int ientry = 0; ientry < nEntries; ++ientry)
  {
    isSavedEntries.at(ientry) = false;

    //distribuzioni forse inutili da togliere?
    if(R91_DA.at(ientry) >= 0.94 && R92_DA.at(ientry) >= 0.94) {
      etaDistrib->Fill(scEta2_DA.at(ientry));
      etaDistrib->Fill(scEta1_DA.at(ientry));
    }

    if(R91_DA.at(ientry) >= 0.94 && R92_DA.at(ientry) >= 0.94) {
      if(fabs(scEta1_DA.at(ientry)) > 1.4442) etaDistribEBcrack->Fill(scEta2_DA.at(ientry));
      if(fabs(scEta2_DA.at(ientry)) > 1.4442) etaDistribEBcrack->Fill(scEta1_DA.at(ientry));
    }

    if(R91_DA.at(ientry) >= 0.94 && R92_DA.at(ientry) >= 0.94) {
      if(fabs(scEta1_DA.at(ientry)) < 1.4442) etaDistribEB->Fill(scEta2_DA.at(ientry));
      if(fabs(scEta2_DA.at(ientry)) < 1.4442) etaDistribEB->Fill(scEta1_DA.at(ientry));
    }

    if( fabs(scEta1_DA.at(ientry)) > 2.5 ||fabs(scEta2_DA.at(ientry)) > 2.5) continue; 

    if(std::string(EBEE) == "EBEB" && (fabs(scEta1_DA.at(ientry)) > 1.4442 )) continue; 
    if(std::string(EBEE) == "EBEB" && (fabs(scEta2_DA.at(ientry)) > 1.4442 )) continue; 
    if(std::string(EBEE) == "notEBEB" && (isEB1_DA.at(ientry) == 1. && isEB2_DA.at(ientry) == 1.) ) continue; 

    if(std::string(LOWHIGH) == "HH" && (R91_DA.at(ientry) < 0.94 || R92_DA.at(ientry) < 0.94) ) continue;
    if(std::string(LOWHIGH) == "notHH" && R91_DA.at(ientry) > 0.94 && R92_DA.at(ientry) > 0.94) continue;

    //      if (R91_DA.at(ientry) < 0.7 || R92_DA.at(ientry) < 0.7 ) continue;

    isSavedEntries.at(ientry) = true;
    
    SorterLC dummy;
    dummy.laserCorr = Ht_DA.at(ientry);
    dummy.entry = ientry;
    sortedEntries.push_back(dummy);
    nSavePts++;   
  }
  std::cout << " this category:: Effective entries = " << nSavePts << std::endl;
  std::cout << " Effective entries sortedEntries.size()= " << sortedEntries.size() << std::endl;


  std::sort(sortedEntries.begin(),sortedEntries.end(),SorterLC());
  std::cout << "DATA sorted in " << EBEE << " - " << nSavePts << " events" << std::endl;


   std::map<int,int> antiMap;
   for(unsigned int iSaved = 0; iSaved < sortedEntries.size(); ++iSaved){
   antiMap[sortedEntries.at(iSaved).entry] = iSaved; 
   //riga sopra mette in antiMap[indice ientry] il numero di entry dopo il sorting   
   //   std::cout << " sortedEntries.at("<< iSaved<<").laserCorr = " << sortedEntries.at(iSaved).laserCorr << std::endl;

   }


  // bins with evtsPerPoint events per bin
   std::cout << " nSavePts = " << nSavePts << std::endl;
   std::cout << " evtsPerPoint = " << evtsPerPoint << std::endl;

   //nBins = numero di eventi tot/eventi per bin ======> quindi numero di bins
  int nBins = std::max(1, int(nSavePts/evtsPerPoint));
   std::cout << " nBins = " << nBins << std::endl;

   //nBinPts = numero di eventi per bin ma calcolato a posteriori
  int nBinPts = int( nSavePts/nBins );
   std::cout << " nBinPts = " << nBinPts << std::endl;
  int nBinTempPts = 0;

  std::cout << "nBins = " << nBins << std::endl;
  
  std::cout << " sortedEntries.at(iSaved).laserCorr = " << sortedEntries.at(0).laserCorr << std::endl;

  std::vector<int> binEntryMax;
  binEntryMax.push_back(0);
  for(int iSaved = 0; iSaved < nSavePts; ++iSaved)
  {
    ++nBinTempPts;
    
    if( nBinTempPts == nBinPts )
    {
      //salvo quanti eventi per bin
      binEntryMax.push_back( iSaved );      
      nBinTempPts = 0;
//        std::cout << "binEntryMax.size() = " << binEntryMax.size() << std::endl;
//        std::cout << " sortedEntries.at(iSaved).laserCorr = " << sortedEntries.at(iSaved).laserCorr << std::endl;
//        std::cout << " sortedEntries.at(iSaved).laserCorr = " << sortedEntries.at(iSaved+1).laserCorr << std::endl;
    }
  }
  binEntryMax.at(nBins) = nSavePts;

  std::cout << " fine : nBins = " << nBins << std::endl;
 

  TVirtualFitter::SetDefaultFitter("Fumili2");
  
  // histogram definition
  
  TH1F** h_EoP_DA = new TH1F*[nBins];
  TH1F** h_EoP_MC = new TH1F*[nBins];
  TH1F** h_Ht = new TH1F*[nBins];
  TH1F** h_Ht_MC = new TH1F*[nBins];
  TH1F* h_Ht_allDA = new TH1F("h_Ht_allDA", "", 5000, 0., 1000.);
  TH1F* h_Ht_allMC = new TH1F("h_Ht_allMC", "", 5000, 0., 1000.);

  TH1F* h_Mass_allDA = new TH1F("h_Mass_allDA", "", 400, 0.5, 1.5);
  TH1F* h_Mass_allMC = new TH1F("h_Mass_allMC", "", 400, 0.5, 1.5);

  TH1F* h_Vtx_DA = new TH1F("h_Vtx_DA", "", 200, 0., 200.);
  TH1F* h_Vtx_MC = new TH1F("h_Vtx_MC", "", 200, 0., 200.);


  h_Ht_allDA->Sumw2(); 
  h_Ht_allMC->Sumw2(); 
  h_Mass_allDA->Sumw2();
  h_Mass_allMC->Sumw2();
  h_Vtx_DA->Sumw2(); 
  h_Vtx_MC->Sumw2(); 


  h_Mass_allDA->SetLineColor(kRed+2);
  h_Ht_allDA->SetLineColor(kRed+2);
  h_Vtx_DA->SetLineColor(kRed+2);

  h_Mass_allMC->SetLineColor(kGreen+2);
  h_Ht_allMC->SetLineColor(kGreen+2);
  h_Vtx_MC->SetLineColor(kGreen+2);

  std::vector<float> EtBinEdge;
  std::vector<float> xNorm_single;

  for(int i = 0; i < nBins; ++i)
  {
    char histoName[80];
    
    sprintf(histoName, "EoP_DA_%d", i);
    h_EoP_DA[i] = new TH1F(histoName, histoName, 600, 0.5, 1.5);
    h_EoP_DA[i] -> SetFillColor(kRed+2);
    h_EoP_DA[i] -> SetFillStyle(3004);
    h_EoP_DA[i] -> SetMarkerStyle(7);
    h_EoP_DA[i] -> SetMarkerColor(kRed+2); 
    h_EoP_DA[i] -> SetLineColor(kRed+2); 
    
    sprintf(histoName, "EoP_MC_%d", i);
    h_EoP_MC[i] = new TH1F(histoName, histoName, 600, 0.5, 1.5);
    h_EoP_MC[i] -> SetFillColor(kGreen+2);
    h_EoP_MC[i] -> SetFillStyle(3004);
    h_EoP_MC[i] -> SetMarkerStyle(7);
    h_EoP_MC[i] -> SetMarkerColor(kGreen+2);
    h_EoP_MC[i] -> SetLineColor(kGreen+2);
    
    sprintf(histoName, "Et_%d", i);
    h_Ht[i] = new TH1F(histoName, histoName, 5000, 0, 1000);
    h_Ht[i]->SetLineColor(kRed+2);

    sprintf(histoName, "Et_MC_%d", i);
    h_Ht_MC[i] = new TH1F(histoName, histoName, 5000, 0, 1000);
    h_Ht_MC[i]->SetLineColor(kGreen+2);

    h_EoP_DA[i]->Sumw2();
    h_EoP_MC[i]->Sumw2();

    h_Ht[i]->Sumw2();
    h_Ht_MC[i]->Sumw2();
  }
  
  std::cout << " Ht_MC.size() = " << Ht_MC.size() << std::endl;
  std::cout << " Ht_DA.size() = " << Ht_DA.size() << std::endl;
  std::cout << " Mass_MC.size() = " << Mass_MC.size() << std::endl;
  std::cout << " Mass_DA.size() = " << Mass_DA.size() << std::endl;


  std::vector<float> x;
  std::vector<float> ex;
  std::vector<float> y;
  std::vector<float> ey;
  
  TGraphErrors* finalGraph = new TGraphErrors();
  
  // function definition
  TF1** f_EoP = new TF1*[nBins];
 
  // loop on the saved and sorted events
  std::cout << std::endl;
  std::cout << "***** Fill and fit histograms *****" << std::endl;

  int DAEntries = Mass_DA.size();

  for(unsigned int ientry = 0; ientry < DAEntries; ++ientry)
    {
    if( (ientry%100000 == 0) ) std::cout << "reading entry " << ientry << std::endl;
      
    if( isSavedEntries.at(ientry) == false) continue;
     
    //iSaved e' l'indice della entry dopo il sorting
     int iSaved = antiMap[ientry];
     int bin = -1;
    
    for(bin = 0; bin < nBins; ++bin)
      if( iSaved >= binEntryMax.at(bin) && iSaved < binEntryMax.at(bin+1) )
    break;

    
    h_EoP_DA[bin]->Fill( Mass_DA.at(ientry)/91.18 );
    h_Ht[bin]->Fill(Ht_DA.at(ientry) );
    h_Ht_allDA->Fill(Ht_DA.at(ientry) );
    h_Mass_allDA->Fill( Mass_DA.at(ientry)/91.18 );

    h_Vtx_DA->Fill(PV_DA.at(int(ientry)) );
    }
  
  std::cout << std::endl;
  
  for(int bin = 0; bin < nBins; bin++)
  {
    std::cout << "h_Ht[bin]->GetEntries() =  " << h_Ht[bin]->GetEntries() << std::endl;
    std::cout << "h_EoP_DA[bin]->GetEntries() =  " << h_EoP_DA[bin]->GetEntries() << std::endl;
    for(int i = 1; i < h_Ht[bin]->GetNbinsX()+1; i++)
      {
	if(h_Ht[bin]->GetBinContent(i) > 0) {
	  EtBinEdge.push_back(h_Ht[bin]->GetBinCenter(i)-h_Ht[bin]->GetBinWidth(i) );
// 	  std::cout << " >>> i = " << EtBinEdge.size() - 1 << std::endl;
// 	  std::cout << " >>> EtBinEdge.at(i) = " << EtBinEdge.at(EtBinEdge.size()-1) << std::endl;
	  break;
	}  
      }
  }

  int MCEntries = Mass_MC.size();

  for(unsigned int ientry = 0; ientry < MCEntries; ++ientry)
    {   
      if( (ientry%100000 == 0) ) std::cout << "reading entry " << ientry << std::endl;

      if( fabs(scEta1_MC.at(ientry)) > 2.5 ||fabs(scEta2_MC.at(ientry)) > 2.5) continue;

      if(std::string(EBEE) == "EBEB" && (fabs(scEta1_MC.at(ientry)) > 1.4442 )) continue;
      if(std::string(EBEE) == "EBEB" && (fabs(scEta2_MC.at(ientry)) > 1.4442 )) continue;
      if(std::string(EBEE) == "notEBEB" && (isEB1_MC.at(ientry) == 1. && isEB2_MC.at(ientry) == 1.) ) continue;
      if(std::string(LOWHIGH) == "HH" && (R91_MC.at(ientry) < 0.94 || R92_MC.at(ientry) < 0.94) ) continue;
      if(std::string(LOWHIGH) == "notHH" && R91_MC.at(ientry) > 0.94 && R92_MC.at(ientry) > 0.94 ) continue;

      //      if (R91_MC.at(ientry) < 0.7 || R92_MC.at(ientry) < 0.7 ) continue;

      for(unsigned int bin = 0; bin < EtBinEdge.size(); ++bin)
	{	  
	  float referenceEt_MC = Ht_MC.at(ientry);
	  if( (bin != EtBinEdge.size()-1 && referenceEt_MC > EtBinEdge.at(bin) && referenceEt_MC < EtBinEdge.at(bin+1)) ||
	      (bin == EtBinEdge.size()-1 && referenceEt_MC > EtBinEdge.at(bin) ) )
	    {
		  h_EoP_MC[(int)bin]->Fill( (Mass_MC.at(ientry))/91.18, puRe.at(ientry));
		  h_Ht_MC[int(bin)]->Fill(referenceEt_MC, puRe.at(ientry) );
		  break;
	    }
	}
      
      h_Ht_allMC->Fill(Ht_MC.at(ientry), puRe.at(ientry));
      h_Vtx_MC->Fill(PV_MC.at(ientry), puRe.at(ientry) );

      h_Mass_allMC->Fill(Mass_MC.at(ientry)/91.18, puRe.at(ientry));
    }
  

  std::cout << " fino a qui ci sono " <<  std::endl;

  for(int i = 0; i < nBins; ++i)
  {
    //------------------------------------
    // Fill the graph for uncorrected data
    // define the fitting function
    // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )

    float xNorm_all = h_Ht_allDA->Integral()/h_Ht_allMC->Integral(); //* h_Ht_allDA->GetBinWidth()/h_Ht_allMC->GetBinWidth();
    h_Ht_allMC->Scale(xNorm_all);
    h_Vtx_MC->Scale(xNorm_all);

    h_EoP_DA[i]->Rebin(2);
    h_EoP_MC[i]->Rebin(2);
    
    h_EoP_MC[i]->Smooth(6);

    float xNorm = h_EoP_DA[i]->Integral()/h_EoP_MC[i]->Integral()*h_EoP_DA[i]->GetBinWidth(1)/h_EoP_MC[i]->GetBinWidth(1);  
    float xNormEt = h_Ht[i]->Integral()/h_Ht_MC[i]->Integral(); //*h_Ht[i]->GetBinWidth()/h_Ht_MC[i]->GetBinWidth();  
    h_EoP_MC[i]->Scale(xNorm);
    h_Ht_MC[i]->Scale(xNormEt);

//     std::cout << " i = " << i << " h_EoP_DA[i]->Integral() = " << h_EoP_DA[i]->Integral() << std::endl;
//     std::cout << " i = " << i << " h_Ht[i]->Integral() = " << h_Ht[i]->Integral() << std::endl;
//     std::cout << " i = " << i << " h_EoP_MC[i]->Integral() = " << h_EoP_MC[i]->Integral() << std::endl;
//     std::cout << " i = " << i << " h_Ht_MC[i]->Integral() = " << h_Ht_MC[i]->Integral() << std::endl;




//     std::cout << " i = " << i << " h_EoP_DA[i]->Integral() = " << h_EoP_DA[i]->Integral() << std::endl;
//     std::cout << " i = " << i << " h_EoP_MC[i]->Integral() = " << h_EoP_MC[i]->Integral() << std::endl;
//     std::cout << " xNorm = " << xNorm << std::endl;
//     std::cout << " xNormEt = " << xNormEt << std::endl;


    
    //    xNorm_single.push_back(xNormEt);
//     if(h_R9[i]->GetMean() < 0.94 && h_R9_MC[i]->GetMean() < 0.94) {
//     h_EoP_DA[i]->Rebin(2);
//     h_EoP_MC[i]->Rebin(2);
//     }

//     float mid = h_EoP_DA[i]->GetBinCenter(h_EoP_DA[i] -> GetMaximumBin());
//     float x_min = mid*0.65;
//     float x_max = mid*1.7;


//     float x_min = 1.;
//     float x_max = 1.;
//     for(int nBinsHisto = 1; nBinsHisto < h_EoP_DA[i]->GetNbinsX(); ++nBinsHisto){
//       if(h_EoP_DA[i]->GetBinContent(nBinsHisto) > 10. && x_min == 1.) x_min = h_EoP_DA[i]->GetBinCenter(nBinsHisto);
//       if(h_EoP_DA[i]->GetBinContent(nBinsHisto) < 10. && x_min != 1. && x_max == 1.) x_max = h_EoP_DA[i]->GetBinCenter(nBinsHisto);
//     }


    histoFunc* templateHistoFunc = new histoFunc(h_EoP_MC[i]);
    char funcName[50];
    sprintf(funcName,"f_EoP_%d",i);
    //    f_EoP[i] = new TF1(funcName, templateHistoFunc, x_min, x_max, 3, "histoFunc");
    f_EoP[i] = new TF1(funcName, templateHistoFunc, 0.7, 1.3, 3, "histoFunc");
    if(std::string(EBEE) != "EBEB")    f_EoP[i] = new TF1(funcName, templateHistoFunc, 0.7, 1.2, 3, "histoFunc");

//     histoFunc* templateHistoFunc = new histoFunc(h_EoP_MC[i]);
//     char funcName[50];
//     sprintf(funcName,"f_EoP_%d",i);
//    f_EoP[i] = new TF1(funcName, templateHistoFunc, 0.8, 1.3, 3, "histoFunc");
    
    f_EoP[i] -> SetParName(0,"Norm"); 
    f_EoP[i] -> SetParName(1,"Scale factor"); 
    f_EoP[i] -> SetLineWidth(1); 
    f_EoP[i] -> SetNpx(10000);

    xNorm = 1.;    

    f_EoP[i] -> FixParameter(0, xNorm);
    //    f_EoP[i] -> SetParameter(1, gRandom->Gaus(1.,0.005));
    f_EoP[i] -> SetParameter(1, 0.99);
    f_EoP[i] -> FixParameter(2, 0.);
    f_EoP[i] -> SetLineColor(kRed+2); 

    TFitResultPtr rp = h_EoP_DA[i] -> Fit(funcName, "QERLS+");
    int fStatus = rp;
    int nTrials = 0;
    while( (fStatus != 0) && (nTrials < 100) )
    {
      rp = h_EoP_DA[i] -> Fit(funcName, "QERLS+");
      fStatus = rp;
      if(fStatus == 0) break;
      ++nTrials;
    }

    double eee = f_EoP[i]->GetParError(1); 
    double k = 1./f_EoP[i]->GetParameter(1);


    // Fill the graph      
    if (fStatus == 0 && eee*k > 0.1*h_EoP_DA[i]->GetRMS()/sqrt(evtsPerPoint))
    {
      x.push_back(h_Ht[i]->GetMean());
      ex.push_back((h_Ht[i]->GetRMS())/sqrt(h_Ht[i]->GetEntries()));
      y.push_back(k-1);
      ey.push_back(eee * k * k);
      
    }
    else
    std::cout << "Fitting uncorrected Et bin: " << i << "   Fail status: " << fStatus << "   sigma: " << eee << std::endl;

  }
   
  for(unsigned int i = 0; i < x.size(); ++i)
  {
    finalGraph->SetPoint(i,  x.at(i) , y.at(i));
    finalGraph->SetPointError(i, ex.at(i), ey.at(i));
  }

  if(year == 2012) finalGraph->SetMarkerColor(kBlue);
  if(year == 2011) finalGraph->SetMarkerColor(kCyan);
  if(strcmp(LOWHIGH,"HIGH")==0 )  finalGraph->GetYaxis()->SetRangeUser(-0.004, 0.014);
  if(strcmp(LOWHIGH,"LOW")==0 )  finalGraph->GetYaxis()->SetRangeUser(-0.03, 0.03);
  finalGraph->GetYaxis()->SetTitle("(M_{ee}/M_{Z})^{data} - (M_{ee}/M_{Z})^{mc}");
  finalGraph->GetXaxis()->SetRangeUser(0., 130.);
  finalGraph->GetXaxis()->SetTitle("Ht");


  std::string plotFolderName = "PLOTS_MZ";
  if(doVsEach == "true") plotFolderName = "PLOTS_MZ"; 

  TFile pippo((plotFolderName+"/results_"+folderName+"_"+string_year+".root").c_str(),"recreate");
  finalGraph->Write("finalGraph");
  h_Ht_allMC->Write();
  h_Ht_allDA->Write();
  h_Vtx_DA->Write();
  h_Vtx_MC->Write();
  etaDistrib->Write();
  etaDistribEB->Write();
  etaDistribEBcrack->Write();
  h_Mass_allDA->Write();
  h_Mass_allMC->Write();

  for(int i = 0; i < nBins; ++i){
    h_EoP_DA[i]->Write();
    h_EoP_MC[i]->Write();
    h_Ht[i]->Write();
    h_Ht_MC[i]->Write();
  }
  pippo.Close();
  

//   // Drawings
//   TPaveStats** s_EoP = new TPaveStats*[nBins];
  
//   TCanvas *c1[100]; 
//   for(int i = 0; i < nBins; ++i)
//   {    
//     char canvasName[50];
//     sprintf(canvasName, "Fits-%0d", i); 
//     c1[i] = new TCanvas(canvasName, canvasName);
//     c1[i]->cd();
//     h_EoP_DA[i] -> GetXaxis() -> SetTitle("M_{ee}/M_{Z}");  
//     h_EoP_DA[i] -> GetYaxis() -> SetRangeUser(0., std::max(h_EoP_DA[i]->GetMaximum(), h_EoP_MC[i]->GetMaximum()) + 10.); 
//     h_EoP_DA[i] -> GetXaxis() -> SetRangeUser(0.5,1.5); 
//     //    h_EoP_DA[i] -> Draw("e");
//     h_EoP_DA[i] -> Draw();
//     gPad->Update();
//     s_EoP[i]= (TPaveStats*)(h_EoP_DA[i]->GetListOfFunctions()->FindObject("stats"));
//     s_EoP[i]->SetTextColor(kRed+2);
//     f_EoP[i]->Draw("same");
    
//     h_EoP_MC[i] -> Draw("same");
    
//     char Name[100];
//     if(PU == 0)      sprintf(Name, (plotFolderName+"/"+folderName+"/noPU_fit_%d_"+string_year+".png").c_str(),i);
//     if(PU == 1)      sprintf(Name, (plotFolderName+"/"+folderName+"/fit_%d_"+string_year+".png").c_str(),i);
//     c1[i] -> Print(Name,".png");     
//   }
  
//   TCanvas *c2[100]; 
//   for(int i = 0; i < nBins; ++i)
//   {
//     char canvasName[50];
//     sprintf(canvasName, "Ht_DA-%0d", i); 
//     c2[i] = new TCanvas(canvasName, canvasName);
//     c2[i]->cd();

//     h_Ht[i]->GetXaxis() -> SetTitle("Ht");
//     h_Ht[i]->GetYaxis()->SetRangeUser(0, std::max(h_Ht[i]->GetMaximum(), h_Ht_MC[i]->GetMaximum()) + 10. ); 
//     if(i<nBins-1) h_Ht[i]->GetXaxis()->SetRangeUser(EtBinEdge.at(i), EtBinEdge.at(i+1)); 
//     else h_Ht[i]->GetXaxis()->SetRangeUser(EtBinEdge.at(i), 150.); 
//     //    h_Ht[i] -> Draw("e");
//     h_Ht[i]->Draw();
//     h_Ht_MC[i]->Draw("same");
//     /*gPad->Update();
//     s_Las[i]= (TPaveStats*)(h_Ht[i]->GetListOfFunctions()->FindObject("stats"));
//     s_Las[i]->SetTextColor(kBlack);*/
    
//     char Name[100];
//     if(PU == 0)      sprintf(Name, (plotFolderName+"/"+folderName+"/noPU_Ht_%d_"+string_year+".png").c_str(),i);
//     if(PU == 1)      sprintf(Name, (plotFolderName+"/"+folderName+"/Ht_%d_"+string_year+".png").c_str(),i);
//     c2[i]->Print(Name,".png");      
//   }



//   TCanvas* Et_spectrum = new TCanvas;
//   gPad->SetLogy();
//   //  h_Ht_allDA->GetYaxis()->SetRangeUser(0.1, 10000.);
//   h_Ht_allDA->GetXaxis()->SetRangeUser(0., 150.);
//   h_Ht_allDA->GetXaxis()->SetTitle("Ht ");
//   h_Ht_allDA->SetMarkerColor(kRed+2);
//   h_Ht_allDA->SetMarkerStyle(7);
//   h_Ht_allDA->Draw("e");
//    for(int jj = 0; jj < nBins; ++jj){
//      h_Ht_MC[jj]->GetXaxis()->SetRangeUser(0., 150.);
//      h_Ht_MC[jj]->Draw("same");
//    }
//   TLegend *tspec = new TLegend(0.64,0.80,0.99,0.99);
//   tspec->SetFillColor(0);
//   tspec->SetTextFont(42); 
//   tspec->AddEntry(h_Ht_allDA,"DATA","PL");
//   tspec->AddEntry(h_Ht_MC[0],"MC ","PL");
//   tspec->Draw(); 
//   Et_spectrum->Print((plotFolderName+"/"+folderName+"/Ht_spectrum_"+string_year+".png").c_str(), ".png");


//   TCanvas* cVtx = new TCanvas();
//   h_Vtx_DA->Draw();
//   h_Vtx_MC->SetLineColor(kGreen+2);
//   h_Vtx_MC->Draw("same");
//   cVtx->Print((plotFolderName+"/"+folderName+"/Vtx_"+string_year+".png").c_str(),".png");


//   std::sort(y.begin(), y.end());
//   std::sort(ey.begin(), ey.end());
//   TCanvas* cplot = new TCanvas("gplot", "gplot",100,100,725,500);
//   cplot->cd();

//   std::cout << " sortato range " << std::endl;

//   TPad *cLeft  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
//   cLeft->SetLeftMargin(0.17); 
//   cLeft->SetRightMargin(0.025); 
//   cLeft->SetBottomMargin(0.17); 
  
//   cLeft->Draw();
 
//   float tYoffset = 1.75; 
//   float tXoffset = 1.6; 
//   float labSize = 0.04;
//   float labSize2 = 0.07;

//   cLeft->cd(); 
//   cLeft->SetGridx();
//   cLeft->SetGridy();
  
//   float x_min = x.at(0)-ex.at(ex.size()-1)-10;
//   float x_max = x.at(x.size()-1)+ex.at(ex.size()-1)+10;
// //   float y_min = y.at(0)-ey.at(ey.size()-1)-0.002;
// //   float y_max = y.at(y.size()-1)+ey.at(ey.size()-1)+0.002;

// //   float y_min = y.at(0)-ey.at(ey.size()-1)-0.005;
// //   float y_max = y.at(y.size()-1)+ey.at(ey.size()-1)+0.005;

//   float y_min = -0.004;
//   float y_max = 0.014;


//   // pad settings
//   TH1F *hPad = (TH1F*)gPad->DrawFrame(0,y_min,130,y_max);
//   hPad->GetXaxis()->SetTitle("H_{T}");
//   hPad->GetYaxis()->SetTitle("(M_{ee}/M_{Z})^{data} - (M_{ee}/M_{Z})^{mc}");
//   hPad->GetYaxis()->SetTitleOffset(tYoffset);
//   hPad->GetXaxis()->SetTitleOffset(tXoffset);
//   hPad->GetXaxis()->SetLabelSize(labSize);
//   hPad->GetXaxis()->SetTitleSize(labSize);
//   hPad->GetYaxis()->SetLabelSize(labSize);
//   hPad->GetYaxis()->SetTitleSize(labSize);
//   finalGraph->Draw("P");
//   cplot->Print((plotFolderName+"/"+folderName+"/EoP_vs_Et_"+string_year+".png").c_str(),".png");
  
  std::cout << " plottato tutto " << std::endl;
  //std::cout << "CREATI I FILES" << std::endl;
  return (0);
}
