// g++ -Wall -o trisCheckScale_AllR9 `root-config --cflags --glibs` ../Utils/setTDRStyle.cc ../Utils/ntupleUtils.cc ../Utils/stabilityUtils.cc ../Utils/ConvoluteTemplate.cc ../Utils/histoFunc.h ../Utils/TPileupReweighting.h trisCheckScale_AllR9.cpp

#include "../Utils/setTDRStyle.h"
#include "../Utils/histoFunc.h"
#include "../Utils/ConvoluteTemplate.h"
#include "../Utils/ntupleUtils.h"
#include "../Utils/stabilityUtils.h"
#include "../Utils/TPileupReweighting.h"
#include "../Utils/GetShervinCorrections.h"


#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TProfile.h"
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


  /////////////// vs R9                                                                                          
  TF1* R9_low_2011 = new TF1("R9_low_2011", "[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3)", 0., 0.94);
  R9_low_2011->SetParameters(-0.018477, -0.0168714, 0.0699123, -0.0278308);
//   TF1* R9_low_2012 = new TF1("R9_low_2012", "[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3)", 0., 0.94);
//   R9_low_2012->SetParameters(-0.222858, 0.931626, -1.29396, 0.598472);
  TF1* R9_low_2012 = new TF1("R9_low_2012", "[0] ", 0., 0.94);
  R9_low_2012->SetParameter(0, 3.05197369770460669e-03);
  TF1* R9_hig_2011 = new TF1("R9_hig_2011", "[0] + [1]*x + [2]*pow(x,2) ", 0.94, 1.02);
  R9_hig_2011->SetParameters(-0.405236, 0.634916, -0.212128);
  TF1* R9_hig_2012 = new TF1("R9_hig_2012", "[0] + [1]*x + [2]*pow(x,2) ", 0.94, 1.02);
  R9_hig_2012->SetParameters(-2.30712976989725455e-01, 2.92312432577749526e-01, -4.51976365389429174e-02);

  TF1* R9_2012 = new TF1("R9_2012", "[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4)", 0., 1.02);
  R9_2012->SetParameters(3.96626752947944194e-01, -2.50481639287510305e+00, 5.83795261983192226e+00, -5.92036026163647744e+00, 2.20583006466489717e+00);


  ////////////// vs Et                                                                                           
  TF1* Et_highR9_2011 = new TF1("Et_highR9_2011", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
  Et_highR9_2011->SetParameters(1.09510022615926274e-02, 2.73732516468546162e-02, -7.40226068649908579e-03);
  TF1* Et_highR9_2012 = new TF1("Et_highR9_2012", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
  Et_highR9_2011->SetParameters(1.68849355223906934e-02, 1.81070156122320122e-02, -8.84331692745581217e-03);
  TF1* Et_lowR9_2011 = new TF1("Et_lowR9_2011", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
  Et_lowR9_2011->SetParameters(1.53853975155959805e-01, 1.40817340033172700e-01, -1.52418979619344919e-01);
  TF1* Et_lowR9_2012 = new TF1("Et_lowR9_2012", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
  Et_lowR9_2012->SetParameters(3.22375044062429561e-02, 6.32870898354409052e-03, -7.86093987242658492e-03);


  //  bool UsePhotonRegression = false;
  bool UsePhotonRegression = true;

  //  bool applyCorrections = true;
  bool applyCorrections = false;

  bool  useShCorr = true;
  //bool  useShCorr = false;

  //-----------------
  // Input parameters
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  char* EBEE = argv[1];

  char* ENE = argv[2];
  int PU = atoi(argv[3]);
  int evtsPerPoint = atoi(argv[4]);
  std::string string_year = argv[5];
  int year = atoi(argv[5]);
  //  std::string doVsEach = argv[6];
  float absEtaMin = -1.;
  float absEtaMax = -1.;
  int IetaMin = -1;
  int IetaMax = -1;
  int IphiMin = -1;
  int IphiMax = -1;

//   if(year == 2012) applyCorrections = false;
//   if(EBEE  == "EE") applyCorrections = false;

  /*
  if(argc >= 4)
  {
    absEtaMin = atof(argv[3]);
    absEtaMax = atof(argv[4]);
  }
  if(argc >= 5)
  {
    IetaMin = atoi(argv[5]);
    IetaMax = atoi(argv[6]);
    IphiMin = atoi(argv[7]);
    IphiMax = atoi(argv[8]);
  }
  */

  std::cout << "EBEE:         " << EBEE         << std::endl;
  //std::cout << "LOWHIGH:      " << LOWHIGH       << std::endl;
  std::cout << "ENE:          " << ENE           << std::endl;
  std::cout << "PU:           " << PU            << std::endl;
  std::cout << "evtsPerPoint: " << evtsPerPoint  << std::endl;
  std::cout << "year:      " << year       << std::endl;
  //  std::cout << "doVsEach:      " << doVsEach       << std::endl;
  std::cout << "IetaMin: "      << IetaMin       << std::endl;
  std::cout << "IetaMax: "      << IetaMax       << std::endl;
  std::cout << "IphiMin: "      << IphiMin       << std::endl;
  std::cout << "IphiMax: "      << IphiMax       << std::endl;
  

  TPileupReweighting* puReweighting;
  //2012 prompt                                                                                                                                     
  if(year == 2012) puReweighting =
    new TPileupReweighting("/afs/cern.ch/work/a/amartell/public/weights/PUweights_DYJetsToLL_Summer12_ABC_TrueNumInteractions.root","pileup");
    //    new TPileupReweighting("/afs/cern.ch/work/a/amartell/public/weights/PUweights_DYJetsToLL_Summer12_Prompt_TrueNumInteractions.root","hweights");
  //2011                                                                                                                                            
  if(year == 2011) puReweighting =
    new TPileupReweighting("/afs/cern.ch/work/a/amartell/public/weights/PUweights_2011_DYJetsToLL_Fall2011_TrueNumInteractions.root", "hweights");

  std::string R9MOD = std::string("HIGH");
  std::string ENERGY = std::string(ENE);
  
  //-------------------
  // Define in/outfiles
  std::string folderName;
  if(PU == 0)
  	folderName = std::string(EBEE)+"_"+std::string(ENE)+"_noPU";
  if(PU == 1)
        folderName = std::string(EBEE)+"_"+std::string(ENE);
  //if( strcmp(LOWHIGH,"")==0 ) folderName = std::string(EBEE);
  //if( strcmp(EBEE,"")==0 ) folderName = std::string(LOWHIGH);
  if( (absEtaMin != -1.) && (absEtaMax != -1.) )
  {
    char absEtaBuffer[50];
    sprintf(absEtaBuffer,"_%.2f-%.2f",absEtaMin,absEtaMax);
    folderName += std::string(absEtaBuffer);
  } 
  
  if( (IetaMin != -1.) && (IetaMax != -1.) )
  {
    char absEtaBuffer[50];
    sprintf(absEtaBuffer,"_Ieta_%d-%d",IetaMin,IetaMax);
    folderName += std::string(absEtaBuffer);
  } 
   
  if( (IphiMin != -1.) && (IphiMax != -1.) )
  {
    char absEtaBuffer[50];
    sprintf(absEtaBuffer,"_Iphi_%d-%d",IphiMin,IphiMax);
    folderName += std::string(absEtaBuffer);
  }   
  
  // Get trees
  std::cout << std::endl;

  std::string nameNtuples = "simpleNtupleEoverP/SimpleNtupleEoverP";
  //  if(year == 2011) nameNtuples = "ntu";
  TChain* ntu_MC = new TChain(nameNtuples.c_str());
  TChain* ntu_DA = new TChain(nameNtuples.c_str());

  //2012                                                                                                  
  if(year == 2012){
    ntu_MC->Add("/tmp/amartell/DYJToLL_M50_TuneZ2S_8TeV-mad_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
    //    ntu_MC->Add("/tmp/amartell/DYJetsToLL_Summer12_START50_V15_noLLR_All_simpleNtuple_mc.root");
    ntu_DA->Add("/tmp/amartell/DoubleElectron_Run2012AB_All_simpleNtuple.root");
    ntu_DA->Add("/tmp/amartell/DoubleElectron_Run2012C-PromptReco-v1-2_All_simpleNtuple.root");
  }
  //2011                                                                                                  
  if(year == 2011){
    ntu_DA->Add("/tmp/amartell/DoubleElectron-RUN2011_AB.root");
    ntu_MC->Add("/tmp/amartell/DYJetsToLL_Fall11_START44_V9B_ok_All_simpleNtuple_mc.root");
  }

  std::cout << "     REFERENCE: " << std::setw(8) << ntu_MC->GetEntries() << " entries" << std::endl;
  std::cout << "     DATA: " << std::setw(8) << ntu_DA->GetEntries() << " entries" << std::endl;


  if (ntu_DA->GetEntries() == 0 || ntu_MC->GetEntries() == 0 )
  {
    std::cout << "Error: At least one file is empty" << std::endl; 
    return -1;
  }
  
  std::vector<int> run_DA, time_DA, Z_DA, PV_DA;
  std::vector<int> run_MC, time_MC, Z_MC, PV_MC;
  std::vector<float> scE_DA, scEt_reg_DA,scE_reg_DA, R9_DA, P_DA, EoP_DA, Et_DA, scEta_DA, ES_DA, isEB_DA, e3x3_DA,e5x5_DA, scERaw_DA;
  std::vector<float> scE_MC, scEt_reg_MC, scE_reg_MC, R9_MC, P_MC, EoP_MC, Et_MC, scEta_MC, ES_MC, isEB_MC, puRe, e3x3_MC, e5x5_MC, scERaw_MC;
  std::vector<int> charge_DA, charge_MC;

  // Set branch addresses
  int isZ,runId,timeStamp,nVtx;
  float npu;
  
  ntu_DA->SetBranchStatus("*",0);
  ntu_DA->SetBranchStatus("runId",1);
  ntu_DA->SetBranchStatus("timeStampHigh",1);
  ntu_DA->SetBranchStatus("isZ",1);
  ntu_DA->SetBranchStatus("PV_n",1);
  ntu_DA->SetBranchAddress("runId", &runId);  
  ntu_DA->SetBranchAddress("timeStampHigh", &timeStamp);  
  ntu_DA->SetBranchAddress("isZ", &isZ);
  ntu_DA->SetBranchAddress("PV_n",&nVtx);
  
  ntu_MC->SetBranchStatus("*",0);
  ntu_MC->SetBranchStatus("PUit_TrueNumInteractions", 1);
  ntu_MC->SetBranchStatus("runId",1);
  ntu_MC->SetBranchStatus("timeStampHigh",1);
  ntu_MC->SetBranchStatus("isZ",1);
  ntu_MC->SetBranchStatus("PV_n",1);
  ntu_MC->SetBranchAddress("PUit_TrueNumInteractions", &npu);
  ntu_MC->SetBranchAddress("runId", &runId);  
  ntu_MC->SetBranchAddress("timeStampHigh", &timeStamp);  
  ntu_MC->SetBranchAddress("isZ", &isZ);
  ntu_MC->SetBranchAddress("PV_n",&nVtx);
  
  // Electron data
  float scEne1, scEneReg1, R9_pho1, EoP1, scEt1, scEta1, elePhi1, ES1, P1, scERaw1, e3x31, e5x51;
  float scEne2, scEneReg2, R9_pho2, EoP2, scEt2, scEta2, elePhi2, ES2, P2, scERaw2, e3x32, e5x52;
  int isEB1,isEB2;
  int ele1_charge, ele2_charge;
 

  ntu_DA->SetBranchStatus("ele1_scE", 1);       ntu_DA->SetBranchAddress("ele1_scE", &scEne1);
  ntu_DA->SetBranchStatus("ele1_scEt", 1);      ntu_DA->SetBranchAddress("ele1_scEt", &scEt1);
  ntu_DA->SetBranchStatus("ele1_scEta", 1);     ntu_DA->SetBranchAddress("ele1_scEta", &scEta1);
  if(!UsePhotonRegression) {
    ntu_DA->SetBranchStatus("ele1_scE_regression", 1);      ntu_DA->SetBranchAddress("ele1_scE_regression", &scEneReg1);
    ntu_DA->SetBranchStatus("ele2_scE_regression",1);       ntu_DA->SetBranchAddress("ele2_scE_regression", &scEneReg2);
  }
  else {
    ntu_DA->SetBranchStatus("ele1_scE_regression_PhotonTuned", 1);    ntu_DA->SetBranchAddress("ele1_scE_regression_PhotonTuned", &scEneReg1);
    ntu_DA->SetBranchStatus("ele2_scE_regression_PhotonTuned",1);     ntu_DA->SetBranchAddress("ele2_scE_regression_PhotonTuned", &scEneReg2);
  }
  ntu_DA->SetBranchStatus("ele1_scERaw",1);      ntu_DA->SetBranchAddress("ele1_scERaw",&scERaw1);
  ntu_DA->SetBranchStatus("ele1_e3x3",1);        ntu_DA->SetBranchAddress("ele1_e3x3", &e3x31);
  ntu_DA->SetBranchStatus("ele1_e5x5",1);        ntu_DA->SetBranchAddress("ele1_e5x5", &e5x51);
  ntu_DA->SetBranchStatus("ele1_EOverP",1);      ntu_DA->SetBranchAddress("ele1_EOverP",&EoP1);
  ntu_DA->SetBranchStatus("ele1_isEB",1);        ntu_DA->SetBranchAddress("ele1_isEB",&isEB1);
  ntu_DA->SetBranchStatus("ele1_es", 1);         ntu_DA->SetBranchAddress("ele1_es", &ES1);
  ntu_DA->SetBranchStatus("ele1_tkP",1);         ntu_DA->SetBranchAddress("ele1_tkP", &P1);
  ntu_DA->SetBranchStatus("ele1_charge",1);      ntu_DA->SetBranchAddress("ele1_charge", &ele1_charge);
  ntu_DA->SetBranchStatus("ele2_scE", 1);        ntu_DA->SetBranchAddress("ele2_scE", &scEne2);
  ntu_DA->SetBranchStatus("ele2_scEta", 1);      ntu_DA->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_DA->SetBranchStatus("ele2_scEt", 1);       ntu_DA->SetBranchAddress("ele2_scEt", &scEt2);
  ntu_DA->SetBranchStatus("ele2_e3x3",1);        ntu_DA->SetBranchAddress("ele2_e3x3", &e3x32);
  ntu_DA->SetBranchStatus("ele2_e5x5",1);        ntu_DA->SetBranchAddress("ele2_e5x5", &e5x52);
  ntu_DA->SetBranchStatus("ele2_scERaw",1);      ntu_DA->SetBranchAddress("ele2_scERaw",&scERaw2);
  ntu_DA->SetBranchStatus("ele2_EOverP",1);      ntu_DA->SetBranchAddress("ele2_EOverP",&EoP2);
  ntu_DA->SetBranchStatus("ele2_isEB",1);        ntu_DA->SetBranchAddress("ele2_isEB",&isEB2);
  ntu_DA->SetBranchStatus("ele2_es", 1);         ntu_DA->SetBranchAddress("ele2_es", &ES2);
  ntu_DA->SetBranchStatus("ele2_tkP",1);         ntu_DA->SetBranchAddress("ele2_tkP", &P2);
  ntu_DA->SetBranchStatus("ele2_charge",1);      ntu_DA->SetBranchAddress("ele2_charge", &ele2_charge);
  ntu_DA->SetBranchStatus("ele1_phi", 1);   ntu_DA->SetBranchAddress("ele1_phi", &elePhi1);
  ntu_DA->SetBranchStatus("ele2_phi", 1);   ntu_DA->SetBranchAddress("ele2_phi", &elePhi2);


  ///////////////////////                                                                                                                     
  ntu_MC->SetBranchStatus("ele1_scE", 1);       ntu_MC->SetBranchAddress("ele1_scE", &scEne1);
  ntu_MC->SetBranchStatus("ele1_scEt", 1);      ntu_MC->SetBranchAddress("ele1_scEt", &scEt1);
  ntu_MC->SetBranchStatus("ele1_scEta", 1);     ntu_MC->SetBranchAddress("ele1_scEta", &scEta1);
  if(!UsePhotonRegression) {
    ntu_MC->SetBranchStatus("ele1_scE_regression", 1);      ntu_MC->SetBranchAddress("ele1_scE_regression", &scEneReg1);
    ntu_MC->SetBranchStatus("ele2_scE_regression",1);       ntu_MC->SetBranchAddress("ele2_scE_regression", &scEneReg2);
  }
  else {
    ntu_MC->SetBranchStatus("ele1_scE_regression_PhotonTuned", 1);    ntu_MC->SetBranchAddress("ele1_scE_regression_PhotonTuned", &scEneReg1);
    ntu_MC->SetBranchStatus("ele2_scE_regression_PhotonTuned",1);     ntu_MC->SetBranchAddress("ele2_scE_regression_PhotonTuned", &scEneReg2);
  }
  ntu_MC->SetBranchStatus("ele1_scERaw",1);      ntu_MC->SetBranchAddress("ele1_scERaw",&scERaw1);
  ntu_MC->SetBranchStatus("ele1_e3x3",1);        ntu_MC->SetBranchAddress("ele1_e3x3", &e3x31);
  ntu_MC->SetBranchStatus("ele1_e5x5",1);        ntu_MC->SetBranchAddress("ele1_e5x5", &e5x51);
  ntu_MC->SetBranchStatus("ele1_EOverP",1);      ntu_MC->SetBranchAddress("ele1_EOverP",&EoP1);
  ntu_MC->SetBranchStatus("ele1_isEB",1);        ntu_MC->SetBranchAddress("ele1_isEB",&isEB1);
  ntu_MC->SetBranchStatus("ele1_es", 1);         ntu_MC->SetBranchAddress("ele1_es", &ES1);
  ntu_MC->SetBranchStatus("ele1_tkP",1);         ntu_MC->SetBranchAddress("ele1_tkP", &P1);
  ntu_MC->SetBranchStatus("ele1_charge",1);      ntu_MC->SetBranchAddress("ele1_charge", &ele1_charge);
  ntu_MC->SetBranchStatus("ele2_scE", 1);        ntu_MC->SetBranchAddress("ele2_scE", &scEne2);
  ntu_MC->SetBranchStatus("ele2_scEta", 1);      ntu_MC->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_MC->SetBranchStatus("ele2_scEt", 1);       ntu_MC->SetBranchAddress("ele2_scEt", &scEt2);
  ntu_MC->SetBranchStatus("ele2_e3x3",1);        ntu_MC->SetBranchAddress("ele2_e3x3", &e3x32);
  ntu_MC->SetBranchStatus("ele2_e5x5",1);        ntu_MC->SetBranchAddress("ele2_e5x5", &e5x52);
  ntu_MC->SetBranchStatus("ele2_scERaw",1);      ntu_MC->SetBranchAddress("ele2_scERaw",&scERaw2);
  ntu_MC->SetBranchStatus("ele2_EOverP",1);      ntu_MC->SetBranchAddress("ele2_EOverP",&EoP2);
  ntu_MC->SetBranchStatus("ele2_isEB",1);        ntu_MC->SetBranchAddress("ele2_isEB",&isEB2);
  ntu_MC->SetBranchStatus("ele2_es", 1);         ntu_MC->SetBranchAddress("ele2_es", &ES2);
  ntu_MC->SetBranchStatus("ele2_tkP",1);         ntu_MC->SetBranchAddress("ele2_tkP", &P2);
  ntu_MC->SetBranchStatus("ele2_charge",1);      ntu_MC->SetBranchAddress("ele2_charge", &ele2_charge);
  ntu_MC->SetBranchStatus("ele1_phi", 1);   ntu_MC->SetBranchAddress("ele1_phi", &elePhi1);
  ntu_MC->SetBranchStatus("ele2_phi", 1);   ntu_MC->SetBranchAddress("ele2_phi", &elePhi2);
  //////////////////////                                                                               
  for(int ientry = 0; ientry < ntu_DA -> GetEntries(); ientry++)
  {
  	if( (ientry%100000 == 0) ) std::cout << "reading DATA entry " << ientry << "\r" << std::flush;
        ntu_DA->GetEntry(ientry);  

	if(isZ == 0) continue;
	//	if(nVtx < 10 || nVtx > 20) continue;

        run_DA.push_back(runId);

	float corrEtR9_1 = 1.;
	float corrEtR9_2 = 1.;


	if(useShCorr == true){
	  corrEtR9_1 = corrEtR9_1 * GetShervingCorrections(scEta1, e3x31/scERaw1, runId);
	  corrEtR9_2 = corrEtR9_2 * GetShervingCorrections(scEta2, e3x32/scERaw2, runId);
	}

// 	std::cout << " corrEtR9_1 = " << corrEtR9_1 << std::endl;
// 	std::cout << " corrEtR9_2 = " << corrEtR9_2 << std::endl;


	charge_DA.push_back(ele1_charge);
	charge_DA.push_back(ele2_charge);

        time_DA.push_back(timeStamp);
        Z_DA.push_back(isZ);
	PV_DA.push_back(nVtx);    
	scE_DA.push_back(scEne1);
	scE_DA.push_back(scEne2);
	scE_reg_DA.push_back(scEneReg1*corrEtR9_1);
	scE_reg_DA.push_back(scEneReg2*corrEtR9_2);
	scEt_reg_DA.push_back(scEneReg1/scEne1*scEt1*corrEtR9_1);
 	scEt_reg_DA.push_back(scEneReg2/scEne2*scEt2*corrEtR9_2);

	R9_DA.push_back(e3x31/scERaw1);
	R9_DA.push_back(e3x32/scERaw2);

	P_DA.push_back(P1);
	P_DA.push_back(P2);

	EoP_DA.push_back(EoP1);
	EoP_DA.push_back(EoP2);
	Et_DA.push_back(scEt1);
        Et_DA.push_back(scEt2);
	scEta_DA.push_back(scEta1);
	scEta_DA.push_back(scEta2);
        ES_DA.push_back(ES1);
	ES_DA.push_back(ES2);
	isEB_DA.push_back(isEB1);
	isEB_DA.push_back(isEB2);
	e5x5_DA.push_back(e5x51);
	e5x5_DA.push_back(e5x52);
	e3x3_DA.push_back(e3x31);
	e3x3_DA.push_back(e3x32);
	scERaw_DA.push_back(scERaw1);
	scERaw_DA.push_back(scERaw2);
  }
  
  std::cout << std::endl;
  float ww = 1.;
  for(int ientry = 0; ientry < ntu_MC -> GetEntries(); ientry++)
  {

        if( (ientry%100000 == 0) ) std::cout << "reading MC entry " << ientry << "\r" << std::flush;
        ntu_MC->GetEntry(ientry);  
 
	if(isZ == 0) continue;

    float R9_ele1 = e3x31/scERaw1;
    if(year == 2012 && isEB1 == 1) R9_ele1 = 0.0010 + 1.0045 * e3x31/scERaw1;
    if(year == 2012 && isEB1 == 0) R9_ele1 = -0.0007 + 1.0086 * e3x31/scERaw1;
    if(year == 2011) R9_ele1 = 1.0035 * e3x31/scERaw1;
    float R9_ele2 = e3x32/scERaw2;
    if(year == 2012 && isEB2 == 1) R9_ele2 = 0.0010 + 1.0045 * e3x32/scERaw2;
    if(year == 2012 && isEB2 == 0) R9_ele2 = -0.0007 + 1.0086 * e3x32/scERaw2;
    if(year == 2011) R9_ele2 = 1.0035 * e3x32/scERaw2;
  

    float energySmearing1 = gRandom->Gaus(1.,0.0075);
    float energySmearing2 = gRandom->Gaus(1.,0.0075);

    if(isEB1 == 1 && year == 2011 && fabs(scEta1) < 1. && R9_ele1 > 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0089);
    if(isEB2 == 1 && year == 2011 && fabs(scEta2) < 1. && R9_ele2 > 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0089);
    if(isEB1 == 1 && year == 2012 && fabs(scEta1) < 1. && R9_ele1 > 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0099);
    if(isEB2 == 1 && year == 2012 && fabs(scEta2) < 1. && R9_ele2 > 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0099);

    if(isEB1 == 1 && year == 2011 && fabs(scEta1) < 1. && R9_ele1 < 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0109);
    if(isEB2 == 1 && year == 2011 && fabs(scEta2) < 1. && R9_ele2 < 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0109);
    if(isEB1 == 1 && year == 2012 && fabs(scEta1) < 1. && R9_ele1 < 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0109);
    if(isEB2 == 1 && year == 2012 && fabs(scEta2) < 1. && R9_ele2 < 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0109);

    if(isEB1 == 1 && year == 2011 && fabs(scEta1) > 1. && R9_ele1 > 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0156);
    if(isEB2 == 1 && year == 2011 && fabs(scEta2) > 1. && R9_ele2 > 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0156);
    if(isEB1 == 1 && year == 2012 && fabs(scEta1) > 1. && R9_ele1 > 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0195);
    if(isEB2 == 1 && year == 2012 && fabs(scEta2) > 1. && R9_ele2 > 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0195);

    if(isEB1 == 1 && year == 2011 && fabs(scEta1) > 1. && R9_ele1 < 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0203);
    if(isEB2 == 1 && year == 2011 && fabs(scEta2) > 1. && R9_ele2 < 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0203);
    if(isEB1 == 1 && year == 2012 && fabs(scEta1) > 1. && R9_ele1 < 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0198);
    if(isEB2 == 1 && year == 2012 && fabs(scEta2) > 1. && R9_ele2 < 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0198);

    //EE 

    if(isEB1 == 0 && year == 2011 && fabs(scEta1) < 2. && R9_ele1 > 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0303);
    if(isEB2 == 0 && year == 2011 && fabs(scEta2) < 2. && R9_ele2 > 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0303);
    if(isEB1 == 0 && year == 2012 && fabs(scEta1) < 2. && R9_ele1 > 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0291);
    if(isEB2 == 0 && year == 2012 && fabs(scEta2) < 2. && R9_ele2 > 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0291);

    if(isEB1 == 0 && year == 2011 && fabs(scEta1) < 2. && R9_ele1 < 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0326);
    if(isEB2 == 0 && year == 2011 && fabs(scEta2) < 2. && R9_ele2 < 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0326);
    if(isEB1 == 0 && year == 2012 && fabs(scEta1) < 2. && R9_ele1 < 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0275);
    if(isEB2 == 0 && year == 2012 && fabs(scEta2) < 2. && R9_ele2 < 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0275);

    if(isEB1 == 0 && year == 2011 && fabs(scEta1) > 2. && R9_ele1 > 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0318);
    if(isEB2 == 0 && year == 2011 && fabs(scEta2) > 2. && R9_ele2 > 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0318);
    if(isEB1 == 0 && year == 2012 && fabs(scEta1) > 2. && R9_ele1 > 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0325);
    if(isEB2 == 0 && year == 2012 && fabs(scEta2) > 2. && R9_ele2 > 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0325);

    if(isEB1 == 0 && year == 2011 && fabs(scEta1) > 2. && R9_ele1 < 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0331);
    if(isEB2 == 0 && year == 2011 && fabs(scEta2) > 2. && R9_ele2 < 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0331);
    if(isEB1 == 0 && year == 2012 && fabs(scEta1) > 2. && R9_ele1 < 0.94 )       energySmearing1 = gRandom->Gaus(1.,0.0356);
    if(isEB2 == 0 && year == 2012 && fabs(scEta2) > 2. && R9_ele2 < 0.94 )       energySmearing2 = gRandom->Gaus(1.,0.0356);


    //    if(reweightZtoH == false && reweightEta == false && reweightR9 == false){
      if(isEB1 == 1 )       energySmearing1 = gRandom->Gaus(1.,0.00078);
      if(isEB2 == 1 )       energySmearing2 = gRandom->Gaus(1.,0.00078);
      if(isEB1 == 1 && year == 2012 && R9_ele1 > 0.94)       energySmearing1 = gRandom->Gaus(1.,0.0008);
      if(isEB2 == 1 && year == 2012 && R9_ele2 > 0.94)       energySmearing2 = gRandom->Gaus(1.,0.0008);
      if(isEB1 == 1 && year == 2011 && R9_ele1 < 0.94)       energySmearing1 = gRandom->Gaus(1.,0.0008);
      if(isEB2 == 1 && year == 2011 && R9_ele2 < 0.94)       energySmearing2 = gRandom->Gaus(1.,0.0008);

      //EE                                                                                             
      if(year == 2011 && isEB1 == 0 )       energySmearing1 = gRandom->Gaus(1.,0.06);
      if(year == 2011 && isEB2 == 0 )       energySmearing2 = gRandom->Gaus(1.,0.06);
      if(year == 2012 && isEB1 == 0 )       energySmearing1 = gRandom->Gaus(1.,0.043);
      if(year == 2012 && isEB2 == 0 )       energySmearing2 = gRandom->Gaus(1.,0.043);
      if(year == 2011 && isEB1 == 0 && R9_ele1 > 0.94)       energySmearing1 = gRandom->Gaus(1.,0.066);
      if(year == 2011 && isEB2 == 0 && R9_ele2 > 0.94)       energySmearing2 = gRandom->Gaus(1.,0.066);
      if(year == 2012 && isEB1 == 0 && R9_ele1 > 0.94)       energySmearing1 = gRandom->Gaus(1.,0.048);
      if(year == 2012 && isEB2 == 0 && R9_ele2 > 0.94)       energySmearing2 = gRandom->Gaus(1.,0.048);

      //    }



	ww = puReweighting->GetWeight((int)npu);
        puRe.push_back(ww);
        run_MC.push_back(runId);
        time_MC.push_back(timeStamp);
        Z_MC.push_back(isZ);
	PV_MC.push_back(nVtx);    
	scE_MC.push_back(scEne1);
	scE_MC.push_back(scEne2);

	scE_reg_MC.push_back(scEneReg1*energySmearing1);
	scE_reg_MC.push_back(scEneReg2*energySmearing2);

	scEt_reg_MC.push_back(scEneReg1/scEne1*scEt1*energySmearing1);
 	scEt_reg_MC.push_back(scEneReg2/scEne2*scEt2*energySmearing2);

	charge_MC.push_back(ele1_charge);
	charge_MC.push_back(ele2_charge);

	P_MC.push_back(P1);
	P_MC.push_back(P2);

	EoP_MC.push_back(EoP1);
	EoP_MC.push_back(EoP2);
	Et_MC.push_back(scEt1);
        Et_MC.push_back(scEt2);
	scEta_MC.push_back(scEta1);
	scEta_MC.push_back(scEta2);
        ES_MC.push_back(ES1);
	ES_MC.push_back(ES2);
	isEB_MC.push_back(isEB1);
	isEB_MC.push_back(isEB2);
	e5x5_MC.push_back(e5x51);
	e5x5_MC.push_back(e5x52);
	e3x3_MC.push_back(e3x31);
	e3x3_MC.push_back(e3x32);
	scERaw_MC.push_back(scERaw1);
	scERaw_MC.push_back(scERaw2);


	R9_MC.push_back(R9_ele1);
	R9_MC.push_back(R9_ele2);
  }
  
  // Loop and sort events
  std::cout << std::endl;
  std::cout << "***** Sort events and define bins *****" << std::endl;
  
  int nEntries = R9_DA.size();
  int nSavePts = 0;
  std::vector<bool> isSavedEntries(nEntries);
  std::vector<SorterLC> sortedEntries;
  
  for(int ientry = 0; ientry < nEntries; ++ientry)
  {
    isSavedEntries.at(ientry) = false;

    if(std::string(ENE) == "scE_reg" && scEt_reg_DA.at(ientry) <  25.) continue; 
    if(std::string(ENE) == "scE" && Et_DA.at(ientry) <  25.) continue; 
    if(std::string(ENE) == "scERaw" && Et_DA.at(ientry) <  25.) continue; 

    if (strcmp(EBEE,"EE")==0 && (fabs(scEta_DA.at(ientry)) < 1.566 || fabs(scEta_DA.at(ientry)) > 2.7 )) continue; // endcap
    if (strcmp(EBEE,"EB")==0 && (fabs(scEta_DA.at(ientry)) > 1.4442 )) continue; // endcap

    if (R9_DA.at(ientry) < 0.7) continue;

    if( (absEtaMin != -1.) && (absEtaMax != -1.) )
    {
	if( (fabs(scEta_DA.at(ientry)) < absEtaMin) || (fabs(scEta_DA.at(ientry)) > absEtaMax) ) continue;
    }
    
    isSavedEntries.at(ientry) = true;

    SorterLC dummy;
    dummy.laserCorr = R9_DA.at(ientry);
    dummy.entry = ientry;
    sortedEntries.push_back(dummy);
    nSavePts++;   
  }
  std::cout << " Effective entries = " << nSavePts << std::endl;
  std::cout << " Effective entries sortedEntries.size()= " << sortedEntries.size() << std::endl;

  std::sort(sortedEntries.begin(),sortedEntries.end(),SorterLC());
  std::cout << "DATA sorted in " << EBEE << " - " << nSavePts << " events" << std::endl;
  
  
   std::map<int,int> antiMap;
   for(unsigned int iSaved = 0; iSaved < sortedEntries.size(); ++iSaved)
   antiMap[sortedEntries.at(iSaved).entry] = iSaved; 
   
  // bins with evtsPerPoint events per bin
   std::cout << " nSavePts = " << nSavePts << std::endl;
   std::cout << " evtsPerPoint = " << evtsPerPoint << std::endl;
  int nBins = std::max(1, int(nSavePts/evtsPerPoint));
   std::cout << " nBins = " << nBins << std::endl;
  int nBinPts = int( nSavePts/nBins );
   std::cout << " nBinPts = " << nBinPts << std::endl;
  int nBinTempPts = 0;

  std::cout << "nBins = " << nBins << std::endl;
  
  std::vector<int> binEntryMax;
  binEntryMax.push_back(0);
  for(int iSaved = 0; iSaved < nSavePts; ++iSaved)
  {
    ++nBinTempPts;
    
    if( nBinTempPts == nBinPts )
    {
      binEntryMax.push_back( iSaved );      
      nBinTempPts = 0;
      //      std::cout << "binEntryMax.size() = " << binEntryMax.size() << std::endl;
    }
  }
  binEntryMax.at(nBins) = nSavePts;
  
  std::cout << " fine : nBins = " << nBins << std::endl;
 
  TVirtualFitter::SetDefaultFitter("Fumili2");
  
  // histogram definition
  
  TH1F** h_EoP_DA = new TH1F*[nBins];
  TH1F** h_EoP_MC = new TH1F*[nBins];
  TH1F** h_R9 = new TH1F*[nBins];
  TH1F** h_R9_MC = new TH1F*[nBins];
  TH1F* h_R9_allDA = new TH1F("h_R9_allDA", "", 2200, 0., 1.1);
  TH1F* h_R9_allMC = new TH1F("h_R9_allMC", "", 2200, 0., 1.1);

  TH2F* h_R9_vsET_MC = new TH2F("h_R9_vsET_MC", "", 200, 0., 200., 220, 0., 1.1);
  TH2F* h_R9_vsET_DA = new TH2F("h_R9_vsET_DA", "", 200, 0., 200., 220, 0., 1.1);
  TProfile* p_R9_vsET_MC = new TProfile("p_R9_vsET_MC", "", 200, 0., 200.);
  TProfile* p_R9_vsET_DA = new TProfile("p_R9_vsET_DA", "", 200, 0., 200.);


  h_R9_allDA->Sumw2();
  h_R9_allMC->Sumw2();
  h_R9_vsET_MC->Sumw2();
  h_R9_vsET_DA->Sumw2();
  p_R9_vsET_MC->Sumw2();
  p_R9_vsET_DA->Sumw2();


  std::vector<float> EtBinEdge;
  EtBinEdge.clear();
  std::vector<float> xNorm_single;
  float xNorm_all;

  for(int i = 0; i < nBins; ++i)
  {
    char histoName[80];
    
    sprintf(histoName, "EoP_DA_%d", i);
    h_EoP_DA[i] = new TH1F(histoName, histoName, 400, 0., 2.);
    h_EoP_DA[i] -> SetFillColor(kRed+2);
    h_EoP_DA[i] -> SetFillStyle(3004);
    h_EoP_DA[i] -> SetMarkerStyle(7);
    h_EoP_DA[i] -> SetMarkerColor(kRed+2); 
    h_EoP_DA[i] -> SetLineColor(kRed+2); 
    
    sprintf(histoName, "EoP_MC_%d", i);
    h_EoP_MC[i] = new TH1F(histoName, histoName, 400, 0., 2.);
    h_EoP_MC[i] -> SetFillColor(kGreen+2);
    h_EoP_MC[i] -> SetFillStyle(3004);
    h_EoP_MC[i] -> SetMarkerStyle(7);
    h_EoP_MC[i] -> SetMarkerColor(kGreen+2);
    h_EoP_MC[i] -> SetLineColor(kGreen+2);
    
    sprintf(histoName, "R9_%d", i);
    h_R9[i] = new TH1F(histoName, histoName, 2200, 0.7, 1.1);
    h_R9[i]->SetLineColor(kRed+2);

    sprintf(histoName, "R9_MC_%d", i);
    h_R9_MC[i] = new TH1F(histoName, histoName, 2200, 0.7, 1.1);
    h_R9_MC[i]->SetLineColor(kGreen+2);

    h_EoP_DA[i]->Sumw2(); 
    h_EoP_MC[i]->Sumw2();
    h_R9[i]->Sumw2();
    h_R9_MC[i]->Sumw2();

  }
  
  std::cout << " R9_MC.size() = " << R9_MC.size() << std::endl;
  std::cout << " R9_DA.size() = " << R9_DA.size() << std::endl;


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
  for(int unsigned ientry = 0; ientry < R9_DA.size(); ++ientry)
  {
    if( (ientry%100000 == 0) ) std::cout << "reading entry " << ientry << std::endl;
      
    if( isSavedEntries.at(ientry) == false) continue;
     
     int iSaved = antiMap[ientry];
     int bin = -1;
    
    for(bin = 0; bin < nBins; ++bin)
      if( iSaved >= binEntryMax.at(bin) && iSaved < binEntryMax.at(bin+1) )
    break;


    if(strcmp(ENE,"scE_reg")==0)
      h_EoP_DA[bin]->Fill((scE_reg_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)));
    if(strcmp(ENE,"scE")==0)
    	h_EoP_DA[bin]->Fill((scE_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)));
    if(strcmp(ENE,"scERaw")==0)
    	h_EoP_DA[bin]->Fill((scERaw_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)));
    if(strcmp(ENE,"e5x5")==0)
    	h_EoP_DA[bin]->Fill((e5x5_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)));
    if(strcmp(ENE,"e3x3")==0)
    	h_EoP_DA[bin]->Fill((e3x3_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)));
    
    h_R9[bin]->Fill(R9_DA.at(ientry));
    h_R9_allDA->Fill(R9_DA.at(ientry));

    if(std::string(ENE) == "scE_reg") {
      h_R9_vsET_DA->Fill(scEt_reg_DA.at(ientry), R9_DA.at(ientry));
      p_R9_vsET_DA->Fill(scEt_reg_DA.at(ientry), R9_DA.at(ientry));
    }
    if(std::string(ENE) == "scE") {
      h_R9_vsET_DA->Fill(Et_DA.at(ientry), R9_DA.at(ientry));
      p_R9_vsET_DA->Fill(Et_DA.at(ientry), R9_DA.at(ientry));
    }
    if(std::string(ENE) == "scERaw") {
      h_R9_vsET_DA->Fill(Et_DA.at(ientry), R9_DA.at(ientry));
      p_R9_vsET_DA->Fill(Et_DA.at(ientry), R9_DA.at(ientry));
    }
  }
  
  std::cout << " dati fillati " << std::endl;
  
  for(int bin = 0; bin < nBins; bin++)
  {
    std::cout << "h_R9[bin]->GetEntries() =  " << h_R9[bin]->GetEntries() << std::endl;
    std::cout << "h_EoP_DA[bin]->GetEntries() =  " << h_EoP_DA[bin]->GetEntries() << std::endl;

    for(int i = 1; i < h_R9[bin]->GetNbinsX()+1; i++)
    {
	if(h_R9[bin]->GetBinContent(i) > 0) {
// 	  std::cout << " if > 0  bin = " << bin << std::endl;
// 	  std::cout << " h_R9[bin]->GetBinCenter(i) = " << h_R9[bin]->GetBinCenter(i) << std::endl;
// 	  std::cout << " h_R9[bin]->GetBinWidth(i) = " << h_R9[bin]->GetBinWidth(i) << std::endl;

	  float valore = h_R9[bin]->GetBinCenter(i)-h_R9[bin]->GetBinWidth(i);

// 	  std::cout << " valore = " << valore << std::endl;
// 	  std::cout << " EtBinEdge.size = " << EtBinEdge.size() << std::endl;

	  EtBinEdge.push_back( valore);
	  //	  std::cout << " EtBinEdge.at(i) = " << EtBinEdge.size() << std::endl;
	  break;
	}  
    }
  }

  std::cout << " EtBinEdge ok " << std::endl;

  for(unsigned int ientry = 0; ientry < R9_MC.size(); ++ientry)
  {
    if( (ientry%100000 == 0) ) std::cout << "reading entry " << ientry << std::endl;

    if (strcmp(EBEE,"EE")==0 && (fabs(scEta_MC.at(ientry)) < 1.566 || fabs(scEta_MC.at(ientry)) > 2.7 )) continue; // endcap
    if (strcmp(EBEE,"EB")==0 && (fabs(scEta_MC.at(ientry)) > 1.4442 )) continue; // endcap

    if(std::string(ENE) == "scE_reg" && scEt_reg_MC.at(ientry) <  25.) continue; 
    if(std::string(ENE) == "scE" && Et_MC.at(ientry) <  25.) continue; 
    if(std::string(ENE) == "scERaw" && Et_MC.at(ientry) <  25.) continue; 
    
    if(R9_MC.at(ientry) < 0.7) continue; 

    if( (absEtaMin != -1.) && (absEtaMax != -1.) )
    {
	if( (fabs(scEta_MC.at(ientry)) < absEtaMin) || (fabs(scEta_MC.at(ientry)) > absEtaMax) ) continue;
    }


    for(unsigned int bin = 0; bin < EtBinEdge.size(); ++bin)
    {
      if( (bin != EtBinEdge.size()-1 && (R9_MC.at(ientry)) > EtBinEdge.at(bin) && (R9_MC.at(ientry)) < EtBinEdge.at(bin+1)) ||
	  (bin == EtBinEdge.size()-1 && (R9_MC.at(ientry)) > EtBinEdge.at(bin) ) )
	{
	  if(PU == 0)
          {
		if(strcmp(ENE,"scE_reg")==0)
    			h_EoP_MC[(int)bin]->Fill((scE_reg_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)));
    		if(strcmp(ENE,"scE")==0)
    			h_EoP_MC[(int)bin]->Fill((scE_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)));
    		if(strcmp(ENE,"scERaw")==0)
    			h_EoP_MC[(int)bin]->Fill((scERaw_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)));
    		if(strcmp(ENE,"e5x5")==0)
    			h_EoP_MC[(int)bin]->Fill((e5x5_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)));
		if(strcmp(ENE,"e3x3")==0)
    			h_EoP_MC[(int)bin]->Fill((e3x3_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)));
	  	h_R9_MC[int(bin)]->Fill(R9_MC.at(ientry));
                break;
	  }
          if(PU == 1)
          {
		if(strcmp(ENE,"scE_reg")==0)
		  h_EoP_MC[(int)bin]->Fill((scE_reg_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)), puRe.at(int(ientry/2)));
    		if(strcmp(ENE,"scE")==0)
    			h_EoP_MC[(int)bin]->Fill((scE_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)),puRe.at(int(ientry/2)));
    		if(strcmp(ENE,"scERaw")==0)
    			h_EoP_MC[(int)bin]->Fill((scERaw_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)),puRe.at(int(ientry/2)));
    		if(strcmp(ENE,"e5x5")==0)
    			h_EoP_MC[(int)bin]->Fill((e5x5_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)),puRe.at(int(ientry/2)));
		if(strcmp(ENE,"e3x3")==0)
    			h_EoP_MC[(int)bin]->Fill((e3x3_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)),puRe.at(int(ientry/2)));
	  	h_R9_MC[int(bin)]->Fill(R9_MC.at(ientry), puRe.at(int(ientry/2)));
		//		if(bin == 34) std::cout << " scE_reg_MC.at(ientry) = " << scE_reg_MC.at(ientry) << std::endl;
                break;
	  }
	} 
    }
      
     if(PU == 0)
       h_R9_allMC->Fill(R9_MC.at(ientry));
     if(PU == 1){
       h_R9_allMC->Fill(R9_MC.at(ientry), puRe.at(int(ientry/2)));

       if(std::string(ENE) == "scE_reg") {
	 h_R9_vsET_MC->Fill(scEt_reg_MC.at(ientry), R9_MC.at(ientry), puRe.at(int(ientry/2)));
	 p_R9_vsET_MC->Fill(scEt_reg_MC.at(ientry), R9_MC.at(ientry), puRe.at(int(ientry/2)));
       }
       if(std::string(ENE) == "scE") {
	 h_R9_vsET_MC->Fill(Et_MC.at(ientry), R9_MC.at(ientry), puRe.at(int(ientry/2)));
	 p_R9_vsET_MC->Fill(Et_MC.at(ientry), R9_MC.at(ientry), puRe.at(int(ientry/2)));
       }
       if(std::string(ENE) == "scERaw") {
	 h_R9_vsET_MC->Fill(Et_MC.at(ientry), R9_MC.at(ientry), puRe.at(int(ientry/2)));
	 p_R9_vsET_MC->Fill(Et_MC.at(ientry), R9_MC.at(ientry), puRe.at(int(ientry/2)));
       }
     }
  }
  std::cout << " ok fino a qui ci sono " << std::endl;

  for(int i = 0; i < nBins; ++i)
  {
    //------------------------------------
    // Fill the graph for uncorrected data
    // define the fitting function
    // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )



    float xNorm = h_EoP_DA[i]->Integral()/h_EoP_MC[i]->Integral();  
    h_EoP_MC[i]->Scale(xNorm);

//     if(i==0)
//     {
//     	std::cout << " h_EoP_DA[i]->Integral() = " << h_EoP_DA[i]->Integral() << std::endl;
//     	std::cout << " h_EoP_MC[i]->Integral() = " << h_EoP_MC[i]->Integral() << std::endl;
//     	std::cout << " xNorm = " << xNorm << std::endl;
//     }



    float xNormEt = h_R9[i]->Integral()/h_R9_MC[i]->Integral();
    h_R9_MC[i]->Scale(xNormEt);
    xNorm_single.push_back(xNormEt);

    xNorm_all = h_R9_allDA->Integral()/h_R9_allMC->Integral();
    h_R9_allMC->Scale(xNorm_all);


    std::string LOWHIGH = "LOW";
    if((h_R9[i]->GetMean() > 0.94 && h_R9_MC[i]->GetMean() > 0.94)) LOWHIGH = "HIGH";

    //    h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);

    //    if(reweightZtoH  == false && reweightEta == false && reweightR9 == false){
    if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "HIGH"){h_EoP_MC[i]->Smooth(1); }
      if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "LOW"){ h_EoP_MC[i]->Smooth(1); }
      if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "HIGH"){ h_EoP_MC[i]->Smooth(2); }
      if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "LOW"){h_EoP_MC[i]->Smooth(2); }

      if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "HIGH")
	{ h_EoP_MC[i]->Smooth(150); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
      if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "LOW")
	{ h_EoP_MC[i]->Smooth(150); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
      if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "HIGH")
	{   h_EoP_MC[i]->Smooth(150); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
      if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "LOW")
	{h_EoP_MC[i]->Smooth(150); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
      //    }



    histoFunc* templateHistoFunc = new histoFunc(h_EoP_MC[i]);
    char funcName[50];
    sprintf(funcName,"f_EoP_%d",i);
    f_EoP[i] = new TF1(funcName, templateHistoFunc, 0.7, 1.3, 3, "histoFunc");
    if(std::string(EBEE) == "EE")  f_EoP[i] = new TF1(funcName, templateHistoFunc, 0.7, 1.4, 3, "histoFunc");

    xNorm = 1.;

    f_EoP[i] -> SetParName(0,"Norm"); 
    f_EoP[i] -> SetParName(1,"Scale factor"); 
    f_EoP[i] -> SetLineWidth(1); 
    f_EoP[i] -> SetNpx(10000);
    
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
      
    if (fStatus==0 && eee*k > 0.1*h_EoP_DA[i]->GetRMS()/sqrt(evtsPerPoint))
    {
      x.push_back(h_R9[i]->GetMean());
      ex.push_back((h_R9[i]->GetRMS())/sqrt(h_R9[i]->GetEntries()));
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
//   if(strcmp(LOWHIGH,"HIGH")==0 )  finalGraph->GetYaxis()->SetRangeUser(-0.004, 0.014);
//   if(strcmp(LOWHIGH,"LOW")==0 )  finalGraph->GetYaxis()->SetRangeUser(-0.03, 0.03);
  finalGraph->GetYaxis()->SetTitle("E/p_{data} - E/p_{mc}");
  finalGraph->GetXaxis()->SetRangeUser(0., 1.1);
  finalGraph->GetXaxis()->SetTitle("R9");

  std::string PLOTS_folderName = "PLOTS_R9";

  TFile pippo( (PLOTS_folderName+"/results_"+folderName+"_"+string_year+"R9.root").c_str(),"recreate");
  finalGraph->Write("finalGraph");
   h_R9_allMC->Write();
   h_R9_allDA->Write();
  h_R9_vsET_MC->Write();
  h_R9_vsET_DA->Write();
  p_R9_vsET_MC->Write();
  p_R9_vsET_DA->Write();
   for(int i = 0; i < nBins; ++i){
     h_EoP_DA[i]->Write();
     h_EoP_MC[i]->Write();
     h_R9[i]->Write();
     h_R9_MC[i]->Write();
   }
  pippo.Close();
  
//   /*
//   // Drawings
//   TPaveStats** s_EoP = new TPaveStats*[nBins];
  
//   TCanvas *c1[100]; 
//   for(int i = 0; i < nBins; ++i)
//   {    
//     char canvasName[50];
//     sprintf(canvasName, "Fits-%0d", i); 
//     c1[i] = new TCanvas(canvasName, canvasName);
//     c1[i] -> cd ();
//     h_EoP_DA[i] -> GetXaxis() -> SetTitle("E/p");  
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
//     if(PU == 0)       sprintf(Name, (PLOTS_folderName+"/"+folderName+"/noPU_fit_%d_"+string_year+".png").c_str(), i);
//     if(PU == 1)      sprintf(Name, (PLOTS_folderName+"/"+folderName+"/PU_fit_%d_"+string_year+".png").c_str(), i);
//     c1[i] -> Print(Name,".png");  
//   }
  
//   TCanvas *c2[100]; 
//   for(int i = 0; i < nBins; ++i)
//   {
//     char canvasName[50];
//     sprintf(canvasName, "R9_DA-%0d", i); 
//     c2[i] = new TCanvas(canvasName, canvasName);
//     c2[i] -> cd ();

//     h_R9[i] -> GetXaxis() -> SetTitle("R9");
//     h_R9[i]->GetYaxis()->SetRangeUser(0, std::max(h_R9[i]->GetMaximum(), h_R9_MC[i]->GetMaximum()) + 10. ); 
//     if(i<nBins-1) h_R9[i]->GetXaxis()->SetRangeUser(EtBinEdge.at(i), EtBinEdge.at(i+1)); 
//     else h_R9[i]->GetXaxis()->SetRangeUser(EtBinEdge.at(i), 1.1); 
//     h_R9[i] -> Draw();
//     h_R9_MC[i]->Draw("same");
    
//     char Name[100];
//     if(PU == 0)      sprintf(Name, (PLOTS_folderName+"/"+folderName+"/noPU_R9_%d_"+string_year+".png").c_str(), i);
//     if(PU == 1)       sprintf(Name, (PLOTS_folderName+"/"+folderName+"/PU_R9_%d_"+string_year+".png").c_str(), i);
//     c2[i] -> Print(Name,".png");  
//   }



//   TCanvas* R9_spectrum = new TCanvas;
//   gPad->SetLogy();
//   h_R9_allDA->GetXaxis()->SetRangeUser(0., 1.1);
//   h_R9_allDA->GetXaxis()->SetTitle("R9 ");
//   h_R9_allDA->SetLineColor(kRed+2);
//   h_R9_allDA->SetMarkerColor(kRed+2);
//   h_R9_allDA->SetMarkerStyle(7);
//   h_R9_allDA->Draw("e");
//    for(int jj = 0; jj < nBins; ++jj){
//      h_R9_MC[jj]->GetXaxis()->SetRangeUser(0., 1.1);
//      h_R9_MC[jj]->Draw("same");
//    }
//   TLegend *tspec = new TLegend(0.01,0.80,0.36,0.99);
//   tspec->SetFillColor(0);
//   tspec->SetTextFont(42); 
//   tspec->AddEntry(h_R9_allDA,"DATA","PL");
//   tspec->AddEntry(h_R9_MC[0],"MC ","PL");
//   tspec->Draw(); 
//   R9_spectrum->Print((PLOTS_folderName+"/"+folderName+"/R9_spectrum_"+string_year+".png").c_str(), ".png");

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

//   float y_min = y.at(0)-ey.at(ey.size()-1)-0.005;
//   float y_max = y.at(y.size()-1)+ey.at(ey.size()-1)+0.005;

//   // pad settings
//   TH1F *hPad = (TH1F*)gPad->DrawFrame(0,y_min,1.1,y_max);
//   hPad->GetXaxis()->SetTitle("R9");
//   hPad->GetYaxis()->SetTitle("E/p_{data}-E/p_{mc}");
//   hPad->GetYaxis()->SetTitleOffset(tYoffset);
//   hPad->GetXaxis()->SetTitleOffset(tXoffset);
//   hPad->GetXaxis()->SetLabelSize(labSize);
//   hPad->GetXaxis()->SetTitleSize(labSize);
//   hPad->GetYaxis()->SetLabelSize(labSize);
//   hPad->GetYaxis()->SetTitleSize(labSize);
//   finalGraph->Draw("P");
//   cplot->Print((PLOTS_folderName+"/"+folderName+"/EoP_vs_R9_"+string_year+".png").c_str(),".png");
//   */  
  
  std::cout << " plottato tutto " << std::endl;
  //std::cout << "CREATI I FILES" << std::endl;
  return (0);
}
