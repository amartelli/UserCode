// g++ -Wall -o trisCheckScale_Allcategories `root-config --cflags --glibs` ../Utils/setTDRStyle.cc ../Utils/ntupleUtils.cc ../Utils/stabilityUtils.cc ../Utils/ConvoluteTemplate.cc ../Utils/histoFunc.h ../Utils/TPileupReweighting.h trisCheckScale_Allcategories.cpp

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


  float totDAevts = 0;
  float DAevtsHIHI = 0;
  
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


  //MC 52X
  //stimate senza PU
//   TF1* Et_highR9_2012 = new TF1("Et_highR9_2012", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
//   Et_highR9_2012->SetParameters(1.71373322900473177e-02, 1.55744254105185699e-02,  -2.11477940336727904e-03);
//   TF1* Et_lowR9_2012 = new TF1("Et_lowR9_2012", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
//   Et_lowR9_2012->SetParameters(2.63075655765558566e-02, 4.57322846169432515e-02, -2.09413281975727485e-02);

  //stimate con PU
//   TF1* Et_highR9_2012 = new TF1("Et_highR9_2012", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
//   Et_highR9_2012->SetParameters(1.71373322900473177e-02, 1.55744254105185699e-02,  -2.11477940336727904e-03);
//   TF1* Et_lowR9_2012 = new TF1("Et_lowR9_2012", "[0] * (1 - exp(-[1] * x) ) +[2] ",0., 100.);
//   Et_lowR9_2012->SetParameters(1.69896128648113487e-02, 1.20797862827948261e-02, -5.86630884749932049e-03);


  //  bool UsePhotonRegression = false;
  bool UsePhotonRegression = true;

  //bool correctEt = false;
    bool correctEt = false;

    // bool reweightZtoH = true;
    bool reweightZtoH = false;

    //bool reweightEta = true;
    bool reweightEta = false;

    //bool reweightR9 = true;
    bool reweightR9 = false;

    //bool useShCorr = true;
    bool useShCorr = false;

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
  float absEtaMin = -1.;
  float absEtaMax = -1.;
  int IetaMin = -1;
  int IetaMax = -1;
  int IphiMin = -1;
  int IphiMax = -1;


  std::cout << "EBEE:         " << EBEE         << std::endl;
  std::cout << "LOWHIGH:      " << LOWHIGH       << std::endl;
  std::cout << "ENE:          " << ENE           << std::endl;
  std::cout << "PU:           " << PU            << std::endl;
  std::cout << "evtsPerPoint: " << evtsPerPoint  << std::endl;
  std::cout << "year:      " << year       << std::endl;
  std::cout << "doVsEach:      " << doVsEach       << std::endl;
  std::cout << "IetaMin: "      << IetaMin       << std::endl;
  std::cout << "IetaMax: "      << IetaMax       << std::endl;
  std::cout << "IphiMin: "      << IphiMin       << std::endl;
  std::cout << "IphiMax: "      << IphiMax       << std::endl;
  

  TPileupReweighting* puReweighting;

//   //2012 prompt           
   if(year == 2012) puReweighting =
 new TPileupReweighting("/afs/cern.ch/work/a/amartell/public/weights/PUweights_DYJetsToLL_Summer12_53X_ShSkim_ABC_TrueNumInteractions.root","hweights");
   //     new TPileupReweighting("/afs/cern.ch/work/a/amartell/public/weights/PUweights_DYJetsToLL_Summer12_Prompt_TrueNumInteractions.root","hweights");

//   //2011                                                                                                                                          
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
  std::string nameNtuplesMC = "simpleNtupleEoverP/SimpleNtupleEoverP";
  //  if(year == 2011) nameNtuples = "ntu";                                                                  
  //  if(year == 2011) nameNtuplesMC = "ntu";                                                                
  //  if(year == 2012) nameNtuplesMC = "simpleNtupleEoverPSh/SimpleNtupleEoverP";
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
  std::vector<float> scE_DA, scEt_reg_DA,scE_reg_DA, R9_DA, P_DA, EoP_DA, Et_DA, scEta_DA, elePhi_DA, ES_DA, isEB_DA, e3x3_DA,e5x5_DA, scERaw_DA;
  std::vector<float> scE_MC, scEt_reg_MC, scE_reg_MC, R9_MC, P_MC, EoP_MC, Et_MC, scEta_MC, elePhi_MC, ES_MC, isEB_MC, puRe, e3x3_MC, e5x5_MC, scERaw_MC;
  std::vector<float> scEtRaw_DA, scEt_3x3_DA, scEt_5x5_DA;
  std::vector<float> scEtRaw_MC, scEt_3x3_MC, scEt_5x5_MC;
  std::vector<float> ele1ele2_scM_DA, ele1ele2_scM_MC;
  std::vector<int> charge_DA, charge_MC;
  // Set branch addresses
  int isZ,runId,timeStamp,nVtx;
  float npu;
  
  ntu_DA->SetBranchStatus("*",0);
  ntu_DA->SetBranchStatus("runId",1);            ntu_DA->SetBranchAddress("runId", &runId);  
  ntu_DA->SetBranchStatus("timeStampHigh",1);    ntu_DA->SetBranchAddress("timeStampHigh", &timeStamp);  
  ntu_DA->SetBranchStatus("isZ",1);              ntu_DA->SetBranchAddress("isZ", &isZ);
  ntu_DA->SetBranchStatus("PV_n",1);             ntu_DA->SetBranchAddress("PV_n",&nVtx);

 
  ntu_MC->SetBranchStatus("*",0);
  ntu_MC->SetBranchStatus("PUit_TrueNumInteractions", 1);      ntu_MC->SetBranchAddress("PUit_TrueNumInteractions", &npu);
  ntu_MC->SetBranchStatus("runId",1);                          ntu_MC->SetBranchAddress("runId", &runId);  
  ntu_MC->SetBranchStatus("timeStampHigh",1);                  ntu_MC->SetBranchAddress("timeStampHigh", &timeStamp);  
  ntu_MC->SetBranchStatus("isZ",1);                            ntu_MC->SetBranchAddress("isZ", &isZ);
  ntu_MC->SetBranchStatus("PV_n",1);                           ntu_MC->SetBranchAddress("PV_n",&nVtx);


  
  // Electron data
  float scEne1, scEneReg1, R9_pho1, EoP1, scEt1, scEta1, elePhi1, ES1, P1, scERaw1, e3x31, e5x51, ele1ele2_scM;
  float scEne2, scEneReg2, R9_pho2, EoP2, scEt2, scEta2, elePhi2, ES2, P2, scERaw2, e3x32, e5x52;
  int isEB1,isEB2;
  int ele1_charge, ele2_charge; 

  ntu_DA->SetBranchStatus("ele1_scE", 1);       ntu_DA->SetBranchAddress("ele1_scE", &scEne1);
  ntu_DA->SetBranchStatus("ele1_scEt", 1);      ntu_DA->SetBranchAddress("ele1_scEt", &scEt1);
  ntu_DA->SetBranchStatus("ele1_scEta", 1);     ntu_DA->SetBranchAddress("ele1_scEta", &scEta1);
  ntu_DA->SetBranchStatus("ele1ele2_scM", 1);   ntu_DA->SetBranchAddress("ele1ele2_scM", &ele1ele2_scM);
  if(!UsePhotonRegression)  {
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
  ntu_MC->SetBranchStatus("ele1ele2_scM", 1);   ntu_MC->SetBranchAddress("ele1ele2_scM", &ele1ele2_scM);
  if(!UsePhotonRegression)  {
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

	++totDAevts;
	if(e3x31/scERaw1 > 0.94 && e3x32/scERaw2 > 0.94) ++DAevtsHIHI;

        run_DA.push_back(runId);

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


	charge_DA.push_back(ele1_charge);
	charge_DA.push_back(ele2_charge);

        time_DA.push_back(timeStamp);
        Z_DA.push_back(isZ);
	PV_DA.push_back(nVtx);    

	scE_DA.push_back(scEne1);
	scE_DA.push_back(scEne2);

	scE_reg_DA.push_back(scEneReg1*corrEtR9_1);
	scE_reg_DA.push_back(scEneReg2*corrEtR9_2);

	float Rt1 = sin(2*atan(exp(-scEta1)) );
	float Rt2 = sin(2*atan(exp(-scEta2)) );

	scEt_reg_DA.push_back(scEneReg1*Rt1*corrEtR9_1);
 	scEt_reg_DA.push_back(scEneReg2*Rt2*corrEtR9_2);

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

	scEtRaw_DA.push_back(scERaw1*scEt1/scEne1);
	scEtRaw_DA.push_back(scERaw2*scEt2/scEne2);
	scEt_3x3_DA.push_back(e3x31*scEt1/scEne1);
	scEt_3x3_DA.push_back(e3x32*scEt2/scEne2);
	scEt_5x5_DA.push_back(e5x51*scEt1/scEne1);
	scEt_5x5_DA.push_back(e5x52*scEt2/scEne2);
  }
  
  std::cout << std::endl;
  float ww = 1.;
  for(int ientry = 0; ientry < ntu_MC -> GetEntries(); ientry++)
  {
  	if( (ientry%100000 == 0) ) std::cout << "reading MC entry " << ientry << "\r" << std::flush;
        ntu_MC->GetEntry(ientry);  
  
	if(isZ == 0) continue;
	//	if(nVtx > 20) continue;


	float R9_ele1 = e3x31/scERaw1;
	if(year == 2012 && isEB1 == 1) R9_ele1 = 0.0010 + 1.0045 * e3x31/scERaw1;
	if(year == 2012 && isEB1 == 0) R9_ele1 = -0.0007 + 1.0086 * e3x31/scERaw1;
	if(year == 2011) R9_ele1 = 1.0035 * e3x31/scERaw1;
	float R9_ele2 = e3x32/scERaw2;
	if(year == 2012 && isEB2 == 1) R9_ele2 = 0.0010 + 1.0045 * e3x32/scERaw2;
	if(year == 2012 && isEB2 == 0) R9_ele2 = -0.0007 + 1.0086 * e3x32/scERaw2;
	if(year == 2011) R9_ele2 = 1.0035 * e3x32/scERaw2;



// 	float energySmearing1 = 1.;
// 	float energySmearing2 = 1.;

/*
	float energySmearing1 = gRandom->Gaus(1.,0.0075);
	float energySmearing2 = gRandom->Gaus(1.,0.0075);
	if(year == 2011 && std::string(LOWHIGH) == "LOW"){	
	  energySmearing1 = gRandom->Gaus(1.,0.008);
	  energySmearing2 = gRandom->Gaus(1.,0.008);
	}
	if(year == 2012 && std::string(LOWHIGH) == "LOW"){
	  energySmearing1 = gRandom->Gaus(1.,0.0075);
	  energySmearing2 = gRandom->Gaus(1.,0.0075);
	}
	if(year == 2012 && std::string(LOWHIGH) == "HIGH"){
	  energySmearing1 = gRandom->Gaus(1.,0.0075);
	  energySmearing2 = gRandom->Gaus(1.,0.0075);
	}
*/

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


 if(reweightZtoH == true && reweightEta == true && reweightR9 == true){
 if(year == 2011 && isEB1 == 0 )       energySmearing1 = gRandom->Gaus(1.,0.055); //energySmearing1 * 1.002;
 if(year == 2011 && isEB2 == 0 )       energySmearing2 = gRandom->Gaus(1.,0.055); //energySmearing2 * 1.002;
 if(year == 2011 && isEB1 == 0 && R9_ele1 > 0.94)       energySmearing1 = gRandom->Gaus(1.,0.06); //energySmearing1 * 1.002;
 if(year == 2011 && isEB2 == 0 && R9_ele2 > 0.94)       energySmearing2 = gRandom->Gaus(1.,0.06); //energySmearing2 * 1.002;

 if(year == 2012 && isEB1 == 0 )       energySmearing1 = gRandom->Gaus(1.,0.045); //energySmearing1 * 1.002;
 if(year == 2012 && isEB2 == 0 )       energySmearing2 = gRandom->Gaus(1.,0.045); //energySmearing2 * 1.002;

 if(year == 2012 && isEB1 == 1 && R9_ele1 > 0.94)       energySmearing1 = gRandom->Gaus(1.,0.0008); //energySmearing1 * 1.002;
 if(year == 2012 && isEB2 == 1 && R9_ele2 > 0.94)       energySmearing2 = gRandom->Gaus(1.,0.0008); //energySmearing2 * 1.002;
 if(year == 2012 && isEB1 == 1 && R9_ele1 < 0.94)       energySmearing1 = gRandom->Gaus(1.,0.00075); //energySmearing1 * 1.002;
 if(year == 2012 && isEB2 == 1 && R9_ele2 < 0.94)       energySmearing2 = gRandom->Gaus(1.,0.00075); //energySmearing2 * 1.002;
 }


 if(reweightZtoH == false && reweightEta == false && reweightR9 == false){
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
 }



//  if(isEB1 == 1 )       energySmearing1 = energySmearing1 * 0.996;
//  if(isEB2 == 1 )       energySmearing2 = energySmearing2 * 0.996;



	charge_MC.push_back(ele1_charge);
	charge_MC.push_back(ele2_charge);

	ww = puReweighting->GetWeight((int)npu);
        puRe.push_back(ww);
        run_MC.push_back(runId);
        time_MC.push_back(timeStamp);
        Z_MC.push_back(isZ);
	PV_MC.push_back(nVtx);    

	scE_MC.push_back(scEne1);
	scE_MC.push_back(scEne2);

	scE_reg_MC.push_back(scEneReg1 * energySmearing1);
	scE_reg_MC.push_back(scEneReg2 * energySmearing2);

	scEt_reg_MC.push_back(scEneReg1/scEne1*scEt1*energySmearing1);
 	scEt_reg_MC.push_back(scEneReg2/scEne2*scEt2*energySmearing2);


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

	scERaw_MC.push_back(scERaw1);
	scERaw_MC.push_back(scERaw2);

	R9_MC.push_back(R9_ele1);
	R9_MC.push_back(R9_ele2);

	scEtRaw_MC.push_back(scERaw1*scEt1/scEne1);
	scEtRaw_MC.push_back(scERaw2*scEt2/scEne2);

	scEt_5x5_MC.push_back(e5x51*scEt1/scEne1);
	scEt_5x5_MC.push_back(e5x52*scEt2/scEne2);

 	e3x3_MC.push_back(R9_ele1*scERaw1);
 	e3x3_MC.push_back(R9_ele2*scERaw2);
 	scEt_3x3_MC.push_back(R9_ele1*scERaw1*scEt1/scEne1);
 	scEt_3x3_MC.push_back(R9_ele2*scERaw2*scEt2/scEne2);
  }
  

  std::cout << "   totDAevts = " << totDAevts << std::endl;
  std::cout << "   DAevtsHIHI = " << DAevtsHIHI << std::endl;
  //  return 200;


  TH1F* h_Et_allDA = new TH1F("h_Et_allDA", "", 1000, 0., 1000.);
  TH1F* h_Et_allMC = new TH1F("h_Et_allMC", "", 1000, 0., 1000.);

  TH1F* h_Eta_allDA = new TH1F("h_Eta_allDA", "", 1000, -3., 3.);
  TH1F* h_Eta_allMC = new TH1F("h_Eta_allMC", "", 1000, -3., 3.);

  TH1F* h_R9_allDA = new TH1F("h_R9_allDA", "", 1000, 0., 1.2);
  TH1F* h_R9_allMC = new TH1F("h_R9_allMC", "", 1000, 0., 1.2);

   h_Et_allDA->Sumw2();
   h_Et_allMC->Sumw2();

   h_Eta_allDA->Sumw2();
   h_Eta_allMC->Sumw2();

   h_R9_allDA->Sumw2();
   h_R9_allMC->Sumw2();

   h_Et_allDA->SetLineColor(kRed+2);
   h_Et_allMC->SetLineColor(kGreen+2);

   h_Eta_allDA->SetLineColor(kRed+2);
   h_Eta_allMC->SetLineColor(kGreen+2);

   h_R9_allDA->SetLineColor(kRed+2);
   h_R9_allMC->SetLineColor(kGreen+2);


  TH1F* h_LHR9_Eta_allMC = new TH1F("h_LHR9_Eta_allMC", "", 1000, -3., 3.);
  TH1F* h_LHR9_Eta_allDA = new TH1F("h_LHR9_Eta_allDA", "", 1000, -3., 3.);

  h_LHR9_Eta_allMC->Sumw2();
  h_LHR9_Eta_allDA->Sumw2();
  h_LHR9_Eta_allDA->SetLineColor(kRed+2);
  h_LHR9_Eta_allMC->SetLineColor(kGreen+2);

  TH1F* h_allR9_Eta_allMC = new TH1F("h_allR9_Eta_allMC", "", 1000, -3., 3.);
  TH1F* h_allR9_Eta_allDA = new TH1F("h_allR9_Eta_allDA", "", 1000, -3., 3.);

  h_allR9_Eta_allMC->Sumw2();
  h_allR9_Eta_allDA->Sumw2();
  h_allR9_Eta_allDA->SetLineColor(kRed+2);
  h_allR9_Eta_allMC->SetLineColor(kGreen+2);


  std::string spectrumFileName = "HggGlobe_Stectra_weights_"+string_year+".root";
  TFile spectrumFile(spectrumFileName.c_str(), "read");  

  TH1F* hweightsMC_Et = (TH1F*)spectrumFile.Get("hweightsMC_Et");
  TH1F* hweightsDA_Et = (TH1F*)spectrumFile.Get("hweightsDa_Et");

   TH1F* hweightsMC_Eta = (TH1F*)spectrumFile.Get("hweightsMC_Eta");
   TH1F* hweightsDA_Eta = (TH1F*)spectrumFile.Get("hweightsDa_Eta");

//   TH1F* hweightsMC_EtaHR9 = (TH1F*)spectrumFile.Get("hweightsMC_Eta");
//   TH1F* hweightsDA_EtaHR9 = (TH1F*)spectrumFile.Get("hweightsDa_Eta");

//   TH1F* hweightsMC_EtaLR9 = (TH1F*)spectrumFile.Get("hweightsMC_EtaHR9");
//   TH1F* hweightsDA_EtaLR9 = (TH1F*)spectrumFile.Get("hweightsDa_EtaLR9");

  TH1F* hweightsMC_EtaR9fr = (TH1F*)spectrumFile.Get("hweightsMC_EtaR9fr");
  TH1F* hweightsDA_EtaR9fr = (TH1F*)spectrumFile.Get("hweightsDA_EtaR9fr");

  TH1F* hweightsMC_EtaLR9fr = (TH1F*)spectrumFile.Get("hweightsMC_EtaLR9fr");
  TH1F* hweightsDA_EtaLR9fr = (TH1F*)spectrumFile.Get("hweightsDA_EtaLR9fr");

  TH1F* hweightsMC_R9 = (TH1F*)spectrumFile.Get("hweightsMC_R9");
  TH1F* hweightsDA_R9 = (TH1F*)spectrumFile.Get("hweightsDa_R9");


  // Loop and sort events
  std::cout << std::endl;
  std::cout << "***** Sort events and define bins *****" << std::endl;
  
  int nEntries = scEt_reg_DA.size();
  int nSavePts = 0;
  std::vector<bool> isSavedEntries(nEntries);
  std::vector<SorterLC> sortedEntries;
  
  for(int ientry = 0; ientry < nEntries; ++ientry)
  {
    isSavedEntries.at(ientry) = false;

    if ( (fabs(scEta_DA.at(ientry)) > 1.4442 && fabs(scEta_DA.at(ientry)) < 1.566) || fabs(scEta_DA.at(ientry)) > 2.7 ) continue; // no cracks
    if(scEt_reg_DA.at(ientry) <  25.) continue;   
    //    if(R9_DA.at(ientry) < 0.7 ) continue;
    //////////// to reweight the distributions ////////////////////////////
    if(reweightR9 == false && reweightEta == false && reweightZtoH == false){
      h_Et_allDA->Fill(scEt_reg_DA.at(ientry));
      h_Eta_allDA->Fill(scEta_DA.at(ientry));
      h_R9_allDA->Fill(R9_DA.at(ientry));
    }

    float Reweight = 1.;
    float Etweight = 0.;
    float Etaweight = 0.;
    float R9weight = 0.;
    if(reweightZtoH){
      Etweight = hweightsDA_Et->GetBinContent(int(scEt_reg_DA.at(ientry) + 1.));
      //      std::cout << " Etweight = " << Etweight << std::endl;
      Reweight = Reweight * Etweight;
    }
    if(reweightEta){
      Etaweight =  hweightsDA_Eta->GetBinContent(int( (scEta_DA.at(ientry) + 3.)/0.006) + 1.);
      if(R9_DA.at(ientry) > 0.94)      Etaweight = Etaweight * hweightsDA_EtaR9fr->GetBinContent(int( (scEta_DA.at(ientry) + 3.)/0.006) + 1.);
      if(R9_DA.at(ientry) < 0.94)      Etaweight = Etaweight * hweightsDA_EtaLR9fr->GetBinContent(int( (scEta_DA.at(ientry) + 3.)/0.006) + 1.);
      //      std::cout << " Etaweight = " << Etaweight << std::endl;
      Reweight = Reweight * Etaweight;
    }
    if(reweightR9){
      R9weight = hweightsDA_R9->GetBinContent(int( (R9_DA.at(ientry) + 0.)/0.0012) + 1.);
      //      std::cout << " R9weight = " << R9weight << std::endl;
      Reweight = Reweight * R9weight;
    }

    if(reweightR9 == true || reweightEta == true || reweightZtoH == true){
      h_Et_allDA->Fill(scEt_reg_DA.at(ientry), Reweight);
      h_Eta_allDA->Fill(scEta_DA.at(ientry), Reweight);
      h_R9_allDA->Fill(R9_DA.at(ientry), Reweight);
    }

    if(R9_DA.at(ientry) > 0.94) h_LHR9_Eta_allDA->Fill(scEta_DA.at(ientry), Reweight);
    if(R9_DA.at(ientry) < 0.94) h_allR9_Eta_allDA->Fill(scEta_DA.at(ientry), Reweight);
    //////////// to reweight the distributions ////////////////////////////

    // save only what is needed for the analysis!!!
    if (strcmp(EBEE,"EE")==0 && (fabs(scEta_DA.at(ientry)) < 1.566 || fabs(scEta_DA.at(ientry)) > 2.7 )) continue; // endcap
    if (strcmp(EBEE,"EB")==0 && (fabs(scEta_DA.at(ientry)) > 1.4442 )) continue; // endcap

    if (strcmp(LOWHIGH,"LOW")==0 && (R9_DA.at(ientry) >= 0.94 || R9_DA.at(ientry) < 0.7) ) continue;
    if (strcmp(LOWHIGH,"HIGH")==0 && R9_DA.at(ientry) < 0.94 ) continue;

    if(doVsEach == "true" && strcmp(ENE,"scE_reg")==0 && scEt_reg_DA.at(ientry) <  25.) continue; 
    if(doVsEach == "true" && strcmp(ENE,"scE")==0 && Et_DA.at(ientry) <  25.) continue;
    if(doVsEach == "true" && strcmp(ENE,"scERaw")==0 && scEtRaw_DA.at(ientry) <  25.) continue;
    if(doVsEach == "true" && strcmp(ENE,"e5x5")==0 && scEt_5x5_DA.at(ientry) <  25.) continue;
    if(doVsEach == "true" && strcmp(ENE,"e3x3")==0 && scEt_3x3_DA.at(ientry) <  25.) continue;
    //    if(doVsEach == "false" && scEt_reg_DA.at(ientry) <  25.) continue;   
    if(doVsEach == "false" && scEt_reg_DA.at(ientry) <  25.) continue;   

    if( (absEtaMin != -1.) && (absEtaMax != -1.) )
    {
	if( (fabs(scEta_DA.at(ientry)) < absEtaMin) || (fabs(scEta_DA.at(ientry)) > absEtaMax) ) continue;
    }
    
    isSavedEntries.at(ientry) = true;

    SorterLC dummy;
    if(doVsEach == "true"){
      if(strcmp(ENE,"scE_reg")==0) dummy.laserCorr = scEt_reg_DA.at(ientry);
      if(strcmp(ENE,"scE")==0) dummy.laserCorr = Et_DA.at(ientry);
      if(strcmp(ENE,"scERaw")==0) dummy.laserCorr = scEtRaw_DA.at(ientry);
      if(strcmp(ENE,"e5x5")==0) dummy.laserCorr = scEt_5x5_DA.at(ientry);
      if(strcmp(ENE,"e3x3")==0) dummy.laserCorr = scEt_3x3_DA.at(ientry);
    }
    else dummy.laserCorr = scEt_reg_DA.at(ientry);

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
    }
  }
  binEntryMax.at(nBins) = nSavePts;
  
  std::cout << " fine : nBins = " << nBins << std::endl;
 
  TVirtualFitter::SetDefaultFitter("Fumili2");
  
  // histogram definition
  
  TH1F** h_EoP_DA = new TH1F*[nBins];
  TH1F** h_EoP_MC = new TH1F*[nBins];
  TH1F** h_Et = new TH1F*[nBins];
  TH1F** h_Et_MC = new TH1F*[nBins];

  TH1F* h_scE_DA = new TH1F("h_scE_DA", "", 4000, 0., 200.);
  TH1F* h_scReg_DA = new TH1F("h_scReg_DA", "", 4000, 0., 200.);
  TH1F* h_scRaw_DA = new TH1F("h_scRaw_DA", "", 4000, 0., 200.);

  TH1F* h_Vtx_DA = new TH1F("h_Vtx_DA", "", 200, 0., 200.);
  TH1F* h_Vtx_MC = new TH1F("h_Vtx_MC", "", 200, 0., 200.);


   h_scE_DA->Sumw2();
   h_scReg_DA->Sumw2();
   h_scRaw_DA->Sumw2();
   h_Vtx_DA->Sumw2();
   h_Vtx_MC->Sumw2();

   h_Vtx_DA->SetLineColor(kRed+2);
   h_Vtx_MC->SetLineColor(kGreen+2);

  std::vector<float> EtBinEdge;
  std::vector<float> xNorm_single;
  float xNorm_all;

  for(int i = 0; i < nBins; ++i)
  {
    char histoName[80];
    
    sprintf(histoName, "EoP_DA_%d", i);
    h_EoP_DA[i] = new TH1F(histoName, histoName, 900, 0., 3.);
    h_EoP_DA[i] -> SetFillColor(kRed+2);
    h_EoP_DA[i] -> SetFillStyle(3004);
    h_EoP_DA[i] -> SetMarkerStyle(7);
    h_EoP_DA[i] -> SetMarkerColor(kRed+2); 
    h_EoP_DA[i] -> SetLineColor(kRed+2); 
    
    sprintf(histoName, "EoP_MC_%d", i);
    h_EoP_MC[i] = new TH1F(histoName, histoName, 900, 0., 3.);
    h_EoP_MC[i] -> SetFillColor(kGreen+2);
    h_EoP_MC[i] -> SetFillStyle(3004);
    h_EoP_MC[i] -> SetMarkerStyle(7);
    h_EoP_MC[i] -> SetMarkerColor(kGreen+2);
    h_EoP_MC[i] -> SetLineColor(kGreen+2);
    
    sprintf(histoName, "Et_%d", i);
    h_Et[i] = new TH1F(histoName, histoName, 1000, 0, 1000);
    h_Et[i]->SetLineColor(kRed+2);

    sprintf(histoName, "Et_MC_%d", i);
    h_Et_MC[i] = new TH1F(histoName, histoName, 1000, 0, 1000);
    h_Et_MC[i]->SetLineColor(kGreen+2);

     h_EoP_DA[i]->Sumw2();
     h_EoP_MC[i]->Sumw2();
     h_Et[i]->Sumw2();
     h_Et_MC[i]->Sumw2();
  }
  
  std::cout << " scEt_reg_MC.size() = " << scEt_reg_MC.size() << std::endl;
  std::cout << " scEt_reg_DA.size() = " << scEt_reg_DA.size() << std::endl;


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

  int DAEntries = scEt_reg_DA.size();
  if(doVsEach == "true" && std::string(ENE) == "scE") DAEntries = Et_DA.size();
  if(doVsEach == "true" && std::string(ENE) == "scERaw") DAEntries = scEtRaw_DA.size();
  if(doVsEach == "true" && std::string(ENE) == "e5x5") DAEntries = scEt_5x5_DA.size();
  if(doVsEach == "true" && std::string(ENE) == "e3x3") DAEntries = scEt_3x3_DA.size();
  for(unsigned int ientry = 0; ientry < DAEntries; ++ientry)
    {
    if( (ientry%100000 == 0) ) std::cout << "reading entry " << ientry << std::endl;
      
    if( isSavedEntries.at(ientry) == false) continue;
     
     int iSaved = antiMap[ientry];
     int bin = -1;
    
    for(bin = 0; bin < nBins; ++bin)
      if( iSaved >= binEntryMax.at(bin) && iSaved < binEntryMax.at(bin+1) )
    break;

    float Reweight = 1.;
    float Etweight = 0.;
    float Etaweight = 0.;
    float R9weight = 0.;
    if(reweightZtoH){
      Etweight = hweightsDA_Et->GetBinContent(int(scEt_reg_DA.at(ientry) + 1.));
      Reweight = Reweight * Etweight;
    }
    if(reweightEta){
      Etaweight = hweightsDA_Eta->GetBinContent(int( (scEta_DA.at(ientry) + 3.)/0.006) + 1.);
      if(R9_DA.at(ientry) > 0.94)      Etaweight = Etaweight * hweightsDA_EtaR9fr->GetBinContent(int( (scEta_DA.at(ientry) + 3.)/0.006) + 1.);
      if(R9_DA.at(ientry) < 0.94)      Etaweight = Etaweight * hweightsDA_EtaLR9fr->GetBinContent(int( (scEta_DA.at(ientry) + 3.)/0.006) + 1.);
      Reweight = Reweight * Etaweight;
    }
    if(reweightR9){
      R9weight = hweightsDA_R9->GetBinContent(int( (R9_DA.at(ientry) + 0.)/0.0012) + 1.);
      Reweight = Reweight * R9weight;
    }

//     std::cout << " scEt_reg_DA.at(ientry) " << scEt_reg_DA.at(ientry) << std::endl;                                
//     std::cout << " bin =  " << int(scEt_reg_DA.at(ientry) + 1.) << std::endl;                                
//     std::cout << " Eweight DA " << Eweight << std::endl;                                

    if(strcmp(ENE,"scE_reg")==0){
      h_EoP_DA[bin]->Fill((scE_reg_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)), Reweight);
      if(doVsEach == "true"){
	h_Et[bin]->Fill(scEt_reg_DA.at(ientry), Reweight);
	//	if(reweightR9 == true || reweightEta == true || reweightZtoH == true) h_Et_allDA->Fill(scEt_reg_DA.at(ientry), Reweight);
      }
    }
    if(strcmp(ENE,"scE")==0){
    	h_EoP_DA[bin]->Fill((scE_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)));
	if(doVsEach == "true"){
	  h_Et[bin]->Fill(Et_DA.at(ientry));
	  h_Et_allDA->Fill(Et_DA.at(ientry));
	}
    }
    if(strcmp(ENE,"scERaw")==0){
    	h_EoP_DA[bin]->Fill((scERaw_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)));
	if(doVsEach == "true"){
	  h_Et[bin]->Fill(scEtRaw_DA.at(ientry));
	  h_Et_allDA->Fill(scEtRaw_DA.at(ientry));
	}
    }
    if(strcmp(ENE,"e5x5")==0){
    	h_EoP_DA[bin]->Fill((e5x5_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)));
	if(doVsEach == "true"){
	  h_Et[bin]->Fill(scEt_5x5_DA.at(ientry));
	  h_Et_allDA->Fill(scEt_5x5_DA.at(ientry));
	}
    }
    if(strcmp(ENE,"e3x3")==0){
    	h_EoP_DA[bin]->Fill((e3x3_DA.at(ientry)-ES_DA.at(ientry))/(P_DA.at(ientry)-ES_DA.at(ientry)));
	if(doVsEach == "true"){
	  h_Et[bin]->Fill(scEt_3x3_DA.at(ientry));
	  h_Et_allDA->Fill(scEt_3x3_DA.at(ientry));
	}
    }

    if(doVsEach == "false"){
      h_Et[bin]->Fill(scEt_reg_DA.at(ientry), Reweight);
      //      if(reweightR9 == true || reweightEta == true || reweightZtoH == true)  h_Et_allDA->Fill(scEt_reg_DA.at(ientry), Reweight); 
    }

//     if(reweightR9 == true || reweightEta == true || reweightZtoH == true){
//       h_Eta_allDA->Fill(scEta_DA.at(ientry), Reweight);
//       h_R9_allDA->Fill(R9_DA.at(ientry), Reweight);
//     }

    h_Vtx_DA->Fill(PV_DA.at(int(ientry/2)) );
    h_scE_DA->Fill(scE_DA.at(ientry) );
    h_scReg_DA->Fill(scE_reg_DA.at(ientry));
    h_scRaw_DA->Fill(scERaw_DA.at(ientry));
    }
  
  std::cout << "data fillati " << std::endl;
  
  for(int bin = 0; bin < nBins; bin++)
  {
    std::cout << "h_Et[bin]->GetEntries() =  " << h_Et[bin]->GetEntries() << std::endl;
    std::cout << "h_EoP_DA[bin]->GetEntries() =  " << h_EoP_DA[bin]->GetEntries() << std::endl;
    for(int i = 1; i < h_Et[bin]->GetNbinsX()+1; i++)
      {
	if(h_Et[bin]->GetBinContent(i) > 0) {
	  EtBinEdge.push_back(h_Et[bin]->GetBinCenter(i)-h_Et[bin]->GetBinWidth(i) );
	  break;
	}  
      }
  }

  std::cout << "prima di MC " << std::endl;

  int MCEntries = scEt_reg_MC.size();
  if(doVsEach == "true" && std::string(ENE) == "scE") MCEntries = Et_MC.size();
  if(doVsEach == "true" && std::string(ENE) == "scERaw") MCEntries = scEtRaw_MC.size();
  if(doVsEach == "true" && std::string(ENE) == "e5x5") MCEntries = scEt_5x5_MC.size();
  if(doVsEach == "true" && std::string(ENE) == "e3x3") MCEntries = scEt_3x3_MC.size();
  for(unsigned int ientry = 0; ientry < MCEntries; ++ientry)
    {    
    if( (ientry%100000 == 0) ) std::cout << "reading entry " << ientry << std::endl;

    if ( (fabs(scEta_MC.at(ientry)) > 1.4442 && fabs(scEta_MC.at(ientry)) < 1.566) || fabs(scEta_MC.at(ientry)) > 2.7 ) continue; // no cracks
    if(scEt_reg_MC.at(ientry) <  25.) continue;   
    //    if(R9_MC.at(ientry) < 0.7 ) continue;

    //       std::cout << " dopo i continue " << std::endl;

    //////////// to reweight the distributions ////////////////////////////
    if(reweightR9 == false && reweightEta == false && reweightZtoH == false){
      h_Et_allMC->Fill(scEt_reg_MC.at(ientry), puRe.at(int(ientry/2)));
      h_Eta_allMC->Fill(scEta_MC.at(ientry), puRe.at(int(ientry/2)));
      h_R9_allMC->Fill(R9_MC.at(ientry), puRe.at(int(ientry/2)));
    }
    //        std::cout << "prima di Reweight " << std::endl;
    float Reweight = 1.;
    float Etweight = 0.;
    float Etaweight = 0.;
    float R9weight = 0.;
    if(reweightZtoH){
      Etweight = hweightsMC_Et->GetBinContent(int(scEt_reg_MC.at(ientry) + 1.));
      Reweight = Reweight * Etweight;
    }
    //        std::cout << " Etweight =  " << Etweight << std::endl;
    if(reweightEta){
      Etaweight = hweightsMC_Eta->GetBinContent(int( (scEta_MC.at(ientry) + 3.)/0.006) + 1);
      if(R9_MC.at(ientry) > 0.94)      Etaweight = Etaweight * hweightsMC_EtaR9fr->GetBinContent(int( (scEta_MC.at(ientry) + 3.)/0.006) + 1);
      if(R9_MC.at(ientry) < 0.94)      Etaweight = Etaweight * hweightsMC_EtaLR9fr->GetBinContent(int( (scEta_MC.at(ientry) + 3.)/0.006) + 1);
      Reweight = Reweight * Etaweight;
    }
    //        std::cout << " Etaweight =  " << Etaweight << std::endl;
    if(reweightR9){
      //      std::cout << " (R9_MC.at(ientry) + 0.)/0.0012  =  " << (R9_MC.at(ientry) + 0.)/0.0012 << std::endl;
      R9weight = hweightsMC_R9->GetBinContent(int( (R9_MC.at(ientry) + 0.)/0.0012) + 1 );
      Reweight = Reweight * R9weight;
    }
    //        std::cout << " R9weight =  " << R9weight << std::endl;
    if(reweightR9 == true || reweightEta == true || reweightZtoH == true){
      h_Et_allMC->Fill(scEt_reg_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);
      h_Eta_allMC->Fill(scEta_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);
      h_R9_allMC->Fill(R9_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);
    }
    //        std::cout << " Reweight = " << Reweight << std::endl;
    if(R9_MC.at(ientry) > 0.94)    h_LHR9_Eta_allMC->Fill(scEta_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);
    if(R9_MC.at(ientry) < 0.94)    h_allR9_Eta_allMC->Fill(scEta_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);
    //////////// to reweight the distributions ////////////////////////////


    if (strcmp(EBEE,"EE")==0 && (fabs(scEta_MC.at(ientry)) < 1.566 || fabs(scEta_MC.at(ientry)) > 2.7 )) continue; // endcap
    if (strcmp(EBEE,"EB")==0 && (fabs(scEta_MC.at(ientry)) > 1.4442 )) continue; // endcap

    if (strcmp(LOWHIGH,"LOW")==0 && (R9_MC.at(ientry) >= 0.94 || R9_MC.at(ientry) < 0.7) ) continue;
    if (strcmp(LOWHIGH,"HIGH")==0 && R9_MC.at(ientry) < 0.94 ) continue;
    
    //    if(charge_MC.at(ientry) > 0) continue;

    if(doVsEach == "true" && strcmp(ENE,"scE_reg")==0 && scEt_reg_MC.at(ientry) <  25.) continue; 
    if(doVsEach == "true" && strcmp(ENE,"scE")==0 && Et_MC.at(ientry) <  25.) continue;
    if(doVsEach == "true" && strcmp(ENE,"scERaw")==0 && scEtRaw_MC.at(ientry) <  25.) continue;
    if(doVsEach == "true" && strcmp(ENE,"e5x5")==0 && scEt_5x5_MC.at(ientry) <  25.) continue;
    if(doVsEach == "true" && strcmp(ENE,"e3x3")==0 && scEt_3x3_MC.at(ientry) <  25.) continue;
    if(doVsEach == "false" && scEt_reg_MC.at(ientry) <  25.) continue;   

    if( (absEtaMin != -1.) && (absEtaMax != -1.) )
    {
	if( (fabs(scEta_MC.at(ientry)) < absEtaMin) || (fabs(scEta_MC.at(ientry)) > absEtaMax) ) continue;
    }
    

    for(unsigned int bin = 0; bin < EtBinEdge.size(); ++bin)
    {

      float referenceEt_MC = scEt_reg_MC.at(ientry);
      if(doVsEach == "true" && strcmp(ENE,"scE")==0) referenceEt_MC = Et_MC.at(ientry);
      if(doVsEach == "true" && strcmp(ENE,"scERaw")==0) referenceEt_MC = scEtRaw_MC.at(ientry);
      if(doVsEach == "true" && strcmp(ENE,"e3x3")==0) referenceEt_MC = scEt_3x3_MC.at(ientry);
      if(doVsEach == "true" && strcmp(ENE,"e5x5")==0) referenceEt_MC = scEt_5x5_MC.at(ientry);

      if( (bin != EtBinEdge.size()-1 && referenceEt_MC > EtBinEdge.at(bin) && referenceEt_MC < EtBinEdge.at(bin+1)) ||
	  (bin == EtBinEdge.size()-1 && referenceEt_MC > EtBinEdge.at(bin) ) )
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
	    h_Et_MC[int(bin)]->Fill(scEt_reg_MC.at(ientry));
	    break;
	  }
          if(PU == 1)
	    {
	      if(strcmp(ENE,"scE_reg")==0){
		  h_EoP_MC[(int)bin]->Fill((scE_reg_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)), 
					   puRe.at(int(ientry/2))*Reweight);
		if(doVsEach == "true")  h_Et_MC[int(bin)]->Fill(scEt_reg_MC.at(ientry),puRe.at(int(ientry/2))*Reweight);
	      }
	      if(strcmp(ENE,"scE")==0){
		  h_EoP_MC[(int)bin]->Fill((scE_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)),puRe.at(int(ientry/2)));
		if(doVsEach == "true") h_Et_MC[int(bin)]->Fill(Et_MC.at(ientry),puRe.at(int(ientry/2)));
	      }
	    if(strcmp(ENE,"scERaw")==0){
	      h_EoP_MC[(int)bin]->Fill((scERaw_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)),puRe.at(int(ientry/2)));
	      if(doVsEach == "true") h_Et_MC[int(bin)]->Fill(scEtRaw_MC.at(ientry), puRe.at(int(ientry/2)));
	    }
	    if(strcmp(ENE,"e5x5")==0){
	      h_EoP_MC[(int)bin]->Fill((e5x5_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)),puRe.at(int(ientry/2)));
	      if(doVsEach == "true") h_Et_MC[int(bin)]->Fill(scEt_5x5_MC.at(ientry), puRe.at(int(ientry/2)) );
	    }
	    if(strcmp(ENE,"e3x3")==0){
	      h_EoP_MC[(int)bin]->Fill((e3x3_MC.at(ientry)-ES_MC.at(ientry))/(P_MC.at(ientry)-ES_MC.at(ientry)),puRe.at(int(ientry/2)));
	      if(doVsEach == "true") h_Et_MC[int(bin)]->Fill(scEt_3x3_MC.at(ientry) , puRe.at(int(ientry/2)));
	    }
	    if(doVsEach == "false") h_Et_MC[int(bin)]->Fill(scEt_reg_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);
	    break;
	  } 
	}
    }
    if(PU == 0){
    	h_Et_allMC->Fill(scEt_reg_MC.at(ientry));
	h_Vtx_MC->Fill(PV_MC.at(int(ientry/2)));
    }
    if(PU == 1){
      if(doVsEach == "true"){
// 	if(strcmp(ENE,"scE_reg")==0)  {
// 	  if(reweightR9 == true || reweightEta == true || reweightZtoH == true)h_Et_allMC->Fill(scEt_reg_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);}
      if(strcmp(ENE,"scE")==0)     	h_Et_allMC->Fill(Et_MC.at(ientry), puRe.at(int(ientry/2)));
      if(strcmp(ENE,"scERaw")==0)     	h_Et_allMC->Fill(scEtRaw_MC.at(ientry), puRe.at(int(ientry/2)));
      if(strcmp(ENE,"e5x5")==0)     	h_Et_allMC->Fill(scEt_5x5_MC.at(ientry), puRe.at(int(ientry/2)));
      if(strcmp(ENE,"e3x3")==0)     	h_Et_allMC->Fill(scEt_3x3_MC.at(ientry), puRe.at(int(ientry/2)));
      }
      else { 
	//	if(reweightR9 == true || reweightEta == true || reweightZtoH == true)h_Et_allMC->Fill(scEt_reg_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);
      }
//       if(reweightR9 == true || reweightEta == true || reweightZtoH == true){
// 	h_Eta_allMC->Fill(scEta_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);
// 	h_R9_allMC->Fill(R9_MC.at(ientry), puRe.at(int(ientry/2))*Reweight);
//       }

      h_Vtx_MC->Fill(PV_MC.at(int(ientry/2)), puRe.at(int(ientry/2)) );
    }    
  }

  std::cout << " fino a qui ci sono " <<  std::endl;

  for(int i = 0; i < nBins; ++i)
  {
    //------------------------------------
    // Fill the graph for uncorrected data
    // define the fitting function
    // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )

    float xNorm_all = h_Et_allDA->Integral()/h_Et_allMC->Integral(); //* h_Et_allDA->GetBinWidth()/h_Et_allMC->GetBinWidth();
    h_Et_allMC->Scale(xNorm_all);
    h_Eta_allMC->Scale(xNorm_all);
    h_R9_allMC->Scale(xNorm_all);
    h_LHR9_Eta_allMC->Scale(xNorm_all);
    h_allR9_Eta_allMC->Scale(xNorm_all);

    float xNorm_allSelection = h_Vtx_DA->Integral()/h_Vtx_MC->Integral();
    h_Vtx_MC->Scale(xNorm_allSelection);

    
  if(reweightZtoH  == false && reweightEta == false && reweightR9 == false){
    if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "HIGH"){h_EoP_MC[i]->Smooth(1); }
    if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "LOW"){ h_EoP_MC[i]->Smooth(1); }
    if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "HIGH"){ h_EoP_MC[i]->Smooth(2); }
    if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "LOW"){h_EoP_MC[i]->Smooth(2); }

    /*
    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "HIGH")
      { h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "LOW")
      { h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
     if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "HIGH")
       {   h_EoP_MC[i]->Smooth(8); h_EoP_DA[i]->Rebin(3); h_EoP_MC[i]->Rebin(3);}
    if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "LOW")
      {h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
    */
    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "HIGH")
      { h_EoP_MC[i]->Smooth(150); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "LOW")
      { h_EoP_MC[i]->Smooth(150); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
     if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "HIGH")
       {   h_EoP_MC[i]->Smooth(150); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
    if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "LOW")
      {h_EoP_MC[i]->Smooth(150); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
  }

  if(reweightZtoH  == true && reweightEta == true && reweightR9 == true){
    if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "HIGH"){h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "LOW") { h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
    if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "HIGH"){ h_EoP_MC[i]->Smooth(6); }
    if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "LOW") {h_EoP_MC[i]->Smooth(2); h_EoP_DA[i]->Rebin(3); h_EoP_MC[i]->Rebin(3);}

    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "HIGH")
      { h_EoP_MC[i]->Smooth(300); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "LOW")
      { h_EoP_MC[i]->Smooth(300); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4); }
    if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "HIGH")
      { h_EoP_MC[i]->Smooth(300); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4); }
    if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "LOW")
      {h_EoP_MC[i]->Smooth(300); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4); }
  }


  /////////////////////
  /*

  if(reweightZtoH  == false && reweightEta == false && reweightR9 == false){
    if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "HIGH"){h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "LOW"){ h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "HIGH"){ h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "LOW"){h_EoP_MC[i]->Smooth(2); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}

    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "HIGH")
      { h_EoP_MC[i]->Smooth(4); h_EoP_MC[i]->Smooth(6); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "LOW")
      { h_EoP_MC[i]->Smooth(4); h_EoP_MC[i]->Smooth(6); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
    if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "HIGH")
      { h_EoP_MC[i]->Smooth(4); h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
    if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "LOW")
      {h_EoP_MC[i]->Smooth(2); h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}
  }

  if(reweightZtoH  == true && reweightEta == true && reweightR9 == true){
    if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "HIGH"){h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EB" && year == 2011 && std::string(LOWHIGH) == "LOW")
      { h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "HIGH"){ h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EB" && year == 2012 && std::string(LOWHIGH) == "LOW")
      {h_EoP_MC[i]->Smooth(2); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4);}

    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "HIGH")
      { h_EoP_MC[i]->Smooth(4); h_EoP_MC[i]->Smooth(6); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EE" && year == 2011 && std::string(LOWHIGH) == "LOW")
      { h_EoP_MC[i]->Smooth(4); h_EoP_MC[i]->Smooth(6); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "HIGH")
      { h_EoP_MC[i]->Smooth(4); h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
    if(std::string(EBEE) == "EE" && year == 2012 && std::string(LOWHIGH) == "LOW")
      {h_EoP_MC[i]->Smooth(2); h_EoP_MC[i]->Smooth(4); h_EoP_DA[i]->Rebin(4); h_EoP_MC[i]->Rebin(4); h_EoP_DA[i]->Rebin(2); h_EoP_MC[i]->Rebin(2);}
  }

   */
    

    float xNorm = h_EoP_DA[i]->Integral()/h_EoP_MC[i]->Integral() * h_EoP_DA[i]->GetBinWidth(1)/h_EoP_MC[i]->GetBinWidth(1);  
    float xNormEt = h_Et[i]->Integral()/h_Et_MC[i]->Integral(); //*h_Et[i]->GetBinWidth()/h_Et_MC[i]->GetBinWidth();  
    h_EoP_MC[i]->Scale(xNorm);
    h_Et_MC[i]->Scale(xNormEt);

//     std::cout << " i = " << i << " h_EoP_DA[i]->Integral() = " << h_EoP_DA[i]->Integral() << std::endl;
//     std::cout << " i = " << i << " h_Et[i]->Integral() = " << h_Et[i]->Integral() << std::endl;
//     std::cout << " i = " << i << " h_EoP_MC[i]->Integral() = " << h_EoP_MC[i]->Integral() << std::endl;
//     std::cout << " i = " << i << " h_Et_MC[i]->Integral() = " << h_Et_MC[i]->Integral() << std::endl;


//     std::cout << " i = " << i << " h_EoP_DA[i]->Integral() = " << h_EoP_DA[i]->Integral() << std::endl;
//     std::cout << " i = " << i << " h_EoP_MC[i]->Integral() = " << h_EoP_MC[i]->Integral() << std::endl;
//     std::cout << " xNorm = " << xNorm << std::endl;
//     std::cout << " xNormEt = " << xNormEt << std::endl;

   
    //    xNorm_single.push_back(xNormEt);


    histoFunc* templateHistoFunc = new histoFunc(h_EoP_MC[i]);
    char funcName[50];
    sprintf(funcName,"f_EoP_%d",i);
    f_EoP[i] = new TF1(funcName, templateHistoFunc, 0.7, 1.3, 3, "histoFunc");
    if(std::string(EBEE) == "EE")   f_EoP[i] = new TF1(funcName, templateHistoFunc, 0.5, 2.5, 3, "histoFunc");
    
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
      x.push_back(h_Et[i]->GetMean());
      ex.push_back((h_Et[i]->GetRMS())/sqrt(h_Et[i]->GetEntries()));
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
  finalGraph->GetYaxis()->SetTitle("E/p_{data} - E/p_{mc}");
  finalGraph->GetXaxis()->SetRangeUser(0., 130.);
  finalGraph->GetXaxis()->SetTitle("Et");


  std::cout << "   totDAevts = " << totDAevts << std::endl;
  std::cout << "   DAevtsHIHI = " << DAevtsHIHI << std::endl;


  std::string plotFolderName = "PLOTS_false";
  if(doVsEach == "true") plotFolderName = "PLOTS_true"; 

  TFile pippo((plotFolderName+"/results_"+folderName+"_"+string_year+".root").c_str(),"recreate");
  finalGraph->Write("finalGraph");
  h_Et_allMC->Write();
  h_Et_allDA->Write();
  h_Eta_allMC->Write();
  h_Eta_allDA->Write();
  h_R9_allMC->Write();
  h_R9_allDA->Write();
  h_Vtx_DA->Write();
  h_Vtx_MC->Write();
  h_scE_DA->Write();
  h_scReg_DA->Write();
  h_scRaw_DA->Write();
  h_LHR9_Eta_allMC->Write();
  h_LHR9_Eta_allDA->Write();
  h_allR9_Eta_allMC->Write();
  h_allR9_Eta_allDA->Write();

  for(int i = 0; i < nBins; ++i){
    h_EoP_DA[i]->Write();
    h_EoP_MC[i]->Write();
    h_Et[i]->Write();
    h_Et_MC[i]->Write();
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
//     c1[i]->cd();
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
//     if(PU == 0)      sprintf(Name, (plotFolderName+"/"+folderName+"/noPU_fit_%d_"+string_year+".png").c_str(),i);
//     if(PU == 1)      sprintf(Name, (plotFolderName+"/"+folderName+"/fit_%d_"+string_year+".png").c_str(),i);
//     c1[i] -> Print(Name,".png");     
//   }
  
//   TCanvas *c2[100]; 
//   for(int i = 0; i < nBins; ++i)
//   {
//     char canvasName[50];
//     sprintf(canvasName, "Et_DA-%0d", i); 
//     c2[i] = new TCanvas(canvasName, canvasName);
//     c2[i]->cd();

//     h_Et[i]->GetXaxis() -> SetTitle("Et");
//     h_Et[i]->GetYaxis()->SetRangeUser(0, std::max(h_Et[i]->GetMaximum(), h_Et_MC[i]->GetMaximum()) + 10. ); 
//     if(i<nBins-1) h_Et[i]->GetXaxis()->SetRangeUser(EtBinEdge.at(i), EtBinEdge.at(i+1)); 
//     else h_Et[i]->GetXaxis()->SetRangeUser(EtBinEdge.at(i), 150.); 
//     //    h_Et[i] -> Draw("e");
//     h_Et[i]->Draw();
//     h_Et_MC[i]->Draw("same");
//     /*gPad->Update();
//     s_Las[i]= (TPaveStats*)(h_Et[i]->GetListOfFunctions()->FindObject("stats"));
//     s_Las[i]->SetTextColor(kBlack);*/
    
//     char Name[100];
//     if(PU == 0)      sprintf(Name, (plotFolderName+"/"+folderName+"/noPU_Et_%d_"+string_year+".png").c_str(),i);
//     if(PU == 1)      sprintf(Name, (plotFolderName+"/"+folderName+"/Et_%d_"+string_year+".png").c_str(),i);
//     c2[i]->Print(Name,".png");      
//   }



//   TCanvas* Et_spectrum = new TCanvas;
//   gPad->SetLogy();
//   //  h_Et_allDA->GetYaxis()->SetRangeUser(0.1, 10000.);
//   h_Et_allDA->GetXaxis()->SetRangeUser(0., 150.);
//   h_Et_allDA->GetXaxis()->SetTitle("Et ");
//   h_Et_allDA->SetMarkerColor(kRed+2);
//   h_Et_allDA->SetMarkerStyle(7);
//   h_Et_allDA->Draw("e");
//    for(int jj = 0; jj < nBins; ++jj){
//      h_Et_MC[jj]->GetXaxis()->SetRangeUser(0., 150.);
//      h_Et_MC[jj]->Draw("same");
//    }
//   TLegend *tspec = new TLegend(0.64,0.80,0.99,0.99);
//   tspec->SetFillColor(0);
//   tspec->SetTextFont(42); 
//   tspec->AddEntry(h_Et_allDA,"DATA","PL");
//   tspec->AddEntry(h_Et_MC[0],"MC ","PL");
//   tspec->Draw(); 
//   Et_spectrum->Print((plotFolderName+"/"+folderName+"/Et_spectrum_"+string_year+".png").c_str(), ".png");


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
//   hPad->GetXaxis()->SetTitle("E_{T}");
//   hPad->GetYaxis()->SetTitle("E/p_{data}-E/p_{mc}");
//   hPad->GetYaxis()->SetTitleOffset(tYoffset);
//   hPad->GetXaxis()->SetTitleOffset(tXoffset);
//   hPad->GetXaxis()->SetLabelSize(labSize);
//   hPad->GetXaxis()->SetTitleSize(labSize);
//   hPad->GetYaxis()->SetLabelSize(labSize);
//   hPad->GetYaxis()->SetTitleSize(labSize);
//   finalGraph->Draw("P");
//   cplot->Print((plotFolderName+"/"+folderName+"/EoP_vs_Et_"+string_year+".png").c_str(),".png");
// */
  
  std::cout << " plottato tutto " << std::endl;
  //std::cout << "CREATI I FILES" << std::endl;
  return (0);

}
