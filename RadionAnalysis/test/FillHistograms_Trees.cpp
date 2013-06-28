
#include "PUreweightingUtils.h"
#include "ConfigParser.h"
#include "ParserUtils.h"
#include "setTDRStyle.h"
#include "drawPlotsUtils.h"

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
#include "TProfile.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THStack.h"
#include "TPad.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

int main(int argc, char** argv)
{

  // Set style options
  setTDRStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110); 
  //gStyle->SetOptStat(0000); 
  gStyle->SetOptFit(1); 

  // Input parameters

  //  bool printOk = true;
  bool printOk = false;
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  
  //Check if all nedeed arguments to parse are there
  if(argc <= 2)
  {
    std::cerr << ">>>>> Program usage: " << argv[0] << " configFileName btagWP" << std::endl ;
    return 1;
  }
 
  /// Parse the config file
  parseConfigFile (argv[1]) ;
  
  std::string inputList = gConfigParser -> readStringOption("Input::inputList");
  std::string inputTree = gConfigParser -> readStringOption("Input::inputTree");
  
  std::string outputName = gConfigParser -> readStringOption("Output::outputName");
  
  bool isReweighted = gConfigParser -> readBoolOption("Input::isReweighted");
  std::string inputMCPUHisto = gConfigParser -> readStringOption("Input::inputMCPUHisto");
  std::string inputDataPUHisto = gConfigParser -> readStringOption("Input::inputDataPUHisto");

  char* btagWP = argv[2];
   
  std::cout << "btagWP: " << btagWP << std::endl;

  if ( strcmp(btagWP,"No")!=0 && strcmp(btagWP,"Loose")!=0 && strcmp(btagWP,"Medium")!=0 && strcmp(btagWP,"Tight")!=0   )
  {
    std::cout << "CHK-STB Error: unknown  btag working point " << btagWP << std::endl;
    std::cout << "CHK-STB Select among No, Loose, Medium and Tight! " << std::endl;
    return -1;
  }

  float btag_wp = -1;
  
  if ( strcmp(btagWP,"No") == 0 ) btag_wp = 0.;
  if ( strcmp(btagWP,"Loose") == 0 ) btag_wp = 0.244;
  if ( strcmp(btagWP,"Medium") == 0 ) btag_wp = 0.679;
  if ( strcmp(btagWP,"Tight") == 0 ) btag_wp = 0.898;

  std::map<int,TChain*> ntu;
  std::map<int, TFile*> Files;

  std::map<int, std::map<int,float> > weight;
  char trees[500];
  FILE *f_trees;
  f_trees = fopen(inputList.c_str(),"r");
  int pos = 0;
  int pos_total = 0;
  
  while(fscanf(f_trees,"%s \n", trees) !=EOF ){
    std::string TREES = std::string(trees);
    if(TREES.find("#") != std::string::npos) continue;
    std::cout << "\nReading input: " << trees << std::endl;
    ntu[pos] = new TChain(inputTree.c_str());
    ntu[pos] -> Add(trees);   
    pos++;
    if(printOk)    std::cout << " >>>> pos = " << pos << std::endl;
  }
  pos_total = pos;
  
  for(int ii = 0; ii < pos_total; ii++){
      if(ntu[ii]->GetEntries() == 0 )
      {
         std::cout << "Error: input file" << ii << " is empty" << std::endl; 
         return -1;
      }
  }
  
  float         pu_n;
  float         weight_an;
  float         evweight;
  float         pu_weight;

  int           ph1_ciclevel;
  int           ph2_ciclevel;
  float         ph1_pfchargedisogood03;
  float         ph2_pfchargedisogood03;
  float         ph1_pfchargedisobad03;
  float         ph2_pfchargedisobad03;
  float         ph1_pfchargedisobad04;
  float         ph2_pfchargedisobad04;
  float         ph1_ecaliso;
  float         ph2_ecaliso;
  float         ph1_ecalisobad;
  float         ph2_ecalisobad;

  float         ph1_badvtx_Et;
  float         ph2_badvtx_Et;
                
  int           category;
  int           ph1_isEB;
  int           ph2_isEB;
  float         rho;

  float         ph1_e;
  float         ph2_e;
  float         ph1_pt;
  float         ph2_pt;
  float         ph1_phi;
  float         ph2_phi;
  float         ph1_eta;
  float         ph2_eta;
  float         ph1_r9;
  float         ph2_r9;
  float         ph1_hoe;
  float         ph2_hoe;
  float         ph1_sieie;
  float         ph2_sieie;
  float         ph1_sieip;
  float         ph2_sieip;
  float         ph1_sipip;
  float         ph2_sipip;
  float         ph1_isconv;
  float         ph2_isconv;
  float         PhotonsMass;
  float         dipho_E;
  float         dipho_pt;
  float         dipho_eta;
  float         dipho_phi;
  float         dipho_cosThetaStar_CS;
  float         dipho_tanhYStar;
  float         dipho_Y;
  int           vtx_ind;
  float         j1_e;
  float         j1_pt;
  float         j1_phi;
  float         j1_eta;
  float         j1_betaStarClassic;
  float         j1_csvBtag;
  float         j1_csvMvaBtag;
  float         j1_jetProbBtag;
  float         j1_tcheBtag;
  float         j1_radionMatched;
  float         j2_e;
  float         j2_pt;
  float         j2_phi;
  float         j2_eta;
  float         j2_betaStarClassic;
  float         j2_csvBtag;
  float         j2_csvMvaBtag;
  float         j2_jetProbBtag;
  float         j2_tcheBtag;
  float         j2_radionMatched;
  float         j3_e;
  float         j3_pt;
  float         j3_phi;
  float         j3_eta;
  float         j3_betaStarClassic;
  float         j3_csvBtag;
  float         j3_csvMvaBtag;
  float         j3_jetProbBtag;
  float         j3_tcheBtag;
  float         j3_radionMatched;
  float         j4_e;
  float         j4_pt;
  float         j4_phi;
  float         j4_eta;
  float         j4_betaStarClassic;
  float         j4_csvBtag;
  float         j4_csvMvaBtag;
  float         j4_jetProbBtag;
  float         j4_tcheBtag;
  float         j4_radionMatched;
  float         JetsMass;
  float         dijet_E;
//   float         dijet_pt;
//   float         dijet_eta;
//   float         dijet_phi;
  float         RadMass;
  float         radion_E;
//   float         radion_pt;
//   float         radion_eta;
//   float         radion_phi;
  
  for(int ii = 0; ii < pos_total; ii++){
      ntu[ii]->SetBranchStatus("*",0);
      ntu[ii]->SetBranchStatus("pu_n",1);
      ntu[ii]->SetBranchStatus("weight",1);
      ntu[ii]->SetBranchStatus("evweight",1);
      ntu[ii]->SetBranchStatus("pu_weight",1);

      ntu[ii]->SetBranchStatus("ph1_ciclevel",1);
      ntu[ii]->SetBranchStatus("ph2_ciclevel",1);
      ntu[ii]->SetBranchStatus("ph1_pfchargedisogood03",1);
      ntu[ii]->SetBranchStatus("ph2_pfchargedisogood03",1);
      ntu[ii]->SetBranchStatus("ph1_pfchargedisobad03",1);
      ntu[ii]->SetBranchStatus("ph2_pfchargedisobad03",1);
      ntu[ii]->SetBranchStatus("ph1_ecaliso",1);
      ntu[ii]->SetBranchStatus("ph2_ecaliso",1);
      ntu[ii]->SetBranchStatus("ph1_ecalisobad",1);
      ntu[ii]->SetBranchStatus("ph2_ecalisobad",1);

      ntu[ii]->SetBranchStatus("ph1_badvtx_Et",1);
      ntu[ii]->SetBranchStatus("ph2_badvtx_Et",1);
      ntu[ii]->SetBranchStatus("ph1_pfchargedisobad04",1);
      ntu[ii]->SetBranchStatus("ph2_pfchargedisobad04",1);
      ntu[ii]->SetBranchStatus("category",1);
      ntu[ii]->SetBranchStatus("ph1_isEB",1);
      ntu[ii]->SetBranchStatus("ph2_isEB",1);
      ntu[ii]->SetBranchStatus("rho",1);

      ntu[ii]->SetBranchStatus("ph1_e",1);
      ntu[ii]->SetBranchStatus("ph2_e",1);
      ntu[ii]->SetBranchStatus("ph1_pt",1);
      ntu[ii]->SetBranchStatus("ph2_pt",1);
      ntu[ii]->SetBranchStatus("ph1_phi",1);
      ntu[ii]->SetBranchStatus("ph2_phi",1);
      ntu[ii]->SetBranchStatus("ph1_eta",1);
      ntu[ii]->SetBranchStatus("ph2_eta",1);
      ntu[ii]->SetBranchStatus("ph1_r9",1);
      ntu[ii]->SetBranchStatus("ph2_r9",1);
      ntu[ii]->SetBranchStatus("ph1_hoe",1);
      ntu[ii]->SetBranchStatus("ph2_hoe",1);
      ntu[ii]->SetBranchStatus("ph1_sieie",1);
      ntu[ii]->SetBranchStatus("ph2_sieie",1);
      ntu[ii]->SetBranchStatus("ph1_sieip",1);
      ntu[ii]->SetBranchStatus("ph2_sieip",1);
      ntu[ii]->SetBranchStatus("ph1_sipip",1);
      ntu[ii]->SetBranchStatus("ph2_sipip",1);
      ntu[ii]->SetBranchStatus("PhotonsMass",1);
      ntu[ii]->SetBranchStatus("dipho_E",1);
      ntu[ii]->SetBranchStatus("dipho_pt",1);
      ntu[ii]->SetBranchStatus("dipho_eta",1);
      ntu[ii]->SetBranchStatus("dipho_phi",1);
      ntu[ii]->SetBranchStatus("dipho_cosThetaStar_CS",1);
      ntu[ii]->SetBranchStatus("dipho_tanhYStar",1);
      ntu[ii]->SetBranchStatus("dipho_Y",1);
      ntu[ii]->SetBranchStatus("vtx_ind",1);
      ntu[ii]->SetBranchStatus("j1_e",1);
      ntu[ii]->SetBranchStatus("j1_pt",1);
      ntu[ii]->SetBranchStatus("j1_phi",1);
      ntu[ii]->SetBranchStatus("j1_eta",1);
      ntu[ii]->SetBranchStatus("j1_betaStarClassic",1);
      ntu[ii]->SetBranchStatus("j1_csvBtag",1);
      ntu[ii]->SetBranchStatus("j1_csvMvaBtag",1);
      ntu[ii]->SetBranchStatus("j1_jetProbBtag",1);
      ntu[ii]->SetBranchStatus("j1_tcheBtag",1);
      ntu[ii]->SetBranchStatus("j1_radionMatched",1);
      ntu[ii]->SetBranchStatus("j2_e",1);
      ntu[ii]->SetBranchStatus("j2_pt",1);
      ntu[ii]->SetBranchStatus("j2_phi",1);
      ntu[ii]->SetBranchStatus("j2_eta",1);
      ntu[ii]->SetBranchStatus("j2_betaStarClassic",1);
      ntu[ii]->SetBranchStatus("j2_csvBtag",1);
      ntu[ii]->SetBranchStatus("j2_csvMvaBtag",1);
      ntu[ii]->SetBranchStatus("j2_jetProbBtag",1);
      ntu[ii]->SetBranchStatus("j2_tcheBtag",1);
      ntu[ii]->SetBranchStatus("j2_radionMatched",1);
      ntu[ii]->SetBranchStatus("j3_e",1);
      ntu[ii]->SetBranchStatus("j3_pt",1);
      ntu[ii]->SetBranchStatus("j3_phi",1);
      ntu[ii]->SetBranchStatus("j3_eta",1);
      ntu[ii]->SetBranchStatus("j3_betaStarClassic",1);
      ntu[ii]->SetBranchStatus("j3_csvBtag",1);
      ntu[ii]->SetBranchStatus("j3_csvMvaBtag",1);
      ntu[ii]->SetBranchStatus("j3_jetProbBtag",1);
      ntu[ii]->SetBranchStatus("j3_tcheBtag",1);
      ntu[ii]->SetBranchStatus("j3_radionMatched",1);
      ntu[ii]->SetBranchStatus("j4_e",1);
      ntu[ii]->SetBranchStatus("j4_pt",1);
      ntu[ii]->SetBranchStatus("j4_phi",1);
      ntu[ii]->SetBranchStatus("j4_eta",1);
      ntu[ii]->SetBranchStatus("j4_betaStarClassic",1);
      ntu[ii]->SetBranchStatus("j4_csvBtag",1);
      ntu[ii]->SetBranchStatus("j4_csvMvaBtag",1);
      ntu[ii]->SetBranchStatus("j4_jetProbBtag",1);
      ntu[ii]->SetBranchStatus("j4_tcheBtag",1);
      ntu[ii]->SetBranchStatus("j4_radionMatched",1);
      ntu[ii]->SetBranchStatus("JetsMass",1);
      ntu[ii]->SetBranchStatus("dijet_E",1);
//       ntu[ii]->SetBranchStatus("dijet_pt",1);
//       ntu[ii]->SetBranchStatus("dijet_phi",1);
//       ntu[ii]->SetBranchStatus("dijet_eta",1);
      ntu[ii]->SetBranchStatus("radion_E",1);
//       ntu[ii]->SetBranchStatus("radion_pt",1);
//       ntu[ii]->SetBranchStatus("radion_eta",1);
//       ntu[ii]->SetBranchStatus("radion_phi",1);
      ntu[ii]->SetBranchAddress("pu_n",&pu_n);
      ntu[ii]->SetBranchAddress("weight",&weight_an);
      ntu[ii]->SetBranchAddress("evweight",&evweight);
      ntu[ii]->SetBranchAddress("pu_weight",&pu_weight);

      ntu[ii]->SetBranchAddress("ph1_ciclevel",&ph1_ciclevel);
      ntu[ii]->SetBranchAddress("ph2_ciclevel",&ph2_ciclevel);
      ntu[ii]->SetBranchAddress("ph1_pfchargedisogood03",&ph1_pfchargedisogood03);
      ntu[ii]->SetBranchAddress("ph2_pfchargedisogood03",&ph2_pfchargedisogood03);
      ntu[ii]->SetBranchAddress("ph1_pfchargedisobad03",&ph1_pfchargedisobad03);
      ntu[ii]->SetBranchAddress("ph2_pfchargedisobad03",&ph2_pfchargedisobad03);
      ntu[ii]->SetBranchAddress("ph1_ecaliso",&ph1_ecaliso);
      ntu[ii]->SetBranchAddress("ph2_ecaliso",&ph2_ecaliso);
      ntu[ii]->SetBranchAddress("ph1_ecalisobad",&ph1_ecalisobad);
      ntu[ii]->SetBranchAddress("ph2_ecalisobad",&ph2_ecalisobad);

      ntu[ii]->SetBranchAddress("ph1_badvtx_Et",&ph1_badvtx_Et);
      ntu[ii]->SetBranchAddress("ph2_badvtx_Et",&ph2_badvtx_Et);
      ntu[ii]->SetBranchAddress("ph1_pfchargedisobad04",&ph1_pfchargedisobad04);
      ntu[ii]->SetBranchAddress("ph2_pfchargedisobad04",&ph2_pfchargedisobad04);
      ntu[ii]->SetBranchAddress("category",&category);
      ntu[ii]->SetBranchAddress("ph1_isEB",&ph1_isEB);
      ntu[ii]->SetBranchAddress("ph2_isEB",&ph2_isEB);
      ntu[ii]->SetBranchAddress("rho",&rho);

      ntu[ii]->SetBranchAddress("ph1_e",&ph1_e);
      ntu[ii]->SetBranchAddress("ph2_e",&ph2_e);
      ntu[ii]->SetBranchAddress("ph1_pt",&ph1_pt);
      ntu[ii]->SetBranchAddress("ph2_pt",&ph2_pt);
      ntu[ii]->SetBranchAddress("ph1_phi",&ph1_phi);
      ntu[ii]->SetBranchAddress("ph2_phi",&ph2_phi);
      ntu[ii]->SetBranchAddress("ph1_eta",&ph1_eta);
      ntu[ii]->SetBranchAddress("ph2_eta",&ph2_eta);
      ntu[ii]->SetBranchAddress("ph1_r9",&ph1_r9);
      ntu[ii]->SetBranchAddress("ph2_r9",&ph2_r9);
      ntu[ii]->SetBranchAddress("ph1_hoe",&ph1_hoe);
      ntu[ii]->SetBranchAddress("ph2_hoe",&ph2_hoe);
      ntu[ii]->SetBranchAddress("ph1_sieie",&ph1_sieie);
      ntu[ii]->SetBranchAddress("ph2_sieie",&ph2_sieie);
      ntu[ii]->SetBranchAddress("ph1_sieip",&ph1_sieip);
      ntu[ii]->SetBranchAddress("ph2_sieip",&ph2_sieip);
      ntu[ii]->SetBranchAddress("ph1_sipip",&ph1_sipip);
      ntu[ii]->SetBranchAddress("ph2_sipip",&ph2_sipip);
      ntu[ii]->SetBranchAddress("PhotonsMass",&PhotonsMass);
      ntu[ii]->SetBranchAddress("dipho_E",&dipho_E);
      ntu[ii]->SetBranchAddress("dipho_pt",&dipho_pt);
      ntu[ii]->SetBranchAddress("dipho_eta",&dipho_eta);
      ntu[ii]->SetBranchAddress("dipho_phi",&dipho_phi);
      ntu[ii]->SetBranchAddress("dipho_cosThetaStar_CS",&dipho_cosThetaStar_CS);
      ntu[ii]->SetBranchAddress("dipho_tanhYStar",&dipho_tanhYStar);
      ntu[ii]->SetBranchAddress("dipho_Y",&dipho_Y);
      ntu[ii]->SetBranchAddress("vtx_ind",&vtx_ind);
      ntu[ii]->SetBranchAddress("j1_e",&j1_e);
      ntu[ii]->SetBranchAddress("j1_pt",&j1_pt);
      ntu[ii]->SetBranchAddress("j1_phi",&j1_phi);
      ntu[ii]->SetBranchAddress("j1_eta",&j1_eta);
      ntu[ii]->SetBranchAddress("j1_betaStarClassic",&j1_betaStarClassic);
      ntu[ii]->SetBranchAddress("j1_csvBtag",&j1_csvBtag);
      ntu[ii]->SetBranchAddress("j1_csvMvaBtag",&j1_csvMvaBtag);
      ntu[ii]->SetBranchAddress("j1_jetProbBtag",&j1_jetProbBtag);
      ntu[ii]->SetBranchAddress("j1_tcheBtag",&j1_tcheBtag);
      ntu[ii]->SetBranchAddress("j1_radionMatched",&j1_radionMatched);
      ntu[ii]->SetBranchAddress("j2_e",&j2_e);
      ntu[ii]->SetBranchAddress("j2_pt",&j2_pt);
      ntu[ii]->SetBranchAddress("j2_phi",&j2_phi);
      ntu[ii]->SetBranchAddress("j2_eta",&j2_eta);
      ntu[ii]->SetBranchAddress("j2_betaStarClassic",&j2_betaStarClassic);
      ntu[ii]->SetBranchAddress("j2_csvBtag",&j2_csvBtag);
      ntu[ii]->SetBranchAddress("j2_csvMvaBtag",&j2_csvMvaBtag);
      ntu[ii]->SetBranchAddress("j2_jetProbBtag",&j2_jetProbBtag);
      ntu[ii]->SetBranchAddress("j2_tcheBtag",&j2_tcheBtag);
      ntu[ii]->SetBranchAddress("j2_radionMatched",&j2_radionMatched);
      ntu[ii]->SetBranchAddress("j3_e",&j3_e);
      ntu[ii]->SetBranchAddress("j3_pt",&j3_pt);
      ntu[ii]->SetBranchAddress("j3_phi",&j3_phi);
      ntu[ii]->SetBranchAddress("j3_eta",&j3_eta);
      ntu[ii]->SetBranchAddress("j3_betaStarClassic",&j3_betaStarClassic);
      ntu[ii]->SetBranchAddress("j3_csvBtag",&j3_csvBtag);
      ntu[ii]->SetBranchAddress("j3_csvMvaBtag",&j3_csvMvaBtag);
      ntu[ii]->SetBranchAddress("j3_jetProbBtag",&j3_jetProbBtag);
      ntu[ii]->SetBranchAddress("j3_tcheBtag",&j3_tcheBtag);
      ntu[ii]->SetBranchAddress("j3_radionMatched",&j3_radionMatched);
      ntu[ii]->SetBranchAddress("j4_e",&j4_e);
      ntu[ii]->SetBranchAddress("j4_pt",&j4_pt);
      ntu[ii]->SetBranchAddress("j4_phi",&j4_phi);
      ntu[ii]->SetBranchAddress("j4_eta",&j4_eta);
      ntu[ii]->SetBranchAddress("j4_betaStarClassic",&j4_betaStarClassic);
      ntu[ii]->SetBranchAddress("j4_csvBtag",&j4_csvBtag);
      ntu[ii]->SetBranchAddress("j4_csvMvaBtag",&j4_csvMvaBtag);
      ntu[ii]->SetBranchAddress("j4_jetProbBtag",&j4_jetProbBtag);
      ntu[ii]->SetBranchAddress("j4_tcheBtag",&j4_tcheBtag);
      ntu[ii]->SetBranchAddress("j4_radionMatched",&j4_radionMatched);
      ntu[ii]->SetBranchAddress("JetsMass",&JetsMass);
      ntu[ii]->SetBranchAddress("dijet_E",&dijet_E);
//       ntu[ii]->SetBranchAddress("dijet_pt",&dijet_pt);
//       ntu[ii]->SetBranchAddress("dijet_phi",&dijet_phi);
//       ntu[ii]->SetBranchAddress("dijet_eta",&dijet_eta);
      ntu[ii]->SetBranchAddress("radion_E",&radion_E);
//       ntu[ii]->SetBranchAddress("radion_pt",&radion_pt);
//       ntu[ii]->SetBranchAddress("radion_eta",&radion_eta);
//       ntu[ii]->SetBranchAddress("radion_phi",&radion_phi);
      
  }

  int nSTEP = 13;

  std::vector<TH1F*>  h_vtx_ind;
  std::vector<TH1F*>  h_PhoEcalIso;
  std::vector<TH1F*>  h_PhoChargedIso;
  std::vector<TH1F*>  h_PhoEcalIsoBad;
  std::vector<TH1F*>  h_PhoChargedIsoBad;
  std::vector<TH1F*>  h_ciclevel;

  std::vector<TH1F*>  h_PhoR9;
  std::vector<TH1F*>  h_PhoHoE;
  std::vector<TH1F*>  h_PhoSieie;
  std::vector<TH1F*>  h_PhoSipip;
  std::vector<TH1F*>  h_PhoSieip;
  std::vector<TH1F*>  h_PhoPt;
  std::vector<TH1F*>  h_PhoEta;
  std::vector<TH1F*>  h_PhoPhi;
  std::vector<TH1F*>  h_diPhoDeltaEta;
  std::vector<TH1F*>  h_diPhoDeltaPhi;
  std::vector<TH1F*>  h_diPhoDeltaR;
  std::vector<TH1F*>  h_diPhoInvMass;
  std::vector<TH1F*>  h_JetNum;
  std::vector<TH1F*>  h_JetNum_btagged;
  std::vector<TH1F*>  h_JetcsvBtag;
  std::vector<TH1F*>  h_JetcsvMvaBtag;
  std::vector<TH1F*>  h_JetProbBtag;
  std::vector<TH1F*>  h_JettcheBtag;
  std::vector<TH1F*>  h_JetPt;
  std::vector<TH1F*>  h_JetEta;
  std::vector<TH1F*>  h_JetPhi;
  std::vector<TH1F*>  h_diJetDeltaEta;
  std::vector<TH1F*>  h_diJetDeltaPhi;
  std::vector<TH1F*>  h_diJetDeltaR;
  std::vector<TH1F*>  h_diJetInvMass;
  std::vector<TH1F*>  h_diJetDeltaEta_maxpt;
  std::vector<TH1F*>  h_diJetDeltaPhi_maxpt;
  std::vector<TH1F*>  h_diJetDeltaR_maxpt;
  std::vector<TH1F*>  h_diJetInvMass_maxpt;
  std::vector<TH1F*>  h_PhoJetDeltaEta;
  std::vector<TH1F*>  h_PhoJetDeltaPhi;
  std::vector<TH1F*>  h_PhoJetDeltaR;
  std::vector<TH1F*>  h_PhoJetDeltaEta_maxpt;
  std::vector<TH1F*>  h_PhoJetDeltaPhi_maxpt;
  std::vector<TH1F*>  h_PhoJetDeltaR_maxpt;
  std::vector<TH1F*>  h_PhoJetDeltaEta_maxpt_btag;
  std::vector<TH1F*>  h_PhoJetDeltaPhi_maxpt_btag;
  std::vector<TH1F*>  h_PhoJetDeltaR_maxpt_btag;
  std::vector<TH1F*>  h_diPhodiJetDeltaEta;
  std::vector<TH1F*>  h_diPhodiJetDeltaPhi;
  std::vector<TH1F*>  h_diPhodiJetDeltaR;
  std::vector<TH1F*>  h_diPhodiJetInvMass;
  std::vector<TH1F*>  h_diPhodiJetInvMass_core;
  std::vector<TH1F*>  h_diPhodiJetInvMass_sidebands;
  std::vector<TH1F*>  h_diPhodiJetInvMass_maxpt;
  std::vector<TH1F*>  h_diPhodiJetInvMass_core_maxpt;
  std::vector<TH1F*>  h_diPhodiJetInvMass_sidebands_maxpt;
  std::vector<TH1F*>  h_diPhodiJetInvMass_btag;
  std::vector<TH1F*>  h_diPhodiJetInvMass_core_btag;
  std::vector<TH1F*>  h_diPhodiJetInvMass_sidebands_btag;
  std::vector<TH1F*>  h_diPhodiJetInvMass_bveto;
  std::vector<TH1F*>  h_diPhodiJetInvMass_core_bveto;
  std::vector<TH1F*>  h_diPhodiJetInvMass_sidebands_bveto;

  /////////////////////////////////////////////////////
  for(int nstep = 0; nstep<nSTEP; ++nstep){
    char histoName[50];
    sprintf(histoName, "_step%d", nstep);
    std::string stepName = std::string(histoName);
    h_vtx_ind.push_back(new TH1F(("h_vtx_ind"+stepName).c_str(), "", 40, 0., 40));
    h_PhoEcalIso.push_back(new TH1F(("h_PhoEcalIso"+stepName).c_str(),"",100,0.,20.));
    h_PhoChargedIso.push_back(new TH1F(("h_PhoChargedIso"+stepName).c_str(),"",100,0.,20.));
    h_PhoEcalIsoBad.push_back(new TH1F(("h_PhoEcalIsoBad"+stepName).c_str(),"",100,0.,20.));
    h_PhoChargedIsoBad.push_back(new TH1F(("h_PhoChargedIsoBad"+stepName).c_str(),"",100,0.,20.));
    h_ciclevel.push_back(new TH1F(("h_ciclevel"+stepName).c_str(),"",100,0.,20.));
    
    h_PhoR9.push_back(new TH1F(("h_PhoR9"+stepName).c_str(),"",50,0.2,1.));
    h_PhoHoE.push_back(new TH1F(("h_PhoHoE"+stepName).c_str(),"",50,0.,0.1));
    h_PhoSieie.push_back(new TH1F(("h_PhoSieie"+stepName).c_str(),"",50,0.,0.07));
    h_PhoSipip.push_back(new TH1F(("h_PhoSipip"+stepName).c_str(),"",50,0.,0.0045));
    h_PhoSieip.push_back(new TH1F(("h_PhoSieip"+stepName).c_str(),"",50,0.,0.0004));
    h_PhoPt.push_back(new TH1F(("h_PhoPt"+stepName).c_str(),"", 180,0.,180.));
    h_PhoEta.push_back(new TH1F(("h_PhoEta"+stepName).c_str(),"",700,-3.5,3.5));
    h_PhoPhi.push_back(new TH1F(("h_PhoPhi"+stepName).c_str(),"",700,-3.5,3.5));
    h_diPhoDeltaEta.push_back(new TH1F(("h_diPhoDeltaEta"+stepName).c_str(),"", 1000,-5.,5.));
    h_diPhoDeltaPhi.push_back(new TH1F(("h_diPhoDeltaPhi"+stepName).c_str(),"", 1000,-4.,4.));
    h_diPhoDeltaR.push_back(new TH1F(("h_diPhoDeltaR"+stepName).c_str(),"", 60,0.,6.));
    h_diPhoInvMass.push_back(new TH1F(("h_diPhoInvMass"+stepName).c_str(),"", 40,100.,180.));
    h_JetNum.push_back(new TH1F(("h_JetNum"+stepName).c_str(),"", 10,0.,10.));
    h_JetNum_btagged.push_back(new TH1F(("h_JetNum_btagged"+stepName).c_str(),"", 10,0.,10.));
    h_JetcsvBtag.push_back(new TH1F(("h_JetcsvBtag"+stepName).c_str(),"", 100,0.,1.));
    h_JetcsvMvaBtag.push_back(new TH1F(("h_JetcsvMvaBtag"+stepName).c_str(),"", 100,0.,1.));
    h_JetProbBtag.push_back(new TH1F(("h_JetProbBtag"+stepName).c_str(),"", 100,0.,1.));
    h_JettcheBtag.push_back(new TH1F(("h_JettcheBtag"+stepName).c_str(),"", 100,0.,1.));
    h_JetPt.push_back(new TH1F(("h_JetPt"+stepName).c_str(),"", 180,0.,180.));
    h_JetEta.push_back(new TH1F(("h_JetEta"+stepName).c_str(),"",800,-4,4));
    h_JetPhi.push_back(new TH1F(("h_JetPhi"+stepName).c_str(),"",700,-3.5,3.5));
    h_diJetDeltaEta.push_back(new TH1F(("h_diJetDeltaEta"+stepName).c_str(),"", 100,-7.,7.));
    h_diJetDeltaPhi.push_back(new TH1F(("h_diJetDeltaPhi"+stepName).c_str(),"", 100,-4.,4.));
    h_diJetDeltaR.push_back(new TH1F(("h_diJetDeltaR"+stepName).c_str(),"", 70,0.,7.));
    h_diJetInvMass.push_back(new TH1F(("h_diJetInvMass"+stepName).c_str(),"", 80,0.,800));
    h_diJetDeltaEta_maxpt.push_back(new TH1F(("h_diJetDeltaEta_maxpt"+stepName).c_str(),"", 100,-7.,7.));
    h_diJetDeltaPhi_maxpt.push_back(new TH1F(("h_diJetDeltaPhi_maxpt"+stepName).c_str(),"", 100,-4.,4.));
    h_diJetDeltaR_maxpt.push_back(new TH1F(("h_diJetDeltaR_maxpt"+stepName).c_str(),"",70,0.,7.));
    h_diJetInvMass_maxpt.push_back(new TH1F(("h_diJetInvMass_maxpt"+stepName).c_str(),"", 80,0.,800));
    h_PhoJetDeltaEta.push_back(new TH1F(("h_PhoJetDeltaEta"+stepName).c_str(),"",100,-5.,5.));
    h_PhoJetDeltaPhi.push_back(new TH1F(("h_PhoJetDeltaPhi"+stepName).c_str(),"",100,-4.,4.));
    h_PhoJetDeltaR.push_back(new TH1F(("h_PhoJetDeltaR"+stepName).c_str(),"",60,0.,6.));
    h_PhoJetDeltaEta_maxpt.push_back(new TH1F(("h_PhoJetDeltaEta_maxpt"+stepName).c_str(),"",100,-5.,5.));
    h_PhoJetDeltaPhi_maxpt.push_back(new TH1F(("h_PhoJetDeltaPhi_maxpt"+stepName).c_str(),"",100,-4.,4.));
    h_PhoJetDeltaR_maxpt.push_back(new TH1F(("h_PhoJetDeltaR_maxpt"+stepName).c_str(),"",60,0.,6.));
    h_PhoJetDeltaEta_maxpt_btag.push_back(new TH1F(("h_PhoJetDeltaEta_maxpt_btag"+stepName).c_str(),"",100,-5.,5.));
    h_PhoJetDeltaPhi_maxpt_btag.push_back(new TH1F(("h_PhoJetDeltaPhi_maxpt_btag"+stepName).c_str(),"",100,-4.,4.));
    h_PhoJetDeltaR_maxpt_btag.push_back(new TH1F(("h_PhoJetDeltaR_maxpt_btag"+stepName).c_str(),"",60,0.,6.));
    h_diPhodiJetDeltaEta.push_back(new TH1F(("h_diPhodiJetDeltaEta"+stepName).c_str(),"",120,-9.,9.));
    h_diPhodiJetDeltaPhi.push_back(new TH1F(("h_diPhodiJetDeltaPhi"+stepName).c_str(),"",100,-4.,4.));
    h_diPhodiJetDeltaR.push_back(new TH1F(("h_diPhodiJetDeltaR"+stepName).c_str(),"",90,0.,9.));
    h_diPhodiJetInvMass.push_back(new TH1F(("h_diPhodiJetInvMass"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_core.push_back(new TH1F(("h_diPhodiJetInvMass_core"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_sidebands.push_back(new TH1F(("h_diPhodiJetInvMass_sidebands"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_maxpt.push_back(new TH1F(("h_diPhodiJetInvMass_maxpt"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_core_maxpt.push_back(new TH1F(("h_diPhodiJetInvMass_core_maxpt"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_sidebands_maxpt.push_back(new TH1F(("h_diPhodiJetInvMass_sidebands_maxpt"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_btag.push_back(new TH1F(("h_diPhodiJetInvMass_btag"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_core_btag.push_back(new TH1F(("h_diPhodiJetInvMass_core_btag"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_sidebands_btag.push_back(new TH1F(("h_diPhodiJetInvMass_sidebands_btag"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_bveto.push_back(new TH1F(("h_diPhodiJetInvMass_bveto"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_core_bveto.push_back(new TH1F(("h_diPhodiJetInvMass_core_bveto"+stepName).c_str(),"",120,0.,1200));
    h_diPhodiJetInvMass_sidebands_bveto.push_back(new TH1F(("h_diPhodiJetInvMass_sidebands_bveto"+stepName).c_str(),"",120,0.,1200));
  }
  /////////////////////////////////////////////////////

  
  std::map<int,TLorentzVector*> phoP4;
  std::map<int,TLorentzVector*> jetP4; 
  std::map<int,float> jet_csvBtag;
  std::map<int,float> jet_csvMvaBtag;
  std::map<int,float> jet_jetProbBtag;
  std::map<int,float> jet_tcheBtag;
  TLorentzVector* sumP4;
  TLorentzVector* PhoSumP4;
  std::map<int,TLorentzVector*> JetSumP4;

  std::map<int,bool> bad_jet;
  
  std::map<int,bool> isBTagCouple;
  std::map<int,bool> isBVetoCouple;

  std::vector<float>* pho_pt;
  std::vector<float>* jet_pt;

  std::map<float,int> phoPt_map;
  std::map<float,int> jetPt_map;

  int jet_couple_num; 
  int numPassingEvents[nSTEP];
  for(int nLocStep = 0; nLocStep<nSTEP; ++nLocStep) numPassingEvents[nLocStep] = 0;

  for(int nn = 0; nn < pos_total; nn++){
      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){
          if(ientry%100000==0) std::cout<<"--- Reading file_" << nn << " entry = "<< ientry <<std::endl;
          ntu[nn]->GetEntry(ientry);

	  int Step_Counter = -1;
	  if(printOk) 	  std::cout << " Step_Counter = " << Step_Counter << std::endl;

	  if(printOk)	  std::cout << "--- Reading file_" << nn << " entry = " << ientry << std::endl;
	  //	  if(ientry > 10)return 10;

          float puRe = 1.;
          puRe = evweight;

	  //ALLEVTS (0)
	  ++numPassingEvents[Step_Counter];
	  ++Step_Counter;

          pho_pt = new std::vector<float>;
          jet_pt = new std::vector<float>;


	  //superTight
	  //EB 
	  bool diPhotonSuperTIGHT_passe = false;
	  if(ph1_ciclevel > 3 && ph2_ciclevel > 3) diPhotonSuperTIGHT_passe = true;

	  /*
	  bool Photon1SuperTIGHT_passe = false;
	  bool Photon2SuperTIGHT_passe = false;
	  //	  std::cout << " category = " << category << std::endl;
	  float val_isosumoet_1 = (ph1_pfchargedisogood03 + ph1_ecaliso + 2.5 - rho * 0.09) * 50. / ph1_pt;
	  float val_isosumoetbad_1 = (ph1_pfchargedisobad04 + ph1_ecalisobad + 2.5 - rho * 0.23) * 50. / ph1_badvtx_Et;
	  float val_trkisooet_1 = ph1_pfchargedisogood03 * 50. / ph1_pt;
	  float val_isosumoet_2 = (ph2_pfchargedisogood03 + ph2_ecaliso + 2.5 - rho * 0.09) * 50. / ph2_pt;
	  float val_isosumoetbad_2 = (ph2_pfchargedisobad04 + ph2_ecalisobad + 2.5 - rho * 0.23) * 50. / ph2_badvtx_Et;
	  float val_trkisooet_2 = ph2_pfchargedisogood03 * 50. / ph2_pt;
	  */

	  /*
	  if( category == 0 && // EBEB HR9 HR9
	      (ph1_isEB == 1 && ph1_r9 > 0.94 && 
	       ph1_sieie < 0.0125 && ph1_hoe < 0.141 && //ph1_ciclevel > 3 &&
	       val_isosumoet_1 < 6.3 && val_isosumoetbad_1 < 18.9 && val_trkisooet_1 < 4.5) &&
	      (ph2_isEB == 1 && ph2_r9 > 0.94 && 
	       ph2_sieie < 0.0125 && ph2_hoe < 0.141 && //ph2_ciclevel > 3 &&
	       val_isosumoet_2 < 6.3 && val_isosumoetbad_2 < 18.9 && val_trkisooet_2 < 4.5) ) diPhotonSuperTIGHT_passe = true;
	  else if( category == 1 && // EBEB !HH
		   (ph1_isEB == 1 && ph1_r9 > 0.33 && 
		    ph1_sieie < 0.0103 && ph1_hoe < 0.138 && //ph1_ciclevel > 3 &&
		    val_isosumoet_1 < 5.6 && val_isosumoetbad_1 < 8. && val_trkisooet_1 < 2.8 &&
		    ph2_isEB == 1 && ph2_r9 > 0.33 && 
		    ph2_sieie < 0.0103 && ph2_hoe < 0.138 && //ph2_ciclevel > 3 &&
		    val_isosumoet_2 < 5.6 && val_isosumoetbad_2 < 8. && val_trkisooet_2 < 2.8) ) diPhotonSuperTIGHT_passe = true;
	  else if( category == 2 && // !EBEB HH
		   (ph1_r9 > 0.94 && 
		    ph1_sieie < 0.029 && ph1_hoe < 0.12 && //ph1_ciclevel > 3 &&
		    val_isosumoet_1 < 5.8 && val_isosumoetbad_1 < 10. && val_trkisooet_1 < 4.  &&
		    ph2_r9 > 0.94 && 
		    ph2_sieie < 0.029 && ph2_hoe < 0.12 && //ph2_ciclevel > 3 &&
		    val_isosumoet_2 < 5.8 && val_isosumoetbad_2 < 10. && val_trkisooet_2 < 4.) ) diPhotonSuperTIGHT_passe = true;
	  else if( category == 3 && // !EBEB !HH
		   (ph1_r9 > 0.37  &&
		    ph1_sieie < 0.028 && ph1_hoe < 0.091 && ph1_ciclevel > 3 && 
		    val_isosumoet_1 < 5.1 &&  val_isosumoetbad_1 < 6.2 && val_trkisooet_1 < 1.62 && 
		    ph2_r9 > 0.37  && 
		    ph2_sieie < 0.028 && ph2_hoe < 0.091 && ph2_ciclevel > 3 && 
		    val_isosumoet_2 < 5.1 &&  val_isosumoetbad_2 < 6.2 && val_trkisooet_2 < 1.62) ) diPhotonSuperTIGHT_passe = true;
	  */
	  //	  if(Photon1SuperTIGHT_passe && Photon2SuperTIGHT_passe) diPhotonSuperTIGHT_passe = true;

	  if(printOk)     std::cout << "--- prima di diPhotonSuperTIGHT_passe = " << diPhotonSuperTIGHT_passe << std::endl;
	  if(diPhotonSuperTIGHT_passe == false) continue;

	  //CIC Hgg (1)
	  ++numPassingEvents[Step_Counter];
	  ++Step_Counter;

	  ///////////////////////////////////////////////////////////
	  TLorentzVector diPhoton1;
	  diPhoton1.SetPtEtaPhiE(ph1_pt,ph1_eta,ph1_phi,ph1_e);
	  TLorentzVector diPhoton2;
	  diPhoton2.SetPtEtaPhiE(ph2_pt,ph2_eta,ph2_phi,ph2_e);
	  TLorentzVector diPhoton = diPhoton1 + diPhoton2;	  
	  ///////////////////////////////////////////////////////////

	  if(printOk)     std::cout << "--- diPhoton.M() = " << diPhoton.M() << std::endl;

	  if(ph1_pt > ph2_pt && (ph2_pt < 25 || ph1_pt < diPhoton.M()/3.) ) continue;
	  if(ph2_pt > ph2_pt && (ph1_pt < 25 || ph2_pt < diPhoton.M()/3.) ) continue;

	  if(diPhoton.M() < 100. || diPhoton.M() > 180.) continue;


	  //CIC Hgg + DiMass cut 100-180 (2)
	  ++numPassingEvents[Step_Counter];
	  ++Step_Counter;

	  if(printOk) {
	  std::cout << " h_vtx_ind.size() = " << h_vtx_ind.size() << std::endl;
	  std::cout << " Step_Counter = " << Step_Counter << std::endl;
	  }
	  //	  if(printOk)     std::cout << "--- before new variables " << std::endl;

	  h_vtx_ind.at(Step_Counter)->Fill(vtx_ind, puRe);


	  //save photons
          phoP4[0] = new TLorentzVector();
	  phoP4[0]->SetPtEtaPhiE(ph1_pt,ph1_eta,ph1_phi,ph1_e);
          
          phoP4[1] = new TLorentzVector();
	  phoP4[1]->SetPtEtaPhiE(ph2_pt,ph2_eta,ph2_phi,ph2_e);
      
          for(int ii = 0; ii < 2; ii++){
	    pho_pt->push_back(phoP4[ii]->Pt());
	    phoPt_map[phoP4[ii]->Pt()] = ii;
          }
          
	  std::map<float,int>::reverse_iterator itPhoMap = phoPt_map.rbegin();
	  int itMaxPho,itNextToMaxPho;
	  itMaxPho = itPhoMap->second;
	  ++itPhoMap;
	  itNextToMaxPho = itPhoMap->second;

	  if(printOk){
	  std::cout << " >>> itMaxPho = " << itMaxPho << std::endl;
	  std::cout << " >>> itNextToMaxPho = " << itNextToMaxPho << std::endl;
	  
	  for(std::map<float,int>::iterator it=phoPt_map.begin(); it!=phoPt_map.end(); ++it)
	      std::cout << it->first << " => " << it->second << '\n';
	  }

	  h_PhoEcalIso.at(Step_Counter)->Fill(ph1_ecaliso/ph1_pt, puRe);
	  h_PhoEcalIso.at(Step_Counter)->Fill(ph2_ecaliso/ph2_pt, puRe);

	  h_PhoEcalIsoBad.at(Step_Counter)->Fill(ph1_pfchargedisogood03 * 50. / ph1_pt, puRe);
	  h_PhoEcalIsoBad.at(Step_Counter)->Fill(ph2_pfchargedisogood03 * 50. / ph2_pt, puRe);

	  h_PhoChargedIso.at(Step_Counter)->Fill( (ph1_pfchargedisogood03+ph1_ecaliso+2.5 - rho * 0.09)*50./ph1_pt, puRe);
	  h_PhoChargedIso.at(Step_Counter)->Fill( (ph2_pfchargedisogood03+ph2_ecaliso+2.5 - rho * 0.09)*50./ph2_pt, puRe);

	  h_PhoChargedIsoBad.at(Step_Counter)->Fill( (ph1_pfchargedisobad04+ph1_ecalisobad+2.5 - rho * 0.23) * 50. / ph1_badvtx_Et, puRe);
	  h_PhoChargedIsoBad.at(Step_Counter)->Fill( (ph2_pfchargedisobad04+ph2_ecalisobad+2.5 - rho * 0.23) * 50. / ph2_badvtx_Et, puRe);

	  h_ciclevel.at(Step_Counter)->Fill(ph1_ciclevel, puRe);
	  h_ciclevel.at(Step_Counter)->Fill(ph2_ciclevel, puRe);

	  h_PhoR9.at(Step_Counter)->Fill(ph1_r9,puRe);
          h_PhoR9.at(Step_Counter)->Fill(ph2_r9,puRe);

          h_PhoHoE.at(Step_Counter)->Fill(ph1_hoe,puRe);
          h_PhoHoE.at(Step_Counter)->Fill(ph2_hoe,puRe);
           
          h_PhoSieie.at(Step_Counter)->Fill(ph1_sieie,puRe);
          h_PhoSieie.at(Step_Counter)->Fill(ph2_sieie,puRe);

          h_PhoSipip.at(Step_Counter)->Fill(ph1_sipip,puRe);
          h_PhoSipip.at(Step_Counter)->Fill(ph2_sipip,puRe);

          h_PhoSieip.at(Step_Counter)->Fill(ph1_sieip,puRe);
          h_PhoSieip.at(Step_Counter)->Fill(ph1_sieip,puRe);
         
	  for(int ii = 0; ii < 2; ii++){
	    h_PhoPt.at(Step_Counter)->Fill(phoP4[ii]->Pt(),puRe); 
	    h_PhoEta.at(Step_Counter)->Fill(phoP4[ii]->Eta(),puRe);
	    h_PhoPhi.at(Step_Counter)->Fill(phoP4[ii]->Phi(),puRe);
	  }


          h_diPhoDeltaEta.at(Step_Counter)->Fill(ph1_eta-ph2_eta,puRe);
	  h_diPhoDeltaPhi.at(Step_Counter)->Fill(phoP4[0]->DeltaPhi(*phoP4[1]),puRe);
          h_diPhoDeltaR.at(Step_Counter)->Fill(phoP4[0]->DeltaR(*phoP4[1]),puRe); 
	  h_diPhoInvMass.at(Step_Counter)->Fill(PhotonsMass,puRe);

	  if(printOk)
	  std::cout << " >>>> prima di sumP4 " << std::endl;

          sumP4 = new TLorentzVector();

          *sumP4 = *phoP4[0]+*phoP4[1];

	  if(printOk)
	  std::cout << " >>>> dopo di sumP4 " << std::endl;

          PhoSumP4 = (TLorentzVector*)sumP4->Clone();
          delete sumP4;

	  if(printOk){
	  std::cout << " >>>> PhoSumP4 Pt " << PhoSumP4->Pt() << std::endl;
	  std::cout << " >>>> PhoSumP4 Eta " << PhoSumP4->Eta() << std::endl;
	  std::cout << " >>>> PhoSumP4 Phi " << PhoSumP4->Phi() << std::endl;
	  std::cout << " >>>> PhoSumP4 E " << PhoSumP4->E() << std::endl;
	  }

	  int n_bad_jets = 0;
          for(int ii = 0; ii < 4; ii++)
	    bad_jet[ii] = false;
          
          jetP4[0] = new TLorentzVector();
          if(j1_pt != 0) jetP4[0]->SetPtEtaPhiE(j1_pt,j1_eta,j1_phi,j1_e);
          else{
	    if(jet_csvMvaBtag[0] < 0) std::cout << "PROOOBLEMMMMM :D " << std::endl;
             bad_jet[0] = true;
             n_bad_jets++;
          }

          jetP4[1] = new TLorentzVector();
          if(j2_pt != 0) jetP4[1]->SetPtEtaPhiE(j2_pt,j2_eta,j2_phi,j2_e);
          else{
	    if(jet_csvMvaBtag[1] < 0) std::cout << "PROOOBLEMMMMM :D " << std::endl;
             bad_jet[1] = true;
             n_bad_jets++;
          }

          jetP4[2] = new TLorentzVector();
          if(j3_pt != 0) jetP4[2]->SetPtEtaPhiE(j3_pt,j3_eta,j3_phi,j3_e);
          else{
	    if(jet_csvMvaBtag[2] < 0) std::cout << "PROOOBLEMMMMM :D " << std::endl;
             bad_jet[2] = true;
             n_bad_jets++;
          }

          jetP4[3] = new TLorentzVector();
          if(j4_pt != 0) jetP4[3]->SetPtEtaPhiE(j4_pt,j4_eta,j4_phi,j4_e);
          else{
	    if(jet_csvMvaBtag[3] < 0) std::cout << "PROOOBLEMMMMM :D " << std::endl;
             bad_jet[3] = true;
             n_bad_jets++;
          }

	  //why Badder? tolto Ra
	  if(n_bad_jets > 2) continue;

	  //NJets >= 2 (3)
	  ++numPassingEvents[Step_Counter];
	  ++Step_Counter;

          jet_csvBtag[0] = j1_csvBtag;
          jet_csvMvaBtag[0] = j1_csvMvaBtag;
          jet_jetProbBtag[0] = j1_jetProbBtag;
          jet_tcheBtag[0] = j1_tcheBtag;

          jet_csvBtag[1] = j2_csvBtag;
          jet_csvMvaBtag[1] = j2_csvMvaBtag;
          jet_jetProbBtag[1] = j2_jetProbBtag;
          jet_tcheBtag[1] = j2_tcheBtag;

          jet_csvBtag[2] = j3_csvBtag;
          jet_csvMvaBtag[2] = j3_csvMvaBtag;
          jet_jetProbBtag[2] = j3_jetProbBtag;
          jet_tcheBtag[2] = j3_tcheBtag;

          jet_csvBtag[3] = j4_csvBtag;
          jet_csvMvaBtag[3] = j4_csvMvaBtag;
          jet_jetProbBtag[3] = j4_jetProbBtag;
          jet_tcheBtag[3] = j4_tcheBtag;
      
          for(int ii = 0; ii < 4; ii++)
              if(jet_csvBtag[ii] < 0) bad_jet[ii] = true;

          int jet_num = 0;
          for(int ii = 0; ii < 4; ii++)
              if(bad_jet[ii] == false) jet_num++;
	  //	  if(jet_num != 4 - n_bad_jets) std::cout << "PROOOBLEMMMMM :D " << std::endl;

	  //Why Badder? Tolto Ra
	  //          if(jet_num < 2) continue;

          h_JetNum.at(Step_Counter)->Fill(jet_num, puRe);

          int n_bjets = 0;
	  for(unsigned int ii = 0; ii < jet_csvBtag.size(); ii++){
              if(bad_jet[ii] == true) continue;
              if(jet_csvBtag[ii] > btag_wp ) n_bjets++;
          }

          if(n_bjets < 1) continue;

	  //NBJets >= 1 (4)
	  ++numPassingEvents[Step_Counter];
	  ++Step_Counter;

          h_JetNum_btagged.at(Step_Counter)->Fill(n_bjets, puRe);
                   
          for(int ii = 0; ii < 4; ii++){
	    if(bad_jet[ii] == true) continue;
	    if(printOk){
	  std::cout << " >>>> jetP4[ii] Pt " << jetP4[ii]->Pt() << std::endl;
	  std::cout << " >>>> jetP4[ii] Eta " << jetP4[ii]->Eta() << std::endl;
	  std::cout << " >>>> jetP4[ii] Phi " << jetP4[ii]->Phi() << std::endl;
	  std::cout << " >>>> jetP4[ii] E " << jetP4[ii]->E() << std::endl;
	    }

              h_JetPt.at(Step_Counter)->Fill(jetP4[ii]->Pt(),puRe); 
              h_JetEta.at(Step_Counter)->Fill(jetP4[ii]->Eta(),puRe);
              h_JetPhi.at(Step_Counter)->Fill(jetP4[ii]->Phi(),puRe);
              
              if(jet_csvBtag[ii] > 0) h_JetcsvBtag.at(Step_Counter)->Fill(jet_csvBtag[ii],puRe);
              if(jet_csvMvaBtag[ii] > 0) h_JetcsvMvaBtag.at(Step_Counter)->Fill(jet_csvMvaBtag[ii],puRe);
              if(jet_jetProbBtag[ii] > 0) h_JetProbBtag.at(Step_Counter)->Fill(jet_jetProbBtag[ii],puRe);
              if(jet_tcheBtag[ii] > 0) h_JettcheBtag.at(Step_Counter)->Fill(jet_tcheBtag[ii],puRe);

              jet_pt->push_back(jetP4[ii]->Pt());
              jetPt_map[jetP4[ii]->Pt()] = ii;
          }

          jet_couple_num = 0;
          
          for(int ii = 0; ii < 4; ii++){
              if(bad_jet[ii] == true) continue;
              for(int jj = 0; jj < 4; jj++){

                  if(bad_jet[jj] == true) continue;     
                  if(ii == jj) continue;
                  h_diJetDeltaEta.at(Step_Counter)->Fill(jetP4[ii]->Eta()-jetP4[jj]->Eta(), puRe);
                  h_diJetDeltaPhi.at(Step_Counter)->Fill(jetP4[ii]->DeltaPhi(*jetP4[jj]),puRe);
                  h_diJetDeltaR.at(Step_Counter)->Fill(jetP4[ii]->DeltaR(*jetP4[jj]),puRe); 
                  sumP4 = new TLorentzVector();
                  *sumP4 = *jetP4[ii]+*jetP4[jj];
                  
                  JetSumP4[jet_couple_num] = (TLorentzVector*)sumP4->Clone();
                  h_diJetInvMass.at(Step_Counter)->Fill(sumP4->M(),puRe);

                  isBTagCouple[jet_couple_num] = false;
                  isBVetoCouple[jet_couple_num] = false;

                  if(jet_csvBtag[ii] > btag_wp) isBTagCouple[jet_couple_num] = true;
                  if(jet_csvBtag[jj] > btag_wp) isBTagCouple[jet_couple_num] = true;
                  if(jet_csvBtag[ii] <= 1-btag_wp && jet_csvBtag[jj] <= btag_wp) isBVetoCouple[jet_couple_num] = true;
                  if(jet_csvBtag[ii] <= btag_wp && jet_csvBtag[jj] <= 1-btag_wp) isBVetoCouple[jet_couple_num] = true;
                   
                  delete sumP4;
              
                 jet_couple_num++;
              }
          }

//           std::sort(pho_pt->begin(),pho_pt->end());
//           std::sort(jet_pt->begin(),jet_pt->end());

	  if(printOk)	  std::cout << " >>>> jetPt_map.size() = " << jetPt_map.size() << std::endl;

	  int itMaxJet,itNextToMaxJet;
	  if(jetPt_map.size() > 1.){
	  std::map<float,int>::reverse_iterator itJetMap = jetPt_map.rbegin();
	  itMaxJet = itJetMap->second;
	  ++itJetMap;
	  itNextToMaxJet = itJetMap->second;
 
	  if(printOk){
	  std::cout << " >>> itMaxJet = " << itMaxJet << std::endl;
	  std::cout << " >>> itNextToMaxJet = " << itNextToMaxJet << std::endl;

	  for(std::map<float,int>::iterator it=jetPt_map.begin(); it!=jetPt_map.end(); ++it)
	      std::cout << it->first << " => " << it->second << '\n';
	  }
//           h_diJetDeltaEta_maxpt.at(Step_Counter)->Fill(jetP4[jetPt_map->at(jet_pt->at(jet_pt->size()-1))]->Eta()-jetP4[jetPt_map->at(jet_pt->at(jet_pt->size()-2))]->Eta(), puRe);
//           h_diJetDeltaPhi_maxpt.at(Step_Counter)->Fill(jetP4[jetPt_map->at(jet_pt->at(jet_pt->size()-1))]->DeltaPhi(*jetP4[jetPt_map->at(jet_pt->at(jet_pt->size()-2))]),puRe);
//           h_diJetDeltaR_maxpt.at(Step_Counter)->Fill(jetP4[jetPt_map->at(jet_pt->at(jet_pt->size()-1))]->DeltaR(*jetP4[jetPt_map->at(jet_pt->at(jet_pt->size()-2))]),puRe); 

          h_diJetDeltaEta_maxpt.at(Step_Counter)->Fill(jetP4[itMaxJet]->Eta()-jetP4[itNextToMaxJet]->Eta(), puRe);
          h_diJetDeltaPhi_maxpt.at(Step_Counter)->Fill(jetP4[itMaxJet]->DeltaPhi(*jetP4[itNextToMaxJet]), puRe);
          h_diJetDeltaR_maxpt.at(Step_Counter)->Fill(jetP4[itMaxJet]->DeltaR(*jetP4[itNextToMaxJet]), puRe); 

          sumP4 = new TLorentzVector();

          *sumP4 = *jetP4[itMaxJet] + *jetP4[itNextToMaxJet];
                  
          h_diJetInvMass_maxpt.at(Step_Counter)->Fill(sumP4->M(),puRe);

          delete sumP4;

          for(int ii = 0; ii < 2; ii++)
              for(int jj = 0; jj < 4; jj++){
                   
                if(bad_jet[jj] == true) continue;

                  h_PhoJetDeltaEta.at(Step_Counter)->Fill(phoP4[ii]->Eta()-jetP4[jj]->Eta(),puRe);
                  h_PhoJetDeltaPhi.at(Step_Counter)->Fill(phoP4[ii]->DeltaPhi(*jetP4[jj]),puRe);
                  h_PhoJetDeltaR.at(Step_Counter)->Fill(phoP4[ii]->DeltaR(*jetP4[jj]),puRe);  
              }
           
          h_PhoJetDeltaEta_maxpt.at(Step_Counter)->Fill(phoP4[itMaxPho]->Eta()-jetP4[itMaxJet]->Eta(),puRe);
          h_PhoJetDeltaPhi_maxpt.at(Step_Counter)->Fill(phoP4[itMaxPho]->DeltaPhi(*jetP4[itMaxJet]),puRe);
          h_PhoJetDeltaR_maxpt.at(Step_Counter)->Fill(phoP4[itMaxPho]->DeltaR(*jetP4[itMaxJet]),puRe);
           
          for(int ii = 0; ii < 4; ii++){
               
	    if(bad_jet[ii] == true) continue;
	    if(jet_csvBtag[ii] < btag_wp) continue;

              h_PhoJetDeltaEta_maxpt_btag.at(Step_Counter)->Fill(phoP4[itMaxPho]->Eta() - jetP4[ii]->Eta(),puRe);
              h_PhoJetDeltaPhi_maxpt_btag.at(Step_Counter)->Fill(phoP4[itMaxPho]->DeltaPhi(*jetP4[ii]),puRe);
              h_PhoJetDeltaR_maxpt_btag.at(Step_Counter)->Fill(phoP4[itMaxPho]->DeltaR(*jetP4[ii]),puRe);

          }

          for(unsigned int jj = 0; jj < JetSumP4.size(); jj++){
      
                  h_diPhodiJetDeltaEta.at(Step_Counter)->Fill(PhoSumP4->Eta()-JetSumP4[jj]->Eta(),puRe);
                  h_diPhodiJetDeltaPhi.at(Step_Counter)->Fill(PhoSumP4->DeltaPhi(*JetSumP4[jj]),puRe);
                  h_diPhodiJetDeltaR.at(Step_Counter)->Fill(PhoSumP4->DeltaR(*JetSumP4[jj]),puRe); 

                  sumP4 = new TLorentzVector();
                   
                  *sumP4 = *PhoSumP4+*JetSumP4[jj];

                  h_diPhodiJetInvMass.at(Step_Counter)->Fill(sumP4->M(),puRe);

                  if(PhoSumP4->M() >= 123 && PhoSumP4->M() <= 127 ) h_diPhodiJetInvMass_core.at(Step_Counter)->Fill(sumP4->M(),puRe);
                  if(PhoSumP4->M() < 123 || PhoSumP4->M() > 127 ) h_diPhodiJetInvMass_sidebands.at(Step_Counter)->Fill(sumP4->M(),puRe);

                  if(isBTagCouple[jj] == true){
                     h_diPhodiJetInvMass_btag.at(Step_Counter)->Fill(sumP4->M(),puRe);
                     if(PhoSumP4->M() >= 123 && PhoSumP4->M() <= 127 ) h_diPhodiJetInvMass_core_btag.at(Step_Counter)->Fill(sumP4->M(),puRe);
                     if(PhoSumP4->M() < 123 || PhoSumP4->M() > 127 ) h_diPhodiJetInvMass_sidebands_btag.at(Step_Counter)->Fill(sumP4->M(),puRe);
                  }
                  
                  if(isBVetoCouple[jj] == true){
                     h_diPhodiJetInvMass_bveto.at(Step_Counter)->Fill(sumP4->M(),puRe);
                     if(PhoSumP4->M() >= 123 && PhoSumP4->M() <= 127 ) h_diPhodiJetInvMass_core_bveto.at(Step_Counter)->Fill(sumP4->M(),puRe);
                     if(PhoSumP4->M() < 123 || PhoSumP4->M() > 127 ) h_diPhodiJetInvMass_sidebands_bveto.at(Step_Counter)->Fill(sumP4->M(),puRe);
                  }
                  
                  delete sumP4;
          }

          sumP4 = new TLorentzVector();
           
          *sumP4 = *jetP4[itMaxJet] + *jetP4[itNextToMaxJet] + *PhoSumP4;

          h_diPhodiJetInvMass_maxpt.at(Step_Counter)->Fill(sumP4->M(),puRe);

          if(PhoSumP4->M() >= 123 && PhoSumP4->M() <= 127 ) h_diPhodiJetInvMass_core_maxpt.at(Step_Counter)->Fill(sumP4->M(),puRe);
          if(PhoSumP4->M() < 123 || PhoSumP4->M() > 127 ) h_diPhodiJetInvMass_sidebands_maxpt.at(Step_Counter)->Fill(sumP4->M(),puRe);
	  
          delete sumP4;
          }

          for(int ii = 0; ii < 2; ii++)
              delete phoP4[ii];

          for(int ii = 0; ii < 4; ii++)
              delete jetP4[ii];

          for(int ii = 0; ii < 4; ii++)
              bad_jet[ii] = false;

          delete pho_pt; 
          delete jet_pt; 
//           delete phoPt_map;
//           delete jetPt_map;
         
	  phoPt_map.clear();
	  jetPt_map.clear();              
      }
  }

  for(int nLocStep=0; nLocStep<nSTEP; ++ nLocStep)
    std::cout << " >>> numPassingEvents[" << nLocStep << "] = " << numPassingEvents[nLocStep] << std::endl;
 
  TFile* output = new TFile(outputName.c_str(),"RECREATE");
  output->cd();

  for(int nLocStep=0; nLocStep<nSTEP; ++ nLocStep){

  h_vtx_ind.at(nLocStep)->Write();

  h_PhoEcalIso.at(nLocStep)->Write();
  h_PhoEcalIsoBad.at(nLocStep)->Write();
  h_PhoChargedIso.at(nLocStep)->Write();
  h_PhoChargedIsoBad.at(nLocStep)->Write();
  h_ciclevel.at(nLocStep)->Write();

  h_PhoR9.at(nLocStep)->Write();
  h_PhoHoE.at(nLocStep)->Write();
  h_PhoSieie.at(nLocStep)->Write();
  h_PhoSipip.at(nLocStep)->Write();
  h_PhoSieip.at(nLocStep)->Write();
  h_PhoPt.at(nLocStep)->Write();
  h_PhoEta.at(nLocStep)->Write();
  h_PhoPhi.at(nLocStep)->Write();
  h_diPhoDeltaEta.at(nLocStep)->Write();
  h_diPhoDeltaPhi.at(nLocStep)->Write();
  h_diPhoDeltaR.at(nLocStep)->Write();
  h_diPhoInvMass.at(nLocStep)->Write();
  h_JetNum.at(nLocStep)->Write();
  h_JetNum_btagged.at(nLocStep)->Write();
  h_JetcsvBtag.at(nLocStep)->Write();
  h_JetcsvMvaBtag.at(nLocStep)->Write();
  h_JetProbBtag.at(nLocStep)->Write();
  h_JettcheBtag.at(nLocStep)->Write();
  h_JetPt.at(nLocStep)->Write();
  h_JetEta.at(nLocStep)->Write();
  h_JetPhi.at(nLocStep)->Write();
  h_diJetDeltaEta.at(nLocStep)->Write();
  h_diJetDeltaPhi.at(nLocStep)->Write();
  h_diJetDeltaR.at(nLocStep)->Write();
  h_diJetInvMass.at(nLocStep)->Write();
  h_diJetDeltaEta_maxpt.at(nLocStep)->Write();
  h_diJetDeltaPhi_maxpt.at(nLocStep)->Write();
  h_diJetDeltaR_maxpt.at(nLocStep)->Write();
  h_diJetInvMass_maxpt.at(nLocStep)->Write();
  h_PhoJetDeltaEta.at(nLocStep)->Write();
  h_PhoJetDeltaPhi.at(nLocStep)->Write();
  h_PhoJetDeltaR.at(nLocStep)->Write();
  h_PhoJetDeltaEta_maxpt.at(nLocStep)->Write();
  h_PhoJetDeltaPhi_maxpt.at(nLocStep)->Write();
  h_PhoJetDeltaR_maxpt.at(nLocStep)->Write(); 
  h_PhoJetDeltaEta_maxpt_btag.at(nLocStep)->Write();
  h_PhoJetDeltaPhi_maxpt_btag.at(nLocStep)->Write();
  h_PhoJetDeltaR_maxpt_btag.at(nLocStep)->Write();
  h_diPhodiJetDeltaEta.at(nLocStep)->Write();
  h_diPhodiJetDeltaPhi.at(nLocStep)->Write();
  h_diPhodiJetDeltaR.at(nLocStep)->Write();
  h_diPhodiJetInvMass.at(nLocStep)->Write();
  h_diPhodiJetInvMass_core.at(nLocStep)->Write();
  h_diPhodiJetInvMass_sidebands.at(nLocStep)->Write();
  h_diPhodiJetInvMass_maxpt.at(nLocStep)->Write();
  h_diPhodiJetInvMass_core_maxpt.at(nLocStep)->Write();
  h_diPhodiJetInvMass_sidebands_maxpt.at(nLocStep)->Write();
  h_diPhodiJetInvMass_btag.at(nLocStep)->Write();
  h_diPhodiJetInvMass_core_btag.at(nLocStep)->Write();
  h_diPhodiJetInvMass_sidebands_btag.at(nLocStep)->Write();
  h_diPhodiJetInvMass_bveto.at(nLocStep)->Write();
  h_diPhodiJetInvMass_core_bveto.at(nLocStep)->Write();
  h_diPhodiJetInvMass_sidebands_bveto.at(nLocStep)->Write();
  }
  output->Close();
  
}


















