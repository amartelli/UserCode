//g++ -Wall -o trisCheckScaleMZ_futyan `root-config --cflags --glibs` -lMinuit2 -L/gwteraw/cmssw/slc5_amd64_gcc462/external/gcc/4.6.2/lib64/ ../Utils/setTDRStyle.cc ../Utils/ntupleUtils.cc ../Utils/stabilityUtils.cc ../Utils/ConvoluteTemplate.cc ../Utils/histoFunc.h ../Utils/TPileupReweighting.h trisCheckScaleMZ_futyan.cpp

#include "../Utils/setTDRStyle.h"
#include "../Utils/histoFunc.h"
#include "../Utils/ConvoluteTemplate.h"
#include "../Utils/ntupleUtils.h"
#include "../Utils/stabilityUtils.h"
#include "../Utils/TPileupReweighting.h"
#include "../Utils/GetShervinCorrections.h"
#include "../Utils/GetSmearings.h"

#include "TROOT.h"
#include "TSystem.h"
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
#include "TFitterMinuit.h"
#include "TArrow.h"
#include "TDirectory.h"
#include "TLorentzVector.h"

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

std::vector<double>* EtBinEdges;
std::vector<double>* HtBinEdges;
unsigned int nEtBins;
unsigned int nHtBins;

TH1F** h_Ht_HtBin_MC;

TH1F** h_Et_EtBin_MC;
TH1F** h_Et_EtBin_DA;
TH1F** h_Et_EtBin_fit_DA;
TH1F** h_Et_EtBin_gausFit_DA;
TH1F** h_Et_EtBin_mean_DA;
TH1F** h_Et_EtBin_recursiveMean_DA;
TH1F** h_Et_EtBin_smallestInterval_DA;

TFile* outFile;

TF1* f_scaleVsEt;
TF1* f_invScaleVsEt;

double invScaleVsEt(double* x, double* par);

void MySetBins(TProfile* p, std::vector<double>& bins);

int MyFindBin(const double& val, const std::vector<double>* binEdges);
int MyFindBin(const double& val, const double& min, const double& max, const double& invWidth);

double MyEval(TGraph* g, const double& x);

void MyFindFit(double& scale, double& scaleErr,
               TH1F* h_MC, TH1F* h_DA);
void MyFindGausFit(double& mean, double& meanErr,
                   std::vector<double>& vals, std::vector<double>& weights,
                   TF1** f_gausFit, const std::string& name);
void MyFindMean(double& mean, double& meanErr,
                std::vector<double>& vals, std::vector<double>& weights);
void MyFindRecursiveMean(double& mean, double& meanErr,
                        std::vector<double>& vals, std::vector<double>& weights,
                        const double& window, const double& tolerance);
void MyFindSmallestInterval(double& mean, double& meanErr, double& min, double& max,
                            std::vector<double>& vals, std::vector<double>& weights,
                            const double& fraction);

double deltaPhi(const double& phi1, const double& phi2);


// MCClosure
bool MCClosure = false;

// Settings for corrections
bool applyEnergySmearings = true;
bool applyEnergyScaleCorr = true;
bool applyEtaR9Reweighting = false;

// bin definition
int nBins_mee = 300;
double meeMin = 0.;
double meeMax = 3.;

int nBins_Ht = 5000;
double HtMin = 0.;
double HtMax = 1000.;

int nSteps = 1;





int main(int argc, char** argv)
{
  // Set style options
  setTDRStyle();
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  
  
  // Fitting functions
  //f_scaleVsEt = new TF1("f_scaleVsEt", "1.+0.01",0., 1000.);
  f_scaleVsEt = new TF1("f_scaleVsEt", "1. + [0] * (1 - exp(-[1] * x) ) +[2] ",0., 1000.);
  f_scaleVsEt -> SetParameters(7.76e-02,3.73e-02,-6.50e-02);
  
  f_invScaleVsEt = new TF1("f_invScaleVsEt",invScaleVsEt,0.,1000.,0);
  
  
  
  
  
  
  //-----------------
  // Input parameters
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  
  char* analysis   = argv[1];
  char* category   = argv[2];
  int evtsPerPoint = atoi(argv[3]);
  int year         = atoi(argv[4]);
  float DphiMax    = 3.15;
  if( argc >= 6 ) DphiMax = atof(argv[5]);
  
  std::cout << "analysis:     " << analysis     << std::endl;
  std::cout << "category:     " << category     << std::endl;
  std::cout << "evtsPerPoint: " << evtsPerPoint << std::endl;
  std::cout << "year:         " << year         << std::endl;
  std::cout << "DphiMax:      " << DphiMax      << std::endl;
  
  
  
  
  
  
  //--------------------
  // Define in/out files
  std::cout << std::endl;
  std::cout << ">>> define in/out files" << std::endl;
  
  // Get trees
  std::string nameNtuplesMC = "DYJetsToLL";
  std::string nameNtuplesDA = "Data";
  if( MCClosure == 1 )
    nameNtuplesDA = "DYJetsToLL";
  
  TChain* ntu_MC = new TChain(nameNtuplesMC.c_str());
  TChain* ntu_DA = new TChain(nameNtuplesDA.c_str());
  
  if( year == 2012 )
  {
    std::string dataFolder = "/gwteray/users/benaglia/HGG/";
    ntu_MC -> Add((dataFolder+"/tree_zee_moriond_preapproval_phoPD_MCtriggers_ptreweight_noEcalIso.root").c_str());
    ntu_DA -> Add((dataFolder+"/tree_zee_moriond_preapproval_phoPD_MCtriggers_ptreweight_noEcalIso.root").c_str());
  }
  
  
  std::cout << ">>>>>>   MC ntuple: " << std::setw(8) << ntu_MC->GetEntries() << " entries" << std::endl;
  std::cout << ">>>>>> DATA ntuple: " << std::setw(8) << ntu_DA->GetEntries() << " entries" << std::endl;
  
  if( ntu_DA->GetEntries() == 0 || ntu_MC->GetEntries() == 0 )
  {
    std::cout << "Error: At least one file is empty" << std::endl;
    return -1;
  }
  
  std::string plotFolderName = "/gwpool/users/benaglia/CALIBRATION/Hgg/MassScaleStudies/NonGlobe/PLOTS_MZfutyan/";
  char DphiChar[50];
  sprintf(DphiChar,"Dphi%dp%02d",int(DphiMax),int(DphiMax*100)%100);
  plotFolderName += std::string(DphiChar);
  if( applyEtaR9Reweighting == true ) plotFolderName += "_etaR9Reweighting";
  plotFolderName += "/";
  gSystem->mkdir(plotFolderName.c_str());
  
  std::string label = std::string(analysis) + "_cat" + std::string(category) + "_" + std::string(DphiChar);
  if( applyEtaR9Reweighting == true )
  {
    label += "_etaR9Reweighting";
  }
  if( MCClosure == true )
    outFile = TFile::Open((plotFolderName+"/trisCheckScaleMZ_futyan_MCClosure_"+label+".root").c_str(),"RECREATE");
  else
    outFile = TFile::Open((plotFolderName+"/trisCheckScaleMZ_futyan_data_"+label+".root").c_str(),"RECREATE");
  outFile -> cd();
  
  
  
  /*
  //---------------
  // PU reweighting
  std::cout << std::endl;
  std::cout << ">>> PU reweighting" << std::endl;

  TPileupReweighting* puReweighting = NULL;

  if( year == 2012 ) puReweighting = new TPileupReweighting("../Pileup/pileup_69p3mb_true_Moriond2013__Summer12_DR53X-PU_S10_START53.root","h_PUweights");  
  */
  
  
  
  //-------------------
  // eta/R9 reweighting
  
  TFile* etaR9reweightFile = TFile::Open("/gwpool/users/benaglia/CALIBRATION/Hgg/MassScaleStudies/NonGlobe/zee_etaR9reweight.root");
  
  TH2F* etaR9reweight_lead    = (TH2F*)( etaR9reweightFile->Get("etaR9reweight_lead") );
  TH2F* etaR9reweight_sublead = (TH2F*)( etaR9reweightFile->Get("etaR9reweight_lead") );
  
  
  
  
  
  
  //----------------
  // Define branches
  std::cout << std::endl;
  std::cout << ">>> define branches" << std::endl;
  
  // vectors
  std::vector<double> scE_reg1_MC, scEt_reg1_MC, Rt1_MC, R91_MC, scEta1_MC;
  std::vector<double> scE_reg2_MC, scEt_reg2_MC, Rt2_MC, R92_MC, scEta2_MC;
  std::vector<double> scE_reg1_DA, scEt_reg1_DA, Rt1_DA, R91_DA, scEta1_DA;
  std::vector<double> scE_reg2_DA, scEt_reg2_DA, Rt2_DA, R92_DA, scEta2_DA;
  std::vector<double> weight_MC, Ht_MC;
  std::vector<double> weight_DA, Ht_DA;
  std::vector<double> mee_MC, Et1_MC, Et2_MC;
  std::vector<double> mee_DA, Et1_DA, Et2_DA;
  std::vector<double> mee_fit_DA, Et1_fit_DA, Et2_fit_DA;
  std::vector<double> mee_gausFit_DA, Et1_gausFit_DA, Et2_gausFit_DA;
  std::vector<double> mee_mean_DA, Et1_mean_DA, Et2_mean_DA;
  std::vector<double> mee_recursiveMean_DA, Et1_recursiveMean_DA, Et2_recursiveMean_DA;
  std::vector<double> mee_smallestInterval_DA, Et1_smallestInterval_DA, Et2_smallestInterval_DA;
  
  // global variables
  int runId,nVtx,cat;
  float weight;
  
  ntu_MC -> SetBranchStatus("*",0);                         
  ntu_MC -> SetBranchStatus("run",   1); ntu_MC -> SetBranchAddress("run",   &runId);
  ntu_MC -> SetBranchStatus("nvtx",  1); ntu_MC -> SetBranchAddress("nvtx",  &nVtx);
  ntu_MC -> SetBranchStatus("weight",1); ntu_MC -> SetBranchAddress("weight",&weight);
  if( std::string(analysis) == "CiC")
    ntu_MC -> SetBranchStatus("category_baseline",1); ntu_MC -> SetBranchAddress("category_baseline",&cat);
  if( std::string(analysis) == "MVA")
    ntu_MC -> SetBranchStatus("category",1); ntu_MC -> SetBranchAddress("category",&cat);
  
  ntu_DA -> SetBranchStatus("*",0);                         
  ntu_DA -> SetBranchStatus("run",   1); ntu_DA -> SetBranchAddress("run",   &runId);
  ntu_DA -> SetBranchStatus("nvtx",  1); ntu_DA -> SetBranchAddress("nvtx",  &nVtx);
  ntu_DA -> SetBranchStatus("weight",1); ntu_DA -> SetBranchAddress("weight",&weight);
  if( std::string(analysis) == "CiC")
    ntu_DA -> SetBranchStatus("category_baseline",1); ntu_DA -> SetBranchAddress("category_baseline",&cat);
  if( std::string(analysis) == "MVA")
    ntu_DA -> SetBranchStatus("category",1); ntu_DA -> SetBranchAddress("category",&cat);
  
  // electron variables
  double scEne1;
  double scEne2;
  float scEneReg1, scEta1, scPhi1, R91;
  float scEneReg2, scEta2, scPhi2, R92;
  
  
  ntu_DA->SetBranchStatus("pho1_energy",     1); ntu_DA->SetBranchAddress("pho1_energy",     &scEne1);
  ntu_DA->SetBranchStatus("pho1_energy_regr",1); ntu_DA->SetBranchAddress("pho1_energy_regr",&scEneReg1);
  ntu_DA->SetBranchStatus("pho1_sceta",      1); ntu_DA->SetBranchAddress("pho1_sceta",      &scEta1);
  ntu_DA->SetBranchStatus("pho1_scphi",      1); ntu_DA->SetBranchAddress("pho1_scphi",      &scPhi1);
  ntu_DA->SetBranchStatus("pho1_r9",         1); ntu_DA->SetBranchAddress("pho1_r9",         &R91);
  
  ntu_DA->SetBranchStatus("pho2_energy",     1); ntu_DA->SetBranchAddress("pho2_energy",     &scEne2);
  ntu_DA->SetBranchStatus("pho2_energy_regr",1); ntu_DA->SetBranchAddress("pho2_energy_regr",&scEneReg2);
  ntu_DA->SetBranchStatus("pho2_sceta",      1); ntu_DA->SetBranchAddress("pho2_sceta",      &scEta2);
  ntu_DA->SetBranchStatus("pho2_scphi",      1); ntu_DA->SetBranchAddress("pho2_scphi",      &scPhi2);
  ntu_DA->SetBranchStatus("pho2_r9",         1); ntu_DA->SetBranchAddress("pho2_r9",         &R92);
  
  ntu_MC->SetBranchStatus("pho1_energy",     1); ntu_MC->SetBranchAddress("pho1_energy",     &scEne1);
  ntu_MC->SetBranchStatus("pho1_energy_regr",1); ntu_MC->SetBranchAddress("pho1_energy_regr",&scEneReg1);
  ntu_MC->SetBranchStatus("pho1_sceta",      1); ntu_MC->SetBranchAddress("pho1_sceta",      &scEta1);
  ntu_MC->SetBranchStatus("pho1_scphi",      1); ntu_MC->SetBranchAddress("pho1_scphi",      &scPhi1);
  ntu_MC->SetBranchStatus("pho1_r9",         1); ntu_MC->SetBranchAddress("pho1_r9",         &R91);
  
  ntu_MC->SetBranchStatus("pho2_energy",     1); ntu_MC->SetBranchAddress("pho2_energy",     &scEne2);
  ntu_MC->SetBranchStatus("pho2_energy_regr",1); ntu_MC->SetBranchAddress("pho2_energy_regr",&scEneReg2);
  ntu_MC->SetBranchStatus("pho2_sceta",      1); ntu_MC->SetBranchAddress("pho2_sceta",      &scEta2);
  ntu_MC->SetBranchStatus("pho2_scphi",      1); ntu_MC->SetBranchAddress("pho2_scphi",      &scPhi2);
  ntu_MC->SetBranchStatus("pho2_r9",         1); ntu_MC->SetBranchAddress("pho2_r9",         &R92);
  
  
  
  
  
  
  //-----------------
  // Loop over events
  std::cout << std::endl;
  std::cout << ">>> loop over events" << std::endl;
  
  
  for(int ientry = 0; ientry < ntu_MC -> GetEntries(); ++ientry)
  {
    if( ientry%100000 == 0 ) std::cout << ">>>>>> reading   MC entry " << ientry << " / " << ntu_MC->GetEntries() << "\r" << std::flush;
    ntu_MC->GetEntry(ientry);  
    
    
    // define variables
    float theta1 = 2*atan(exp(-scEta1));
    float theta2 = 2*atan(exp(-scEta2));
    float Rt1 = sin(theta1);
    float Rt2 = sin(theta2);
    
    if( MCClosure == false && applyEnergySmearings == true )
    {
      float energySmearing1 = 1. + gRandom->Gaus(0.,GetSmearings(scEta1,R91,year,scEta1<1.5?1:0));
      float energySmearing2 = 1. + gRandom->Gaus(0.,GetSmearings(scEta2,R92,year,scEta2<1.5?1:0));
      
      scEneReg1 *= energySmearing1;
      scEneReg2 *= energySmearing2;
    }
    
    float leadEta = scEta1;
    float leadR9 = R91;
    float subleadEta = scEta2;
    float subleadR9 = R92;
    if( scEneReg1*Rt1 < scEneReg2*Rt2 )
    {
      leadEta = scEta2;
      leadR9 = R92;
      subleadEta = scEta1;
      subleadR9 = R91;
    }
    if( applyEtaR9Reweighting == true )
    {
      float etaR9Weight_lead    = etaR9reweight_lead    -> GetBinContent(etaR9reweight_lead->FindBin(leadEta,leadR9));
      float etaR9Weight_sublead = etaR9reweight_sublead -> GetBinContent(etaR9reweight_sublead->FindBin(subleadEta,subleadR9));    
      
      weight *= etaR9Weight_lead;
      weight *= etaR9Weight_sublead;
    }
    
    TLorentzVector p1; p1.SetPtEtaPhiE(scEneReg1*Rt1,scEta1,scPhi1,scEneReg1);
    TLorentzVector p2; p2.SetPtEtaPhiE(scEneReg2*Rt2,scEta2,scPhi2,scEneReg2);
    float mee = sqrt( 4. * scEneReg1 * scEneReg2 * pow(sin(0.5*(p1.Vect()).Angle(p2.Vect())),2) ) / 91.18;
    float Dphi = deltaPhi(scPhi1,scPhi2);
        
    
    // apply cuts
    if( MCClosure == true ) if( ientry%2 == 1 ) continue;
    if( mee < 0. || mee > 2.5 ) continue;
    if( (atoi(category) != -1) && (atoi(category) != cat) ) continue;
    if( Dphi > DphiMax ) continue;
        
    
    // fill vectors
    scE_reg1_MC.push_back(scEneReg1);
    scE_reg2_MC.push_back(scEneReg2);
    
    scEt_reg1_MC.push_back(scEneReg1*Rt1);
    scEt_reg2_MC.push_back(scEneReg2*Rt2);
    
    Rt1_MC.push_back(Rt1);
    Rt2_MC.push_back(Rt2);
    
    scEta1_MC.push_back(scEta1);
    scEta2_MC.push_back(scEta2);
    
    R91_MC.push_back(R91);
    R92_MC.push_back(R92);
    
    mee_MC.push_back(mee);
    Ht_MC.push_back(scEneReg1*Rt1 + scEneReg2*Rt2);
    Et1_MC.push_back(scEneReg1*Rt1);
    Et2_MC.push_back(scEneReg2*Rt2);
    weight_MC.push_back(weight);
  }
  std::cout << std::endl;
  
  
  for(int ientry = 0; ientry < ntu_DA -> GetEntries(); ientry++)
  {
    if( ientry%100000 == 0 ) std::cout << ">>>>>> reading DATA entry " << ientry << " / " << ntu_DA->GetEntries() << "\r" << std::flush;
    ntu_DA->GetEntry(ientry);  
    
    
    // define variables
    float theta1 = 2*atan(exp(-scEta1));
    float theta2 = 2*atan(exp(-scEta2));
    float Rt1 = sin(theta1);
    float Rt2 = sin(theta2);
    
    if( MCClosure == true )
    {
      scEneReg1 *= ( f_scaleVsEt -> Eval(scEneReg1*Rt1) );
      scEneReg2 *= ( f_scaleVsEt -> Eval(scEneReg2*Rt2) );
    }
    
    if( applyEnergyScaleCorr == true )
    {
      scEneReg1 *= GetShervingCorrections(scEta1,R91,runId);
      scEneReg2 *= GetShervingCorrections(scEta2,R92,runId);
    }
    
    float leadEta = scEta1;
    float leadR9 = R91;
    float subleadEta = scEta2;
    float subleadR9 = R92;
    if( scEneReg1*Rt1 < scEneReg2*Rt2 )
    {
      leadEta = scEta2;
      leadR9 = R92;
      subleadEta = scEta1;
      subleadR9 = R91;
    }
    if( applyEtaR9Reweighting == true )
    {
      float etaR9Weight_lead    = etaR9reweight_lead    -> GetBinContent(etaR9reweight_lead->FindBin(leadEta,leadR9));
      float etaR9Weight_sublead = etaR9reweight_sublead -> GetBinContent(etaR9reweight_sublead->FindBin(subleadEta,subleadR9));    
      
      weight *= etaR9Weight_lead;
      weight *= etaR9Weight_sublead;
    }

    TLorentzVector p1; p1.SetPtEtaPhiE(scEneReg1*Rt1,scEta1,scPhi1,scEneReg1);
    TLorentzVector p2; p2.SetPtEtaPhiE(scEneReg2*Rt2,scEta2,scPhi2,scEneReg2);
    float mee = sqrt( 4. * scEneReg1 * scEneReg2 * pow(sin(0.5*(p1.Vect()).Angle(p2.Vect())),2) ) / 91.18;    
    float Dphi = deltaPhi(scPhi1,scPhi2);    
    
    
    // apply cuts
    if( MCClosure == true ) if( ientry%2 == 0 ) continue;
    if( mee < 0. || mee > 2.5 ) continue;
    if( (atoi(category) != -1) && (atoi(category) != cat) ) continue;
    if( Dphi > DphiMax ) continue;
    
    
    // fill vectors
    scE_reg1_DA.push_back(scEneReg1);
    scE_reg2_DA.push_back(scEneReg2);
    
    scEt_reg1_DA.push_back(scEneReg1*Rt1);
    scEt_reg2_DA.push_back(scEneReg2*Rt2);
    
    Rt1_DA.push_back(Rt1);
    Rt2_DA.push_back(Rt2);
    
    scEta1_DA.push_back(scEta1);
    scEta2_DA.push_back(scEta2);
    
    R91_DA.push_back(R91);
    R92_DA.push_back(R92);
    
    mee_DA.push_back(mee);
    Ht_DA.push_back(scEneReg1*Rt1 + scEneReg2*Rt2);
    mee_fit_DA.push_back(mee);
    mee_gausFit_DA.push_back(mee);
    mee_mean_DA.push_back(mee);
    mee_recursiveMean_DA.push_back(mee);
    mee_smallestInterval_DA.push_back(mee);
    
    Et1_DA.push_back(scEneReg1*Rt1);
    Et2_DA.push_back(scEneReg2*Rt2);
    Et1_fit_DA.push_back(scEneReg1*Rt1);
    Et2_fit_DA.push_back(scEneReg2*Rt2);
    Et1_gausFit_DA.push_back(scEneReg1*Rt1);
    Et2_gausFit_DA.push_back(scEneReg2*Rt2);
    Et1_mean_DA.push_back(scEneReg1*Rt1);
    Et2_mean_DA.push_back(scEneReg2*Rt2);
    Et1_recursiveMean_DA.push_back(scEneReg1*Rt1);
    Et2_recursiveMean_DA.push_back(scEneReg2*Rt2);
    Et1_smallestInterval_DA.push_back(scEneReg1*Rt1);
    Et2_smallestInterval_DA.push_back(scEneReg2*Rt2);
    
    weight_DA.push_back(weight);
  }
  std::cout << std::endl;
  
  
  
  //------------
  // sort events
  std::cout << std::endl;
  std::cout << ">>> sort MC events vs. Ht" << std::endl;
  
  int nEntries = Ht_MC.size();
  int nSavePts = 0;
  std::vector<SorterLC> sortedEntries;

  for(int ientry = 0; ientry < nEntries; ++ientry)
  {
    SorterLC dummy;
    dummy.laserCorr = Ht_MC.at(ientry);
    dummy.entry = ientry;
    sortedEntries.push_back(dummy);
    nSavePts++;   
  }
  
  std::cout << ">>>>>> Sorting variable " << "Ht" << std::endl;
  std::cout << ">>>>>> Effective entries: " << nSavePts << std::endl;
  std::cout << ">>>>>> sortedEntries.size(): " << sortedEntries.size() << std::endl;
  std::sort(sortedEntries.begin(),sortedEntries.end(),SorterLC());
  
  
  
  //------------
  // define bins
  std::cout << std::endl;
  std::cout << ">>> define bins" << std::endl;
   
  HtBinEdges = new std::vector<double>;
  HtBinEdges -> push_back( Ht_MC.at(sortedEntries.at(0).entry) );

  int nBinTempPts = 0;
  for(int iSaved = 0; iSaved < nSavePts; ++iSaved)
  {
    ++nBinTempPts;
    
    if( nBinTempPts == evtsPerPoint )
    {
      HtBinEdges -> push_back( Ht_MC.at(sortedEntries.at(iSaved).entry) );
      nBinTempPts = 0;
    }
  }
  HtBinEdges -> push_back( Ht_MC.at(sortedEntries.at(nSavePts-1).entry) );
  
  nHtBins = HtBinEdges->size() - 1;
  for(unsigned int i = 0; i < nHtBins; ++i)
    std::cout << ">>> Ht bin " << i << ":   [" << HtBinEdges->at(i) << "," << HtBinEdges->at(i+1) << "]" << std::endl;
  std::cout << std::endl;
  
  
  
  EtBinEdges = new std::vector<double>;
  for(unsigned int HtBinEdgeIt = 0; HtBinEdgeIt < HtBinEdges->size(); ++HtBinEdgeIt)
    EtBinEdges -> push_back( 0.5 * HtBinEdges->at(HtBinEdgeIt) );
  nEtBins = EtBinEdges->size()-1;
  
  for(unsigned int i = 0; i < nEtBins; ++i)
    std::cout << ">>> Et bin " << i << ":   [" << EtBinEdges->at(i) << "," << EtBinEdges->at(i+1) << "]" << std::endl;
  std::cout << std::endl;
  
  
  
  
  
  
  //------------------
  // define histograms
  std::cout << std::endl;
  std::cout << ">>> define histograms" << std::endl;
  
  TH1F* h_scEta_MC = new TH1F("h_scEta_MC","",500,-2.5,2.5);
  h_scEta_MC -> Sumw2();
  TH1F* h_scEta_DA = new TH1F("h_scEta_DA","",500,-2.5,2.5);
  h_scEta_DA -> Sumw2();
  TH1F* h_R9_MC = new TH1F("h_R9_MC","",500,-1.,1.);
  h_R9_MC -> Sumw2();
  TH1F* h_R9_DA = new TH1F("h_R9_DA","",500,-1.,1.);
  h_R9_DA -> Sumw2();
  TH1F* h_Ht_MC = new TH1F("h_Ht_MC","",500,0.,500.);
  h_Ht_MC -> Sumw2();
  TH1F* h_Ht_DA = new TH1F("h_Ht_DA","",500,0.,500.);
  h_Ht_DA -> Sumw2();
  TH1F* h_mee_MC = new TH1F("h_mee_MC","",480,60.,120.);
  h_mee_MC -> Sumw2();
  TH1F* h_mee_DA = new TH1F("h_mee_DA","",480,60.,120.);
  h_mee_DA -> Sumw2();
  
  TGraphAsymmErrors* scale_fit_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_fit_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_fit_DAOverMC = new TGraphAsymmErrors();
  
  TGraphAsymmErrors* scale_gausFit_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_gausFit_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_gausFit_DAOverMC = new TGraphAsymmErrors();
  
  TGraphAsymmErrors* scale_mean_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_mean_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_mean_DAOverMC = new TGraphAsymmErrors();
  
  TGraphAsymmErrors* scale_recursiveMean_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_recursiveMean_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_recursiveMean_DAOverMC = new TGraphAsymmErrors();
  
  TGraphAsymmErrors* scale_smallestInterval_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_smallestInterval_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_smallestInterval_DAOverMC = new TGraphAsymmErrors();
  
  std::vector<double>* mee_HtBin_MC    = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_MC = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_DA                      = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_fit_DA                  = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_gausFit_DA              = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_mean_DA                 = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_recursiveMean_DA        = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_smallestInterval_DA     = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_DA                   = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_fit_DA               = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_gausFit_DA           = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_mean_DA              = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_recursiveMean_DA     = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_smallestInterval_DA  = new std::vector<double>[nHtBins];
    
  TH1F** h_mee_HtBin_MC = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_fit_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_gausFit_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_mean_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_recursiveMean_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_smallestInterval_DA = new TH1F*[nHtBins];
  
  TF1** f_gausFit_HtBin_MC = new TF1*[nHtBins];
  TF1** f_gausFit_HtBin_DA = new TF1*[nHtBins];
  
  h_Ht_HtBin_MC = new TH1F*[nHtBins];
  
  h_Et_EtBin_MC = new TH1F*[nEtBins];
  h_Et_EtBin_DA = new TH1F*[nEtBins];
  h_Et_EtBin_fit_DA = new TH1F*[nEtBins];
  h_Et_EtBin_gausFit_DA = new TH1F*[nEtBins];
  h_Et_EtBin_mean_DA = new TH1F*[nEtBins];
  h_Et_EtBin_recursiveMean_DA = new TH1F*[nEtBins];
  h_Et_EtBin_smallestInterval_DA = new TH1F*[nEtBins];
  
  for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
  {
    char histoName[50];
    
    sprintf(histoName,"h_mee_HtBin%d_MC",HtBin);
    h_mee_HtBin_MC[HtBin] = new TH1F(histoName,"",nBins_mee,meeMin,meeMax);
    //h_mee_HtBin_MC[HtBin] -> SetFillColor(kGreen+2);
    //h_mee_HtBin_MC[HtBin] -> SetFillStyle(3004);
    //h_mee_HtBin_MC[HtBin] -> SetMarkerStyle(7);
    //h_mee_HtBin_MC[HtBin] -> SetMarkerColor(kGreen+2);
    //h_mee_HtBin_MC[HtBin] -> SetLineColor(kGreen+2);
    h_mee_HtBin_MC[HtBin]->Sumw2();
    
    sprintf(histoName, "Ht_HtBin%d_MC",HtBin);
    h_Ht_HtBin_MC[HtBin] = new TH1F(histoName,"",5000,0.,1000.);
    h_Ht_HtBin_MC[HtBin] -> SetLineColor(kGreen+2);
    h_Ht_HtBin_MC[HtBin] -> Sumw2();
  }
  
  for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
  {
    char histoName[50];
    
    sprintf(histoName, "Et_EtBin%d_MC",EtBin);
    h_Et_EtBin_MC[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
    h_Et_EtBin_MC[EtBin] -> SetLineColor(kGreen+2);
    h_Et_EtBin_MC[EtBin] -> Sumw2();
  }
  
  
  
  
  //----------------
  // fill histograms
  std::cout << std::endl;
  std::cout << ">>> fill histograms" << std::endl;
  
  int MCEntries = mee_MC.size();
  for(int ientry = 0; ientry < MCEntries; ++ientry)
  {   
    if( (ientry%100000 == 0) ) std::cout << "reading   MC entry " << ientry << " / " << MCEntries << "\r" << std::flush;
    
    
    int HtBin = MyFindBin(Ht_MC.at(ientry),HtBinEdges);
    if( HtBin == -1 ) continue;
    
    int EtBin1 = MyFindBin(Et1_MC.at(ientry),EtBinEdges);
    int EtBin2 = MyFindBin(Et2_MC.at(ientry),EtBinEdges);
    if( EtBin1 == -1 ) continue;
    if( EtBin2 == -1 ) continue;
    
    //std::cout << "HtBin: " << HtBin << "   EtBin1: " << EtBin1 << "   EtBin2: " << EtBin2 << std::endl;
    mee_HtBin_MC[HtBin].push_back( mee_MC.at(ientry) );
    weight_HtBin_MC[HtBin].push_back( weight_MC.at(ientry) );
    
    h_mee_HtBin_MC[HtBin] -> Fill( mee_MC.at(ientry),weight_MC.at(ientry) );
    h_Ht_HtBin_MC[HtBin]  -> Fill(  Ht_MC.at(ientry),weight_MC.at(ientry) );
    
    h_Et_EtBin_MC[EtBin1] -> Fill( Et1_MC.at(ientry),weight_MC.at(ientry) );
    h_Et_EtBin_MC[EtBin2] -> Fill( Et2_MC.at(ientry),weight_MC.at(ientry) );
    
    h_scEta_MC -> Fill( scEta1_MC.at(ientry),weight_MC.at(ientry) );
    h_scEta_MC -> Fill( scEta2_MC.at(ientry),weight_MC.at(ientry) );
    
    h_R9_MC -> Fill( R91_MC.at(ientry),weight_MC.at(ientry) );
    h_R9_MC -> Fill( R92_MC.at(ientry),weight_MC.at(ientry) );
    
    h_Ht_MC -> Fill( Ht_MC.at(ientry),weight_MC.at(ientry) );
    h_Ht_MC -> Fill( Ht_MC.at(ientry),weight_MC.at(ientry) );
    
    h_mee_MC -> Fill( mee_MC.at(ientry)*91.18,weight_MC.at(ientry) );
    h_mee_MC -> Fill( mee_MC.at(ientry)*91.18,weight_MC.at(ientry) );
  }
  std::cout << std::endl;
  
  
  
  for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
  {
    double x = h_Ht_HtBin_MC[HtBin]->GetMean();
    
    scale_fit_DA -> SetPoint(HtBin,x,1.);
    scale_fit_MC -> SetPoint(HtBin,x,1.);
    scale_fit_DAOverMC -> SetPoint(HtBin,x,1.);
    
    scale_gausFit_DA -> SetPoint(HtBin,x,1.);
    scale_gausFit_MC -> SetPoint(HtBin,x,1.);
    scale_gausFit_DAOverMC -> SetPoint(HtBin,x,1.);
    
    scale_mean_DA -> SetPoint(HtBin,x,1.);
    scale_mean_MC -> SetPoint(HtBin,x,1.);
    scale_mean_DAOverMC -> SetPoint(HtBin,x,1.);
    
    scale_recursiveMean_DA -> SetPoint(HtBin,x,1.);
    scale_recursiveMean_MC -> SetPoint(HtBin,x,1.);
    scale_recursiveMean_DAOverMC -> SetPoint(HtBin,x,1.);
    
    scale_smallestInterval_DA -> SetPoint(HtBin,x,1.);
    scale_smallestInterval_MC -> SetPoint(HtBin,x,1.);
    scale_smallestInterval_DAOverMC -> SetPoint(HtBin,x,1.);
  }
  
  
  
  for(int step = 1; step < nSteps+1; ++step)
  {
    std::cout << std::endl;
    std::cout << "****** step " << step << " ******" << std::endl;
    
    
    TProfile* p_avgEtCorr_fit              = new TProfile("p_avgEtCorr_fit",             "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_gausFit          = new TProfile("p_avgEtCorr_gausFit",         "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_mean             = new TProfile("p_avgEtCorr_mean",            "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_recursiveMean    = new TProfile("p_avgEtCorr_recursiveMean",   "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_smallestInterval = new TProfile("p_avgEtCorr_smallestInterval","",nEtBins,30.,1000.);
    MySetBins(p_avgEtCorr_fit,*EtBinEdges);
    MySetBins(p_avgEtCorr_gausFit,*EtBinEdges);
    MySetBins(p_avgEtCorr_mean,*EtBinEdges);
    MySetBins(p_avgEtCorr_recursiveMean,*EtBinEdges);
    MySetBins(p_avgEtCorr_smallestInterval,*EtBinEdges);
    
    TProfile* p_avgEtCorr_rebin2_fit              = new TProfile("p_avgEtCorr_rebin2_fit",             "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_rebin2_gausFit          = new TProfile("p_avgEtCorr_rebin2_gausFit",         "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_rebin2_mean             = new TProfile("p_avgEtCorr_rebin2_mean",            "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_rebin2_recursiveMean    = new TProfile("p_avgEtCorr_rebin2_recursiveMean",   "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_rebin2_smallestInterval = new TProfile("p_avgEtCorr_rebin2_smallestInterval","",nEtBins,30.,1000.);
    MySetBins(p_avgEtCorr_rebin2_fit,*EtBinEdges);
    MySetBins(p_avgEtCorr_rebin2_gausFit,*EtBinEdges);
    MySetBins(p_avgEtCorr_rebin2_mean,*EtBinEdges);
    MySetBins(p_avgEtCorr_rebin2_recursiveMean,*EtBinEdges);
    MySetBins(p_avgEtCorr_rebin2_smallestInterval,*EtBinEdges);
    p_avgEtCorr_rebin2_fit              -> Rebin(2);
    p_avgEtCorr_rebin2_gausFit          -> Rebin(2);
    p_avgEtCorr_rebin2_mean             -> Rebin(2);
    p_avgEtCorr_rebin2_recursiveMean    -> Rebin(2);
    p_avgEtCorr_rebin2_smallestInterval -> Rebin(2);
    
    for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
    {
        mee_HtBin_DA[HtBin].clear();
        mee_HtBin_fit_DA[HtBin].clear();
        mee_HtBin_gausFit_DA[HtBin].clear();
        mee_HtBin_mean_DA[HtBin].clear();
        mee_HtBin_recursiveMean_DA[HtBin].clear();
        mee_HtBin_smallestInterval_DA[HtBin].clear();
        weight_HtBin_DA[HtBin].clear();
        weight_HtBin_fit_DA[HtBin].clear();
        weight_HtBin_gausFit_DA[HtBin].clear();
        weight_HtBin_mean_DA[HtBin].clear();
        weight_HtBin_recursiveMean_DA[HtBin].clear();
        weight_HtBin_smallestInterval_DA[HtBin].clear();
    }	
    for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
    {
      char histoName[50];
      
      sprintf(histoName,"h_mee_HtBin%d_DA_step%d",HtBin,step);
      h_mee_HtBin_DA[HtBin] = new TH1F(histoName,"",nBins_mee,meeMin,meeMax);
      h_mee_HtBin_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_fit_DA_step%d",HtBin,step);
      h_mee_HtBin_fit_DA[HtBin] = new TH1F(histoName,"",nBins_mee,meeMin,meeMax);
      h_mee_HtBin_fit_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_gausFit_DA_step%d",HtBin,step);
      h_mee_HtBin_gausFit_DA[HtBin] = new TH1F(histoName,"",nBins_mee,meeMin,meeMax);
      h_mee_HtBin_gausFit_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_mean_DA_step%d",HtBin,step);
      h_mee_HtBin_mean_DA[HtBin] = new TH1F(histoName,"",nBins_mee,meeMin,meeMax);
      h_mee_HtBin_mean_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_recursiveMean_DA_step%d",HtBin,step);
      h_mee_HtBin_recursiveMean_DA[HtBin] = new TH1F(histoName,"",nBins_mee,meeMin,meeMax);
      h_mee_HtBin_recursiveMean_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_smallestInterval_DA_step%d",HtBin,step);
      h_mee_HtBin_smallestInterval_DA[HtBin] = new TH1F(histoName,"",nBins_mee,meeMin,meeMax);
      h_mee_HtBin_smallestInterval_DA[HtBin] -> Sumw2();
    }
    for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
    {
      char histoName[50];
      
      sprintf(histoName,"h_Et_EtBin%d_DA_step%d",EtBin,step);
      h_Et_EtBin_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_Et_EtBin%d_fit_DA_step%d",EtBin,step);
      h_Et_EtBin_fit_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_fit_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_fit_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_Et_EtBin%d_gausFit_DA_step%d",EtBin,step);
      h_Et_EtBin_gausFit_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_gausFit_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_gausFit_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_Et_EtBin%d_mean_DA_step%d",EtBin,step);
      h_Et_EtBin_mean_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_mean_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_mean_DA[EtBin] -> Sumw2();
            
      sprintf(histoName,"h_Et_EtBin%d_recursiveMean_DA_step%d",EtBin,step);
      h_Et_EtBin_recursiveMean_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_recursiveMean_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_recursiveMean_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_Et_EtBin%d_smallestInterval_DA_step%d",EtBin,step);
      h_Et_EtBin_smallestInterval_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_smallestInterval_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_smallestInterval_DA[EtBin] -> Sumw2();
    }
    
    
    int DAEntries = mee_DA.size();
    for(int ientry = 0; ientry < DAEntries; ++ientry)
    {   
      if( (ientry%100000 == 0) ) std::cout << "reading DATA entry " << ientry << " / " << DAEntries << "\r" << std::flush;
      
      
      double k1 = 1.;
      double k2 = 1.;
      double Et1 = Et1_DA.at(ientry)/k1;
      double Et2 = Et2_DA.at(ientry)/k2;
      int EtBin1 = MyFindBin(Et1,EtBinEdges);
      int EtBin2 = MyFindBin(Et2,EtBinEdges);
      double Ht = Et1 + Et2;
      int HtBin = MyFindBin(Ht,HtBinEdges);
      double mee = mee_DA.at(ientry)/sqrt(k1*k2);
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_DA[HtBin].push_back( mee );
        weight_HtBin_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        
        h_mee_HtBin_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      
      k1 = MyEval(scale_fit_DAOverMC,2.*Et1_fit_DA.at(ientry));
      k2 = MyEval(scale_fit_DAOverMC,2.*Et2_fit_DA.at(ientry));
      Et1 = Et1_fit_DA.at(ientry)/k1;
      Et2 = Et2_fit_DA.at(ientry)/k2;
      Et1_fit_DA.at(ientry) = Et1;
      Et2_fit_DA.at(ientry) = Et2;
      p_avgEtCorr_fit -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_fit -> Fill(Et2,Et2/Et2_DA.at(ientry));
      p_avgEtCorr_rebin2_fit -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_rebin2_fit -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_fit_DA.at(ientry)/sqrt(k1*k2);
      mee_fit_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_fit_DA[HtBin].push_back( mee );
        weight_HtBin_fit_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_fit_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_fit_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_fit_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      k1 = MyEval(scale_mean_DAOverMC,2.*Et1_mean_DA.at(ientry));
      k2 = MyEval(scale_mean_DAOverMC,2.*Et2_mean_DA.at(ientry));
      Et1 = Et1_mean_DA.at(ientry)/k1;
      Et2 = Et2_mean_DA.at(ientry)/k2;
      Et1_mean_DA.at(ientry) = Et1;
      Et2_mean_DA.at(ientry) = Et2;
      p_avgEtCorr_mean -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_mean -> Fill(Et2,Et2/Et2_DA.at(ientry));
      p_avgEtCorr_rebin2_mean -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_rebin2_mean -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_mean_DA.at(ientry)/sqrt(k1*k2);
      mee_mean_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_mean_DA[HtBin].push_back( mee );
        weight_HtBin_mean_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_mean_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_mean_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_mean_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      k1 = MyEval(scale_gausFit_DAOverMC,2.*Et1_gausFit_DA.at(ientry));
      k2 = MyEval(scale_gausFit_DAOverMC,2.*Et2_gausFit_DA.at(ientry));
      Et1 = Et1_gausFit_DA.at(ientry)/k1;
      Et2 = Et2_gausFit_DA.at(ientry)/k2;
      Et1_gausFit_DA.at(ientry) = Et1;
      Et2_gausFit_DA.at(ientry) = Et2;
      p_avgEtCorr_gausFit -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_gausFit -> Fill(Et2,Et2/Et2_DA.at(ientry));
      p_avgEtCorr_rebin2_gausFit -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_rebin2_gausFit -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_gausFit_DA.at(ientry)/sqrt(k1*k2);
      mee_gausFit_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_gausFit_DA[HtBin].push_back( mee );
        weight_HtBin_gausFit_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_gausFit_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_gausFit_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_gausFit_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      k1 = MyEval(scale_recursiveMean_DAOverMC,2.*Et1_recursiveMean_DA.at(ientry));
      k2 = MyEval(scale_recursiveMean_DAOverMC,2.*Et2_recursiveMean_DA.at(ientry));
      Et1 = Et1_recursiveMean_DA.at(ientry)/k1;
      Et2 = Et2_recursiveMean_DA.at(ientry)/k2;
      Et1_recursiveMean_DA.at(ientry) = Et1;
      Et2_recursiveMean_DA.at(ientry) = Et2;
      p_avgEtCorr_recursiveMean -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_recursiveMean -> Fill(Et2,Et2/Et2_DA.at(ientry));
      p_avgEtCorr_rebin2_recursiveMean -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_rebin2_recursiveMean -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_recursiveMean_DA.at(ientry)/sqrt(k1*k2);
      mee_recursiveMean_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_recursiveMean_DA[HtBin].push_back( mee );
        weight_HtBin_recursiveMean_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_recursiveMean_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_recursiveMean_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_recursiveMean_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      k1 = MyEval(scale_smallestInterval_DAOverMC,2.*Et1_smallestInterval_DA.at(ientry));
      k2 = MyEval(scale_smallestInterval_DAOverMC,2.*Et2_smallestInterval_DA.at(ientry));
      Et1 = Et1_smallestInterval_DA.at(ientry)/k1;
      Et2 = Et2_smallestInterval_DA.at(ientry)/k2;
      Et1_smallestInterval_DA.at(ientry) = Et1;
      Et2_smallestInterval_DA.at(ientry) = Et2;
      p_avgEtCorr_smallestInterval -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_smallestInterval -> Fill(Et2,Et2/Et2_DA.at(ientry));
      p_avgEtCorr_rebin2_smallestInterval -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_rebin2_smallestInterval -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_smallestInterval_DA.at(ientry)/sqrt(k1*k2);
      mee_smallestInterval_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_smallestInterval_DA[HtBin].push_back( mee );
        weight_HtBin_smallestInterval_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_smallestInterval_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_smallestInterval_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_smallestInterval_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      
      if( step == 1 )
      {
        h_scEta_DA -> Fill( scEta1_DA.at(ientry),weight_DA.at(ientry) );
        h_scEta_DA -> Fill( scEta2_DA.at(ientry),weight_DA.at(ientry) );
        
        h_R9_DA -> Fill( R91_DA.at(ientry),weight_DA.at(ientry) );
        h_R9_DA -> Fill( R92_DA.at(ientry),weight_DA.at(ientry) );
        
        h_Ht_DA -> Fill( Ht_DA.at(ientry),weight_DA.at(ientry) );
        h_Ht_DA -> Fill( Ht_DA.at(ientry),weight_DA.at(ientry) );
        
        h_mee_DA -> Fill( mee_DA.at(ientry)*91.18,weight_DA.at(ientry) );
        h_mee_DA -> Fill( mee_DA.at(ientry)*91.18,weight_DA.at(ientry) );
      }
    }
    std::cout << std::endl;
    
    
    
    
    
    
    //---------------------
    // find scale estimator
    std::cout << std::endl;
    std::cout << ">>> find scale estimator" << std::endl;
    
    std::vector<double> smallestIntervalMins_MC;
    std::vector<double> smallestIntervalMaxs_MC;
    std::vector<double> smallestIntervalMins_DA;
    std::vector<double> smallestIntervalMaxs_DA;
    
    for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
    {
      std::cout << ">>> HtBin::" << HtBin << std::endl;
      std::cout << ">>>>>> nEvents:                  " << mee_HtBin_DA[HtBin].size() << std::endl;
      std::cout << ">>>>>> nEvents_fit:              " << mee_HtBin_fit_DA[HtBin].size() << std::endl;
      std::cout << ">>>>>> nEvents_gausFit:          " << mee_HtBin_gausFit_DA[HtBin].size() << std::endl;
      std::cout << ">>>>>> nEvents_mean:             " << mee_HtBin_mean_DA[HtBin].size() << std::endl;
      std::cout << ">>>>>> nEvents_recursiveMean:    " << mee_HtBin_recursiveMean_DA[HtBin].size() << std::endl;
      std::cout << ">>>>>> nEvents_smallestInterval: " << mee_HtBin_smallestInterval_DA[HtBin].size() << std::endl;
      
      double x = h_Ht_HtBin_MC[HtBin]->GetMean();
      double ex = h_Ht_HtBin_MC[HtBin]->GetMeanError();
      double exlow = ex;
      double exhig = ex;
      
      h_mee_HtBin_MC[HtBin] -> Scale(1./h_mee_HtBin_MC[HtBin]->Integral());
      h_mee_HtBin_DA[HtBin] -> Scale(1./h_mee_HtBin_DA[HtBin]->Integral());
      h_mee_HtBin_fit_DA[HtBin] -> Scale(1./h_mee_HtBin_fit_DA[HtBin]->Integral());
      h_mee_HtBin_gausFit_DA[HtBin] -> Scale(1./h_mee_HtBin_gausFit_DA[HtBin]->Integral());
      h_mee_HtBin_mean_DA[HtBin] -> Scale(1./h_mee_HtBin_mean_DA[HtBin]->Integral());
      h_mee_HtBin_recursiveMean_DA[HtBin] -> Scale(1./h_mee_HtBin_recursiveMean_DA[HtBin]->Integral());
      h_mee_HtBin_smallestInterval_DA[HtBin] -> Scale(1./h_mee_HtBin_smallestInterval_DA[HtBin]->Integral());

      
      // fit
      if( mee_HtBin_fit_DA[HtBin].size() > 3 )
      {
        double scale_MC = 1.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        MyFindFit(scale_DA,scaleErr_DA,h_mee_HtBin_MC[HtBin],h_mee_HtBin_fit_DA[HtBin]);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_fit_MC -> SetPoint(HtBin,x,scale_MC);
        scale_fit_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_fit_DA -> SetPoint(HtBin,x,scale_DA);
        scale_fit_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_fit_DAOverMC -> SetPoint(HtBin,x,y);
        scale_fit_DAOverMC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // gausFit
      if( mee_HtBin_gausFit_DA[HtBin].size() > 3 )
      {
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        char funcName_MC[50];
        sprintf(funcName_MC,"f_gausFit_HtBin%d_step%d_MC",HtBin,step);
        char funcName_DA[50];
        sprintf(funcName_DA,"f_gausFit_HtBin%d_step%d_DA",HtBin,step);
        MyFindGausFit(scale_MC,scaleErr_MC,mee_HtBin_MC[HtBin],weight_HtBin_MC[HtBin],&(f_gausFit_HtBin_MC[HtBin]),std::string(funcName_MC));
        MyFindGausFit(scale_DA,scaleErr_DA,mee_HtBin_gausFit_DA[HtBin],weight_HtBin_gausFit_DA[HtBin],&(f_gausFit_HtBin_DA[HtBin]),std::string(funcName_DA));
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_gausFit_MC -> SetPoint(HtBin,x,scale_MC);
        scale_gausFit_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_gausFit_DA -> SetPoint(HtBin,x,scale_DA);
        scale_gausFit_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_gausFit_DAOverMC -> SetPoint(HtBin,x,y);
        scale_gausFit_DAOverMC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // mean
      if( mee_HtBin_mean_DA[HtBin].size() > 3 )
      {
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        MyFindMean(scale_MC,scaleErr_MC,mee_HtBin_MC[HtBin],weight_HtBin_MC[HtBin]);
        MyFindMean(scale_DA,scaleErr_DA,mee_HtBin_mean_DA[HtBin],weight_HtBin_mean_DA[HtBin]);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_mean_MC -> SetPoint(HtBin,x,scale_MC);
        scale_mean_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_mean_DA -> SetPoint(HtBin,x,scale_DA);
        scale_mean_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_mean_DAOverMC -> SetPoint(HtBin,x,y);
        scale_mean_DAOverMC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // recursive mean
      if( mee_HtBin_recursiveMean_DA[HtBin].size() > 3 )
      {
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        MyFindRecursiveMean(scale_MC,scaleErr_MC,mee_HtBin_MC[HtBin],weight_HtBin_MC[HtBin],5./91.18,0.001);
        MyFindRecursiveMean(scale_DA,scaleErr_DA,mee_HtBin_recursiveMean_DA[HtBin],weight_HtBin_recursiveMean_DA[HtBin],5./91.18,0.001);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_recursiveMean_MC -> SetPoint(HtBin,x,scale_MC);
        scale_recursiveMean_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_recursiveMean_DA -> SetPoint(HtBin,x,scale_DA);
        scale_recursiveMean_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_recursiveMean_DAOverMC -> SetPoint(HtBin,x,y);
        scale_recursiveMean_DAOverMC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // smallest interval
      if( mee_HtBin_smallestInterval_DA[HtBin].size() > 3 )
      {
        double min_MC;
        double max_MC;
        double min_DA;
        double max_DA;
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        MyFindSmallestInterval(scale_MC,scaleErr_MC,min_MC,max_MC,mee_HtBin_MC[HtBin],weight_HtBin_MC[HtBin],0.68);
        MyFindSmallestInterval(scale_DA,scaleErr_DA,min_DA,max_DA,mee_HtBin_smallestInterval_DA[HtBin],weight_HtBin_smallestInterval_DA[HtBin],0.68);
        smallestIntervalMins_MC.push_back(min_MC);
        smallestIntervalMaxs_MC.push_back(max_MC);
        smallestIntervalMins_DA.push_back(min_DA);
        smallestIntervalMaxs_DA.push_back(max_DA);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_smallestInterval_MC -> SetPoint(HtBin,x,scale_MC);
        scale_smallestInterval_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_smallestInterval_DA -> SetPoint(HtBin,x,scale_DA);
        scale_smallestInterval_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_smallestInterval_DAOverMC -> SetPoint(HtBin,x,y);
        scale_smallestInterval_DAOverMC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
      }
    }
    
    

    
    outFile -> cd();
    
    char dirName[50];
    sprintf(dirName,"step%d",step);
    TDirectory* baseDir = outFile -> mkdir(dirName);
    TDirectory* subDir = NULL;
    baseDir -> cd();
    
    scale_fit_MC           -> Write("scale_fit_MC");
    scale_fit_DA           -> Write("scale_fit_DA");
    scale_fit_DAOverMC     -> Write("scale_fit_DAOverMC");
    p_avgEtCorr_fit        -> Write();
    p_avgEtCorr_rebin2_fit -> Write();
    
    scale_gausFit_MC           -> Write("scale_gausFit_MC");
    scale_gausFit_DA           -> Write("scale_gausFit_DA");
    scale_gausFit_DAOverMC     -> Write("scale_gausFit_DAOverMC");
    p_avgEtCorr_gausFit        -> Write();
    p_avgEtCorr_rebin2_gausFit -> Write();
    
    scale_mean_MC           -> Write("scale_mean_MC");
    scale_mean_DA           -> Write("scale_mean_DA");
    scale_mean_DAOverMC     -> Write("scale_mean_DAOverMC");
    p_avgEtCorr_mean        -> Write();
    p_avgEtCorr_rebin2_mean -> Write();
    
    scale_recursiveMean_MC           -> Write("scale_recursiveMean_MC");
    scale_recursiveMean_DA           -> Write("scale_recursiveMean_DA");
    scale_recursiveMean_DAOverMC     -> Write("scale_recursiveMean_DAOverMC");
    p_avgEtCorr_recursiveMean        -> Write();
    p_avgEtCorr_rebin2_recursiveMean -> Write();
    
    scale_smallestInterval_MC           -> Write("scale_smallestInterval_MC");
    scale_smallestInterval_DA           -> Write("scale_smallestInterval_DA");
    scale_smallestInterval_DAOverMC     -> Write("scale_smallestInterval_DAOverMC");
    p_avgEtCorr_smallestInterval        -> Write();
    p_avgEtCorr_rebin2_smallestInterval -> Write();
    
    
    delete p_avgEtCorr_fit;
    delete p_avgEtCorr_gausFit;
    delete p_avgEtCorr_mean;
    delete p_avgEtCorr_recursiveMean;
    delete p_avgEtCorr_smallestInterval;
    delete p_avgEtCorr_rebin2_fit;
    delete p_avgEtCorr_rebin2_gausFit;
    delete p_avgEtCorr_rebin2_mean;
    delete p_avgEtCorr_rebin2_recursiveMean;
    delete p_avgEtCorr_rebin2_smallestInterval;
    
    
    baseDir -> cd();
    subDir = baseDir -> mkdir("mee_HtBin");
    subDir -> cd();
    
    std::string outputPdf_MC = plotFolderName + "h_mee_HtBin_MC_"   + std::string(analysis)+"_cat"+ std::string(category) + ".pdf";
    std::string outputPdf_DA = plotFolderName + "h_mee_HtBin_DATA_" + std::string(analysis)+"_cat"+ std::string(category) + ".pdf";
    
    for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
    {
      h_mee_HtBin_MC[HtBin] -> Write();
      h_mee_HtBin_DA[HtBin] -> Write();
      h_mee_HtBin_fit_DA[HtBin] -> Write();
      h_mee_HtBin_gausFit_DA[HtBin] -> Write();
      h_mee_HtBin_mean_DA[HtBin] -> Write();
      h_mee_HtBin_recursiveMean_DA[HtBin] -> Write();
      h_mee_HtBin_smallestInterval_DA[HtBin] -> Write();
      
      f_gausFit_HtBin_MC[HtBin] -> Write();
      f_gausFit_HtBin_DA[HtBin] -> Write();
      

      
      if( step == 1 )
      {
        TCanvas* c_DAOverMC = new TCanvas();
        c_DAOverMC -> cd();
        c_DAOverMC -> SetGridx();
        c_DAOverMC -> SetGridy();
        
        char axisTitle[50];
        sprintf(axisTitle,"m_{ee}/m_{Z}^{PDG}   -   H_{T} #in [%d,%d]",int(HtBinEdges->at(HtBin)),int(HtBinEdges->at(HtBin+1)));
        h_mee_HtBin_MC[HtBin] -> GetXaxis() -> SetTitle(axisTitle);
        sprintf(axisTitle,"event fraction");
        h_mee_HtBin_MC[HtBin] -> GetYaxis() -> SetTitle(axisTitle);
        h_mee_HtBin_MC[HtBin] -> GetXaxis() -> SetLabelSize(0.04);
        h_mee_HtBin_MC[HtBin] -> GetYaxis() -> SetLabelSize(0.04);
        h_mee_HtBin_MC[HtBin] -> GetXaxis() -> SetTitleSize(0.05);
        h_mee_HtBin_MC[HtBin] -> GetYaxis() -> SetTitleSize(0.05);
        h_mee_HtBin_MC[HtBin] -> SetLineColor(kBlack);
        h_mee_HtBin_MC[HtBin] -> SetLineWidth(1);
        h_mee_HtBin_MC[HtBin] -> GetXaxis() -> SetRangeUser(0.65,1.34999);
        float maximum = std::max(h_mee_HtBin_MC[HtBin]->GetMaximum(),h_mee_HtBin_DA[HtBin]->GetMaximum());
        h_mee_HtBin_MC[HtBin] -> SetMaximum( 1.1*maximum );
        
        h_mee_HtBin_MC[HtBin] -> Draw("hist");
        h_mee_HtBin_DA[HtBin] -> Draw("P,same");
        
        char pngName[200];
        sprintf(pngName,"%s/h_mee_cat%s_%s_HtBin%03d_DAOverMC.png",plotFolderName.c_str(),category,analysis,HtBin);
        c_DAOverMC -> Print(pngName,"png");
        
        
        TCanvas* c_MC = new TCanvas();
        c_MC -> cd();
        c_MC -> SetGridx();
        c_MC -> SetGridy();
        
        double x,y;
        
        scale_mean_MC -> GetPoint(HtBin,x,y);
        TArrow* line_mean_MC = new TArrow(y,0.,y,h_mee_HtBin_MC[HtBin]->GetBinContent(h_mee_HtBin_MC[HtBin]->FindBin(y)));
        line_mean_MC -> SetLineColor(kBlack);
        line_mean_MC -> SetLineWidth(2);
        
        scale_recursiveMean_MC -> GetPoint(HtBin,x,y);
        TArrow* line_recursiveMean_MC = new TArrow(y,0.,y,h_mee_HtBin_MC[HtBin]->GetMaximum());
        line_recursiveMean_MC -> SetLineColor(kBlue);
        line_recursiveMean_MC -> SetLineWidth(2);
        TArrow* line_recursiveMean_min_MC = new TArrow(y-5./91.18,0.,y-5./91.18,h_mee_HtBin_MC[HtBin]->GetMaximum());
        line_recursiveMean_min_MC -> SetLineColor(kBlue);
        line_recursiveMean_min_MC -> SetLineWidth(2);
        TArrow* line_recursiveMean_max_MC = new TArrow(y+5./91.18,0.,y+5./91.18,h_mee_HtBin_MC[HtBin]->GetMaximum());
        line_recursiveMean_max_MC -> SetLineColor(kBlue);
        line_recursiveMean_max_MC -> SetLineWidth(2);
        
        scale_smallestInterval_MC -> GetPoint(HtBin,x,y);
        TArrow* line_smallestInterval_MC = new TArrow(y,0.,y,h_mee_HtBin_MC[HtBin]->GetBinContent(h_mee_HtBin_MC[HtBin]->FindBin(y)));
        line_smallestInterval_MC -> SetLineColor(kRed);
        line_smallestInterval_MC -> SetLineWidth(2);
        
        h_mee_HtBin_MC[HtBin] -> Draw("HIST");
        
        TH1F* clone = (TH1F*)( h_mee_HtBin_MC[HtBin]->Clone() );
        for(int bin = 1; bin < clone->GetNbinsX(); ++bin)
          if( (clone->GetBinCenter(bin) < smallestIntervalMins_MC.at(HtBin)) ||
              (clone->GetBinCenter(bin) > smallestIntervalMaxs_MC.at(HtBin)) )
            clone -> SetBinContent(bin,0.);
        clone -> SetFillColor(kYellow);
        clone -> SetLineWidth(0);
        clone -> Draw("HIST,same");
        
        line_mean_MC              -> Draw("same");
        line_recursiveMean_MC     -> Draw("same");
        line_recursiveMean_min_MC -> Draw("same");
        line_recursiveMean_max_MC -> Draw("same");
        line_smallestInterval_MC  -> Draw("same");
        
        
        TCanvas* c_DA = new TCanvas();
        c_DA -> cd();
        c_DA -> SetGridx();
        c_DA -> SetGridy();
        
        sprintf(axisTitle,"m_{ee} (GeV/c^{2})   -   H_{T} #in [%d,%d]",int(HtBinEdges->at(HtBin)),int(HtBinEdges->at(HtBin+1)));
        h_mee_HtBin_DA[HtBin] -> GetXaxis() -> SetTitle(axisTitle);
        sprintf(axisTitle,"event fraction");
        h_mee_HtBin_DA[HtBin] -> GetYaxis() -> SetTitle(axisTitle);
        h_mee_HtBin_DA[HtBin] -> GetXaxis() -> SetLabelSize(0.04);
        h_mee_HtBin_DA[HtBin] -> GetYaxis() -> SetLabelSize(0.04);
        h_mee_HtBin_DA[HtBin] -> GetXaxis() -> SetTitleSize(0.05);
        h_mee_HtBin_DA[HtBin] -> GetYaxis() -> SetTitleSize(0.05);
        h_mee_HtBin_DA[HtBin] -> SetLineColor(kBlack);
        h_mee_HtBin_DA[HtBin] -> SetLineWidth(1);
        h_mee_HtBin_DA[HtBin] -> GetXaxis() -> SetRangeUser(0.65,1.34999);
        h_mee_HtBin_DA[HtBin] -> SetMaximum( 1.1*h_mee_HtBin_DA[HtBin]->GetMaximum() );
        
        
        scale_mean_DA -> GetPoint(HtBin,x,y);
        TArrow* line_mean_DA = new TArrow(y,0.,y,h_mee_HtBin_DA[HtBin]->GetBinContent(h_mee_HtBin_DA[HtBin]->FindBin(y)));
        line_mean_DA -> SetLineColor(kBlack);
        line_mean_DA -> SetLineWidth(2);
        
        scale_recursiveMean_DA -> GetPoint(HtBin,x,y);
        TArrow* line_recursiveMean_DA = new TArrow(y,0.,y,h_mee_HtBin_DA[HtBin]->GetMaximum());
        line_recursiveMean_DA -> SetLineColor(kBlue);
        line_recursiveMean_DA -> SetLineWidth(2);
        TArrow* line_recursiveMean_min_DA = new TArrow(y-5./91.18,0.,y-5./91.18,h_mee_HtBin_DA[HtBin]->GetMaximum());
        line_recursiveMean_min_DA -> SetLineColor(kBlue);
        line_recursiveMean_min_DA -> SetLineWidth(2);
        TArrow* line_recursiveMean_max_DA = new TArrow(y+5./91.18,0.,y+5./91.18,h_mee_HtBin_DA[HtBin]->GetMaximum());
        line_recursiveMean_max_DA -> SetLineColor(kBlue);
        line_recursiveMean_max_DA -> SetLineWidth(2);
        
        scale_smallestInterval_DA -> GetPoint(HtBin,x,y);
        TArrow* line_smallestInterval_DA = new TArrow(y,0.,y,h_mee_HtBin_DA[HtBin]->GetBinContent(h_mee_HtBin_DA[HtBin]->FindBin(y)));
        line_smallestInterval_DA -> SetLineColor(kRed);
        line_smallestInterval_DA -> SetLineWidth(2);
        
        h_mee_HtBin_DA[HtBin] -> Draw("HIST");
        
        clone = (TH1F*)( h_mee_HtBin_DA[HtBin]->Clone() );
        for(int bin = 1; bin < clone->GetNbinsX(); ++bin)
          if( (clone->GetBinCenter(bin) < smallestIntervalMins_DA.at(HtBin)) ||
              (clone->GetBinCenter(bin) > smallestIntervalMaxs_DA.at(HtBin)) )
            clone -> SetBinContent(bin,0.);
        clone -> SetFillColor(kYellow);
        clone -> SetLineWidth(0);
        clone -> Draw("HIST,same");
        
        line_mean_DA              -> Draw("same");
        line_recursiveMean_DA     -> Draw("same");
        line_recursiveMean_min_DA -> Draw("same");
        line_recursiveMean_max_DA -> Draw("same");
        line_smallestInterval_DA  -> Draw("same");
        
        
        if( HtBin == 0 )
        {
          c_MC -> Print((outputPdf_MC+"[").c_str());
          c_DA -> Print((outputPdf_DA+"[").c_str());
        }
        {
          c_MC -> RedrawAxis();
          c_MC -> Print(outputPdf_MC.c_str());
          c_DA -> RedrawAxis();
          c_DA -> Print(outputPdf_DA.c_str());
        }
        if( HtBin == nHtBins-1 )
        {
          c_MC -> RedrawAxis();
          c_MC -> Print((outputPdf_MC+"]").c_str());
          c_DA -> RedrawAxis();
          c_DA -> Print((outputPdf_DA+"]").c_str());
        }
      }
    }
    
    
    
    baseDir -> cd();
    subDir = baseDir -> mkdir("Et_EtBin");
    subDir -> cd();
    
    for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
    {
      h_Et_EtBin_MC[EtBin] -> Scale(1./h_Et_EtBin_MC[EtBin]->Integral());
      h_Et_EtBin_MC[EtBin] -> Write();
      
      h_Et_EtBin_DA[EtBin] -> Scale(1./h_Et_EtBin_DA[EtBin]->Integral());
      h_Et_EtBin_DA[EtBin] -> Write();
      h_Et_EtBin_fit_DA[EtBin] -> Scale(1./h_Et_EtBin_fit_DA[EtBin]->Integral());
      h_Et_EtBin_fit_DA[EtBin] -> Write();
      h_Et_EtBin_gausFit_DA[EtBin] -> Scale(1./h_Et_EtBin_gausFit_DA[EtBin]->Integral());
      h_Et_EtBin_gausFit_DA[EtBin] -> Write();
      h_Et_EtBin_mean_DA[EtBin] -> Scale(1./h_Et_EtBin_mean_DA[EtBin]->Integral());
      h_Et_EtBin_mean_DA[EtBin] -> Write();
      h_Et_EtBin_recursiveMean_DA[EtBin] -> Scale(1./h_Et_EtBin_recursiveMean_DA[EtBin]->Integral());
      h_Et_EtBin_recursiveMean_DA[EtBin] -> Write();
      h_Et_EtBin_smallestInterval_DA[EtBin] -> Scale(1./h_Et_EtBin_smallestInterval_DA[EtBin]->Integral());
      h_Et_EtBin_smallestInterval_DA[EtBin] -> Write();
    }
    
    outFile -> cd();
  }
  
  
  
  outFile -> mkdir("Ht_HtBin");
  outFile -> cd("Ht_HtBin");
  
  for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
  {
    h_Ht_HtBin_MC[HtBin] -> Scale(1./h_Ht_HtBin_MC[HtBin]->Integral());
    h_Ht_HtBin_MC[HtBin] -> Write();
  }
  
  outFile -> cd();
  
  f_scaleVsEt -> Write();
  f_invScaleVsEt -> Write();
  
  h_scEta_MC -> Scale(1./h_scEta_MC->Integral());
  h_scEta_MC -> Write();
  h_scEta_DA -> Scale(1./h_scEta_DA->Integral());
  h_scEta_DA -> Write();
  
  h_R9_MC -> Scale(1./h_R9_MC->Integral());
  h_R9_MC -> Write();
  h_R9_DA -> Scale(1./h_R9_DA->Integral());
  h_R9_DA -> Write();
  
  h_Ht_MC -> Scale(1./h_Ht_MC->Integral());
  h_Ht_MC -> Write();
  h_Ht_DA -> Scale(1./h_Ht_DA->Integral());
  h_Ht_DA -> Write();
  
  h_mee_MC -> Scale(1./h_mee_MC->Integral());
  h_mee_MC -> Write();
  h_mee_DA -> Scale(1./h_mee_DA->Integral());
  h_mee_DA -> Write();
  
  outFile -> Close();
  
  
  
  return 0;
}






void MySetBins(TProfile* p, std::vector<double>& bins)
{
  double* binEdges = new double[bins.size()];
  for(unsigned int i = 0; i < bins.size(); ++i)
    binEdges[i] = bins.at(i);
  
  p -> SetBins(bins.size()-1,binEdges);
}



int MyFindBin(const double& val, const std::vector<double>* binEdges)
{
  for(unsigned int bin = 0; bin < binEdges->size()-1; ++bin)
  {
    if( (val >= binEdges->at(bin)) && (val < binEdges->at(bin+1)) )
      return bin;
  }
  return -1;
}



int MyFindBin(const double& val, const double& min, const double& max, const double& invWidth)
{
  if( val < min || val >= max ) return -1;
  return int( (val - min) * invWidth );
}



double MyEval(TGraph* g, const double& x)
{
  double x1,x2,y1,y2;
  g -> GetPoint(0,x1,y1);
  g -> GetPoint(g->GetN()-1,x2,y2);
  
  if( x  < x1 ) return y1;
  if( x >= x2 ) return y2;
  
  for(int point = 0; point < g->GetN()-1; ++point)
  {
    g -> GetPoint(point,  x1,y1);
    g -> GetPoint(point+1,x2,y2);
    
    if( x >= x1 && x < x2 )
    {
      return y1 + (y2-y1)/(x2-x1)*(x-x1);
    }
  }
  
  return 1.;
}



void MyFindFit(double& scale, double& scaleErr,
               TH1F* h_MC, TH1F* h_DA)
{
  float xNorm = h_DA->Integral() / h_MC->Integral() * h_DA->GetBinWidth(1) / h_MC->GetBinWidth(1);  
  h_MC -> Scale(xNorm);
  
  
  
  histoFunc* templateHistoFunc = new histoFunc(h_MC);
  char funcName[50];
  sprintf(funcName,"f_template");
  
  TF1* f_template = new TF1(funcName,templateHistoFunc,0.9,1.1,3,"histoFunc");
  
  f_template -> SetParName(0,"Norm"); 
  f_template -> SetParName(1,"Scale factor"); 
  f_template -> SetLineWidth(1); 
  f_template -> SetNpx(10000);
  f_template -> SetLineColor(kRed+2); 
  
  f_template->FixParameter(0,1.);
  f_template->SetParameter(1,0.99);
  f_template->FixParameter(2,0.);
  
  TFitResultPtr rp = h_DA -> Fit(funcName,"QERLS+");
  int fStatus = rp;
  int nTrials = 0;
  while( (fStatus != 0) && (nTrials < 10) )
  {
    rp = h_DA -> Fit(funcName,"QERLS+");
    fStatus = rp;
    if( fStatus == 0 ) break;
    ++nTrials;
  }
  
  double k   = f_template->GetParameter(1);
  double eee = f_template->GetParError(1); 
  
  scale = 1/k;
  scaleErr = eee/k/k;
  
  delete f_template;
}



void MyFindMean(double& mean, double& meanErr,
                std::vector<double>& vals, std::vector<double>& weights)
{
  std::cout << ">>>>>> MyFindMean" << std::endl;
  
  TH1F* h_temp = new TH1F("h_temp","",nBins_mee,meeMin,meeMax);
  h_temp -> Sumw2();
  
  for(unsigned int point = 0; point < vals.size(); ++point)
  {
    h_temp -> Fill( vals.at(point),weights.at(point) );
  }
  
  mean    = h_temp -> GetMean();
  meanErr = h_temp -> GetMeanError();
  delete h_temp;
}



void MyFindGausFit(double& mean, double& meanErr,
                   std::vector<double>& vals, std::vector<double>& weights,
                   TF1** f_gausFit, const std::string& name)
{
  std::cout << ">>>>>> MyFindGausFit" << std::endl;
  
  TH1F* h_temp = new TH1F("h_temp","",nBins_mee,meeMin,meeMax);
  h_temp -> Sumw2();
  
  for(unsigned int point = 0; point < vals.size(); ++point)
  {
    h_temp -> Fill( vals.at(point),weights.at(point) );
  }
  h_temp -> Scale(1./h_temp->Integral());
  
  mean    = h_temp -> GetMean();
  meanErr = h_temp -> GetRMS();
  
  TF1* f_temp = new TF1("f_temp","[0]*exp(-1.*(x-[1])*(x-[1])/2/[2]/[2])",0.75,1.25);
  f_temp -> SetParameters(h_temp->GetMaximum(),1.,meanErr);
  f_temp -> FixParameter(0,h_temp->GetMaximum());
  f_temp -> FixParameter(1,1.);
  f_temp -> SetParLimits(2,0.,1.);
  h_temp -> Fit("f_temp","QNR","");
  
  mean    = f_temp -> GetParameter(1);
  meanErr = f_temp -> GetParameter(2);
  
  (*f_gausFit) = new TF1(name.c_str(),"[0]*exp(-1.*(x-[1])*(x-[1])/2/[2]/[2])",mean-1.5*meanErr,mean+1.5*meanErr);
  (*f_gausFit) -> SetParameters(0.1,mean,meanErr);
  (*f_gausFit) -> SetParLimits(0,0.,1.);
  (*f_gausFit) -> SetParLimits(1,0.5,1.5);
  (*f_gausFit) -> SetParLimits(2,0.,1.);
  h_temp -> Fit(name.c_str(),"QNR","");
  
  mean    = (*f_gausFit) -> GetParameter(1);
  meanErr = (*f_gausFit) -> GetParError(1);
  
  delete f_temp;
  delete h_temp;
}



void MyFindRecursiveMean(double& mean, double& meanErr,
                         std::vector<double>& vals, std::vector<double>& weights,
                         const double& window, const double& tolerance)
{
  std::cout << ">>>>>> MyFindRecursiveMean" << std::endl;
  
  int trial = 0;
  double oldMean = -1.;
  double delta = 999999;
  
  while( delta > tolerance )
  {
    TH1F* h_temp = new TH1F("h_temp","",nBins_mee,meeMin,meeMax);
    h_temp -> Sumw2();
    
    for(unsigned int point = 0; point < vals.size(); ++point)
    {
      if( (trial == 0) && (fabs(vals.at(point) -      1.) > 2.*window) ) continue;
      if( (trial  > 0) && (fabs(vals.at(point) - oldMean) >    window) ) continue;
      
      h_temp -> Fill( vals.at(point),weights.at(point) );
    }
    
    mean    = h_temp -> GetMean();
    meanErr = h_temp -> GetMeanError();
    delete h_temp;
    
    if( fabs(mean-oldMean) > tolerance )
    {
      oldMean = mean;
      ++trial;
    }
    else
    {
      return;
    }
  }
  
}



void MyFindSmallestInterval(double& mean, double& meanErr, double& min, double& max,
                            std::vector<double>& vals, std::vector<double>& weights,
                            const double& fraction)
{
  std::cout << ">>>>>> MyFindSmallestInterval" << std::endl;
  std::sort(vals.begin(),vals.end());
  
  unsigned int nPoints = vals.size();
  unsigned int maxPoints = (unsigned int)(fraction * nPoints);
  
  unsigned int minPoint = 0;
  unsigned int maxPoint = 0;
  double delta = 999999.;
  for(unsigned int point = 0; point < nPoints-maxPoints; ++point)
  {
    double tmpMin = vals.at(point);
    double tmpMax = vals.at(point+maxPoints-1);
    if( tmpMax-tmpMin < delta )
    {
      delta = tmpMax - tmpMin;
      min = tmpMin;
      max = tmpMax;
      minPoint = point;
      maxPoint = point + maxPoints - 1;
    }
  }
  
  TH1F* h_temp = new TH1F("h_temp","",nBins_mee,meeMin,meeMax);
  h_temp -> Sumw2();
  
  for(unsigned int point = minPoint; point < maxPoint; ++point)
    h_temp -> Fill( vals.at(point),weights.at(point) );
  
  mean    = h_temp -> GetMean();
  meanErr = h_temp -> GetMeanError();
  delete h_temp;
}



double invScaleVsEt(double* x, double* par)
{
  double xx = x[0];
  
  if( xx  < EtBinEdges->at(0)       ) return 1. / f_scaleVsEt -> Eval( h_Et_EtBin_DA[0]        ->GetMean() );
  if( xx >= EtBinEdges->at(nEtBins) ) return 1. / f_scaleVsEt -> Eval( h_Et_EtBin_DA[nEtBins-1]->GetMean() );
  
  int EtBin = MyFindBin(xx,EtBinEdges);
  return 1. / f_scaleVsEt -> Eval(h_Et_EtBin_DA[EtBin]->GetMean());
  
  return 0.;
}



double deltaPhi(const double& phi1, const double& phi2)
{
  double deltaphi = fabs(phi1 - phi2);
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308 - deltaphi;
  return deltaphi;
}
