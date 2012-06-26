#include "Zutils.h"
#include "setTDRStyle.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TTree.h"
#include "TVirtualFitter.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>

#define a 0.5346
#define b 0.2166
#define FWHMZ 2.4952
 
#define NormalizationReReco  995487.
#define NormalizationPrompt  1.
#define NormalizationMC 1603411.

int main(int argc, char **argv){

   //set the style
  setTDRStyle();
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110); 


 /// Acquisition from cfg file
 
 if(argc != 2){
 std::cerr << " >>>>> analysis.cpp::usage: " << argv[0] << " configFileName" << std::endl ;
 return 1;
 }

 parseConfigFile (argv[1]) ;


 bool debug = gConfigParser->readBoolOption("Input::debug");
 if(debug) std::cout << " Start " << std::endl;
 if(debug) std::cout << " argc = " << argc << std::endl;
 if(debug) std::cout << " argv[0] = " << argv[0] << std::endl;
 if(debug) std::cout << " argv[1] = " << argv[1] << std::endl;
 
 std::string type;
 type = gConfigParser->readStringOption("Input::type");
 std::cout << " >>>>> Input::type = " << type << std::endl;
 std::cout << std::endl;

 float reScale = 1;
 if(type == "MC")reScale = 1.*NormalizationReReco/NormalizationMC;
 
 std::string treeName  = gConfigParser->readStringOption("Input::treeName");
 std::cout << " Input Tree Name = " << treeName << std::endl;

 std::string inputFileReReco  = gConfigParser->readStringOption("Input::inputFileReReco");
 std::string inputFilePrompt  = gConfigParser->readStringOption("Input::inputFilePrompt");
 std::string inputFileMC  = gConfigParser->readStringOption("Input::inputFileMC");

 std::cout << " Input file ReReco = " << inputFileReReco << std::endl;
 std::cout << " Input file Prompt = " << inputFilePrompt << std::endl;
 std::cout << " Input file MC = " << inputFileMC << std::endl;

 std::string WeightforMC =  gConfigParser->readStringOption("Input::WeightforMC");
 std::cout << " Weights for MC = " << WeightforMC << std::endl;

 bool vsRun = gConfigParser->readBoolOption("Input::vsRun");
 std::cout << " vsRun = " << vsRun << std::endl;
 bool allOthers = gConfigParser->readBoolOption("Input::allOthers");
 std::cout << " allOthers = " << allOthers << std::endl;

 unsigned int numTotBin = gConfigParser->readIntOption("Input::numTotBin");
 
 std::string outputFile = gConfigParser->readStringOption("Output::outputFile");
 std::cout << " Output File = " << outputFile << std::endl;
 std::string outputTable = gConfigParser->readStringOption("Output::outputFileTable");
 std::cout << " Output Table File = " << outputTable << std::endl;


 /// Input data infos
  TChain* tree = new TChain(treeName.c_str());

  if(type == "ReReco") FillChain(*tree,inputFileReReco.c_str());
  if(type == "Prompt") FillChain(*tree,inputFilePrompt.c_str());
  if(type == "MC") FillChain(*tree,inputFileMC.c_str());

  std::cout << " type = " << type << ": " << std::setw(8) << tree->GetEntries() << " entries" << std::endl; 
  
  if(tree->GetEntries() == 0 ){
    std::cout << ">>>recalibZ::Error: at least one file is empty" << std::endl; 
    return -1;
  }

  std::vector<std::string> FitCategories;
  FitCategories = gConfigParser->readStringListOption("Input::FitCategories");
  
  std::cout << " >>>>> Input::FitCategories size = " << FitCategories.size() << std::endl;  
  std::cout << " >>>>> >>>>>  "; 
  for (unsigned int iCat = 0; iCat < FitCategories.size(); iCat++){
    std::cout << " " << FitCategories.at(iCat) << ", ";
  }
  std::cout << std::endl; 

  std::vector<std::string> FitCategoriesVsRun;
  FitCategoriesVsRun = gConfigParser->readStringListOption("Input::FitCategoriesVsRun");
  
  std::cout << " >>>>> Input::FitCategoriesVsRun size = " << FitCategoriesVsRun.size() << std::endl;  
  std::cout << " >>>>> >>>>>  "; 
  for (unsigned int iCat = 0; iCat < FitCategoriesVsRun.size(); iCat++){
    std::cout << " " << FitCategoriesVsRun.at(iCat) << ", ";
  }
  std::cout << std::endl; 
 

  //--- weights for MC
  TFile weightsFile (WeightforMC.c_str(),"READ"); 
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();
  

 /// option Infos
  int nbinZ  = gConfigParser->readIntOption("Option::nbinZ");
  double mZ_Max = gConfigParser->readDoubleOption("Option::mZMax");
  double mZ_Min = gConfigParser->readDoubleOption("Option::mZMin");
  double scaleEB = gConfigParser->readDoubleOption("Option::scaleEB");
  double scaleEE = gConfigParser->readDoubleOption("Option::scaleEE");
  int nPoints  = gConfigParser->readIntOption("Option::nPoints");

  if(debug){
  std::cout<<" nbinZ = "<<nbinZ<<std::endl;
  std::cout<<" mZ_Max = "<<mZ_Max<<std::endl;
  std::cout<<" mZ_Min = "<<mZ_Min<<std::endl;
  std::cout<<" scaleEB = "<<scaleEB<<std::endl;
  std::cout<<" scaleEE = "<<scaleEE<<std::endl;
  std::cout<<" nPoints = "<<nPoints<<std::endl;
  }

 ///**** Book histos
 if(debug) std::cout << "Book histos" <<  std::endl;
 std::map<int, std::vector<int> > runIdEventId_DATA; 
 std::vector<int> runIdVec;
 std::map<std::string, std::vector<TH1F*> > Zmass_regression_runId;
 if(debug) std::cout << "Book histos  1 " <<  std::endl;
 std::map<std::string, TH1F* > Zmass;
 std::map<std::string, TH1F* > Zmass_regression;

 if(debug) std::cout << "Book histos  2 " <<  std::endl;

 std::string VtxIstoName = "h_Vtx_"+type;
 TH1F* h_Vtx = new TH1F(VtxIstoName.c_str(),"", 100, 0., 100.);

 TH1F* Zscale = new TH1F("Zscale","", int(FitCategories.size()), 0., FitCategories.size());
 TH1F* ZScb = new TH1F("ZScb","", int(FitCategories.size()), 0., FitCategories.size());
 TH1F* ZSol = new TH1F("ZSol","", int(FitCategories.size()), 0., FitCategories.size());

 if(debug) std::cout << "Book histos  3 " <<  std::endl;
 std::map<std::string, TH1F* > Zscale_runId;
 std::map<std::string, TH1F* > ZScb_runId;
 std::map<std::string, TH1F* > ZSol_runId;

 if(allOthers){
 if(debug) std::cout << "Book histos  allOthers " <<  std::endl;
   for(unsigned int i = 0; i < FitCategories.size(); ++i){
     std::string category = FitCategories.at(i);

     std::string histoName = "h_Zmass"+type+"_"+category;
     Zmass[category] = new TH1F(histoName.c_str(),"", nbinZ, mZ_Min, mZ_Max);
     Zmass[category]->Sumw2();
     
     std::string histoName2 = "h_Zmass"+type+"_regression_"+category;
     Zmass_regression[category] = new TH1F(histoName2.c_str(),"", nbinZ, mZ_Min, mZ_Max);
     Zmass_regression[category]->Sumw2();
   }
 }

 if(vsRun) {
   if(debug) std::cout << " vsRun define histos " << std::endl; 

   char binRun[100];
   std::stringstream dummy;
   std::string Sdummy;

   for(unsigned int i = 0; i < FitCategoriesVsRun.size(); ++i){
     std::string category = FitCategoriesVsRun.at(i);

     for(unsigned int jj=0; jj<numTotBin; ++jj){
  
       if(debug)std::cout << " jj  " << jj << std::endl;
       sprintf(binRun, "runId_bin%d", int(jj));
       if(debug)       printf(binRun);
       dummy << binRun;
       if(debug) std::cout << " dummy = " << dummy <<  std::endl;
       dummy >> Sdummy;
       dummy.clear();
       if(debug) std::cout << " Sdummy " << Sdummy << std::endl;
       std::string histoName = "h_Zmass"+type+"_regression_"+category+"_"+Sdummy;

       if(debug) std::cout << " histoName " << histoName << std::endl;
       Zmass_regression_runId[category].push_back( new TH1F(histoName.c_str(), "", nbinZ, mZ_Min, mZ_Max) );
       if(debug) std::cout << " histoName 2 " << histoName << std::endl;
       Zmass_regression_runId[category].at(jj)->Sumw2();

       if(debug) std::cout << " Sdummy fatto " << std::endl;
     }

     if(debug) std::cout << " vsRun fine categories " << std::endl;

     std::string histoNamevsRun1 = "h_Zscale"+type+"_regression_"+category+"_vsRunId";
     Zscale_runId[category] = new TH1F(histoNamevsRun1.c_str(), "", numTotBin, 0., numTotBin);

     std::string histoNamevsRun2 = "h_ZScb"+type+"_regression_"+category+"_vsRunId";
     ZScb_runId[category] = new TH1F(histoNamevsRun2.c_str(), "", numTotBin, 0., numTotBin);

     std::string histoNamevsRun3 = "h_ZSol"+type+"_regression_"+category+"_vsRunId";
     ZSol_runId[category] = new TH1F(histoNamevsRun3.c_str(), "", numTotBin, 0., numTotBin);
     
   }
 }



 /// Set branch addresses
 if(debug) std::cout << "Set branch addresses" <<  std::endl;
 int isZ;
 float ele1ele2_scM,ele1ele2_scM_regression;
 int ele1_isEB,ele2_isEB;
 float ele1_scEta,ele2_scEta,ele1_scE,ele2_scE,ele1_es,ele2_es,ele1_scERaw,ele2_scERaw, ele1_scE_regression,
       ele2_scE_regression,ele1_e3x3,ele2_e3x3;
 int ele1_seedIeta,ele1_seedIphi,ele2_seedIeta,ele2_seedIphi,ele1_seedIx,ele2_seedIx,ele1_seedIy,ele2_seedIy,ele1_seedZside,ele2_seedZside;
 int PUit_NumInteractions;
 int PV_n;
 int runId;

 tree->SetBranchAddress("isZ", &isZ);
 tree->SetBranchAddress("ele1ele2_scM", &ele1ele2_scM);
 tree->SetBranchAddress("ele1ele2_scM_regression", &ele1ele2_scM_regression);
 tree->SetBranchAddress("ele1_isEB",  &ele1_isEB);
 tree->SetBranchAddress("ele2_isEB",  &ele2_isEB);
 if(type == "MC") tree->SetBranchAddress("PUit_NumInteractions", &PUit_NumInteractions);
 tree->SetBranchAddress("PV_n", &PV_n);
 tree->SetBranchAddress("runId", &runId);

 tree->SetBranchAddress("ele1_scEta", &ele1_scEta);
 tree->SetBranchAddress("ele2_scEta", &ele2_scEta);

 tree->SetBranchAddress("ele1_seedIeta", &ele1_seedIeta);
 tree->SetBranchAddress("ele1_seedIphi", &ele1_seedIphi);
 tree->SetBranchAddress("ele2_seedIeta", &ele2_seedIeta);
 tree->SetBranchAddress("ele2_seedIphi",  &ele2_seedIphi);
 tree->SetBranchAddress("ele1_seedIx",  &ele1_seedIx);
 tree->SetBranchAddress("ele2_seedIx",  &ele2_seedIx);
 tree->SetBranchAddress("ele1_seedIy",  &ele1_seedIy);
 tree->SetBranchAddress("ele2_seedIy",  &ele2_seedIy);
 tree->SetBranchAddress("ele1_seedZside",  &ele1_seedZside);
 tree->SetBranchAddress("ele2_seedZside",  &ele2_seedZside);


 tree->SetBranchAddress("ele1_scE", &ele1_scE);
 tree->SetBranchAddress("ele2_scE", &ele2_scE);
 tree->SetBranchAddress("ele1_scERaw", &ele1_scERaw);
 tree->SetBranchAddress("ele2_scERaw", &ele2_scERaw);
 tree->SetBranchAddress("ele1_e3x3", &ele1_e3x3);
 tree->SetBranchAddress("ele2_e3x3", &ele2_e3x3);
 tree->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression);
 tree->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression);
 tree->SetBranchAddress("ele1_es", &ele1_es);
 tree->SetBranchAddress("ele2_es", &ele2_es);

 //*** Loop on tree **//
 if(debug)  std::cout << " Fill with " << type << " Events " <<std::endl;
 int nEntries = tree->GetEntries();
 std::cout << " Events = " << nEntries << std::endl;
 std::vector<std::string> eventCategory;
 int nEntries_Z = 0;
 for(int iEntry=0; iEntry<nEntries; ++iEntry){

  if( (iEntry%100000 == 0) ) std::cout << " reading saved entry " << iEntry << std::endl;
  tree -> GetEntry(iEntry);
  double weight = 1.;
  if(type == "MC") weight = w[PUit_NumInteractions];


  //only the Z
  if (isZ != 1) continue;
  ++nEntries_Z;
  h_Vtx->Fill(PV_n, weight);


  if(vsRun) runIdEventId_DATA[runId].push_back(iEntry);

   if(!allOthers) continue;

   eventCategory.clear();

   float ele1_R9 = ele1_e3x3/ele1_scERaw; 
   float ele2_R9 = ele2_e3x3/ele2_scERaw; 

   if( (ele1_seedZside == 0) && (ele2_seedZside == 0) )  eventCategory.push_back("EB-EB");
   else if( fabs(ele1_seedZside) == 1 && fabs(ele2_seedZside) == 1 ) eventCategory.push_back("EE-EE");
   else eventCategory.push_back("EB-EE");

   if( (ele1_seedZside == 0) && (ele2_seedZside == 0) && ele1_scEta > 0. && ele2_scEta > 0.) eventCategory.push_back("EBp");
   else if( (ele1_seedZside == 0) && (ele2_seedZside == 0) && ele1_scEta < 0. && ele2_scEta < 0.) eventCategory.push_back("EBm");
   
   if( (ele1_seedZside == 1) && (ele2_seedZside == 1)) eventCategory.push_back("EEp");
   else if( (ele1_seedZside == -1) && (ele2_seedZside == -1)) eventCategory.push_back("EEm");

   if( (ele1_seedZside != 0) && (ele2_seedZside != 0) && ele1_R9 > 0.94 && ele2_R9 > 0.94) eventCategory.push_back("EE_R9_g");
   else if( (ele1_seedZside != 0) && (ele2_seedZside != 0) && ele1_R9 < 0.94 && ele2_R9 < 0.94) eventCategory.push_back("EE_R9_l");

   if( (ele1_seedZside == 0) && (ele2_seedZside == 0) && ele1_R9 > 0.94 && ele2_R9 > 0.94) eventCategory.push_back("EB_R9_g");
   else if( (ele1_seedZside == 0) && (ele2_seedZside == 0) && ele1_R9 < 0.94 && ele2_R9 < 0.94)eventCategory.push_back("EB_R9_l");
   
   if( (ele1_seedZside == 0) && (ele2_seedZside == 0) && fabs(ele1_scEta) < 1. && fabs(ele2_scEta) < 1. ) eventCategory.push_back("EB_Eta_l");
   else if( (ele1_seedZside == 0) && (ele2_seedZside == 0) && fabs(ele1_scEta) > 1. && fabs(ele2_scEta) > 1. ) eventCategory.push_back("EB_Eta_g");

   if( (fabs(ele1_seedZside) == 1) && (fabs(ele2_seedZside) == 1) && fabs(ele1_scEta) < 2. && fabs(ele2_scEta) < 2. ) eventCategory.push_back("EE_Eta_l");
   else if( (fabs(ele1_seedZside) == 1) && (fabs(ele2_seedZside) == 1) && 
	    fabs(ele1_scEta) > 2. && fabs(ele2_scEta) > 2. ) eventCategory.push_back("EE_Eta_g");

   
   for(unsigned int i=0; i<FitCategories.size(); ++i){
     for(unsigned int j=0; j<eventCategory.size(); ++j){
       if(FitCategories.at(i) == eventCategory.at(j)){
	 Zmass[eventCategory.at(j)]->Fill( ele1ele2_scM * sqrt(scaleEB*scaleEB));
	 Zmass_regression[eventCategory.at(j)]->Fill( ele1ele2_scM * sqrt((ele1_scE_regression/ele1_scE)*(ele2_scE_regression/ele2_scE)) * 
							sqrt(scaleEB*scaleEB));
       }
     }
   }

 }


 if(vsRun){
   float numero = (1.*nEntries_Z)/ numTotBin;
   //   int nBin = 0;
   int nEvents = 0;
   std::vector<std::string> eventCategory_vsRun;

   for(std::map<int, std::vector<int> >::iterator mapIt=runIdEventId_DATA.begin(); mapIt!=runIdEventId_DATA.end(); ++mapIt){
     std::vector<int> dummy = mapIt->second;
     for(unsigned int jj=0; jj<dummy.size(); ++jj){
       tree->GetEntry(dummy.at(jj));

       double weight = 1.;
       if(type == "MC") weight = w[PUit_NumInteractions];

       int nBin = static_cast<int>( nEvents/numero );
       //       if(debug) std::cout << " nBin " << nBin << std::endl;
//        if(debug && nBin != 0){
// 	 std::cout << " int(nEvents/nBin)" << int(nEvents/nBin) << std::endl;
// 	 std::cout << " nEvents % int(nEvents/nBin) " << nEvents % int(nEvents/nBin) << std::endl;
//        }
       if(nEvents == 0 || (nBin != 0 && nEvents % int(numero) == 0)  || nEvents == nEntries_Z-1) {
         if(debug){
	   std::cout << " nEvents " << nEvents << std::endl;
	   std::cout << " numero " << numero << std::endl;
	   std::cout << " int nEvents/numero " << int(nEvents/numero) << std::endl;
	   std::cout << " runId " << runId << std::endl;
	   std::cout << " nBin =  " << nBin << std::endl;
         }
         runIdVec.push_back(runId);
       }

       eventCategory_vsRun.clear();

       if( (ele1_seedZside == 0) && (ele2_seedZside == 0) ){ eventCategory_vsRun.push_back("EB-EB");
// 	 Zmass_regression_runId[eventCategory_vsRun].at(nBin)->Fill( ele1ele2_scM * sqrt((ele1_scE_regression/ele1_scE)*(ele2_scE_regression/ele2_scE)) *
// 								     sqrt(scaleEB*scaleEB));
       }
       if( (fabs(ele1_seedZside) == 1) && (fabs(ele2_seedZside) == 1) ){eventCategory_vsRun.push_back("EE-EE");
// 	 Zmass_regression_runId[eventCategory_vsRun].at(nBin)->Fill( ele1ele2_scM * sqrt((ele1_scE_regression/ele1_scE)*(ele2_scE_regression/ele2_scE)) *
// 								    sqrt(scaleEE*scaleEE));
       }

       for(unsigned int i=0; i<FitCategoriesVsRun.size(); ++i){
	 for(unsigned int j=0; j<eventCategory_vsRun.size(); ++j){
	   if(FitCategoriesVsRun.at(i) == eventCategory_vsRun.at(j)){
	     Zmass_regression_runId[eventCategory_vsRun.at(j)].at(nBin)->Fill( ele1ele2_scM * sqrt(scaleEB*scaleEB),weight );
	   }
	 }
       }

       ++nEvents;
       if(debug)       std::cout << " runId = " << runId << std::endl;
     }
   }

   if(debug)   std::cout << " nEvents = " << nEvents << std::endl;
 }
 if(debug) std::cout << " runIdVec.size() = " << runIdVec.size() << std::endl;


 /// Z Lineshape Tool
if(allOthers){
  std::string nameTable = outputTable+"_"+type+"_vsCategories.txt";
   std::ofstream outTableFile (nameTable.c_str(),std::ios::out);

   outTableFile<<"\\begin{table}[!htb]"<<std::endl;
   outTableFile<<"\\begin{center}"<<std::endl;

   outTableFile<<"\\begin{tabular}{c|c|c|c|c|c}"<<std::endl;

   outTableFile<<"\\hline"<<std::endl;
   outTableFile<<"\\hline"<<std::endl;

   outTableFile << " & " << " & " << " & " << " & " << " & " << " \\\\ " <<std::endl;

   outTableFile << " Category " << " & " << " type " << " & " << " Event " << " & " 
		<< "$ \\Delta M_" << type << "}$ " << " & " 
		<< "\\sigma_{cb}^{" << type << "}/" <<  "\\Delta M_{" << type << "}$ " << " & "
		<< "\\sigma_{OL}^{" << type << "}/" <<  "\\Delta M_{" << type << "}$ " << " \\\\ " << std::endl; 


   outTableFile << " & " << " & " << " & " << " & " << " & " << " \\\\ " <<std::endl;


   outTableFile << "\\hline" <<std::endl;
   outTableFile << "\\hline" <<std::endl;


   for(unsigned int i = 0; i < FitCategories.size(); ++i){
     std::string cat = FitCategories.at(i);

     std::string energyType = "Reg";

     if(debug) std::cout << " Fitting: energyType = " << energyType << " category = " << cat << " type  = " << type << std::endl;
     BinnedFitZPeak(cat, 1, Zmass_regression[cat], type, reScale, nPoints, mZ_Min, mZ_Max, energyType);

     
     std::string funcName = "bw_cb_"+cat+"_"+energyType;
     std::pair<double,double> extreme = 
       breitWigner_crystalBallLowFWHM(Zmass_regression[cat]->GetFunction(funcName.c_str()), mZ_Min,mZ_Max);
     double sigma = sqrt(TMath::Power(((extreme.second-extreme.first)-a*FWHMZ),2)-b*FWHMZ*FWHMZ)/(2.*sqrt(2.*log(2.)));

     outTableFile << " & " << " & " << " & " << " & " << " & " << " \\\\ " <<std::endl;

  int Entries = Zmass_regression[cat]->GetEntries(); 
  float DM = Zmass_regression[cat]->GetFunction(funcName.c_str())->GetParameter(3);
  float DM_err = Zmass_regression[cat]->GetFunction(funcName.c_str())->GetParError(3);
  float Scb = Zmass_regression[cat]->GetFunction(funcName.c_str())->GetParameter(4);
  float Scb_err = Zmass_regression[cat]->GetFunction(funcName.c_str())->GetParError(4);
  float Sol_err = (0.001 *(extreme.second-extreme.first - a*FWHMZ*FWHMZ) / 2. /
		      (sqrt( 2.*log(2) * ( pow(extreme.second-extreme.first - a*FWHMZ*FWHMZ, 2) - b*FWHMZ*FWHMZ )) ) );

  outTableFile << cat << "-" << energyType << " & " << type << " & "<< Entries << " & " 
	       << DM << " $\\pm$ " << DM_err << " & "  
	       << Scb / (91.18+DM) << " & " 
	       << sigma / (91.18+DM) <<  " \\\\ " <<std::endl;
     

  Zscale->SetBinContent(i+1, DM);
  Zscale->SetBinError(i+1, DM_err);
  ZScb->SetBinContent(i+1, Scb / (91.18+DM));
  ZScb->SetBinError(i+1, sqrt(Scb_err*Scb_err + DM_err*DM_err) / (91.18+DM) );
  ZSol->SetBinContent(i+1, sigma / (91.18+DM));
  ZSol->SetBinError(i+1, sqrt(Sol_err*Sol_err + DM_err*DM_err) / (91.18+DM));

  energyType = "NoReg";

  if(debug) std::cout << " Fitting: energyType = " << energyType << " category = " << cat << " type  = " << type << std::endl;
  BinnedFitZPeak(cat, 1, Zmass[cat], type, reScale, nPoints, mZ_Min, mZ_Max, energyType);
  
  funcName = "bw_cb_"+cat+"_"+energyType;
  extreme = breitWigner_crystalBallLowFWHM(Zmass[cat]->GetFunction(funcName.c_str()), mZ_Min,mZ_Max);
  sigma = sqrt(TMath::Power(((extreme.second-extreme.first)-a*FWHMZ),2)-b*FWHMZ*FWHMZ)/(2.*sqrt(2.*log(2.)));


  outTableFile << " & " << " & " << " & " << " & " << " & " << " \\\\ " <<std::endl;

  Entries = Zmass[cat]->GetEntries(); 
  DM = Zmass[cat]->GetFunction(funcName.c_str())->GetParameter(3);
  DM_err = Zmass[cat]->GetFunction(funcName.c_str())->GetParError(3);
  Scb = Zmass[cat]->GetFunction(funcName.c_str())->GetParameter(4);

  outTableFile << cat << "-" << energyType << " & " << type << " & "<< Entries << " & "
	       << DM << " $\\pm$ " << DM_err << " & "
               << Scb / (91.18+DM) << " & "
               << sigma / (91.18+DM) <<  " \\\\ " <<std::endl;


  if(debug) std::cout << " END!! Fitting: category " << cat << " energyType = " << energyType << std::endl;
   }

   outTableFile << " \\hline " << std::endl;
   outTableFile << " \\hline " << std::endl;
   
  outTableFile << " \\end{tabular} " << std::endl;
  outTableFile << " \\end{center} " << std::endl;
  outTableFile << " \\end{table} " << std::endl;
 }


//vsRun
 if(vsRun){
   std::string nameTable = outputTable+"_"+type+"_vsRun.txt";
   std::ofstream outTableFile (nameTable.c_str(),std::ios::out);  
 
   outTableFile << " \\begin{table}[!htb] " << std::endl;
   outTableFile << " \\begin{center} " << std::endl;

   outTableFile << "\\begin{tabular}{c|c|c|c|c|c}" << std::endl;

   outTableFile << "\\hline" << std::endl;
   outTableFile << "\\hline" << std::endl;

   outTableFile << " & " << " & " << " & " << " & " << " & " << " & " << " \\\\" << std::endl;

   outTableFile << " Category" << " & " << "Run-bin " << " & " << " type " << " & " << " Event " << " & " 
		<< " $\\Delta M_{" << type << "}$ "  << " & " 
		<< " $\\sigma_{cb}^{" << type << "}/ (M + \\Delta M_{" << type << "}$ " << " & "
		<< " $\\sigma_{ol}^{" << type << "}/ (M + \\Delta M_{" << type << "}$ " << " \\\\ " << std::endl;

   outTableFile << "\\hline" << std::endl;
   outTableFile << "\\hline" << std::endl;

   outTableFile << " & " << " & " << " & " << " & " << " & " << " & " << " \\\\ " << std::endl;

   for(unsigned int jj=0; jj<FitCategoriesVsRun.size(); ++jj){
     std::string cat = FitCategoriesVsRun.at(jj);
     for(unsigned int ii=0; ii<numTotBin; ++ii){
     
     std::string energyType = "Reg";  

     if(debug) std::cout << " Fitting vs RUN: energyType = " << energyType << " category = " << cat << " type  = " << type << std::endl;
     BinnedFitZPeak(cat, 1, Zmass_regression_runId[cat].at(ii), type, reScale, nPoints, mZ_Min, mZ_Max, energyType);    


     std::string funcName = "bw_cb_"+cat+"_"+energyType;
     std::pair<double,double> extreme = 
       breitWigner_crystalBallLowFWHM(Zmass_regression_runId[cat].at(ii)->GetFunction(funcName.c_str()),mZ_Min,mZ_Max);
     double sigma = sqrt(TMath::Power((extreme.second-extreme.first)-a*FWHMZ,2)-b*FWHMZ*FWHMZ)/(2.*sqrt(2.*log(2.)));

     outTableFile << " & " << " & " << " & " << " & " << " & " << " & " << " \\\\ " << std::endl;

     int Entries = Zmass_regression_runId[cat].at(ii)->GetEntries(); 
     float DM = Zmass_regression_runId[cat].at(ii)->GetFunction(funcName.c_str())->GetParameter(3);
     float DM_err = Zmass_regression_runId[cat].at(ii)->GetFunction(funcName.c_str())->GetParError(3);
     float Scb = Zmass_regression_runId[cat].at(ii)->GetFunction(funcName.c_str())->GetParameter(4);
     float Scb_err = Zmass_regression_runId[cat].at(ii)->GetFunction(funcName.c_str())->GetParError(4);
     float Sol_err = (0.001 *(extreme.second-extreme.first - a*FWHMZ*FWHMZ) / 2. /
		      (sqrt( 2.*log(2) * ( pow(extreme.second-extreme.first - a*FWHMZ*FWHMZ, 2) - b*FWHMZ*FWHMZ )) ) );

     outTableFile << cat << " " << energyType << " & " << runIdVec.at(ii) << "-" << runIdVec.at(ii+1) << " & " << type << " & " << Entries << " & " 
	       << DM << " $\\pm$ " << DM_err << " & " 
	       << Scb / (91.18+DM) << " & " 
	       << sigma / (91.18+DM) << " \\\\ " <<std::endl;


  Zscale_runId[cat]->SetBinContent(ii+1, DM);
  Zscale_runId[cat]->SetBinError(ii+1, DM_err);
  ZScb_runId[cat]->SetBinContent(ii+1, Scb / (91.18+DM));
  ZScb_runId[cat]->SetBinError(ii+1, sqrt(Scb_err*Scb_err + DM_err*DM_err) / (91.18+DM) );
  ZSol_runId[cat]->SetBinContent(ii+1, sigma / (91.18+DM));
  ZSol_runId[cat]->SetBinError(ii+1, sqrt(Sol_err*Sol_err + DM_err*DM_err) / (91.18+DM));

     }
   }
 }


//

 /// Output save 
 std::cout << "recalibZ::Saving and Closing" << std::endl;

 /// Output infos
 std::string outputFileName = outputFile+"_"+type+".root";
 TFile* outputTFile = new TFile(outputFileName.c_str(),"RECREATE");

 if(allOthers){
 Zscale->Write();
 ZScb->Write();
 ZSol->Write();

   for(unsigned int i = 0; i < FitCategories.size(); ++i){
     std::string cat = FitCategories.at(i);

     Zmass[cat]->Write();
     Zmass_regression[cat]->Write();
   }

   h_Vtx->Write();
 }

 if(vsRun){

   for(unsigned int jj=0; jj<FitCategoriesVsRun.size(); ++jj){
     std::string cat = FitCategoriesVsRun.at(jj);
     for(unsigned int ii=0; ii<numTotBin; ++ii){
       Zmass_regression_runId[cat].at(ii)->Write();
     }

     Zscale_runId[cat]->Write();
     ZScb_runId[cat]->Write();
     ZSol_runId[cat]->Write();
   }
 }

 outputTFile -> Close();

 if(debug) std::cout << "recalibZ::Closed now Printing" << std::endl;

   
  std::string outDir = "summaryPlots_"+type;
  Color_t iStyle = kGreen+2;

  if(type == "MC") iStyle = kRed+1;

   if(allOthers){

     Zscale->SetLineColor(iStyle);
     ZScb->SetLineColor(iStyle);
     ZSol->SetLineColor(iStyle);

     Zscale->SetMarkerColor(iStyle);
     ZScb->SetMarkerColor(iStyle);
     ZSol->SetMarkerColor(iStyle);

     Zscale->SetLineWidth(2);
     ZScb->SetLineWidth(2);
     ZSol->SetLineWidth(2);
     
     Zscale->SetMarkerStyle(7);
     ZScb->SetMarkerStyle(7);
     ZSol->SetMarkerStyle(7);
     

//      for(unsigned int i = 0; i < FitCategories.size(); ++i)
//        std::cout << FitCategories.at(i).c_str() << std::endl;

     TCanvas* cScale = new TCanvas();
     gPad->SetGrid();
     Zscale->GetYaxis()->SetRangeUser(-1., 2.);
     Zscale->GetYaxis()->SetTitle("#Delta M");
     Zscale->Draw();
     for(unsigned int i = 0; i < FitCategories.size(); ++i)
       Zscale->GetXaxis()->SetBinLabel(i+1, FitCategories.at(i).c_str());
     cScale->Print((outDir+"/Zscale_"+type+".png").c_str(),"png");
     
     TCanvas* cScb = new TCanvas();
     gPad->SetGrid();
     ZScb->GetYaxis()->SetRangeUser(0.005, 0.05);
     ZScb->GetYaxis()->SetTitle("#sigma_{cb} / (M + #Delta M)");
     ZScb->Draw();
     for(unsigned int i = 0; i < FitCategories.size(); ++i)
       ZScb->GetXaxis()->SetBinLabel(i+1, FitCategories.at(i).c_str());
     cScb->Print((outDir+"/ZScb_"+type+".png").c_str(),"png");
     
     TCanvas* cSol = new TCanvas();
     gPad->SetGrid();
     ZSol->GetYaxis()->SetRangeUser(0.005, 0.05);
     ZSol->GetYaxis()->SetTitle("#sigma_{ol} / (M + #Delta M)");
     ZSol->Draw();
     for(unsigned int i = 0; i < FitCategories.size(); ++i)
       ZSol->GetXaxis()->SetBinLabel(i+1, FitCategories.at(i).c_str());
     cSol->Print((outDir+"/ZSol_"+type+".png").c_str(),"png");
   }


 if(vsRun){

   TCanvas* cScale_runId[FitCategoriesVsRun.size()];
   TCanvas* cScb_runId[FitCategoriesVsRun.size()];
   TCanvas* cSol_runId[FitCategoriesVsRun.size()];
   char labelName[100];
     
   for(unsigned int ii=0; ii<FitCategoriesVsRun.size(); ++ii){
     std::string cat = FitCategoriesVsRun.at(ii);
     Zscale_runId[cat]->SetLineColor(iStyle);
     ZScb_runId[cat]->SetLineColor(iStyle);
     ZSol_runId[cat]->SetLineColor(iStyle);
   
     Zscale_runId[cat]->SetMarkerColor(iStyle);
     ZScb_runId[cat]->SetMarkerColor(iStyle);
     ZSol_runId[cat]->SetMarkerColor(iStyle);

     Zscale_runId[cat]->SetLineWidth(2);
     ZScb_runId[cat]->SetLineWidth(2);
     ZSol_runId[cat]->SetLineWidth(2);

     Zscale_runId[cat]->SetMarkerStyle(7);
     ZScb_runId[cat]->SetMarkerStyle(7);
     ZSol_runId[cat]->SetMarkerStyle(7);


     cScale_runId[ii] = new TCanvas;
     gPad->SetGrid();
     //   ZscaleDATA_EB_runId->GetYaxis()->SetRangeUser(-1., 2.);
     Zscale_runId[cat]->GetYaxis()->SetTitle("#Delta M");
     Zscale_runId[cat]->Draw();
     for(unsigned int yy=0; yy<numTotBin; ++yy){
       sprintf(labelName, "%d-%d", runIdVec.at(yy), runIdVec.at(yy+1));
       Zscale_runId[cat]->GetXaxis()->SetBinLabel(yy+1, labelName);
     }
     cScale_runId[ii]->Print((outDir+"/Zscale_"+type+"_"+cat+"_vsRunId.png").c_str(),"png");

     cScb_runId[ii] = new TCanvas;
     gPad->SetGrid();
     //   ZscaleDATA_EB_runId->GetYaxis()->SetRangeUser(-1., 2.);
     ZScb_runId[cat]->GetYaxis()->SetTitle("#sigma_{cb} / (M + #Delta M)");
     ZScb_runId[cat]->Draw();
     for(unsigned int yy=0; yy<numTotBin; ++yy){
       sprintf(labelName, "%d-%d", runIdVec.at(yy), runIdVec.at(yy+1));
       ZScb_runId[cat]->GetXaxis()->SetBinLabel(yy+1, labelName);
     }
     cScb_runId[ii]->Print((outDir+"/ZScb_"+type+"_"+cat+"_vsRunId.png").c_str(),"png");

     cSol_runId[ii] = new TCanvas;
     gPad->SetGrid();
     //   ZscaleDATA_EB_runId->GetYaxis()->SetRangeUser(-1., 2.);
     ZSol_runId[cat]->GetYaxis()->SetTitle("#sigma_{ol} / (M + #Delta M)");
     ZSol_runId[cat]->Draw();
     for(unsigned int yy=0; yy<numTotBin; ++yy){
       sprintf(labelName, "%d-%d", runIdVec.at(yy), runIdVec.at(yy+1));
       ZSol_runId[cat]->GetXaxis()->SetBinLabel(yy+1, labelName);
     }
     cSol_runId[ii]->Print((outDir+"/ZSol_"+type+"_"+cat+"_vsRunId.png").c_str(),"png");
     
   }
 }
   
 return 0;
}
