//****** simple macro to compute PU weights ******
{
  /*
  //*** mc file - 2012
  TChain* ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
  ntu->Add("/media/DATA/ALCARECO/DYSummer12/DYJets-Summer12.root");
  TH1F* hmc = new TH1F("hmc","hmc",100,0.,100.);
  ntu->Draw("PUit_TrueNumInteractions >> hmc","","goff");
  */

  //*** mc file - 2011

  //   TChain* ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
   //   ntu->Add("/tmp/amartell/DYJToLL_M50_TuneZ2S_8TeV-mad_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
   TChain* ntu = new TChain("simpleNtupleEoverPSh/SimpleNtupleEoverP");
   ntu->Add("/tmp/amartell/DYJets-Summer12-START53-noSkim.root");
   TH1F* hmc = new TH1F("hmc","hmc",100,0.,100.);
   //  ntu->Draw("PUit_NumInteractions >> hmc","","goff");
   ntu->Draw("PUit_TrueNumInteractions >> hmc","","goff");


  /*
  //*** data file 2012
  //TFile* fda = TFile::Open("True_PilePU_PostRerecoICHEP2012.root");
  //  TFile* fda = TFile::Open("True_PilePU_Prompt2012.root");
  TFile* fda = TFile::Open("Observed_PilePU_Prompt2012.root");
  */

  //*** data file 2011
  TFile* fda = TFile::Open("PUweights_DYJetsToLL_Summer12_ABC_TrueNumInteractions.root");

  TH1F* pileup = (TH1F*)fda->Get("pileup");

  TH1F* hdata = new TH1F("hdata","hdata",100,0,100);
  //TH1F* hdata = new TH1F("hdata","hdata",100,0,100);
  for(int ibin = 1; ibin < 101; ibin++){
    hdata->SetBinContent(ibin, pileup->GetBinContent(ibin));
  }
 
  //*** compute weights
  TH1F* hweights = (TH1F*)hdata->Clone("hweights");
  hweights->Reset();
  hweights->Divide(hdata, hmc, 1./hdata->GetSumOfWeights(), 1./hmc->GetSumOfWeights());
  //hweights->Divide(hdata, hmc, 1., 1.);
  TFile* fout = new TFile("PUweights_DYJetsToLL_Summer12_53X_ShSkim_ABC_TrueNumInteractions.root","recreate");
  hweights->Write("hweights");
  fout->Close();

}
