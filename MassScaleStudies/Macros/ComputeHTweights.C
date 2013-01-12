//****** simple macro to compute PU weights ******


void ComputeHTweights()
{
  // histograms
  TH1F* h_mc = new TH1F("h_mc","",1000,0.,1000.);
  TH1F* h_da = new TH1F("h_da","",1000,0.,1000.);
  TH1F* temp = NULL;
  
  // files
  //   TFile* f_mc = TFile::Open("../NonGlobe/PLOTS_MZ/results_EBEB_HH_scE_reg_2012.root");
  //   TFile* f_da = TFile::Open("../NonGlobe/PLOTS_MZ/results_EBEB_HH_scE_reg_2012.root");

//   TFile* f_mc = TFile::Open("../NonGlobe/PLOTS_MZ/results_EBEB_notHH_scE_reg_2012.root");
//   TFile* f_da = TFile::Open("../NonGlobe/PLOTS_MZ/results_EBEB_notHH_scE_reg_2012.root");

//    TFile* f_mc = TFile::Open("../NonGlobe/PLOTS_MZ/results_notEBEB_HH_scE_reg_2012.root");
//    TFile* f_da = TFile::Open("../NonGlobe/PLOTS_MZ/results_notEBEB_HH_scE_reg_2012.root");

   TFile* f_mc = TFile::Open("../NonGlobe/PLOTS_MZ/results_notEBEB_notHH_scE_reg_2012.root");
   TFile* f_da = TFile::Open("../NonGlobe/PLOTS_MZ/results_notEBEB_notHH_scE_reg_2012.root");
  
  
  // fill mc histogram
  std::cout << "\n\n\n>>> MC pileup histogram" << std::endl;
  
  temp = (TH1F*)( f_mc->Get("h_Ht_allMC") );
  temp->Rebin(5);
  for(int bin = 1; bin <= temp->GetNbinsX(); ++bin)
  {
    float binCenter  = temp -> GetBinCenter(bin) ;
    float binContent = temp -> GetBinContent(bin);
    std::cout << "bin: " << bin << "   binCenter: " << binCenter << "   binContent: " << binContent << std::endl;
    
    h_mc -> Fill(binCenter,binContent);
  }
  h_mc -> Scale(1./h_mc->Integral());
  std::cout << "Integral: " << h_mc -> Integral() << std::endl;
  
  
  // fill da histogram
  std::cout << "\n\n\n>>> DA pileup histogram" << std::endl;
  
  temp = (TH1F*)( f_da->Get("h_Ht_allMC") );
  temp->Rebin(5);
  for(int bin = 1; bin <= temp->GetNbinsX(); ++bin)
  {
    float binCenter  = temp -> GetBinCenter(bin) ;
    float binContent = temp -> GetBinContent(bin);
    std::cout << "bin: " << bin << "   binCenter: " << binCenter << "   binContent: " << binContent << std::endl;
    
    h_da -> Fill(binCenter,binContent);
  }
  h_da -> Scale(1./h_da->Integral());
  std::cout << "Integral: " << h_da -> Integral() << std::endl;
  
  
  // save PU weights
  //  TFile* f = new TFile("../HT/EBEB_HH_2012.root","RECREATE");
  //  TFile* f = new TFile("../HT/EBEB_notHH_2012.root","RECREATE");
  //   TFile* f = new TFile("../HT/notEBEB_HH_2012.root","RECREATE");
   TFile* f = new TFile("../HT/notEBEB_notHH_2012.root","RECREATE");
  f -> cd();
  
  TH1F* h_PUweights = (TH1F*)(h_da->Clone("h_HTweights"));
  h_PUweights -> Divide(h_mc);
  
  h_da -> Write();
  h_mc -> Write();
  h_PUweights -> Write();
  
  f -> Close();
}
