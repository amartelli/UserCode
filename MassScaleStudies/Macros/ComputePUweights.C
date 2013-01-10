//****** simple macro to compute PU weights ******


void ComputePUweights()
{
  // histograms
  TH1F* h_mc = new TH1F("h_mc","",100,-0.5,99.5);
  TH1F* h_da = new TH1F("h_da","",100,-0.5,99.5);
  TH1F* temp = NULL;
  
  // files
  TFile* f_mc = TFile::Open("../Pileup/pileup_DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM.root");
  TFile* f_da = TFile::Open("../Pileup/pileup_69p3mb_true_Moriond2013.root");
  
  
  
  // fill mc histogram
  std::cout << "\n\n\n>>> MC pileup histogram" << std::endl;
  
  temp = (TH1F*)( f_mc->Get("PUDumper/nPUtrue") );
  for(int bin = 1; bin <= temp->GetNbinsX(); ++bin)
  {
    float binCenter  = temp -> GetBinCenter(bin);
    float binContent = temp -> GetBinContent(bin);
    std::cout << "bin: " << bin << "   binCenter: " << binCenter << "   binContent: " << binContent << std::endl;
    
    h_mc -> Fill(binCenter,binContent);
  }
  h_mc -> Scale(1./h_mc->Integral());
  std::cout << "Integral: " << h_mc -> Integral() << std::endl;
  
  
  // fill da histogram
  std::cout << "\n\n\n>>> DA pileup histogram" << std::endl;
  
  temp = (TH1F*)( f_da->Get("pileup") );
  for(int bin = 1; bin <= temp->GetNbinsX(); ++bin)
  {
    float binCenter  = temp -> GetBinCenter(bin) - 0.5;
    float binContent = temp -> GetBinContent(bin);
    std::cout << "bin: " << bin << "   binCenter: " << binCenter << "   binContent: " << binContent << std::endl;
    
    h_da -> Fill(binCenter,binContent);
  }
  h_da -> Scale(1./h_da->Integral());
  std::cout << "Integral: " << h_da -> Integral() << std::endl;
  
  
  // save PU weights
  TFile* f = new TFile("../Pileup/pileup_69p3mb_true_Moriond2013__Summer12_DR53X-PU_S10_START53.root","RECREATE");
  f -> cd();
  
  TH1F* h_PUweights = (TH1F*)(h_da->Clone("h_PUweights"));
  h_PUweights -> Divide(h_mc);
  
  h_da -> Write();
  h_mc -> Write();
  h_PUweights -> Write();
  
  f -> Close();
}
