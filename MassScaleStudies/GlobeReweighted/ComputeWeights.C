//****** simple macro to compute  weights ******
{
          std::string year = "2011";
  //        std::string year = "2012";
  std::string HIGHLOW = "HIGH";
  //    std::string HIGHLOW = "LOW";

  TFile File2012( ("results_Globe_"+year+".root").c_str(), "read");

  TH1F* h_Et_allHgg_scE_reg = (TH1F*)File2012.Get("h_Et_allggh_scE_reg");
  h_Et_allHgg_scE_reg->Add((TH1F*)File2012.Get("h_Et_allvbf_scE_reg"));
  h_Et_allHgg_scE_reg->Add((TH1F*)File2012.Get("h_Et_allwzh_scE_reg"));
  h_Et_allHgg_scE_reg->Add((TH1F*)File2012.Get("h_Et_alltth_scE_reg"));

  TH1F* h_Eta_allHgg_scE_reg = (TH1F*)File2012.Get("h_Eta_allggh_scE_reg");
  h_Eta_allHgg_scE_reg->Add((TH1F*)File2012.Get("h_Eta_allvbf_scE_reg"));
  h_Eta_allHgg_scE_reg->Add((TH1F*)File2012.Get("h_Eta_allwzh_scE_reg"));
  h_Eta_allHgg_scE_reg->Add((TH1F*)File2012.Get("h_Eta_alltth_scE_reg"));

  ///////////
  TH1F* h_Eta_allHgg_HR9 = (TH1F*)File2012.Get("h_Eta_allggh_HR9");
  h_Eta_allHgg_HR9->Add((TH1F*)File2012.Get("h_Eta_allvbf_HR9"));
  h_Eta_allHgg_HR9->Add((TH1F*)File2012.Get("h_Eta_allwzh_HR9"));
  h_Eta_allHgg_HR9->Add((TH1F*)File2012.Get("h_Eta_alltth_HR9"));

  TH1F* h_Eta_allHgg_LR9 = (TH1F*)File2012.Get("h_Eta_allggh_LR9");
  h_Eta_allHgg_LR9->Add((TH1F*)File2012.Get("h_Eta_allvbf_LR9"));
  h_Eta_allHgg_LR9->Add((TH1F*)File2012.Get("h_Eta_allwzh_LR9"));
  h_Eta_allHgg_LR9->Add((TH1F*)File2012.Get("h_Eta_alltth_LR9"));
  ///////////

  TH1F* h_R9_allHgg_scE_reg = (TH1F*)File2012.Get("h_R9_allggh_scE_reg");
  h_R9_allHgg_scE_reg->Add((TH1F*)File2012.Get("h_R9_allvbf_scE_reg"));
  h_R9_allHgg_scE_reg->Add((TH1F*)File2012.Get("h_R9_allwzh_scE_reg"));
  h_R9_allHgg_scE_reg->Add((TH1F*)File2012.Get("h_R9_alltth_scE_reg"));

  TFile File_SN( ("results_Ntuples_"+year+".root").c_str(), "read");

  TH1F* h_Et_allZee_scE_reg = (TH1F*)File_SN->Get("h_Et_allMC");
  TH1F* h_Et_allDa_scE_reg = (TH1F*)File_SN->Get("h_Et_allDA");

  TH1F* h_Eta_allZee_scE_reg = (TH1F*)File_SN->Get("h_Eta_allMC");
  TH1F* h_Eta_allDa_scE_reg = (TH1F*)File_SN->Get("h_Eta_allDA");

  //////////
  TH1F* h_HR9_Eta_allMC = (TH1F*)File_SN->Get("h_LHR9_Eta_allMC");
  TH1F* h_HR9_Eta_allDA = (TH1F*)File_SN->Get("h_LHR9_Eta_allDA");

  TH1F* h_LR9_Eta_allMC = (TH1F*)File_SN->Get("h_allR9_Eta_allMC");
  TH1F* h_LR9_Eta_allDA = (TH1F*)File_SN->Get("h_allR9_Eta_allDA");
  //////////
  TH1F* h_R9_allZee_scE_reg = (TH1F*)File_SN->Get("h_R9_allMC");
  TH1F* h_R9_allDa_scE_reg = (TH1F*)File_SN->Get("h_R9_allDA");


  ////////////////////////   Scale

  h_Et_allZee_scE_reg->Scale(1./h_Et_allZee_scE_reg->Integral());
  h_Et_allDa_scE_reg->Scale(1./h_Et_allDa_scE_reg->Integral());

  h_Eta_allZee_scE_reg->Scale(1./h_Eta_allZee_scE_reg->Integral());
  h_Eta_allDa_scE_reg->Scale(1./h_Eta_allDa_scE_reg->Integral());

  h_HR9_Eta_allMC->Scale(1./h_HR9_Eta_allMC->Integral());
  h_HR9_Eta_allDA->Scale(1./h_HR9_Eta_allDA->Integral());

  h_LR9_Eta_allMC->Scale(1./h_LR9_Eta_allMC->Integral());
  h_LR9_Eta_allDA->Scale(1./h_LR9_Eta_allDA->Integral());

  h_R9_allZee_scE_reg->Scale(1./h_R9_allZee_scE_reg->Integral());
  h_R9_allDa_scE_reg->Scale(1./h_R9_allDa_scE_reg->Integral());
  /// Hgg
  h_Et_allHgg_scE_reg->Scale(1./h_Et_allHgg_scE_reg->Integral());
  h_Eta_allHgg_scE_reg->Scale(1./h_Eta_allHgg_scE_reg->Integral());
  h_Eta_allHgg_HR9->Scale(1./h_Eta_allHgg_HR9->Integral());
  h_Eta_allHgg_LR9->Scale(1./h_Eta_allHgg_LR9->Integral());
  h_R9_allHgg_scE_reg->Scale(1./h_R9_allHgg_scE_reg->Integral());

  //*** compute weights

TH1F* hweightsMC_Et = (TH1F*)h_Et_allHgg_scE_reg->Clone("hweightsMC_Et");
hweightsMC_Et->Reset();
hweightsMC_Et->Divide(h_Et_allHgg_scE_reg, h_Et_allZee_scE_reg, 1./h_Et_allHgg_scE_reg->GetSumOfWeights(), 1./h_Et_allZee_scE_reg->GetSumOfWeights());

TH1F* hweightsDa_Et = (TH1F*)h_Et_allHgg_scE_reg->Clone("hweightsDa_Et");
hweightsDa_Et->Reset();
hweightsDa_Et->Divide(h_Et_allHgg_scE_reg, h_Et_allDa_scE_reg, 1./h_Et_allHgg_scE_reg->GetSumOfWeights(), 1./h_Et_allDa_scE_reg->GetSumOfWeights());

  std::cout << " hweightsMC_Et->GetNbinsX() = " << hweightsMC_Et->GetNbinsX() << std::endl;
  std::cout << " hweightsDa_Et->GetNbinsX() = " << hweightsDa_Et->GetNbinsX() << std::endl;
/////////////////
TH1F* hweightsMC_Eta = (TH1F*)h_Eta_allHgg_scE_reg->Clone("hweightsMC_Eta");
hweightsMC_Eta->Reset();
hweightsMC_Eta->Divide(h_Eta_allHgg_scE_reg, h_Eta_allZee_scE_reg, 1./h_Eta_allHgg_scE_reg->GetSumOfWeights(), 1./h_Eta_allZee_scE_reg->GetSumOfWeights());

TH1F* hweightsDa_Eta = (TH1F*)h_Eta_allHgg_scE_reg->Clone("hweightsDa_Eta");
hweightsDa_Eta->Reset();
hweightsDa_Eta->Divide(h_Eta_allHgg_scE_reg, h_Eta_allDa_scE_reg, 1./h_Eta_allHgg_scE_reg->GetSumOfWeights(), 1./h_Eta_allDa_scE_reg->GetSumOfWeights());

  std::cout << " hweightsMC_Eta->GetNbinsX() = " << hweightsMC_Eta->GetNbinsX() << std::endl;
  std::cout << " hweightsDA_Eta->GetNbinsX() = " << hweightsDa_Eta->GetNbinsX() << std::endl;

/////////////////  Eta highR9 fraction
TH1F* Hgg_EtaR9fraction = (TH1F*)h_Eta_allHgg_HR9->Clone("Hgg_EtaR9fraction");
Hgg_EtaR9fraction->Reset();
Hgg_EtaR9fraction->Divide(h_Eta_allHgg_HR9, h_Eta_allHgg_scE_reg, 1./h_Eta_allHgg_HR9->GetSumOfWeights(), 1./h_Eta_allHgg_scE_reg->GetSumOfWeights());

TH1F* Zee_EtaR9fraction_DA = (TH1F*)h_HR9_Eta_allDA->Clone("Zee_EtaR9fraction_DA");
Zee_EtaR9fraction_DA->Reset();
Zee_EtaR9fraction_DA->Divide(h_HR9_Eta_allDA, h_Eta_allDa_scE_reg, 1./h_HR9_Eta_allDA->GetSumOfWeights(), 1./h_Eta_allDa_scE_reg->GetSumOfWeights());

TH1F* Zee_EtaR9fraction_MC = (TH1F*)h_HR9_Eta_allMC->Clone("Zee_EtaR9fraction_MC");
Zee_EtaR9fraction_MC->Reset();
Zee_EtaR9fraction_MC->Divide(h_HR9_Eta_allMC, h_Eta_allZee_scE_reg, 1./h_HR9_Eta_allMC->GetSumOfWeights(), 1./h_Eta_allZee_scE_reg->GetSumOfWeights());
///////////////   weights
TH1F* hweightsMC_EtaR9fr = (TH1F*)Hgg_EtaR9fraction->Clone("hweightsMC_EtaR9fr");
hweightsMC_EtaR9fr->Reset();
hweightsMC_EtaR9fr->Divide(Hgg_EtaR9fraction, Zee_EtaR9fraction_MC, 1./Hgg_EtaR9fraction->GetSumOfWeights(), 1./Zee_EtaR9fraction_MC->GetSumOfWeights());

TH1F* hweightsDA_EtaR9fr = (TH1F*)Hgg_EtaR9fraction->Clone("hweightsDA_EtaR9fr");
hweightsDA_EtaR9fr->Reset();
hweightsDA_EtaR9fr->Divide(Hgg_EtaR9fraction, Zee_EtaR9fraction_DA, 1./Hgg_EtaR9fraction->GetSumOfWeights(), 1./Zee_EtaR9fraction_DA->GetSumOfWeights());

  std::cout << " hweightsMC_EtaR9fr->GetNbinsX() = " << hweightsMC_EtaR9fr->GetNbinsX() << std::endl;
  std::cout << " hweightsDA_EtaR9fr->GetNbinsX() = " << hweightsDA_EtaR9fr->GetNbinsX() << std::endl;

/////////////////  Eta lowR9 fraction
TH1F* Hgg_EtaLR9fraction = (TH1F*)h_Eta_allHgg_LR9->Clone("Hgg_EtaLR9fraction");
Hgg_EtaLR9fraction->Reset();
Hgg_EtaLR9fraction->Divide(h_Eta_allHgg_LR9, h_Eta_allHgg_scE_reg, 1./h_Eta_allHgg_LR9->GetSumOfWeights(), 1./h_Eta_allHgg_scE_reg->GetSumOfWeights());

TH1F* Zee_EtaLR9fraction_DA = (TH1F*)h_LR9_Eta_allDA->Clone("Zee_EtaLR9fraction_DA");
Zee_EtaLR9fraction_DA->Reset();
Zee_EtaLR9fraction_DA->Divide(h_LR9_Eta_allDA, h_Eta_allDa_scE_reg, 1./h_LR9_Eta_allDA->GetSumOfWeights(), 1./h_Eta_allDa_scE_reg->GetSumOfWeights());

TH1F* Zee_EtaLR9fraction_MC = (TH1F*)h_LR9_Eta_allMC->Clone("Zee_EtaLR9fraction_MC");
Zee_EtaLR9fraction_MC->Reset();
Zee_EtaLR9fraction_MC->Divide(h_LR9_Eta_allMC, h_Eta_allZee_scE_reg, 1./h_LR9_Eta_allMC->GetSumOfWeights(), 1./h_Eta_allZee_scE_reg->GetSumOfWeights());
///////////////   weights
TH1F* hweightsMC_EtaLR9fr = (TH1F*)Hgg_EtaLR9fraction->Clone("hweightsMC_EtaLR9fr");
hweightsMC_EtaLR9fr->Reset();
hweightsMC_EtaLR9fr->Divide(Hgg_EtaLR9fraction, Zee_EtaLR9fraction_MC, 
			    1./Hgg_EtaLR9fraction->GetSumOfWeights(), 1./Zee_EtaLR9fraction_MC->GetSumOfWeights());

TH1F* hweightsDA_EtaLR9fr = (TH1F*)Hgg_EtaLR9fraction->Clone("hweightsDA_EtaLR9fr");
hweightsDA_EtaLR9fr->Reset();
hweightsDA_EtaLR9fr->Divide(Hgg_EtaLR9fraction, Zee_EtaLR9fraction_DA, 
			    1./Hgg_EtaLR9fraction->GetSumOfWeights(), 1./Zee_EtaLR9fraction_DA->GetSumOfWeights());

  std::cout << " hweightsMC_EtaLR9fr->GetNbinsX() = " << hweightsMC_EtaLR9fr->GetNbinsX() << std::endl;
  std::cout << " hweightsDA_EtaLR9fr->GetNbinsX() = " << hweightsDA_EtaLR9fr->GetNbinsX() << std::endl;


/////////////////  HR9
TH1F* hweightsMC_EtaHR9 = (TH1F*)h_Eta_allHgg_HR9->Clone("hweightsMC_EtaHR9");
hweightsMC_EtaHR9->Reset();
hweightsMC_EtaHR9->Divide(h_Eta_allHgg_HR9, h_HR9_Eta_allMC, 1./h_Eta_allHgg_HR9->GetSumOfWeights(), 1./h_HR9_Eta_allMC->GetSumOfWeights());

TH1F* hweightsDA_EtaHR9 = (TH1F*)h_Eta_allHgg_HR9->Clone("hweightsDA_EtaHR9");
hweightsDA_EtaHR9->Reset();
hweightsDA_EtaHR9->Divide(h_Eta_allHgg_HR9, h_HR9_Eta_allDA, 1./h_Eta_allHgg_HR9->GetSumOfWeights(), 1./h_HR9_Eta_allDA->GetSumOfWeights());

  std::cout << " hweightsMC_EtaHR9->GetNbinsX() = " << hweightsMC_EtaHR9->GetNbinsX() << std::endl;
  std::cout << " hweightsDA_EtaHR9->GetNbinsX() = " << hweightsDA_EtaHR9->GetNbinsX() << std::endl;

/////////////////  LR9
TH1F* hweightsMC_EtaLR9 = (TH1F*)h_Eta_allHgg_LR9->Clone("hweightsMC_EtaLR9");
hweightsMC_EtaLR9->Reset();
hweightsMC_EtaLR9->Divide(h_Eta_allHgg_LR9, h_LR9_Eta_allMC, 1./h_Eta_allHgg_LR9->GetSumOfWeights(), 1./h_LR9_Eta_allMC->GetSumOfWeights());

TH1F* hweightsDA_EtaLR9 = (TH1F*)h_Eta_allHgg_LR9->Clone("hweightsDA_EtaLR9");
hweightsDA_EtaLR9->Reset();
hweightsDA_EtaLR9->Divide(h_Eta_allHgg_LR9, h_LR9_Eta_allDA, 1./h_Eta_allHgg_LR9->GetSumOfWeights(), 1./h_LR9_Eta_allDA->GetSumOfWeights());

  std::cout << " hweightsMC_EtaLR9->GetNbinsX() = " << hweightsMC_EtaLR9->GetNbinsX() << std::endl;
  std::cout << " hweightsDA_EtaLR9->GetNbinsX() = " << hweightsDA_EtaLR9->GetNbinsX() << std::endl;

/////////////////
TH1F* hweightsMC_R9 = (TH1F*)h_R9_allHgg_scE_reg->Clone("hweightsMC_R9");
hweightsMC_R9->Reset();
hweightsMC_R9->Divide(h_R9_allHgg_scE_reg, h_R9_allZee_scE_reg, 1./h_R9_allHgg_scE_reg->GetSumOfWeights(), 1./h_R9_allZee_scE_reg->GetSumOfWeights());

TH1F* hweightsDa_R9 = (TH1F*)h_R9_allHgg_scE_reg->Clone("hweightsDa_R9");
hweightsDa_R9->Reset();
hweightsDa_R9->Divide(h_R9_allHgg_scE_reg, h_R9_allDa_scE_reg, 1./h_R9_allHgg_scE_reg->GetSumOfWeights(), 1./h_R9_allDa_scE_reg->GetSumOfWeights());

  std::cout << " hweightsMC_R9->GetNbinsX() = " << hweightsMC_R9->GetNbinsX() << std::endl;
  std::cout << " hweightsDA_R9->GetNbinsX() = " << hweightsDa_R9->GetNbinsX() << std::endl;

  ///////// write
  TFile* fout = new TFile(("HggGlobe_Stectra_weights_"+year+".root").c_str(),"recreate");
  hweightsMC_Et->Write();
  hweightsDa_Et->Write();
  hweightsMC_Eta->Write();
  hweightsMC_EtaLR9->Write();
  hweightsDA_EtaLR9->Write();
  hweightsMC_EtaHR9->Write();
  hweightsDA_EtaHR9->Write();
  hweightsMC_EtaR9fr->Write();
  hweightsDA_EtaR9fr->Write();
  hweightsMC_EtaLR9fr->Write();
  hweightsDA_EtaLR9fr->Write();

  hweightsDa_Eta->Write();
  hweightsMC_R9->Write();
  hweightsDa_R9->Write();
  fout->Close();
}
