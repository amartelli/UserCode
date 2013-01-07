#include "TGraphErrors.h"
#include "TPad.h"
#include "TFile.h"
#include <iostream>
#include <algorithm>

void doComparisonPlotsEE(){

  //  gROOT->ProcessLine(".x /Users/Arabella/Public/style.C");


  ////////////////////////////  da qui forse e' utile
   std::string plotDir1 = "PLOTS_false_EBEE";
   std::string plotDir2 = "PLOTS_false_EBEE_corrected_noSh";
   std::string plotDirOut = "2011vs2012_false_noSh_EE";
   
    std::string f11_h1 = plotDir1+"/results_EE_HIGH_scE_reg_2011.root";
    std::string f12_h1 = plotDir1+"/results_EE_HIGH_scE_reg_2012.root";
    std::string f11_l1 = plotDir1+"/results_EE_LOW_scE_reg_2011.root";
    std::string f12_l1 = plotDir1+"/results_EE_LOW_scE_reg_2012.root";

    std::string f11_h2 = plotDir2+"/results_EE_HIGH_scE_reg_2011.root";
    std::string f12_h2 = plotDir2+"/results_EE_HIGH_scE_reg_2012.root";
    std::string f11_l2 = plotDir2+"/results_EE_LOW_scE_reg_2011.root";
    std::string f12_l2 = plotDir2+"/results_EE_LOW_scE_reg_2012.root";

    TFile F12_h1(f12_h1.c_str(),"read");
    TFile F11_h1(f11_h1.c_str(),"read");

    TFile F12_l1(f12_l1.c_str(),"read");
    TFile F11_l1(f11_l1.c_str(),"read");

    TFile F12_h2(f12_h2.c_str(),"read");
    TFile F11_h2(f11_h2.c_str(),"read");

    TFile F12_l2(f12_l2.c_str(),"read");
    TFile F11_l2(f11_l2.c_str(),"read");


    TGraphErrors* g2012_h1 = (TGraphErrors*)F12_h1.Get("finalGraph");
    TGraphErrors* g2011_h1 = (TGraphErrors*)F11_h1.Get("finalGraph");
    TGraphErrors* g2012_l1 = (TGraphErrors*)F12_l1.Get("finalGraph");
    TGraphErrors* g2011_l1 = (TGraphErrors*)F11_l1.Get("finalGraph");

    TGraphErrors* g2012_h2 = (TGraphErrors*)F12_h2.Get("finalGraph");
    TGraphErrors* g2011_h2 = (TGraphErrors*)F11_h2.Get("finalGraph");
    TGraphErrors* g2012_l2 = (TGraphErrors*)F12_l2.Get("finalGraph");
    TGraphErrors* g2011_l2 = (TGraphErrors*)F11_l2.Get("finalGraph");

    g2012_h1->SetMarkerColor(kGreen+2);
    g2012_l1->SetMarkerColor(kGreen+2);
    g2011_h1->SetMarkerColor(kRed+2);
    g2011_l1->SetMarkerColor(kRed+2);

    g2012_h1->SetMarkerStyle(20);
    g2011_h1->SetMarkerStyle(20);
    g2012_l1->SetMarkerStyle(21);
    g2011_l1->SetMarkerStyle(21);

    g2012_h2->SetMarkerColor(kGreen+2);
    g2012_l2->SetMarkerColor(kGreen+2);
    g2011_h2->SetMarkerColor(kRed+2);
    g2011_l2->SetMarkerColor(kRed+2);

    g2012_h2->SetMarkerStyle(20);
    g2011_h2->SetMarkerStyle(20);
    g2012_l2->SetMarkerStyle(21);
    g2011_l2->SetMarkerStyle(21);

    TCanvas* cA = new TCanvas("gA", "gA",100,100,725,500);
    cA->cd();
    TPad *cLeft  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
    cLeft->SetLeftMargin(0.17);
    cLeft->SetRightMargin(0.025);
    cLeft->SetBottomMargin(0.17);
    cLeft->Draw();
    float tYoffset = 1.75;
    float tXoffset = 1.6;
    float labSize = 0.04;
    float labSize2 = 0.07;
    cLeft->cd();
    cLeft->SetGridx();
    cLeft->SetGridy();
    
    // pad settings                                                                     
    TH1F *hPad = (TH1F*)gPad->DrawFrame(20,-0.015,100,0.015);
    hPad->GetXaxis()->SetTitle("E_{T}");
    hPad->GetYaxis()->SetTitle("E/p_{data}-E/p_{mc}");
    hPad->GetYaxis()->SetTitleOffset(tYoffset);
    hPad->GetXaxis()->SetTitleOffset(tXoffset);
    hPad->GetXaxis()->SetLabelSize(labSize);
    hPad->GetXaxis()->SetTitleSize(labSize);
    hPad->GetYaxis()->SetLabelSize(labSize);
    hPad->GetYaxis()->SetTitleSize(labSize);

    TLegend *tspec = new TLegend(0.64,0.80,0.99,0.99);
    tspec->SetFillColor(0);
    tspec->SetTextFont(42);

    tspec->AddEntry(g2012_h1,"2012 EE hR9","PL");
    tspec->AddEntry(g2011_h1,"2011 EE hR9","PL");
    tspec->AddEntry(g2012_l1,"2012 EE lR9","PL");
    tspec->AddEntry(g2011_l1,"2011 EE lR9","PL");

    g2012_h1->Draw("p");
    g2011_h1->Draw("p, same");
    g2012_l1->Draw("p, same");
    g2011_l1->Draw("p, same");

    g2012_h2->Draw("p, same");
    g2011_h2->Draw("p, same");
    g2012_l2->Draw("p, same");
    g2011_l2->Draw("p, same");
    
    tspec->Draw("same");
    cA->Print( (plotDirOut+"/EoP_vs_Et_Hgg.png").c_str(),".png");

    /////////////////////////////////   Spettri

    TH1F* Eta_2012_h1 = (TH1F*)F12_h1.Get("h_Eta_allDA");
    TH1F* Eta_2011_h1 = (TH1F*)F11_h1.Get("h_Eta_allDA");

    TH1F* Eta_2012_h2 = (TH1F*)F12_h2.Get("h_Eta_allDA");
    TH1F* Eta_2011_h2 = (TH1F*)F11_h2.Get("h_Eta_allDA");

    Eta_2012_h1->SetLineColor(kGreen+2);
    Eta_2011_h1->SetLineColor(kRed+2);

    Eta_2012_h1->SetLineWidth(2);  
    Eta_2011_h1->SetLineWidth(2);  

    Eta_2012_h2->SetLineColor(kGreen+2);
    Eta_2011_h2->SetLineColor(kRed+2);

    Eta_2012_h2->SetLineWidth(2);  
    Eta_2011_h2->SetLineWidth(2);  

    Eta_2012_h2->SetLineStyle(2);  
    Eta_2011_h2->SetLineStyle(2);  

    Eta_2012_h1->Rebin(10);
    Eta_2011_h1->Rebin(10);
    Eta_2012_h2->Rebin(10);
    Eta_2011_h2->Rebin(10);

    TLegend *tspec1 = new TLegend(0.64,0.80,0.99,0.99);
    tspec1->SetFillColor(0);
    tspec1->SetTextFont(42);

    tspec1->AddEntry(Eta_2012_h1,"2012","L");
    tspec1->AddEntry(Eta_2011_h1,"2011","L");
    tspec1->AddEntry(Eta_2012_h2,"2012 - Hgg (Et, #eta rew.)","L");
    tspec1->AddEntry(Eta_2011_h2,"2011 - Hgg (Et, #eta rew.)","L");

    TCanvas* c2 = new TCanvas("g2A", "g2A",100,100,725,500);
    Eta_2012_h1->GetXaxis()->SetTitle("Eta");
    Eta_2012_h1->GetXaxis()->SetRangeUser(-3., 3.);
    Eta_2012_h1->GetYaxis()->SetRangeUser(0.1, 100000.);
    Eta_2012_h1->Draw("hist");
    Eta_2011_h1->Draw("hist, same");
    Eta_2012_h2->Draw("hist, same");
    Eta_2011_h2->Draw("hist, same");
    tspec1->Draw("same");
    c2->Print( (plotDirOut+"/Eta_HggEta.png").c_str(),".png");


    //    fine spettri ripesati
    TH1F* Eta_HR9_2012_1 = (TH1F*)F12_h1.Get("h_LHR9_Eta_allDA");
    TH1F* Eta_HR9_2011_1 = (TH1F*)F11_h1.Get("h_LHR9_Eta_allDA");

    TH1F* Eta_HR9_2012_2 = (TH1F*)F12_h2.Get("h_LHR9_Eta_allDA");
    TH1F* Eta_HR9_2011_2 = (TH1F*)F11_h2.Get("h_LHR9_Eta_allDA");

    Eta_HR9_2012_1->SetLineColor(kGreen+2);
    Eta_HR9_2011_1->SetLineColor(kRed+2);
    Eta_HR9_2012_1->SetLineWidth(2);
    Eta_HR9_2011_1->SetLineWidth(2);

    Eta_HR9_2012_2->SetLineColor(kGreen+2);
    Eta_HR9_2011_2->SetLineColor(kRed+2);
    Eta_HR9_2012_2->SetLineWidth(2);
    Eta_HR9_2011_2->SetLineWidth(2);

    Eta_HR9_2012_2->SetLineStyle(2);
    Eta_HR9_2011_2->SetLineStyle(2);

    Eta_HR9_2012_1->Rebin(10);
    Eta_HR9_2011_1->Rebin(10);
    Eta_HR9_2012_2->Rebin(10);
    Eta_HR9_2011_2->Rebin(10);

    ////////////////////////////
    //frac 2012

    std::cout << " Eta_HR9_2012_1->GetNbinsX() = " << Eta_HR9_2012_1->GetNbinsX() << std::endl;
    std::cout << " Eta_2012_h1->GetNbinsX() = " << Eta_2012_h1->GetNbinsX() << std::endl;

    TH1F* LH_fraction_2012 = (TH1F*)Eta_HR9_2012_1->Clone("LH_fraction_2012");
    LH_fraction_2012->Reset();
    LH_fraction_2012->Divide(Eta_HR9_2012_1, Eta_2012_h1, 1./Eta_HR9_2012_1->GetSumOfWeights(), 1./Eta_2012_h1->GetSumOfWeights());
    //LH_fraction_2012->Divide(Eta_2012_h1, LH_2012, 1., 1.);

    //reweighted
    TH1F* LH_fraction_2012_rew = (TH1F*)Eta_2012_h2->Clone("LH_fraction_2012_rew");
    LH_fraction_2012_rew->Reset();
    LH_fraction_2012_rew->Divide(Eta_HR9_2012_2, Eta_2012_h2, 1./Eta_HR9_2012_2->GetSumOfWeights(), 1./Eta_2012_h2->GetSumOfWeights());
    //LH_fraction_2012_rew->Divide(Eta_2012_h2, LH_2012_rew, 1., 1.);

    TLegend *tspec3 = new TLegend(0.64,0.80,0.99,0.99);
    tspec3->SetFillColor(0);
    tspec3->SetTextFont(42);

    tspec3->AddEntry(LH_fraction_2012,"2012 fraction r9>0.94","L");
    tspec3->AddEntry(LH_fraction_2012_rew,"2012 fraction r9>0.94 - Hgg (Et, #eta rew.)","L");

    TCanvas* c_LH_fraction_2012 = new TCanvas();
    //    LH_fraction_2012->GetXaxis()->SetRangeUser(0., 3.);
     LH_fraction_2012->Draw("hist");
     LH_fraction_2012_rew->Draw("hist, same");
     tspec3->Draw("same");
//     Eta_HR9_2012_1->Draw("hist");
//     Eta_HR9_2012_2->Draw("hist, same");
//     Eta_2012_h1->Draw("hist, same");
//     Eta_2012_h2->Draw("hist, same");
    c_LH_fraction_2012->Print( (plotDirOut+"/LH_fraction_2012.png").c_str(),".png");


    /////////////////////////////////
    TH1F* LH_fraction_2011 = (TH1F*)Eta_HR9_2011_1->Clone("LH_fraction_2011");
    LH_fraction_2011->Reset();
    LH_fraction_2011->Divide(Eta_HR9_2011_1, Eta_2011_h1, 1./Eta_HR9_2011_1->GetSumOfWeights(), 1./Eta_2011_h1->GetSumOfWeights());
    //LH_fraction_2011->Divide(Eta_2011_h1, LH_2011, 1., 1.);

    //reweighted
    TH1F* LH_fraction_2011_rew = (TH1F*)Eta_HR9_2011_2->Clone("LH_fraction_2011_rew");
    LH_fraction_2011_rew->Reset();
    LH_fraction_2011_rew->Divide(Eta_HR9_2011_2, Eta_2011_h2, 1./Eta_HR9_2011_2->GetSumOfWeights(), 1./Eta_2011_h2->GetSumOfWeights());
    //LH_fraction_2011_rew->Divide(Eta_2011_h2, LH_2011_rew, 1., 1.);

    TLegend *tspec4 = new TLegend(0.64,0.80,0.99,0.99);
    tspec4->SetFillColor(0);
    tspec4->SetTextFont(42);
    tspec4->AddEntry(LH_fraction_2011,"2011 fraction r9>0.94","L");
    tspec4->AddEntry(LH_fraction_2011_rew,"2011 fraction r9>0.94 - Hgg (Et, #eta rew.)","L");

    TCanvas* c_LH_fraction_2011 = new TCanvas();
    //    LH_fraction_2011->GetXaxis()->SetRangeUser(0., 3.);
    LH_fraction_2011->Draw("hist");
    LH_fraction_2011_rew->Draw("hist, same");
    tspec4->Draw("same");
    c_LH_fraction_2011->Print( (plotDirOut+"/LH_fraction_2011.png").c_str(),".png");


    // fino a qui 2

    //////////////from Hgg

    TFile File2012( "fromHgg/results_Globe_2012.root", "read");
    TFile File2011( "fromHgg/results_Globe_2011.root", "read");

    TH1F* h_Eta_allHgg_HR92012 = (TH1F*)File2012.Get("h_Eta_allggh_HR9");
    h_Eta_allHgg_HR92012->Add((TH1F*)File2012.Get("h_Eta_allvbf_HR9"));
    h_Eta_allHgg_HR92012->Add((TH1F*)File2012.Get("h_Eta_allwzh_HR9"));
    h_Eta_allHgg_HR92012->Add((TH1F*)File2012.Get("h_Eta_alltth_HR9"));

    TH1F* h_Eta_allHgg_HR92011 = (TH1F*)File2011.Get("h_Eta_allggh_HR9");
    h_Eta_allHgg_HR92011->Add((TH1F*)File2011.Get("h_Eta_allvbf_HR9"));
    h_Eta_allHgg_HR92011->Add((TH1F*)File2011.Get("h_Eta_allwzh_HR9"));
    h_Eta_allHgg_HR92011->Add((TH1F*)File2011.Get("h_Eta_alltth_HR9"));

    TH1F* h_Eta_allHgg_LH2012 = (TH1F*)File2012.Get("h_Eta_allggh_scE_reg");
    h_Eta_allHgg_LH2012->Add((TH1F*)File2012.Get("h_Eta_allvbf_scE_reg"));
    h_Eta_allHgg_LH2012->Add((TH1F*)File2012.Get("h_Eta_allwzh_scE_reg"));
    h_Eta_allHgg_LH2012->Add((TH1F*)File2012.Get("h_Eta_alltth_scE_reg"));

    TH1F* h_Eta_allHgg_LH2011 = (TH1F*)File2011.Get("h_Eta_allggh_scE_reg");
    h_Eta_allHgg_LH2011->Add((TH1F*)File2011.Get("h_Eta_allvbf_scE_reg"));
    h_Eta_allHgg_LH2011->Add((TH1F*)File2011.Get("h_Eta_allwzh_scE_reg"));
    h_Eta_allHgg_LH2011->Add((TH1F*)File2011.Get("h_Eta_alltth_scE_reg"));

    h_Eta_allHgg_HR92012->Rebin(10);
    h_Eta_allHgg_HR92011->Rebin(10);
    h_Eta_allHgg_LH2012->Rebin(10);
    h_Eta_allHgg_LH2011->Rebin(10);

    h_Eta_allHgg_HR92012->SetLineWidth(2);
    h_Eta_allHgg_HR92011->SetLineWidth(2);
    h_Eta_allHgg_LH2012->SetLineWidth(2);
    h_Eta_allHgg_LH2011->SetLineWidth(2);

    //LH from Hgg
    TH1F* LH_fraction_2011_Hgg = (TH1F*)h_Eta_allHgg_HR92011->Clone("LH_fraction_2011_Hgg");
    LH_fraction_2011_Hgg->Reset();
    LH_fraction_2011_Hgg->Divide(h_Eta_allHgg_HR92011, h_Eta_allHgg_LH2011, 
				 1./h_Eta_allHgg_HR92011->GetSumOfWeights(), 1./h_Eta_allHgg_LH2011->GetSumOfWeights());
    //LH_fraction_2011_Hgg->Divide(h_Eta_allHgg_H2011, LH_2011_Hgg, 1., 1.);

    TH1F* LH_fraction_2012_Hgg = (TH1F*)h_Eta_allHgg_HR92012->Clone("LH_fraction_2012_Hgg");
    LH_fraction_2012_Hgg->Reset();
    LH_fraction_2012_Hgg->Divide(h_Eta_allHgg_HR92012, h_Eta_allHgg_LH2012, 
				 1./h_Eta_allHgg_HR92012->GetSumOfWeights(), 1./h_Eta_allHgg_LH2012->GetSumOfWeights());
    //LH_fraction_2012_Hgg->Divide(h_Eta_allHgg_H2012, LH_2012_Hgg, 1., 1.);

    TLegend *tspec5 = new TLegend(0.64,0.80,0.99,0.99);
    tspec5->SetFillColor(0);
    tspec5->SetTextFont(42);
    tspec5->AddEntry(LH_fraction_2011,"2011 fraction r9>0.94","L");
    tspec5->AddEntry(LH_fraction_2011_Hgg,"2011 fraction r9>0.94 - HggGlobe","L");

    TCanvas* c_LH_fraction_2011_Hgg = new TCanvas();
    //    LH_fraction_2011_Hgg->GetXaxis()->SetRangeUser(0., 3.);
    LH_fraction_2011_Hgg->Draw("hist");
    LH_fraction_2011->Draw("hist, same");
    tspec5->Draw("same");
    c_LH_fraction_2011_Hgg->Print( (plotDirOut+"/LH_fraction_2011_Hgg.png").c_str(),".png");

    TLegend *tspec6 = new TLegend(0.64,0.80,0.99,0.99);
    tspec6->SetFillColor(0);
    tspec6->SetTextFont(42);
    tspec6->AddEntry(LH_fraction_2012,"2012 fraction r9>0.94","L");
    tspec6->AddEntry(LH_fraction_2012_Hgg,"2012 fraction r9>0.94 - HggGlobe","L");

    TCanvas* c_LH_fraction_2012_Hgg = new TCanvas();
    //    LH_fraction_2012_Hgg->GetXaxis()->SetRangeUser(0., 3.);
    LH_fraction_2012_Hgg->Draw("hist");
    LH_fraction_2012->Draw("hist, same");
    tspec6->Draw("same");
    c_LH_fraction_2012_Hgg->Print( (plotDirOut+"/LH_fraction_2012_Hgg.png").c_str(),".png");


    ///////////////////// normalized
//     LH_fraction_2011_Hgg->Scale(1./LH_fraction_2011_Hgg->Integral());
//     LH_fraction_2011->Scale(1./LH_fraction_2011->Integral());
//     LH_fraction_2011_rew->Scale(1./LH_fraction_2011_rew->Integral());

    //LH_fraction_2011_Hgg->Scale(1.*LH_fraction_2011->Integral()/LH_fraction_2011_Hgg->Integral()); 
    // //  LH_fraction_2011->Scale(1./LH_fraction_2011->Integral());         
    //LH_fraction_2011_rew->Scale(1.*LH_fraction_2011->Integral()/LH_fraction_2011_rew->Integral()); 

    //     LH_fraction_2011_Hgg->Scale(1./LH_fraction_2011_Hgg->GetEntries());
    //     LH_fraction_2011->Scale(1./LH_fraction_2011->GetEntries());
    //     LH_fraction_2011_rew->Scale(1./LH_fraction_2011_rew->GetEntries());


    TLegend *tspec7 = new TLegend(0.64,0.80,0.99,0.99);
    tspec7->SetFillColor(0);
    tspec7->SetTextFont(42);
    tspec7->AddEntry(LH_fraction_2011,"2011 fraction r9>0.94","L");
    tspec7->AddEntry(LH_fraction_2011_rew,"2011 fraction r9>0.94 - Hgg (Et, #eta rew.)","L");
    tspec7->AddEntry(LH_fraction_2011_Hgg,"2011 fraction r9>0.94 - HggGlobe","L");
    TCanvas* c_LH_fraction_2011_all = new TCanvas();
    //    LH_fraction_2011_Hgg->GetXaxis()->SetRangeUser(0., 3.);
    LH_fraction_2011->GetXaxis()->SetTitle("#eta");
    LH_fraction_2011->Draw("hist");
    LH_fraction_2011_Hgg->Draw("hist, same");
    LH_fraction_2011_rew->Draw("hist, same");
    tspec7->Draw("same");
    c_LH_fraction_2011_all->Print( (plotDirOut+"/LH_fraction_2011_all.png").c_str(),".png");



//     LH_fraction_2012_Hgg->Scale(1./LH_fraction_2012_Hgg->Integral());
//     LH_fraction_2012->Scale(1./LH_fraction_2012->Integral());
//     LH_fraction_2012_rew->Scale(1./LH_fraction_2012_rew->Integral());

    //    LH_fraction_2012_Hgg->Scale(1.*LH_fraction_2012->Integral()/LH_fraction_2012_Hgg->Integral());
    ////    LH_fraction_2012->Scale(1./LH_fraction_2012->Integral());        
    //    LH_fraction_2012_rew->Scale(1.*LH_fraction_2012->Integral()/LH_fraction_2012_rew->Integral());

    //     LH_fraction_2012_Hgg->Scale(1./LH_fraction_2012_Hgg->GetEntries());
    //     LH_fraction_2012->Scale(1./LH_fraction_2012->GetEntries());
    //     LH_fraction_2012_rew->Scale(1./LH_fraction_2012_rew->GetEntries());

    TLegend *tspec8 = new TLegend(0.64,0.80,0.99,0.99);
    tspec8->SetFillColor(0);
    tspec8->SetTextFont(42);
    tspec8->AddEntry(LH_fraction_2012,"2012 fraction r9>0.94","L");
    tspec8->AddEntry(LH_fraction_2012_rew,"2012 fraction r9>0.94 - Hgg (Et, #eta rew.)","L");
    tspec8->AddEntry(LH_fraction_2012_Hgg,"2012 fraction r9>0.94 - HggGlobe","L");
    TCanvas* c_LH_fraction_2012_all = new TCanvas();
    //    LH_fraction_2012_Hgg->GetXaxis()->SetRangeUser(0., 3.);
    LH_fraction_2012->GetXaxis()->SetTitle("#eta");
    LH_fraction_2012->Draw("hist");
    LH_fraction_2012_Hgg->Draw("hist, same");
    LH_fraction_2012_rew->Draw("hist, same");
    tspec8->Draw("same");
    c_LH_fraction_2012_all->Print( (plotDirOut+"/LH_fraction_2012_all.png").c_str(),".png");

}

