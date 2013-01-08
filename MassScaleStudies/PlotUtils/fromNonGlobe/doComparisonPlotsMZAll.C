#include "TGraphErrors.h"
#include "TPad.h"
#include "TFile.h"
#include <iostream>
#include <algorithm>

void doComparisonPlotsMZAll(){

  //  gROOT->ProcessLine(".x /Users/Arabella/Public/style.C");

  std::string plotDir = "../../NonGlobe/PLOTS_MZAll";
  std::string plotDirOut = "2011vs2012_MZ";

//    std::string plotDir = "PLOTS_MZ_Sh";
//    std::string plotDirOut = "2011vs2012_MZ_Sh";

//       std::string plotDir = "PLOTS_MZmom";
//       std::string plotDirOut = "2011vs2012_MZmom";
     
  for(int ii=0; ii<4; ++ii){
    std::string category;
    if(ii == 0)  category = "EB_HIGH_scE_reg";
    if(ii == 1)  category = "EB_LOW_scE_reg";
    if(ii == 2)  category = "notEBEB_HIGH_scE_reg";
    if(ii == 3)  category = "notEBEB_LOW_scE_reg";

    std::string file2011 = plotDir+"/results_"+category+"_2011.root";
    std::string file2012 = plotDir+"/results_"+category+"_2012.root";

    TFile f2012(file2012.c_str(),"read");
    TFile f2011(file2011.c_str(),"read");

   TGraphErrors* graph2012 = (TGraphErrors*)f2012.Get("finalGraph");
   TGraphErrors* graph2011 = (TGraphErrors*)f2011.Get("finalGraph");

   if(ii == 0 || ii == 2){
     graph2012->SetMarkerColor(kGreen+2);
     graph2011->SetMarkerColor(kRed+2);
   }
   if(ii == 1 || ii == 3){
     graph2012->SetMarkerColor(kOrange-3);
     graph2011->SetMarkerColor(kAzure+8);
   }

   if(ii == 0 || ii == 2){
     graph2012->SetMarkerStyle(20);
     graph2011->SetMarkerStyle(20);
   }
   if(ii == 1 || ii == 3){
     graph2012->SetMarkerStyle(21);
     graph2011->SetMarkerStyle(21);
   }

  TCanvas* cplot = new TCanvas("gplot", "gplot",100,100,725,500);
  cplot->cd();
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
  TH1F *hPad = (TH1F*)gPad->DrawFrame(20.,-0.005,150.,0.01);
  hPad->GetXaxis()->SetTitle("H_{T}");
  hPad->GetYaxis()->SetTitle("(Mee/MZ)_{data}-(Mee/MZ)_{mc}");
  hPad->GetYaxis()->SetTitleOffset(tYoffset);
  hPad->GetXaxis()->SetTitleOffset(tXoffset);
  hPad->GetXaxis()->SetLabelSize(labSize);
  hPad->GetXaxis()->SetTitleSize(labSize);
  hPad->GetYaxis()->SetLabelSize(labSize);
  hPad->GetYaxis()->SetTitleSize(labSize);
  TLegend *tspec = new TLegend(0.64,0.80,0.99,0.99);
  tspec->SetFillColor(0);
  tspec->SetTextFont(42);
  if(ii == 0){
  tspec->AddEntry(graph2012,"2012 EB hR9","PL");
  tspec->AddEntry(graph2011,"2011 EB hR9","PL");
  }
  if(ii == 1){
  tspec->AddEntry(graph2012,"2012 EB lR9","PL");
  tspec->AddEntry(graph2011,"2011 EB lR9","PL");
  }
  if(ii == 2){
  tspec->AddEntry(graph2012,"2012 !EBEB hR9","PL");
  tspec->AddEntry(graph2011,"2011 !EBEB hR9","PL");
  }
  if(ii == 3){
  tspec->AddEntry(graph2012,"2012 !EBEB lR9","PL");
  tspec->AddEntry(graph2011,"2011 !EBEB lR9","PL");
  }

  graph2012->Draw("P");
  graph2011->Draw("P,same");

  tspec->Draw("same");
  cplot->Print( (plotDirOut+"/EoP_vs_Et_"+category+".png").c_str(),".png");

  //  cplot->Print( (plotDirOut+"/"+histoName+".png").c_str(),".png");

  }

  // tutto EB in 1 plot
  std::string f11_h = plotDir+"/results_EB_HIGH_scE_reg_2011.root";
  std::string f12_h = plotDir+"/results_EB_HIGH_scE_reg_2012.root";
  std::string f11_l = plotDir+"/results_EB_LOW_scE_reg_2011.root";
  std::string f12_l = plotDir+"/results_EB_LOW_scE_reg_2012.root";

  TFile F12_h(f12_h.c_str(),"read");
  TFile F11_h(f11_h.c_str(),"read");

  TFile F12_l(f12_l.c_str(),"read");
  TFile F11_l(f11_l.c_str(),"read");

  TGraphErrors* g2012_h = (TGraphErrors*)F12_h.Get("finalGraph");
  TGraphErrors* g2011_h = (TGraphErrors*)F11_h.Get("finalGraph");
  TGraphErrors* g2012_l = (TGraphErrors*)F12_l.Get("finalGraph");
  TGraphErrors* g2011_l = (TGraphErrors*)F11_l.Get("finalGraph");


  g2012_h->SetMarkerColor(kGreen+2);
  g2012_l->SetMarkerColor(kOrange-3);
  g2011_h->SetMarkerColor(kRed+2);
  g2011_l->SetMarkerColor(kAzure+8);


    g2012_h->SetMarkerStyle(20);
    g2011_h->SetMarkerStyle(20);
    g2012_l->SetMarkerStyle(21);
    g2011_l->SetMarkerStyle(21);

  TCanvas* c = new TCanvas("g", "g",100,100,725,500);
  c->cd();
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
  TH1F *hPad = (TH1F*)gPad->DrawFrame(20,-0.005,150,0.01);
  hPad->GetXaxis()->SetTitle("H_{T}");
  hPad->GetYaxis()->SetTitle("(Mee/MZ)_{data}-(Mee/MZ)_{mc}");
  hPad->GetYaxis()->SetTitleOffset(tYoffset);
  hPad->GetXaxis()->SetTitleOffset(tXoffset);
  hPad->GetXaxis()->SetLabelSize(labSize);
  hPad->GetXaxis()->SetTitleSize(labSize);
  hPad->GetYaxis()->SetLabelSize(labSize);
  hPad->GetYaxis()->SetTitleSize(labSize);

  TLegend *tspec = new TLegend(0.64,0.80,0.99,0.99);
  tspec->SetFillColor(0);
  tspec->SetTextFont(42);

  tspec->AddEntry(g2012_h,"2012 EB hR9","PL");
  tspec->AddEntry(g2011_h,"2011 EB hR9","PL");
  tspec->AddEntry(g2012_l,"2012 EB lR9","PL");
  tspec->AddEntry(g2011_l,"2011 EB lR9","PL");

  g2012_h->Draw("p");
  g2011_h->Draw("p, same");
  g2012_l->Draw("p, same");
  g2011_l->Draw("p, same");
    
  tspec->Draw("same");
  c->Print( (plotDirOut+"/EoP_vs_Et_EB.png").c_str(),".png");

  ///////////////////////////////////////////////////////
  // tutto EE in 1 plot
  std::string f11_hE = plotDir+"/results_notEBEB_HIGH_scE_reg_2011.root";
  std::string f12_hE = plotDir+"/results_notEBEB_HIGH_scE_reg_2012.root";
  std::string f11_lE = plotDir+"/results_notEBEB_LOW_scE_reg_2011.root";
  std::string f12_lE = plotDir+"/results_notEBEB_LOW_scE_reg_2012.root";

  TFile F12_hE(f12_hE.c_str(),"read");
  TFile F11_hE(f11_hE.c_str(),"read");

  TFile F12_lE(f12_lE.c_str(),"read");
  TFile F11_lE(f11_lE.c_str(),"read");

  TGraphErrors* g2012_hE = (TGraphErrors*)F12_hE.Get("finalGraph");
  TGraphErrors* g2011_hE = (TGraphErrors*)F11_hE.Get("finalGraph");
  TGraphErrors* g2012_lE = (TGraphErrors*)F12_lE.Get("finalGraph");
  TGraphErrors* g2011_lE = (TGraphErrors*)F11_lE.Get("finalGraph");


  g2012_hE->SetMarkerColor(kGreen+2);
  g2012_lE->SetMarkerColor(kOrange-3);
  g2011_hE->SetMarkerColor(kRed+2);
  g2011_lE->SetMarkerColor(kAzure+8);

    g2012_hE->SetMarkerStyle(20);
    g2011_hE->SetMarkerStyle(20);
    g2012_lE->SetMarkerStyle(21);
    g2011_lE->SetMarkerStyle(21);

  TCanvas* cE = new TCanvas("gE", "gE",100,100,725,500);
  cE->cd();
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
  TH1F *hPad = (TH1F*)gPad->DrawFrame(20,-0.005,150,0.01);
  hPad->GetXaxis()->SetTitle("H_{T}");
  hPad->GetYaxis()->SetTitle("(Mee/MZ)_{data}-(Mee/MZ)_{mc}");
  hPad->GetYaxis()->SetTitleOffset(tYoffset);
  hPad->GetXaxis()->SetTitleOffset(tXoffset);
  hPad->GetXaxis()->SetLabelSize(labSize);
  hPad->GetXaxis()->SetTitleSize(labSize);
  hPad->GetYaxis()->SetLabelSize(labSize);
  hPad->GetYaxis()->SetTitleSize(labSize);

  TLegend *tspecE = new TLegend(0.64,0.80,0.99,0.99);
  tspecE->SetFillColor(0);
  tspecE->SetTextFont(42);

  tspecE->AddEntry(g2012_hE,"2012 !EBEB hR9","PL");
  tspecE->AddEntry(g2011_hE,"2011 !EBEB hR9","PL");
  tspecE->AddEntry(g2012_lE,"2012 !EBEB lR9","PL");
  tspecE->AddEntry(g2011_lE,"2011 !EBEBE lR9","PL");

  g2012_hE->Draw("p");
  g2011_hE->Draw("p, same");
  g2012_lE->Draw("p, same");
  g2011_lE->Draw("p, same");
    
  tspecE->Draw("same");
  cE->Print( (plotDirOut+"/EoP_vs_Et_notEBEB.png").c_str(),".png");

}

