#include "TGraphErrors.h"
#include "TPad.h"
#include "TFile.h"
#include <iostream>
#include <algorithm>

void doComparisonPlotsMZAll(){

  gROOT->ProcessLine(".x /Users/Arabella/Public/style.C");

//     std::string plotDir = "PLOTS_MZAll";
//     std::string plotDirOut = "2011vs2012_MZAll";
//   float xmin = -0.005;
//   float xmax = 0.015;

     std::string plotDir = "PLOTS_MZAll_Sh";
     std::string plotDirOut = "2011vs2012_MZAll_Sh";
  float xmin = -0.005;
  float xmax = 0.015;

//     std::string plotDir = "PLOTS_MZP";
//     std::string plotDirOut = "2011vs2012_MZP";

//         std::string plotDir = "PLOTS_MZP_Sh";
//         std::string plotDirOut = "2011vs2012_MZP_Sh";
     
  for(int ii=0; ii<4; ++ii){
    if(ii == 0)  std::string category = "EB_HIGH_scE_reg";
    if(ii == 1)  std::string category = "EB_LOW_scE_reg";
    if(ii == 2)  std::string category = "notEBEB_HIGH_scE_reg";
    if(ii == 3)  std::string category = "notEBEB_LOW_scE_reg";

    std::string file1112 = plotDir+"/results_"+category+"_1112.root";

    TFile f1112(file1112.c_str(),"read");

   TGraphErrors* graph1112 = (TGraphErrors*)f1112.Get("finalGraph");

   if(ii == 0 || ii == 2)   graph1112->SetMarkerColor(kGreen+2);
   if(ii == 1 || ii == 3)   graph1112->SetMarkerColor(kRed+2);


   if(ii == 0 || ii == 2){
     graph1112->SetMarkerStyle(20);
   }
   if(ii == 1 || ii == 3){
     graph1112->SetMarkerStyle(21);
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
  TH1F *hPad = (TH1F*)gPad->DrawFrame(20.,xmin,150.,xmax);
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
  tspec->AddEntry(graph1112,"2011+2012 EB hR9","PL");
  }
  if(ii == 1){
  tspec->AddEntry(graph1112,"2011+2012 EB lR9","PL");
  }
  if(ii == 2){
  tspec->AddEntry(graph1112,"2011+2012 !EBEB hR9","PL");
  }
  if(ii == 3){
  tspec->AddEntry(graph1112,"2011+2012 !EBEB lR9","PL");
  }

  graph1112->Draw("P");
  //  graph2011->Draw("P,same");

  tspec->Draw("same");
  cplot->Print( (plotDirOut+"/EoP_vs_Et_"+category+".png").c_str(),".png");

  //  cplot->Print( (plotDirOut+"/"+histoName+".png").c_str(),".png");

  }

  // tutto EB in 1 plot
  std::string f1112_h = plotDir+"/results_EB_HIGH_scE_reg_1112.root";
  std::string f1112_l = plotDir+"/results_EB_LOW_scE_reg_1112.root";

  TFile F1112_h(f1112_h.c_str(),"read");
  TFile F1112_l(f1112_l.c_str(),"read");

  TGraphErrors* g1112_h = (TGraphErrors*)F1112_h.Get("finalGraph");
  TGraphErrors* g1112_l = (TGraphErrors*)F1112_l.Get("finalGraph");

    g1112_h->SetMarkerColor(kGreen+2);
    g1112_l->SetMarkerColor(kRed+2);

    g1112_h->SetMarkerStyle(20);
    g1112_l->SetMarkerStyle(21);

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
  TH1F *hPad = (TH1F*)gPad->DrawFrame(20,xmin,150,xmax);
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

  tspec->AddEntry(g1112_h,"2011+2012 EB hR9","PL");
  tspec->AddEntry(g1112_l,"2011+2012 EB lR9","PL");

  g1112_h->Draw("p");
  g1112_l->Draw("p, same");
    
  tspec->Draw("same");
  c->Print( (plotDirOut+"/EoP_vs_Et_EB.png").c_str(),".png");

  ///////////////////////////////////////////////////////
  // tutto EE in 1 plot
  std::string f1112E_h = plotDir+"/results_notEBEB_HIGH_scE_reg_1112.root";
  std::string f1112E_l = plotDir+"/results_notEBEB_LOW_scE_reg_1112.root";

  TFile F1112E_h(f1112E_h.c_str(),"read");
  TFile F1112E_l(f1112E_l.c_str(),"read");

  TGraphErrors* g1112E_h = (TGraphErrors*)F1112E_h.Get("finalGraph");
  TGraphErrors* g1112E_l = (TGraphErrors*)F1112E_l.Get("finalGraph");

    g1112E_h->SetMarkerColor(kGreen+2);
    g1112E_l->SetMarkerColor(kRed+2);

    g1112E_h->SetMarkerStyle(20);
    g1112E_l->SetMarkerStyle(21);

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
  TH1F *hPad = (TH1F*)gPad->DrawFrame(20,xmin,150,xmax);
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

  tspec->AddEntry(g1112E_h,"2011+2012 !EBEB hR9","PL");
  tspec->AddEntry(g1112E_l,"2011+2012 !EBEB lR9","PL");

  g1112E_h->Draw("p");
  g1112E_l->Draw("p, same");
    
  tspec->Draw("same");
  c->Print( (plotDirOut+"/EoP_vs_Et_notEBEB.png").c_str(),".png");


}

