#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include <iostream>
#include <cstring>
#include <cctype>
using namespace std;

// #include "AliHFInvMassFitter.h"
// #include "AliHFMassFitter.h"
#endif


//modify the invariant mass histogram for better visualization; accepts a TH1F object and returns a TH1F object
TH1F* PlotHistos(TH1F* mc){
  mc->SetMarkerColor(kBlack);
  mc->SetMarkerSize(0.8);
  mc->SetMarkerStyle(8);
  mc->SetStats(0);
  mc->SetTitleSize(0.05);
  mc->SetYTitle("Counts");
  mc->GetYaxis()->CenterTitle();
  mc->GetYaxis()->SetTitleOffset(1.2);
  mc->GetYaxis()->SetLabelSize(0.05);
  mc->GetYaxis()->SetTitleSize(0.06);
  mc->GetXaxis()->SetTitleSize(0.05);
  //mc->GetXaxis()->SetRangeUser(1.105, 1.13);
  mc->GetXaxis()->SetTitleOffset(1.2);
  mc->GetXaxis()->SetLabelSize(0.05);
  mc->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  mc->GetXaxis()->CenterTitle();
  
//  mc->GetYaxis()->SetRangeUser(0, 3000);
  return mc;
}
//modify signal function for better visualization; accepts a TF1 object and returns a TF1 object
TF1* PlotSignal(TF1* sig_func)
  {
  sig_func->SetNpx(1000);
  sig_func->SetLineColor(kGreen);
  sig_func->SetLineStyle(2);
  sig_func->SetLineWidth(4);
  return sig_func;
  }
//modify background function for better visualization; accepts a TF1 object and returns a TF1 object
TF1* PlotBackground(TF1* bac_func)
  {  
  bac_func->SetNpx(100);
  //bac_func->SetLineStyle(3);
  bac_func->SetLineColor(kBlue);
  bac_func->SetLineWidth(4);
  return bac_func;
  }
 
//modify total mass fit function for better visualization; accepts a TF1 object and returns a TF1 object
TF1* PlotMass(TF1* M_func)
  { 
  M_func->SetNpx(10000);
  //M_func->SetLineStyle(3);
  M_func->SetLineWidth(4);
  M_func->SetLineColor(kRed);
  return M_func;
  }

//creates a canvas
TCanvas* CreateCanvas()
  {
    auto* canvas = new TCanvas(" canvas ","", 500,500);
    canvas->SetBottomMargin(0.15);
    canvas->SetLeftMargin (0.15);
    return canvas;
  }
  
TPad* CreatePad1()
  { TPad* pad1 = new TPad (" pad1 "," pad1 " ,0 ,0.3 ,1 ,1);
    pad1 -> SetBottomMargin (0);
    pad1 -> SetLeftMargin (0.15);
    //pad1 . SetLogy();
    pad1 -> Draw ();
    return pad1;
  }

 TPad* CreatePad2()
 {
   TPad* pad2 = new TPad (" pad2 "," pad2 " ,0 ,0.05 ,1 ,0.3);
   pad2 -> SetGrid();
   pad2 -> SetBottomMargin (0.5);
   pad2 -> SetLeftMargin (0.15);
   //pad1 . SetLogy();
   pad2 -> Draw ();
   return pad2;
 }

//plots legend; accepts a TH1F and 3 TF1 objects, returns a TLegend object
TLegend* legend_plot(TH1F* h1, TF1* M_func, TF1* sig_func, TF1* bac_func)
{
  TLegend* legend = new TLegend(0.5,0.4,0.8,0.7);
  legend->AddEntry(h1,"#Lambda hyperon","pe,X0");
  legend->AddEntry(M_func,"A[#frac{(1-B)}{#sqrt{2#pi}/#sigma_{1}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{(#sigma_{1})^{2}}}+#frac{(B)}{#sqrt{2#pi}/#sigma_{2}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{#sigma_{2}^{2}}}]+pol2","l");

  legend->AddEntry(sig_func,"A[#frac{(1-B)}{#sqrt{2#pi}/#sigma_{1}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{(#sigma_{1})^{2}}}+#frac{(B)}{#sqrt{2#pi}/#sigma_{2}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{#sigma_{2}^{2}}}]","l");
  //legend->AddEntry(bac_func,"C+Dx+Ex^{2}","l");
  //legend->AddEntry(bac_func,"C","l");
  legend->AddEntry(bac_func,"#frac{C}{m_{f}-m_{i}}+D(x-0.5(m_{f}-m_{i}))","l");
  legend -> SetLineWidth (0);
  legend->SetTextSize(0.023);
  legend->Draw();
  return legend;
  }
//creates a TH2F object 
TH2F* hist2d_pt_y()
{
  TH2F* h4 = new TH2F("", "", 15,0,3,15,0,3);
  //h4->GetZaxis()->SetRangeUser (0.7 ,1.2);
  h4->GetYaxis()->SetRangeUser (0 ,3);
  h4->GetXaxis()->SetRangeUser (0 ,3);
  h4->GetXaxis()->SetTitle("y_{Lab}");
  h4->GetYaxis()->SetTitle("p_{T} (GeV/#it{c}");
  h4->SetStats(0);
  //h4->Draw("colz");
  //gStyle->SetPaintTextFormat("4.4f");//4.4 for yield
  //h4->Draw("TEXT SAME");
  return h4;
}

//Plots a histogram which is filled with data-fit / delta data
TH1F* dif_hist(TH1F* h3) 
{
    h3 -> SetLineWidth(3);
    h3 -> SetStats (0);
    h3 -> GetXaxis() -> SetTitle("Mass (GeV/c^{2})");
    h3 -> SetTitle ("");
    h3 -> GetXaxis () -> SetLabelSize (0.15);
    h3 -> GetXaxis() -> CenterTitle();
    h3 -> GetXaxis () -> SetTitleSize (0.15);
    h3 -> GetYaxis () -> SetLabelSize (0.15);
    h3 -> GetYaxis () -> SetTitleSize (0.15);
    h3 -> GetXaxis () -> SetTitleOffset (1.1);
    h3 -> GetYaxis () -> SetTitleOffset (0.4);
    //ratio . GetYaxis ()-> SetTitle (" Data /MC");
    //207,512 divisions
    h3 -> GetYaxis ()-> SetNdivisions (104);
    h3->SetLineColor(kBlack);
    h3->SetYTitle("d-f/#Deltad");
    //h3 -> GetXaxis () -> SetRangeUser(1.105,1.13);
    return h3;
}

//creates efficiency files
void eff_creator(vector<Double_t> &y_bin_low1, vector<Double_t> &pt_bin_low1, vector<Double_t> &yield, vector<Double_t> &yield_error, vector<Double_t> &sim_mc_yield){

  //TTree* t = new TTree("t1","tree");
  
  
  TFile *f1 = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/lambda_qa_dcm.root");
  TDirectory* dir = (TDirectory*)f1->Get("SimParticles_McLambda");
  TH2D *mc_spectra = (TH2D*)dir->Get("SimParticles_rapidity_SimParticles_pT_McLambda");
  
  TFile *f3 = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/lambda_qa_urqmd.root");
  TDirectory* dir1 = (TDirectory*)f3->Get("SimParticles_McLambda");
  TH2D *Mc = (TH2D*)dir1->Get("SimParticles_rapidity_SimParticles_pT_McLambda");
  Mc->SetName("Mc"); Mc->SetTitle("Mc"); 

  TH2F* h4 = new TH2F("recons_urqmd", "recons_urqmd", 15,0,3,15,0,3);
  //h4->GetZaxis()->SetRangeUser (0.7 ,1.2);
  h4->GetYaxis()->SetRangeUser (0 ,3);
  h4->GetXaxis()->SetRangeUser (0 ,3);
  h4->GetXaxis()->SetTitle("y_{Lab}");
  h4->GetYaxis()->SetTitle("p_{T} (GeV/#it{c}");
  h4->SetStats(0);
  
  TH2F* h6 = new TH2F("Mc_urqmd", "Mc_urqmd", 15,0,3,15,0,3);
  //h4->GetZaxis()->SetRangeUser (0.7 ,1.2);
  h6->GetYaxis()->SetRangeUser (0 ,3);
  h6->GetXaxis()->SetRangeUser (0 ,3);
  h6->GetXaxis()->SetTitle("y_{Lab}");
  h6->GetYaxis()->SetTitle("p_{T} (GeV/#it{c}");
  h6->SetStats(0);
  //use clone method for h5
  TH2F* h5 = new TH2F("urqmd_Efficiency", "urqmd_Efficiency", 15,0,3,15,0,3);
  TH2F* recons = new TH2F("recons", "recons", 15,0,3,15,0,3);
    
  for (int i=0;i<225;i++)
  {
    Double_t y= y_bin_low1[i];
    Double_t pT=pt_bin_low1[i];
    Double_t y_bin = int((y+0.1)/0.2 + 1);
    Double_t pT_bin = int((pT+0.1)/0.2 + 1);
    h4->SetBinContent(y_bin, pT_bin, yield[i]);
    h4->SetBinError(y_bin, pT_bin,yield_error[i]);
    h5->SetBinContent(y_bin, pT_bin, yield[i]);
    h5->SetBinError(y_bin, pT_bin,yield_error[i]);
    auto a = mc_spectra->GetBinContent(y_bin, pT_bin);
    auto da = mc_spectra->GetBinError(y_bin, pT_bin);
    h6->SetBinContent(y_bin, pT_bin, a);
    h6->SetBinError(y_bin, pT_bin,da);
    if (sim_mc_yield[i]>0){recons->SetBinContent(y_bin,pT_bin,sim_mc_yield[i]);}
    
  }
  /*h4->Rebin(3);
  h5->Rebin(3);
  h6->Rebin(3);
  Mc->Rebin(3);
  recons->Rebin(3);*/
  
  //URQMD file creator
  auto ratio = h5->Divide(mc_spectra);
  TFile* f = new TFile("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/new_urqmd_efficiency_pt_y_yield_bdt_cut_0.8.root","recreate");
  f->cd();
  h4->Write();h5->Write();h6->Write();
  f->Write();
  f->Close();

  //DCM file creator
  TFile* ff = new TFile("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/new_dcm_100_efficiency_pt_y_yield_bdt_cut_0.8.root","recreate");
  ff->cd();
  Mc->Write(); recons->Write();
  ff->Write();
  ff->Close();
}
//A double gaussian 
   Double_t double_gaus(Double_t *x, Double_t *par)
{
   Double_t f1 = par[0]*((1.-par[3])/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]));
   Double_t f2 = par[0]*(par[3]/TMath::Sqrt(2.*TMath::Pi())/par[4]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[4]/par[4]));
   Double_t f3=(f1+f2);
   return f3;
}

   Double_t gaus1_func(Double_t *x, Double_t *par)
{
   Double_t f1 = par[0]*((1.-par[3])/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]));
   return f1;
}

   Double_t gaus2_func(Double_t *x, Double_t *par)
{
   Double_t f2 = par[0]*(par[3]/TMath::Sqrt(2.*TMath::Pi())/par[4]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[4]/par[4]));
   return f2;
}

