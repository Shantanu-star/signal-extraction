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
#include <TString.h>
using namespace std;
#endif


template<typename h>
auto PlotHistos(h hist)
{
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerSize(0.4);
  hist->SetMarkerStyle(8);
  hist->SetStats(0);
  hist->SetTitleSize(0.05);
  hist->SetYTitle("counts");
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetMaxDigits(3);
  hist->GetYaxis()->SetTitleOffset(0.7);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleSize(0.05);
  //hist->GetXaxis()->SetRangeUser(1.105, 1.13);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  hist->GetXaxis()->CenterTitle();
  
//  hist->GetYaxis()->SetRangeUser(0, 3000);
  return hist;
}

template<typename h>
auto dif_hist(h hist) 
{
    hist -> SetLineWidth(1);
    hist->SetMarkerColor(kBlack);
    hist->SetMarkerStyle(8);
    hist -> SetMarkerSize(0.1);
    hist -> SetStats (0);
    hist -> GetXaxis() -> SetTitle("Mass (GeV/c^{2})");
    hist -> SetTitle ("");
    hist -> GetXaxis () -> SetLabelSize (0.1);
    hist -> GetXaxis()  -> CenterTitle();
    hist -> GetXaxis () -> SetTitleSize (0.15);
    hist -> GetYaxis () -> SetLabelSize (0.1);
    hist -> GetYaxis () -> SetTitleSize (0.1);
    hist -> GetXaxis () -> SetTitleOffset (1.1);
    hist -> GetYaxis () -> SetTitleOffset (0.3);
    //ratio . GetYaxis ()-> SetTitle (" Data /MC");
    //207,512 divisions
    hist -> GetYaxis ()-> SetNdivisions (104);
    hist->SetLineColor(kBlack);
    hist->SetYTitle("(d-f)/#Deltad");
    //hist -> GetXaxis () -> SetRangeUser(1.105,1.13);
    return hist;
}

TF1* PlotSignal(TF1* sig_func)
  {
  sig_func->SetNpx(1000);
  sig_func->SetLineColor(kGreen);
  sig_func->SetLineStyle(2);
  sig_func->SetLineWidth(4);
  sig_func->SetLineWidth(2);
  return sig_func;
  }

TF1* PlotBackground(TF1* bac_func)
  {  
  bac_func->SetNpx(100);
  //bac_func->SetLineStyle(3);
  bac_func->SetLineColor(kBlue);
  bac_func->SetLineWidth(4);
  return bac_func;
  }
 
TF1* PlotMass(TF1* M_func)
  { 
  M_func->SetNpx(10000);
  M_func->SetLineStyle(3);
  M_func->SetLineWidth(4);
  M_func->SetLineColor(kRed);
  return M_func;
  }
  
TCanvas* CreateCanvas(TString name)
  {
    auto* canvas = new TCanvas(name,"", 500,500);
    canvas->SetBottomMargin(0.15);
    canvas->SetLeftMargin (0.15);
    return canvas;
  }
  
TPad* CreatePad1()
  { TPad* pad1 = new TPad (" pad1 "," pad1 " ,0 ,0.3 ,1 ,1);
    pad1 -> SetBottomMargin (0.01);
    pad1 -> SetLeftMargin (0.15);
    //pad1 -> SetLogy();
    pad1 -> Draw ();
    return pad1;
  }

 TPad* CreatePad2()
 {
   TPad* pad2 = new TPad (" pad2 "," pad2 " ,0 ,0.05 ,1 ,0.3);
   //pad2 -> SetGrid();
   pad2 -> SetBottomMargin (0.5);
   pad2 -> SetLeftMargin (0.15);
   pad2 -> Draw ();
   return pad2;
 }

TLegend* legend_plot(TH1F* h1, TF1* M_func, TF1* sig_func, TF1* bac_func)
{
  TLegend* legend = new TLegend(0.5,0.1,0.8,0.3);
  legend->AddEntry(h1,"#Lambda hyperon","pe,X0");
  //legend->AddEntry(M_func,"A[#frac{(1-B)}{#sqrt{2#pi}/#sigma_{1}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{(#sigma_{1})^{2}}}+#frac{(B)}{#sqrt{2#pi}/#sigma_{2}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{#sigma_{2}^{2}}}]+pol2","l");
  legend->AddEntry(M_func,"DSCB+pol2","l");
  legend->AddEntry(sig_func,"DSCB","l");

  //legend->AddEntry(sig_func,"A[#frac{(1-B)}{#sqrt{2#pi}/#sigma_{1}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{(#sigma_{1})^{2}}}+#frac{(B)}{#sqrt{2#pi}/#sigma_{2}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{#sigma_{2}^{2}}}]","l");
  legend->AddEntry(bac_func,"C+Dx+Ex^{2}","l");
  //legend->AddEntry(bac_func,"C","l");
  //legend->AddEntry(bac_func,"#frac{C}{m_{f}-m_{i}}+D(x-0.5(m_{f}-m_{i}))","l");
  legend -> SetLineWidth (0);
  legend->SetTextSize(0.023);
  legend->Draw();
  return legend;
  }
  
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



void ratio_plot_hists(TH1F* h1, TF1* func, Int_t &bins, Double_t min_mass, Double_t max_mass)
{
  TH1F *func_hist = new TH1F("func_hist", "", bins, min_mass, max_mass);
  TH1F* diff_hist = new TH1F("diff_hist", "",bins, min_mass, max_mass);
  for (int i =h1->FindBin(min_mass); i<h1->FindBin(max_mass); i++)
  {
    Double_t f_value= func->Eval(h1->GetBinCenter(i));
    Double_t t_value = h1->GetBinContent(i);
    func_hist->SetBinContent(i,f_value);
    if (h1->GetBinError(i) > 0){
      diff_hist->SetBinContent(i,(t_value-f_value)/h1->GetBinError(i));}
  }
  TH1F* clone =(TH1F*) h1->Clone(); clone->Sumw2();
  clone->Divide(func_hist);
  clone->GetYaxis()->SetRangeUser(0.8,1.2);
  clone = dif_hist(clone);
  clone->SetYTitle("ratio = data/fit");
/*diff_hist = dif_hist(diff_hist);
  diff_hist->Draw();
  TLine* line = new TLine (min_mass,0 ,max_mass ,0);
*/
  clone->Draw();
  TLine* line = new TLine (min_mass,1 ,max_mass ,1);
  line -> SetLineColor ( kRed );
  line->Draw("SAME");
  delete func_hist, diff_hist;
}


void eff_creator(vector<Double_t> &y_bin_low1, vector<Double_t> &pt_bin_low1, vector<Double_t> &yield, vector<Double_t> &yield_error, vector<Double_t> &sim_mc_yield,const char* data, const char* MC)
{  
  
  TFile *file_urqmd = TFile::Open(data);
  TH2D *mc_spectra = (TH2D*)file_urqmd->Get("SimParticles_McLambda/SimParticles_rapidity_SimParticles_pT_McLambda");
  
  TFile *file_dcm = TFile::Open(MC);
  TH2D *Mc = (TH2D*)file_dcm->Get("SimParticles_McLambda/SimParticles_rapidity_SimParticles_pT_McLambda");
  Mc->SetName("Mc"); Mc->SetTitle("Mc"); 

  TH2F* h4 = new TH2F("recons_urqmd", "recons_urqmd", 15,0,3,15,0,3);  h4->SetStats(0);
  //h4->GetZaxis()->SetRangeUser (0.7 ,1.2);
  h4->GetYaxis()->SetRangeUser (0 ,3);  h4->GetXaxis()->SetRangeUser (0 ,3);  h4->GetXaxis()->SetTitle("y_{Lab}");  h4->GetYaxis()->SetTitle("p_{T} (GeV/#it{c}");
  
  
  TH2F* h6 = new TH2F("Mc_urqmd", "Mc_urqmd", 15,0,3,15,0,3);h6->SetStats(0);
  //h4->GetZaxis()->SetRangeUser (0.7 ,1.2);
  h6->GetYaxis()->SetRangeUser (0 ,3);  h6->GetXaxis()->SetRangeUser (0 ,3);  h6->GetXaxis()->SetTitle("y_{Lab}");  h6->GetYaxis()->SetTitle("p_{T} (GeV/#it{c}");
  
  //use clone method for h5
  TH2F* h5 = new TH2F("urqmd_Efficiency", "urqmd_Efficiency", 15,0,3,15,0,3);
  TH2F* recons = new TH2F("recons", "recons", 15,0,3,15,0,3);
    
  for (int i=0;i<225;i++)
  {
    Double_t y= y_bin_low1[i];
    Double_t pT=pt_bin_low1[i];
    Double_t y_bin = int((y+0.1)/0.2 + 1);
    Double_t pT_bin = int((pT+0.1)/0.2 + 1);
    h4->SetBinContent(y_bin, pT_bin, yield[i]);           h4->SetBinError(y_bin, pT_bin,yield_error[i]);
    h5->SetBinContent(y_bin, pT_bin, yield[i]);           h5->SetBinError(y_bin, pT_bin,yield_error[i]);
    auto a = mc_spectra->GetBinContent(y_bin, pT_bin);    auto da = mc_spectra->GetBinError(y_bin, pT_bin);
    h6->SetBinContent(y_bin, pT_bin, a);                  h6->SetBinError(y_bin, pT_bin,da);
    if (sim_mc_yield[i]>0){recons->SetBinContent(y_bin,pT_bin,sim_mc_yield[i]);}
  }

  
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



TH2F* Divide2DHisto(const TH2F* h1, const TH2F* h2, const TString& name) {

  if(h1->GetNbinsX() != h2->GetNbinsX() || h1->GetNbinsY() != h2->GetNbinsY()) {
    throw std::runtime_error("Histograms should have the same number of bins!");
  }

  TH2F* res = dynamic_cast<TH2F*>(h1->Clone(name));
  res->Divide(h2);

  for(int i = 1; i <= h1->GetNbinsX(); ++i){
    for(int j = 1; j <= h1->GetNbinsY(); ++j) {
      const auto a = h1->GetBinContent(i, j);
      const auto b = h2->GetBinContent(i, j);
      const auto da = h1->GetBinError(i, j);
      const auto db = h2->GetBinError(i, j);
      if(a == 0 || b == 0) {
        continue;
      }
      double error = res->GetBinContent(i, j) * sqrt( (da * da) / (a * a) + (db * db) / (b * b) );
      res->SetBinError(i, j, error);
    }
  }
  return res;
}

TH2F* eff_creator_new(vector<Double_t> &y_bin_low1, vector<Double_t> &pt_bin_low1, vector<Double_t> &yield, vector<Double_t> &yield_error, vector<Double_t> &sim_mc_yield,TH2F* Mc_sim, Float_t &BDT, TH2F* Mc_sim_data, vector<Double_t> &pt_y_yield_bdt_0,Int_t check, vector<Double_t> mc_yield_data)
{  
  
  TH2F* recons_data = new TH2F("recons_data", "recons_data", 10,0,3,10,0,3); 
  recons_data->GetYaxis()->SetRangeUser (0 ,3);  recons_data->GetXaxis()->SetRangeUser (0 ,3);  recons_data->GetXaxis()->SetTitle("y_{Lab}");  recons_data->GetYaxis()->SetTitle("p_{T} (GeV/#it{c}");recons_data->SetStats(0);
  
  TH2F* recons_data_mc = new TH2F("recons_data_mc", "recons_data_mc", 10,0,3,10,0,3); 
  
  TH2F* recons_sim = new TH2F("recons_sim", "recons_sim", 10,0,3,10,0,3);
  
  //TH2F* recons_sim_bdt0 = new TH2F("recons_sim_bdt0", "recons_sim_bdt0", 10,0,3,10,0,3);
    
  for (int i=0;i<100;i++)
  {
    Double_t y= y_bin_low1[i];
    Double_t pT=pt_bin_low1[i];
    Double_t y_bin = int((y+0.1)/0.3 + 1);
    Double_t pT_bin = int((pT+0.1)/0.3 + 1);
    recons_data->SetBinContent(y_bin, pT_bin, yield[i]);           recons_data->SetBinError(y_bin, pT_bin,yield_error[i]);
    recons_sim->SetBinContent(y_bin,pT_bin,sim_mc_yield[i]);
    //recons_sim_bdt0->SetBinContent(y_bin,pT_bin,pt_y_yield_bdt_0[i]);//was added if we want ML efficiency correction
    recons_data_mc->SetBinContent(y_bin,pT_bin,mc_yield_data[i]);
  }
  //auto eff = Divide2DHisto(recons_sim, Mc_sim, "eff");
  //auto eff_corr_yield = Divide2DHisto(recons_data, eff, "eff_corr");//change recons_data_mc back to recons_data
  //auto ratio = Divide2DHisto(eff_corr_yield,Mc_sim_data,"ratio");
  //return eff_corr_yield;
  
  TH2F* recons_clone =(TH2F*) recons_data->Clone(); recons_clone->Sumw2();//changed recons_data to recons_data_mc; change back
  TH2F* recons_sim_clone = (TH2F*) recons_sim->Clone(); recons_sim_clone->Sumw2();
  recons_sim_clone->Divide(Mc_sim);
  recons_clone->Divide(recons_sim_clone);
  //recons_clone->Sumw2();Mc_sim_data->Sumw2();
  if (check==1)
  {
  recons_clone->Divide(Mc_sim_data);
  }
  delete recons_sim, recons_sim_clone, recons_data, recons_data_mc;
  return recons_clone;  
  
  
}  


TH2F* eff_creator_sim(vector<Double_t> &y_bin_low1, vector<Double_t> &pt_bin_low1, vector<Double_t> &yield, vector<Double_t> &yield_error, vector<Double_t> &sim_mc_yield,TH2F* Mc_sim, Float_t &BDT, TH2F* Mc_sim_data, vector<Double_t> &pt_y_yield_bdt_0,Int_t check, vector<Double_t> mc_yield_data)
{  
  
  
  TH2F* recons_data_mc = new TH2F("recons_data_mc1", "recons_data_mc1", 10,0,3,10,0,3); 
  
  TH2F* recons_sim = new TH2F("recons_sim1", "recons_sim1", 10,0,3,10,0,3);
  
    
  for (int i=0;i<100;i++)
  {
    Double_t y= y_bin_low1[i];
    Double_t pT=pt_bin_low1[i];
    Double_t y_bin = int((y+0.1)/0.3 + 1);
    Double_t pT_bin = int((pT+0.1)/0.3 + 1);
    recons_sim->SetBinContent(y_bin,pT_bin,sim_mc_yield[i]);
    recons_data_mc->SetBinContent(y_bin,pT_bin,mc_yield_data[i]);
  }
  //auto eff = Divide2DHisto(recons_sim, Mc_sim, "eff");
  //auto eff_corr_yield = Divide2DHisto(recons_data, eff, "eff_corr");//change recons_data_mc back to recons_data
  //auto ratio = Divide2DHisto(eff_corr_yield,Mc_sim_data,"ratio");
  //return eff_corr_yield;
  
  TH2F* recons_clone =(TH2F*) recons_data_mc->Clone(); recons_clone->Sumw2();//changed recons_data to recons_data_mc; change back
  TH2F* recons_sim_clone = (TH2F*) recons_sim->Clone(); recons_sim_clone->Sumw2();
  recons_sim_clone->Divide(Mc_sim);
  recons_clone->Divide(recons_sim_clone);
  //recons_clone->Sumw2();Mc_sim_data->Sumw2();
  if (check==1)
  {
  recons_clone->Divide(Mc_sim_data);
  }
  delete recons_sim, recons_sim_clone, recons_data_mc;
  return recons_clone;  
  
  
} 




TH1F* corrected_yield_hist(vector<Float_t> &bin_content_vector, double range_min,double range_max )
{
    TH1F* hist_corr_yield = new TH1F("hist_corr_yield","",bin_content_vector.size(),range_min,range_max);
      //hist_corr_yield->GetXaxis()->SetTitle("Corrected yield");
    hist_corr_yield->SetYTitle("Counts");
      for (int i=0;i<bin_content_vector.size();i++)
      {
        hist_corr_yield->Fill(bin_content_vector[i]);
      }
      return hist_corr_yield;
}



TGraphMultiErrors* bdt_vs_corrected_yield
(vector<Float_t> &bin_content_vector, vector<Float_t> &cut_value,vector<Float_t> &bin_content_error_vector, Float_t &opt_bdt_min, Float_t &opt_bdt_max)
{
  Double_t x[bin_content_vector.size()];             Double_t y[bin_content_vector.size()];          Double_t exl[bin_content_vector.size()];          Double_t exh[bin_content_vector.size()];   Double_t eylstat[bin_content_vector.size()];       Double_t eyhstat[bin_content_vector.size()]; Double_t eylsys[bin_content_vector.size()];Double_t eyhsys[bin_content_vector.size()];   
   for (int i=0;i<=bin_content_vector.size();i++)
   {
     x[i]=cut_value[i];  y[i]=bin_content_vector[i];  exl[i]=0; exh[i]=0;  eylstat[i]=bin_content_error_vector[i]; eyhstat[i]=bin_content_error_vector[i];
     
     if( (cut_value[i]>opt_bdt_min) && (cut_value[i]<opt_bdt_max) ){   eylsys[i] = bin_content_error_vector[i]; eyhsys[i] = bin_content_error_vector[i]; }
     else{eylsys[i]=0;     eyhsys[i]=0;}
   }
   TGraphMultiErrors* gme = new TGraphMultiErrors("gme", "", bin_content_vector.size(), x, y, exl, exh, eylstat, eyhstat);
   gme->AddYError(bin_content_vector.size(), eylsys, eyhsys);
   gme->GetXaxis()->SetTitle("BDT selection");
   gme->GetYaxis()->SetTitle("Corrected yield");
   gme->GetAttLine(0)->SetLineColor(kRed);//change to blue for reconstructed and green to simulated
   gme->SetMarkerColor( kRed );
   gme->GetAttLine(1)->SetLineColor(kBlue);//change to red for reconstructed
   gme->GetAttFill(1)->SetFillStyle(0);
   gme->SetMarkerStyle(20);
   return gme;
}

template<typename T>
void Print_temp(string words, T value)
{
  std::cout<<words<< value << std::endl;
}



template<typename T, typename U>
class latex_legend_sim {   
  public:          
    
    TLatex* latex_sim_data_draw(T Fitter, Int_t sigma_for_yield, U sum_mc) 
    { // Constructor with parameters
      Double_t signal_s, signal_errs;
      Fitter->Signal(sigma_for_yield,signal_s,signal_errs);  
      TLatex* latex_sim_data = new TLatex ();
      latex_sim_data -> SetNDC ();                latex_sim_data -> SetTextSize (0.03);
      latex_sim_data -> DrawLatex (0.5 ,0.55, Form(  "#color[3]{yield = in %0.1d #sigma %0.1f#pm %0.1f ; MC %0.2f}",sigma_for_yield ,signal_s, signal_errs, sum_mc   ));
      latex_sim_data -> DrawLatex (0.5 ,0.52, Form(  "#color[3]{Yield/MC = %0.3f#pm%0.3f}", (signal_s)/(sum_mc) ,(  (signal_s)/ (sum_mc)   )*(TMath::Sqrt( (signal_errs/signal_s )*(signal_errs/signal_s ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc) ))    ));
      latex_sim_data -> DrawLatex (0.5 ,0.45, Form("#color[3]{#sigma = %0.4f #pm %0.5f;} ",Fitter->GetSigma(),Fitter->GetSigmaUncertainty()));
      latex_sim_data -> DrawLatex (0.5 ,0.4, Form("#color[3]{#mu = %0.4f #pm %0.5f;}", Fitter->GetMean(),Fitter->GetMeanUncertainty()));
      latex_sim_data -> DrawLatex (0.5 ,0.35, Form("#color[3]{#chi^{2}_{red} = %0.3f }", Fitter->GetReducedChiSquare()));
      return  latex_sim_data;
    }
    
    
   TLegend* legend_signal_only_draw(T Fitter)
    {
      TLegend* legend_signal_only = new TLegend(0.5,0.7,0.8,0.85);
      legend_signal_only->SetTextSize(0.023);
      legend_signal_only->AddEntry(Fitter->GetHistoClone(),"#Lambda hyperon signal only","pe,X0");
      legend_signal_only->AddEntry(Fitter->GetSignalFunc(),"DSCB","l");
      return legend_signal_only;  
    }
    
    
};




template<typename T, typename U>
  TLatex* latex_sim_data_draw(T Fitter, Int_t sigma_for_yield, U sum_mc)
  {
  Double_t signal_s, signal_errs;
  Fitter->Signal(sigma_for_yield,signal_s,signal_errs);  
  
  TLatex* latex_sim_data = new TLatex ();
  latex_sim_data -> SetNDC ();                latex_sim_data -> SetTextSize (0.03);
  latex_sim_data -> DrawLatex (0.5 ,0.55, Form(  "#color[3]{yield = in %0.1i #sigma %0.1f#pm %0.1f ; MC %0.2f}",sigma_for_yield ,signal_s, signal_errs, sum_mc   ));
  latex_sim_data -> DrawLatex (0.5 ,0.52, Form(  "#color[3]{Yield/MC = %0.3f#pm%0.3f}", (signal_s)/(sum_mc) ,(  (signal_s)/ (sum_mc)   )*(TMath::Sqrt( (signal_errs/signal_s )*(signal_errs/signal_s ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc) ))    ));
  latex_sim_data -> DrawLatex (0.5 ,0.45, Form("#color[3]{#sigma = %0.4f #pm %0.5f;} ",Fitter->GetSigma(),Fitter->GetSigmaUncertainty()));
  latex_sim_data -> DrawLatex (0.5 ,0.4, Form("#color[3]{#mu = %0.4f #pm %0.5f;}", Fitter->GetMean(),Fitter->GetMeanUncertainty()));
  latex_sim_data -> DrawLatex (0.5 ,0.35, Form("#color[3]{#chi^{2}_{red} = %0.3f }", Fitter->GetReducedChiSquare()));
  return latex_sim_data;
  }
  

template<typename T>
TLegend* legend_signal_only_draw(T Fitter)
{
  TLegend* legend_signal_only = new TLegend(0.5,0.7,0.8,0.85);
  legend_signal_only->SetTextSize(0.023);
  legend_signal_only->AddEntry(Fitter->GetHistoClone(),"#Lambda hyperon signal only","pe,X0");
  legend_signal_only->AddEntry(Fitter->GetSignalFunc(),"DSCB","l");
  return legend_signal_only;  
}


template<typename T>
TLegend* legend_data_draw(T fitter)
{
  TLegend* legend = new TLegend(0.5,0.7,0.8,0.85);                  legend->SetTextSize(0.023);
  legend->AddEntry(fitter->GetHistoClone(),"#Lambda hyperon","pe,X0");
  legend->AddEntry(fitter->GetMassFunc(),Form("DSCB+pol%i",fitter->GetBackgroundRecalcFunc()->GetNpar()-1),"l");
  legend->AddEntry(fitter->GetSignalFunc(),"DSCB","l");
  legend->AddEntry(fitter->GetBackgroundRecalcFunc(),Form("pol%i",fitter->GetBackgroundRecalcFunc()->GetNpar()-1),"l");
  legend->AddEntry(fitter->GetInitialBackgroundFit(),"Initial background fit","l");
  legend->AddEntry(fitter->GetBackgroundFullRangeFunc(),"Initial background extended","l");
  return legend; 
}

