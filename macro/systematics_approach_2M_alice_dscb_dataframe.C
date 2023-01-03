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
#include "plot_tools.C"
#include <TString.h>
#include <vector>
#include <numeric>
#endif

enum { kGaus = 0,
  k2Gaus=1, kPol2=2, kNoBk=3, kPol3=6};
  

  Int_t typeb = kPol3;
  Int_t types = 3;
  
  using std::cout;
  using std::endl;
  
 
  
  void systematics_approach_2M_alice_dscb_dataframe(TString data,TString sim, TString fNameOut = "") 
    {
      ROOT::EnableImplicitMT(4);
      ROOT::RDataFrame df_sim("t_BDT", sim.Data());      ROOT::RDataFrame df_data("t_BDT", data.Data());
      vector<Int_t> bins={250};
      vector<Float_t> initial_mass={17}; vector<Float_t> final_mass={20};Double_t sigma_range = 3;Double_t sigma_for_yield = 11;Int_t sigma_sides=8; Int_t back_pol_degree=4; Int_t sigma_for_yield0 = 3;
      vector<Float_t> pt_min{0.9,1.2};
      vector<Float_t> y_min{1.5,1.8};
      vector<Float_t> cut{0.53};
      
      std::string name="canvas";  TCanvas* canvas_data = CreateCanvas(name);    std::string name0="canvas_sim";  TCanvas* canvas_sim = CreateCanvas(name0);
      
      TH1F* hist_sim_data =  (TH1F*)(df_sim.Filter(Form("(pT >= %.2f) && (pT < %.2f) && (rapidity >= %0.2f) && (rapidity < %0.2f) && (BDT_score > %.3f) && (issignal == 1)",pt_min[0],pt_min[0+1],y_min[0],y_min[0+1],cut[0])).Histo1D({"hist_inv_mass", "hist_inv_mass",bins[0],1.077,1.2}, "mass"))->Clone();
      hist_sim_data = PlotHistos(hist_sim_data);
      
      TH1F* hist_data =  (TH1F*)(df_data.Filter(Form("(pT >= %.2f) && (pT < %.2f) && (rapidity >= %0.2f) && (rapidity < %0.2f)&& (BDT_score > %.3f) && (issignal < 2)",pt_min[0],pt_min[0+1],y_min[0],y_min[0+1],cut[0])).Histo1D({"hist_inv_mass", "hist_inv_mass",bins[0],1.077,1.2}, "mass"))->Clone(); 
      hist_data = PlotHistos(hist_data);

      
      Double_t minMassForFit_sim = hist_sim_data->GetMean()-sigma_range*(hist_sim_data->GetRMS());         Double_t maxMassForFit_sim = hist_sim_data->GetMean()+sigma_range*(hist_sim_data->GetRMS());
      AliHFInvMassFitter* fitter_mc_signal_DSCB = new AliHFInvMassFitter(hist_sim_data, minMassForFit_sim, maxMassForFit_sim,    AliHFInvMassFitter::ETypeOfBkg::kNoBk, types);
      fitter_mc_signal_DSCB->SetParticlePdgMass(1.115683);                         fitter_mc_signal_DSCB->SetUseLikelihoodFit();
      fitter_mc_signal_DSCB->SetBoundGaussianMean(1.11567,1.113,1.119);            fitter_mc_signal_DSCB->SetBoundGaussianSigma(0.0012,2);
      fitter_mc_signal_DSCB->SetBoundDSCBa1(1,0,10);                               fitter_mc_signal_DSCB->SetBoundDSCBn1(1,0,100);
      fitter_mc_signal_DSCB->SetBoundDSCBa2(1,0,10);                               fitter_mc_signal_DSCB->SetBoundDSCBn2(1,0,100);
      
      Bool_t ok = fitter_mc_signal_DSCB->MassFitter(kFALSE);
      TF1 *sig_func_sim_data = fitter_mc_signal_DSCB->GetSignalFunc(); sig_func_sim_data = PlotSignal(sig_func_sim_data);//signal func only
      
      canvas_sim->Draw();
      TPad* pad0_sim = CreatePad1();                pad0_sim -> Clear();                pad0_sim -> SetLogy();      
      hist_sim_data->Draw("pe,X0");
      sig_func_sim_data->Draw("same");
      

      TLatex* latex_sim_data = (TLatex*) latex_sim_data_draw(fitter_mc_signal_DSCB,sigma_for_yield, hist_sim_data->GetEntries());
      latex_sim_data -> DrawLatex (0.5 ,0.3, Form("DCM signal only "));
      latex_sim_data->Draw("same");
      TLegend* legend_sim_data = (TLegend*) legend_signal_only_draw(fitter_mc_signal_DSCB);
      legend_sim_data->Draw("same");
      
      canvas_sim->cd();
      TPad* pad1_sim = CreatePad2();                pad1_sim -> Clear();                pad1_sim->cd();
      ratio_plot_hists(hist_sim_data, sig_func_sim_data, bins[0], minMassForFit_sim, maxMassForFit_sim); 
      
      
      
      Double_t minMassForFit = fitter_mc_signal_DSCB->GetMean()- initial_mass[0] * fitter_mc_signal_DSCB->GetSigma();
      Double_t maxMassForFit = fitter_mc_signal_DSCB->GetMean()+ final_mass[0]   * fitter_mc_signal_DSCB->GetSigma();
      //total mass fit
      AliHFInvMassFitter* fitter = new AliHFInvMassFitter(hist_data, minMassForFit, maxMassForFit,  typeb, types);
      fitter->SetInitialGaussianMean(sig_func_sim_data->GetParameter(1));                
      fitter->SetFixDSCBa1(sig_func_sim_data->GetParameter(3));                          fitter->SetFixDSCBn1(sig_func_sim_data->GetParameter(4));
      fitter->SetFixDSCBa2(sig_func_sim_data->GetParameter(5));                          fitter->SetFixDSCBn2(sig_func_sim_data->GetParameter(6));
      fitter->SetFixGaussianSigma(sig_func_sim_data->GetParameter(2));                   fitter->SetPolDegreeForBackgroundFit(back_pol_degree); 
      fitter->SetParticlePdgMass(1.115683);                                             fitter->SetNSigma4SideBands(sigma_sides); 
      fitter->SetUseLikelihoodFit();
      Bool_t ok1 = fitter->MassFitter(kFALSE);
      
      TF1 *sig_func_data = fitter->GetSignalFunc();                      sig_func_data = PlotSignal(sig_func_data);//signal func only
      TF1 *bac_func_data = fitter->GetBackgroundRecalcFunc();            bac_func_data = PlotBackground(bac_func_data);
      TF1 *total_func_data = fitter->GetMassFunc();                      total_func_data = PlotMass(total_func_data);
      TF1 *initial_bac = fitter->GetInitialBackgroundFit();  
      TF1 *initial_bac_extend = fitter->GetBackgroundFullRangeFunc(); initial_bac_extend->SetLineStyle(3); initial_bac_extend->SetLineColor(kMagenta);
      
      canvas_data->Clear();
      canvas_data->Draw();
      TPad* pad0_data = CreatePad1();                pad0_data -> Clear();                pad0_data -> SetLogy();      
      hist_data->Draw("pe,X0");
      bac_func_data->Draw("SAME");
      total_func_data->Draw("SAME");//ESAME
      sig_func_data->Draw("same");
      initial_bac->Draw("same");  initial_bac_extend->Draw("same");
      auto legend =  legend_data_draw(fitter);
      legend->Draw("same"); 
               
      
      canvas_data->cd();
      TPad* pad1_data = CreatePad2();                pad1_data -> Clear();                pad1_data->cd();
      ratio_plot_hists(hist_data, total_func_data, bins[0], minMassForFit, maxMassForFit); 

      
    
  }
  
  int main(int argc,char** argv) {
    TString data_path = argv[1];
    TString sim_path = argv[2];
    systematics_approach_2M_alice_dscb_dataframe(data_path, sim_path);
    return 0;
  }
          

                

  
  
