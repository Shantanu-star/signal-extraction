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
//using namespace std;

// #include "AliHFInvMassFitter.h"
// #include "AliHFMassFitter.h"
#endif

enum { kGaus = 0,
  k2Gaus=1, kPol2=2, kNoBk=3, kPol3=6};
  

  Int_t typeb = kPol3;
  Int_t types = 3;
  
  using std::cout;
  using std::endl;
  
  
/*0.45 ,0.452 ,0.454 ,0.456 ,0.458 ,0.46 ,0.462 ,0.464 ,0.466 ,0.468 ,0.47 ,0.472 ,0.474 ,0.476 ,0.478 ,0.48 ,0.482 ,0.484 ,0.486 ,0.488 ,0.49 ,0.492 ,0.494 ,0.496 ,0.498 ,0.5 ,0.502 ,0.504 ,0.506 ,0.508 ,0.51 ,0.512 ,0.514 ,0.516 ,0.518 ,0.52 ,0.522 ,0.524 ,0.526 ,0.528 ,0.53 ,0.532 ,0.534 ,0.536 ,0.538 ,0.54 ,0.542 ,0.544 ,0.546 ,0.548 ,0.55 ,0.552 ,0.554 ,0.556 ,0.558 ,0.56 ,0.562 ,0.564 ,0.566 ,0.568 ,0.57 ,0.572 ,0.574 ,0.576 ,0.578 ,0.58 ,0.582 ,0.584 ,0.586 ,0.588 ,0.59 ,0.592 ,0.594 ,0.596 ,0.598 ,0.6 ,0.602 ,0.604 ,0.606 ,0.608 ,0.61 ,0.612 ,0.614 ,0.616 ,0.618 ,0.62 ,0.622 ,0.624 ,0.626 ,0.628 ,0.63 ,0.632 ,0.634 ,0.636 ,0.638 ,0.64 ,0.642 ,0.644 ,0.646 ,0.648 ,0.65 ,0.652 ,0.654 ,0.656 ,0.658 ,0.66 ,0.662 ,0.664 ,0.666 ,0.668 ,0.67 ,0.672 ,0.674 ,0.676 ,0.678 ,0.68 ,0.682 ,0.684 ,0.686 ,0.688 ,0.69 ,0.692 ,0.694 ,0.696 ,0.698 ,0.7 ,0.702 ,0.704 ,0.706 ,0.708 ,0.71 ,0.712 ,0.714 ,0.716 ,0.718 ,0.72 ,0.722 ,0.724 ,0.726 ,0.728 ,0.73 ,0.732 ,0.734 ,0.736 ,0.738 ,0.74 ,0.742 ,0.744 ,0.746 ,0.748 ,0.75 ,0.752 ,0.754 ,0.756 ,0.758 ,0.76 ,0.762 ,0.764 ,0.766 ,0.768 ,0.77 ,0.772 ,0.774 ,0.776 ,0.778 ,0.78 ,0.782 ,0.784 ,0.786 ,0.788 ,0.79 ,0.792 ,0.794 ,0.796 ,0.798 ,0.8 ,0.802 ,0.804 ,0.806 ,0.808 ,0.81 ,0.812 ,0.814 ,0.816 ,0.818 ,0.82 ,0.822 ,0.824 ,0.826 ,0.828 ,0.83 ,0.832 ,0.834 ,0.836 ,0.838 ,0.84 ,0.842 ,0.844 ,0.846 ,0.848 ,0.85 ,0.852 ,0.854 ,0.856 ,0.858 ,0.86 ,0.862 ,0.864 ,0.866 ,0.868 ,0.87 ,0.872 ,0.874 ,0.876 ,0.878 ,0.88 ,0.882 ,0.884 ,0.886 ,0.888 ,0.89 ,0.892 ,0.894 ,0.896 ,0.898 ,0.9 ,0.902 ,0.904 ,0.906 ,0.908 ,0.91 ,0.912 ,0.914 ,0.916 ,0.918 ,0.92 ,0.922 ,0.924 ,0.926 ,0.928 ,0.93 ,0.932 ,0.934 ,0.936 ,0.938 ,0.94 ,0.942 ,0.944 ,0.946 ,0.948 ,0.95 ,0.952 ,0.954 ,0.956 ,0.958 ,0.96 ,0.962 ,0.964 ,0.966 ,0.968 ,0.97 ,0.972 ,0.974 ,0.976 ,0.978 ,0.98
 */  
//opt is 0.936 
 
  
  void systematics_approach_2M_alice_dscb(TString data,TString sim, TString fNameOut = "") 
    {
    TGaxis::SetMaxDigits(4);
    
    TFile* f = TFile::Open(data.Data(), "read");    
    TTree* t = f->Get<TTree>("t_BDT"); 
    
    TFile* file_sim = TFile::Open(sim.Data(), "read"); 
    TTree* t1 = file_sim->Get<TTree>("t_BDT");
    
    
    Float_t mass, pT, rapidity, MCmass, MCpT, MCrapidity, xgb_preds, MCxgb_preds;Int_t issignal, MCissignal; vector<Int_t> bins={500};  vector<Float_t> mu, mu_error, sigma1, sigma1_error, n1, n1_error, a1, a1_error,n2, n2_error, a2, a2_error, back0, back0_error, back1, back1_error,back2, back2_error, back3, back3_error, back4,back4_error,back5,back5_error, bin_content_vector, bin_content_error_vector, bin_content_vector_sim, bin_content_vector_sim_error, bin_content_vector_sigma_method, bin_content_error_vector_sigma_method,bin_content_vector_mis_match_removed,bin_content_error_vector_mis_match_removed, significance,significance_error, raw_yield_only,raw_yield_only_error,cut_value;
    vector<Float_t> cut={0.45 ,0.452 ,0.454 ,0.456 ,0.458 ,0.46 ,0.462 ,0.464 ,0.466 ,0.468 ,0.47 ,0.472 ,0.474 ,0.476 ,0.478 ,0.48 ,0.482 ,0.484 ,0.486 ,0.488 ,0.49 ,0.492 ,0.494 ,0.496 ,0.498 ,0.5 ,0.502 ,0.504 ,0.506 ,0.508 ,0.51 ,0.512 ,0.514 ,0.516 ,0.518 ,0.52 ,0.522 ,0.524 ,0.526 ,0.528 ,0.53 ,0.532 ,0.534 ,0.536 ,0.538 ,0.54 ,0.542 ,0.544 ,0.546 ,0.548 ,0.55 ,0.552 ,0.554 ,0.556 ,0.558 ,0.56 ,0.562 ,0.564 ,0.566 ,0.568 ,0.57 ,0.572 ,0.574 ,0.576 ,0.578 ,0.58 ,0.582 ,0.584 ,0.586 ,0.588 ,0.59 ,0.592 ,0.594 ,0.596 ,0.598 ,0.6 ,0.602 ,0.604 ,0.606 ,0.608 ,0.61 ,0.612 ,0.614 ,0.616 ,0.618 ,0.62 ,0.622 ,0.624 ,0.626 ,0.628 ,0.63 ,0.632 ,0.634 ,0.636 ,0.638 ,0.64 ,0.642 ,0.644 ,0.646 ,0.648 ,0.65 ,0.652 ,0.654 ,0.656 ,0.658 ,0.66 ,0.662 ,0.664 ,0.666 ,0.668 ,0.67 ,0.672 ,0.674 ,0.676 ,0.678 ,0.68 ,0.682 ,0.684 ,0.686 ,0.688 ,0.69 ,0.692 ,0.694 ,0.696 ,0.698 ,0.7 ,0.702 ,0.704 ,0.706 ,0.708 ,0.71 ,0.712 ,0.714 ,0.716 ,0.718 ,0.72 ,0.722 ,0.724 ,0.726 ,0.728 ,0.73 ,0.732 ,0.734 ,0.736 ,0.738 ,0.74 ,0.742 ,0.744 ,0.746 ,0.748 ,0.75 ,0.752 ,0.754 ,0.756 ,0.758 ,0.76 ,0.762 ,0.764 ,0.766 ,0.768 ,0.77 ,0.772 ,0.774 ,0.776 ,0.778 ,0.78 ,0.782 ,0.784 ,0.786 ,0.788 ,0.79 ,0.792 ,0.794 ,0.796 ,0.798 ,0.8 ,0.802 ,0.804 ,0.806 ,0.808 ,0.81 ,0.812 ,0.814 ,0.816 ,0.818 ,0.82 ,0.822 ,0.824 ,0.826 ,0.828 ,0.83 ,0.832 ,0.834 ,0.836 ,0.838 ,0.84 ,0.842 ,0.844 ,0.846 ,0.848 ,0.85 ,0.852 ,0.854 ,0.856 ,0.858 ,0.86 ,0.862 ,0.864 ,0.866 ,0.868 ,0.87 ,0.872 ,0.874 ,0.876 ,0.878 ,0.88 ,0.882 ,0.884 ,0.886 ,0.888 ,0.89 ,0.892 ,0.894 ,0.896 ,0.898 ,0.9 ,0.902 ,0.904 ,0.906 ,0.908 ,0.91 ,0.912 ,0.914 ,0.916 ,0.918 ,0.92 ,0.922 ,0.924 ,0.926 ,0.928 ,0.93 ,0.932 ,0.934 ,0.936 ,0.938 ,0.94 ,0.942 ,0.944 ,0.946 ,0.948 ,0.95 ,0.952 ,0.954 ,0.956 ,0.958 ,0.96 ,0.962 ,0.964 ,0.966 ,0.968 ,0.97 ,0.972 ,0.974 ,0.976 ,0.978 ,0.98  };
    vector<Float_t> initial_mass={25}; vector<Float_t> final_mass={67};Double_t sigma_range = 4.5;Double_t sigma_for_yield = 11;Int_t sigma_sides=11; Int_t back_pol_degree=4;
    
    t->SetBranchAddress("mass", &mass);t->SetBranchAddress("pT",&pT);t->SetBranchAddress("rapidity",&rapidity);t->SetBranchAddress("BDT_score",&xgb_preds);t->SetBranchAddress("issignal",&issignal);
    t1->SetBranchAddress("pT",&MCpT);t1->SetBranchAddress("rapidity",&MCrapidity);t1->SetBranchAddress("mass", &MCmass);t1->SetBranchAddress("BDT_score", &MCxgb_preds);t1->SetBranchAddress("issignal", &MCissignal);
    
    std::string name="canvas";  TCanvas* canvas = CreateCanvas(name);
    std::string name0="canvas0";  TCanvas* canvas0 = CreateCanvas(name0);
    std::string name01="canvas01";  TCanvas* canvas01 = CreateCanvas(name01);
    
    for (int l=0;l<cut.size();l++){
      
      vector<Double_t> y_bin_low1, pt_bin_low1, yield, yield_error,  mc_yield, sim_mc_yield, mc_yield_for_bin_vec, pt_y_yield_bdt_0, yield_sigma_method, yield_error_sigma_method, yield_minus_mismatch, yield_minus_mismatch_err;
    
    Double_t bin_size = 0.3;
    float y_bin_low=0.9, y_bin_up =1.2;
    for (int i = 0; i<1; i++)
      {
      y_bin_low = y_bin_low+bin_size;      y_bin_up = y_bin_up+bin_size;
      float pt_bin_low =0.,         pt_bin_up =0.3;
      
      for (int ii = 0; ii<1; ii++)
        {
        pt_bin_low = pt_bin_low+bin_size;  pt_bin_up = pt_bin_up+bin_size;
      
        vector<Double_t> yield_for_sigma, yield_for_sigma_error;
        Double_t mc_yield_for_bin, middle_corrected_yield;
        for (int k=0;k<bins.size();k++){ 
        TH1F* real_data = new TH1F("real_data", Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
        real_data -> SetStats (0);
        TH1F* real_data_signal = new TH1F("real_data_signal",Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
        TH1F* real_data_bac = new TH1F("real_data_bac",Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
        TH1F* MC_data_signal = new TH1F("MC_data_signal",Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
        Float_t sum = 0; Float_t sum_mc = 0; Float_t sum_bdt_0=0;
        for (int it =0; it < t->GetEntries(); it++)
          {
          t->GetEntry(it);
        //MC original cut less than 2 is applied remove this
          if ( (pT<pt_bin_up) && (pT>pt_bin_low) && (rapidity >y_bin_low) && (rapidity <y_bin_up) && (xgb_preds>cut[l]) &&(issignal<2) ){
            real_data->Fill(mass);
            if(issignal==1){
            sum_mc+=1;     
            real_data_signal->Fill(mass);
            }
            if(issignal==0){real_data_bac->Fill(mass);}
            }
          }
          
          for (int it1 =0; it1 < t1->GetEntries(); it1++)
          {
          t1->GetEntry(it1);
          //MC original cut less than 2 is applied remove this
          if ( (MCpT<pt_bin_up) && (MCpT>pt_bin_low) && (MCrapidity >y_bin_low) && (MCrapidity <y_bin_up) && (MCissignal==1) && (MCxgb_preds>cut[l]))
          {
            MC_data_signal->Fill(MCmass);
            sum +=1; 
          }
            
           }
          
              int j=0;
            
              if (sum>1){
                j++;
                
                Double_t minMassForFit0 = MC_data_signal->GetMean()-sigma_range*(MC_data_signal->GetRMS());         Double_t maxMassForFit0 = MC_data_signal->GetMean()+sigma_range*(MC_data_signal->GetRMS());              

                //DSCB fit on MC data signal only
                AliHFInvMassFitter* fitter_mc_signal_DSCB = new AliHFInvMassFitter(MC_data_signal, minMassForFit0, maxMassForFit0,   AliHFInvMassFitter::ETypeOfBkg::kNoBk, types);
                fitter_mc_signal_DSCB->SetParticlePdgMass(1.115683);                   fitter_mc_signal_DSCB->SetUseLikelihoodFit();
                fitter_mc_signal_DSCB->SetBoundGaussianMean(1.11567,1.113,1.119);            fitter_mc_signal_DSCB->SetBoundGaussianSigma(0.0012,2);
                fitter_mc_signal_DSCB->SetBoundDSCBa1(1,0,10);                               fitter_mc_signal_DSCB->SetBoundDSCBn1(1,1,100);
                fitter_mc_signal_DSCB->SetBoundDSCBa2(1,0,10);                               fitter_mc_signal_DSCB->SetBoundDSCBn2(1,1,100);
                Bool_t ok = fitter_mc_signal_DSCB->MassFitter(kFALSE);
                
                //DG fit on MC data signal only
                AliHFInvMassFitter* fitter_mc_signal_DG = new AliHFInvMassFitter(h01, minMassForFit0, maxMassForFit0,   AliHFInvMassFitter::ETypeOfBkg::kNoBk, AliHFInvMassFitter::ETypeOfSgn::k2Gaus);
                fitter_mc_signal_DG->SetInitialGaussianMean(1.1156);                                fitter_mc_signal_DG->SetInitialGaussianSigma(0.05);
                fitter_mc_signal_DG->SetInitialFrac2Gaus(0.2);                                      fitter_mc_signal_DG->SetInitialSecondGaussianSigma(0.0009);
                Bool_t ok_DG = fitter_mc_signal_DG->MassFitter(kFALSE);
                if (ok) {
                 
                TF1 *sig_func_MC_DSCB = fitter_mc_signal_DSCB->GetSignalFunc(); sig_func_MC_DSCB = PlotSignal(sig_func_MC_DSCB);
                TF1 *sig_func_MC_DG = fitter_mc_signal_DG->GetSignalFunc(); sig_func_MC_DG = PlotSignal(sig_func_MC_DG);
               
                //Real data signal only fit with DSCB
                AliHFInvMassFitter* fitter0 = new AliHFInvMassFitter(real_data_signal, minMassForFit0, maxMassForFit0,   AliHFInvMassFitter::ETypeOfBkg::kNoBk, types);               
                fitter0->SetInitialGaussianMean(sig_func_MC_DSCB->GetParameter(1));            fitter0->SetInitialGaussianSigma(sig_func_MC_DSCB->GetParameter(2));
                fitter0->SetInitiala1(sig_func_MC_DSCB->GetParameter(3));                     fitter0->SetInitialn1(sig_func_MC_DSCB->GetParameter(4));
                fitter0->SetInitiala2(sig_func_MC_DSCB->GetParameter(5));                     fitter0->SetInitialn2(sig_func_MC_DSCB->GetParameter(6));
               
                fitter0->SetParticlePdgMass(1.115683); fitter0->SetUseLikelihoodFit();
               
                Bool_t ok11 = fitter0->MassFitter(kFALSE);
                TF1 *sig_func_real_data = fitter0->GetSignalFunc(); sig_func_real_data = PlotSignal(sig_func_real_data);
                 
                 Double_t signal_sig, signal_errsig, signal_s, signal_errs, signal_b, signal_errb;
                  fitter0->Signal(sigma_for_yield,signal_s,signal_errs);
                  fitter0->Background(sigma_for_yield,signal_b,signal_errb);
                  fitter0->Significance(sigma_for_yield,signal_sig,signal_errsig);
///*               
                canvas0->Draw();                canvas0->Clear();
                TPad* pad0 = CreatePad1();                pad0 -> Clear();                pad0 -> SetLogy();
                real_data_signal=PlotHistos(real_data_signal);
                real_data_signal->Draw("pe,X0");
                sig_func_real_data->Draw("same");
                auto latex0 = new TLatex ();
                latex0 -> SetNDC ();                latex0 -> SetTextSize (0.03);
                latex0 -> DrawLatex (0.5 ,0.55, Form(  "yield = in %0.1f #sigma %0.1f#pm %0.1f ; MC %0.2f",sigma_for_yield ,signal_s, signal_errs, sum_mc   ));
                latex0 -> DrawLatex (0.5 ,0.52, Form(  "Yield/MC = %0.3f#pm%0.3f", (signal_s)/(sum_mc) ,(  (signal_s)/ (sum_mc)   )*(TMath::Sqrt( (signal_errs/signal_s )*(signal_errs/signal_s ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc) ))    ));
                latex0 -> DrawLatex (0.5 ,0.45, Form("#sigma = %0.4f #pm %0.5f; ",sig_func_real_data->GetParameter(2),sig_func_real_data->GetParError(2)));
                latex0 -> DrawLatex (0.5 ,0.4, Form("#mu = %0.4f #pm %0.5f;", sig_func_real_data->GetParameter(1),sig_func_real_data->GetParError(1)));
                latex0 -> DrawLatex (0.5 ,0.35, Form("#chi^{2}_{red} = %0.3f ", fitter0->GetReducedChiSquare()));
                latex0 -> DrawLatex (0.5 ,0.3, Form("URQMD signal only DSCB fit "));
                latex0->Draw("same");
                TLine* line_signal_select_low = new TLine (sig_func_real_data->GetParameter(1)-sigma_for_yield* fitter0->GetSigma(),0 ,sig_func_real_data->GetParameter(1)-sigma_for_yield* fitter0->GetSigma() ,60000);  line_signal_select_low -> SetLineColor ( kYellow );    line_signal_select_low -> SetLineWidth(2);
                  TLine* line_signal_select_high = new TLine (sig_func_real_data->GetParameter(1)+sigma_for_yield* fitter0->GetSigma(),0 ,sig_func_real_data->GetParameter(1)+sigma_for_yield* fitter0->GetSigma() ,60000);  line_signal_select_high -> SetLineColor ( kYellow );    line_signal_select_high -> SetLineWidth(2);
                line_signal_select_low->Draw("same");
                line_signal_select_high->Draw("same");
                TLegend* legend_signal_only = new TLegend(0.5,0.7,0.8,0.85);
                legend_signal_only->SetTextSize(0.023);
                legend_signal_only->AddEntry(real_data_signal,"#Lambda hyperon signal only","pe,X0");
                legend_signal_only->AddEntry(sig_func_real_data,"DSCB","l");
                legend_signal_only->Draw("same");  
                
                canvas0->cd();
                TPad* pad02 = CreatePad2();                pad02 -> Clear();                pad02->cd();
                ratio_plot_hists(real_data_signal, sig_func_real_data, bins[k], minMassForFit0, maxMassForFit0);
//*/            
                  

                Double_t minMassForFit = sig_func_MC_DSCB->GetParameter(1)- initial_mass[0] * ( sig_func_MC_DSCB->GetParameter(2));
                Double_t maxMassForFit = sig_func_MC_DSCB->GetParameter(1)+ final_mass[0]   * ( sig_func_MC_DSCB->GetParameter(2));
                cout<<"max fit range "<<maxMassForFit<<endl;
                
                
                
                
                
                
                //total mass fit
                AliHFInvMassFitter* fitter = new AliHFInvMassFitter(real_data, minMassForFit, maxMassForFit,  typeb, types);
                
                fitter->SetInitialGaussianMean(sig_func_MC_DSCB->GetParameter(1));                fitter->SetInitialGaussianSigma(sig_func_MC_DSCB->GetParameter(2));
                fitter->SetInitiala1(sig_func_MC_DSCB->GetParameter(3));                          fitter->SetInitialn1(sig_func_MC_DSCB->GetParameter(4));
                fitter->SetInitiala2(sig_func_MC_DSCB->GetParameter(5));                          fitter->SetInitialn2(sig_func_MC_DSCB->GetParameter(6));
                
                fitter->SetPolDegreeForBackgroundFit(back_pol_degree); 
                fitter->SetParticlePdgMass(1.115683);
                fitter->SetNSigma4SideBands(sigma_sides); fitter->SetFitOption("L,E");
                
                
               
                Bool_t ok1 = fitter->MassFitter(kFALSE);//kTRUE to draw
                if (ok1) {
                                                      
                  TF1 *sig_func = fitter->GetSignalFunc(); sig_func = PlotSignal(sig_func);
                  TF1 *bac_func= fitter->GetBackgroundRecalcFunc(); bac_func = PlotBackground(bac_func);
                  TF1 *M_func = fitter->GetMassFunc();M_func = PlotMass(M_func);
                  TF1 *initial_bac = fitter->GetInitialBackgroundFit();  
                  TF1 *initial_bac_extend = fitter->GetBackgroundFullRangeFunc(); initial_bac_extend->SetLineStyle(3); initial_bac_extend->SetLineColor(kMagenta);
                  
                  
                  
                  
                //mis-matched signal
                
                AliHFInvMassFitter* fitter_back = new AliHFInvMassFitter(real_data_bac, minMassForFit, maxMassForFit,  typeb,AliHFInvMassFitter::ETypeOfSgn::k2Gaus);
                
                
                fitter_back->SetInitialGaussianMean(sig_func->GetParameter(1));                   fitter_back->SetInitialGaussianSigma(sig_func->GetParameter(2));
                fitter_back->SetInitialFrac2Gaus(0.2);                                            fitter_back->SetInitialSecondGaussianSigma(0.1);
                //fitter_back->SetInitiala1(sig_func0->GetParameter(3));                          fitter_back->SetInitialn1(sig_func0->GetParameter(4));
                //fitter_back->SetInitiala2(sig_func0->GetParameter(5));                          fitter_back->SetInitialn2(sig_func0->GetParameter(6));
                fitter_back->SetPolDegreeForBackgroundFit(back_pol_degree);                 fitter_back->SetParticlePdgMass(1.115683);
                fitter_back->SetNSigma4SideBands(sigma_sides); //fitter_back->SetFitOption("L,E");*/
                Bool_t ok_back = fitter_back->MassFitter(kFALSE);
                
                TF1 *sig_back = fitter_back->GetSignalFunc(); sig_back = PlotSignal(sig_back);
                TF1 *back_bac_func= fitter_back->GetBackgroundRecalcFunc(); back_bac_func = PlotBackground(back_bac_func);
                TF1 *back_M_func = fitter_back->GetMassFunc();back_M_func = PlotMass(back_M_func);
                Int_t sigma_for_yield0 = 3;
                Double_t  back_s, back_errs;
                fitter_back->Signal(sigma_for_yield0,back_s,back_errs);
                real_data_bac->GetYaxis()->SetRangeUser(0, 1.06*real_data_bac->GetBinContent(real_data_bac->FindBin(1.11567)));  
                canvas01->Draw();                canvas01->Clear();
                TPad* pad001 = CreatePad1();                pad001 -> Clear();                //pad001 -> SetLogy();
                real_data_bac=PlotHistos(real_data_bac);
                
                real_data_bac->SetYTitle("Counts");
                real_data_bac->Draw("pe,X0");
                back_bac_func->Draw("same");
                back_M_func->Draw("same");
                sig_back->Draw("same");
                auto latex00 = new TLatex ();
                latex00 -> SetNDC ();                latex00 -> SetTextSize (0.03);
                latex00 -> DrawLatex (0.5 ,0.55, Form(  "yield in %d #sigma = %0.1f#pm %0.1f ;",sigma_for_yield0, back_s, back_errs  ));
                latex00 -> DrawLatex (0.5 ,0.45, Form("#sigma = %0.4f #pm %0.5f; ",sig_back->GetParameter(2),sig_back->GetParError(2)));
                latex00 -> DrawLatex (0.5 ,0.4, Form("#sigma_{2} = %0.4f #pm %0.5f; ",sig_back->GetParameter(4),sig_back->GetParError(4)));
                latex00 -> DrawLatex (0.5 ,0.35, Form("#mu = %0.4f #pm %0.5f;", sig_back->GetParameter(1),sig_back->GetParError(1)));
                latex00 -> DrawLatex (0.5 ,0.3, Form("#chi^{2}_{red} = %0.3f ", fitter_back->GetReducedChiSquare()));
                latex00->Draw("same");
                TLine* line_mis_match_low = new TLine (sig_func->GetParameter(1)-sigma_sides* fitter0->GetSigma(),0 ,sig_func->GetParameter(1)-sigma_sides* fitter0->GetSigma() ,60000);  line_mis_match_low -> SetLineColor ( kRed );    line_mis_match_low -> SetLineWidth(2);
                  TLine* line_mis_match_high = new TLine (sig_func->GetParameter(1)+sigma_sides* fitter0->GetSigma(),0 ,sig_func->GetParameter(1)+sigma_sides* fitter0->GetSigma() ,60000);  line_mis_match_high -> SetLineColor ( kRed );    line_mis_match_high -> SetLineWidth(2);
                  line_mis_match_low->Draw("same");
                  line_mis_match_high->Draw("same");
                  TLine* line_mis_match_low1 = new TLine (sig_back->GetParameter(1)-sigma_for_yield0* fitter_back->GetSigma(),0 ,sig_back->GetParameter(1)-sigma_for_yield0* fitter_back->GetSigma() ,60000);  line_mis_match_low1 -> SetLineColor ( kGreen );    line_mis_match_low -> SetLineWidth(2);
                  TLine* line_mis_match_high1 = new TLine (sig_back->GetParameter(1)+sigma_for_yield0* fitter_back->GetSigma(),0 ,sig_back->GetParameter(1)+sigma_for_yield0* fitter_back->GetSigma() ,60000);  line_mis_match_high1 -> SetLineColor ( kGreen );    line_mis_match_high -> SetLineWidth(2);
                  line_mis_match_low1->Draw("same");
                  line_mis_match_high1->Draw("same");
                  TLegend* legend_back_signal = new TLegend(0.5,0.7,0.8,0.85);
                  legend_back_signal->SetTextSize(0.023);
                  legend_back_signal->AddEntry(real_data_bac,"#Lambda hyperon mis-match","pe,X0");
                  legend_back_signal->AddEntry(back_M_func,"DG+pol4","l");
                  legend_back_signal->AddEntry(sig_back,"DG","l");
                  legend_back_signal->AddEntry(back_bac_func,"pol4","l");
                  legend_back_signal->Draw("same");   
                
                canvas01->cd();
                TPad* pad002 = CreatePad2();                pad002 -> Clear();                pad002->cd();
                ratio_plot_hists(real_data_bac, back_M_func, bins[k], minMassForFit, maxMassForFit);
                  
                  
                  
                  //Double_t cov_mat = M_func->GetCovarianceMatrix();
                  real_data = PlotHistos(real_data);
///*                  
                  canvas ->Draw();
                  canvas -> Clear();
                  
                  //fitter->DrawHere(canvas, 3, 2);
                  TPad* pad1 = CreatePad1();
                  pad1 -> Clear();
                  pad1->SetLogy();
                  pad1->cd();
                  real_data->Draw("pe,X0");
                  bac_func->Draw("SAME");
                  M_func->Draw("SAME");//ESAME
                  sig_func->Draw("SAME");
                  
                  initial_bac->Draw("same");  initial_bac_extend->Draw("same");
                  
                  TLine* line = new TLine (sig_func_MC_DSCB->GetParameter(1)-sigma_sides* fitter0->GetSigma(),0 ,sig_func_MC_DSCB->GetParameter(1)-sigma_sides* fitter0->GetSigma() ,60000);  line -> SetLineColor ( kRed );    line -> SetLineWidth(2);
                  TLine* line1 = new TLine (sig_func_MC_DSCB->GetParameter(1)+sigma_sides* fitter0->GetSigma(),0 ,sig_func_MC_DSCB->GetParameter(1)+sigma_sides* fitter0->GetSigma() ,60000);  line1 -> SetLineColor ( kRed );    line1 -> SetLineWidth(2);
                  line->Draw("same");
                  line1->Draw("same");
                  TLine* line0 = new TLine (sig_func->GetParameter(1)-sigma_for_yield* fitter->GetSigma(),0 ,sig_func->GetParameter(1)-sigma_for_yield* fitter->GetSigma() ,60000);  line0 -> SetLineColor ( kGreen );    line0 -> SetLineWidth(2);
                  TLine* line01 = new TLine (sig_func->GetParameter(1)+sigma_for_yield* fitter->GetSigma(),0 ,sig_func->GetParameter(1)+sigma_for_yield* fitter->GetSigma() ,60000);  line01 -> SetLineColor ( kGreen );    line01 -> SetLineWidth(2);
                  line0->Draw("same");
                  line01->Draw("same");
//*/                  
                  
                  
                  Double_t sig, errsig, s, errs, b, errb;
                  fitter->Signal(sigma_for_yield,s,errs);
                  fitter->Background(sigma_for_yield,b,errb);
                  fitter->Significance(sigma_for_yield,sig,errsig);
                  
                  Double_t mis_match_removed = s-back_s;  Double_t mis_match_removed_err = TMath::Sqrt( (errs*errs) + (back_errs*back_errs) -2*errs*back_errs    );
                  
///*
                  TLegend* legend = new TLegend(0.5,0.7,0.8,0.85);
                  legend->SetTextSize(0.023);
                  legend->AddEntry(real_data,"#Lambda hyperon","pe,X0");
                  legend->AddEntry(M_func,"DSCB+pol4","l");
                  legend->AddEntry(sig_func,"DSCB","l");
                  legend->AddEntry(bac_func,"pol4","l");
                  legend->AddEntry(initial_bac,"Initial background fit","l");
                  legend->AddEntry(initial_bac_extend,"Initial background extended","l");
                  legend->Draw("same");   
                  
                  auto latex = new TLatex ();
                  latex -> SetNDC ();
                  latex -> SetTextSize (0.03);
                  latex -> DrawLatex (0.5 ,0.65, Form("#mu = %0.4f; #sigma = %0.4f", sig_func->GetParameter(1), sig_func->GetParameter(2) ));
                  //latex -> DrawLatex (0.5 ,0.25, Form("yield = %0.1f#pm %0.1f ; MC %0.2f", fitter->GetRawYield(), fitter->GetRawYieldError(), sum_mc ));
                  //latex -> DrawLatex (0.5 ,0.2, Form("#color[2]{Yield/MC = %0.3f#pm%0.3f}", (fitter->GetRawYield())/(sum_mc),(  (fitter->GetRawYield())/ (sum_mc)   )*(TMath::Sqrt( (fitter->GetRawYieldError()/fitter->GetRawYield() )*(fitter->GetRawYieldError()/fitter->GetRawYield() ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc) ))  ));
                  latex -> DrawLatex (0.5 ,0.25, Form("yield from signal in %0.1f #sigma = %0.1f#pm %0.1f",sigma_for_yield, s, errs ));
                  latex -> DrawLatex (0.5 ,0.2, Form("#color[2]{Yield/MC = %0.3f#pm%0.3f}", (s)/(sum_mc),(  (s)/(sum_mc)   )*(TMath::Sqrt( (errs/s )*(errs/s ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc)-2*( (TMath::Sqrt(sum_mc) * errs)  / (sum_mc*s)) ))  ));
                  latex -> DrawLatex (0.5 ,0.15, Form("#color[2]{(Yield- mis-match yield) /MC = %0.3f#pm%0.5f}", (mis_match_removed)/(sum_mc),(  (mis_match_removed)/(sum_mc)   )*(TMath::Sqrt( (mis_match_removed_err/mis_match_removed )*(mis_match_removed_err/mis_match_removed ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc)-2*( (TMath::Sqrt(sum_mc) * mis_match_removed_err)  / (sum_mc*mis_match_removed)) ))  ));
                  latex -> DrawLatex (0.5 ,0.08, Form("#chi^{2}_{red} = %0.3f ", fitter->GetReducedChiSquare()));
                  latex->Draw("same");
                  
                  
                  canvas->cd();
                  TPad* pad2 = CreatePad2();
                  pad2 -> Clear();
                  pad2->cd();

                  ratio_plot_hists(real_data, M_func, bins[k], minMassForFit, maxMassForFit);
                  
                  
               
                  
                  
                  
                  if (j==1){
                    canvas0->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf(",cut[l]));
                    canvas01->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                    canvas->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                            }
                  else{
                    canvas0->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                    canvas01->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                      canvas->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                      }
                      
//*/                  
                  ///
                  //Create a dataframe
                  //object and store
                  //all these variables in it
                  //write this data frame
                  y_bin_low1.push_back(y_bin_low);
                  pt_bin_low1.push_back(pt_bin_low);
                  sim_mc_yield.push_back(sum);
                  

                  yield_sigma_method.push_back(s);
                  yield_error_sigma_method.push_back(errs);
                  yield_minus_mismatch.push_back(mis_match_removed);
                  yield_minus_mismatch_err.push_back(mis_match_removed_err);
                  yield.push_back(fitter->GetRawYield());
                  yield_error.push_back(fitter->GetRawYieldError());
                  significance.push_back(sig);
                  significance_error.push_back(errsig);
                  mc_yield.push_back(sum_mc);
                  pt_y_yield_bdt_0.push_back(sum_bdt_0);
                  raw_yield_only.push_back(s);
                  raw_yield_only_error.push_back(errs);
                  
                  mu.push_back(sig_func->GetParameter(1));
                  mu_error.push_back(sig_func->GetParError(1));
                  sigma1.push_back(sig_func->GetParameter(2));
                  sigma1_error.push_back(sig_func->GetParError(2));
                  a1.push_back(sig_func->GetParameter(3));
                  a1_error.push_back(sig_func->GetParError(3));
                  n1.push_back(sig_func->GetParameter(4));
                  n1_error.push_back(sig_func->GetParError(4));
                  a2.push_back(sig_func->GetParameter(5));
                  a2_error.push_back(sig_func->GetParError(5));
                  n2.push_back(sig_func->GetParameter(6));
                  n2_error.push_back(sig_func->GetParError(6));
                  
                  
                  back0.push_back(bac_func->GetParameter(0));
                  back0_error.push_back(bac_func->GetParError(0));
                  back1.push_back(bac_func->GetParameter(1));
                  back1_error.push_back(bac_func->GetParError(1));
                  back2.push_back(bac_func->GetParameter(2));
                  back2_error.push_back(bac_func->GetParError(2));
                  back3.push_back(bac_func->GetParameter(3));
                  back3_error.push_back(bac_func->GetParError(3));
                  back4.push_back(bac_func->GetParameter(4));
                  back4_error.push_back(bac_func->GetParError(4));
                  back5.push_back(bac_func->GetParameter(5));
                  back5_error.push_back(bac_func->GetParError(5));
                  
                  
                        }
                        }
                else if (!ok) {
                  cout << ".........I am sorry" << endl;
                  //fitter0->GetHistoClone()->Draw();
                  y_bin_low1.push_back(y_bin_low);
                  pt_bin_low1.push_back(pt_bin_low);
                  yield.push_back(0.0);
                  yield_error.push_back(0.0);
                  mc_yield.push_back(0.0);
                  sim_mc_yield.push_back(0.0);
                  pt_y_yield_bdt_0.push_back(0);
                              }
                //delete fitter0;
      
                         }           
        
        delete real_data; 
          delete MC_data_signal,real_data_signal;
        }
        
        }
                    
      }
      
      
      
    //correction of yield
    TFile *file_sim = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/2M/simulated_files/dcm_prim_signal_m_200_400_bins_10.root");
    TH2F *Mc_sim = (TH2F*)file_sim->Get("SimParticles_McLambda200_400/SimParticles_rapidity_SimParticles_pT_McLambda200_400");
    
    TFile *file_data = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/2M/simulated_files/urqmd_prim_signal_m_200_400_bins_10.root");
    TH2F *Mc_sim_data = (TH2F*)file_data->Get("SimParticles_McLambda200_400/SimParticles_rapidity_SimParticles_pT_McLambda200_400");
    
    TH2F* h2d = eff_creator_new(y_bin_low1, pt_bin_low1, yield, yield_error, sim_mc_yield, Mc_sim, cut[l], Mc_sim_data,pt_y_yield_bdt_0 ,0, mc_yield);
    TH2F* h2d1 = eff_creator_sim(y_bin_low1, pt_bin_low1, yield, yield_error, sim_mc_yield, Mc_sim, cut[l], Mc_sim_data,pt_y_yield_bdt_0 ,0, mc_yield);
    TH2F* h2d_sigma_method = eff_creator_new(y_bin_low1, pt_bin_low1, yield_sigma_method, yield_error_sigma_method, sim_mc_yield, Mc_sim, cut[l], Mc_sim_data,pt_y_yield_bdt_0 ,0, mc_yield);
    TH2F* h2d_mismatch_removed = eff_creator_new(y_bin_low1, pt_bin_low1, yield_minus_mismatch, yield_minus_mismatch_err, sim_mc_yield, Mc_sim, cut[l], Mc_sim_data,pt_y_yield_bdt_0 ,0, mc_yield);

    if (h2d->GetBinContent( h2d_sigma_method->FindBin(1.4,0.5))>0){bin_content_vector.push_back(h2d->GetBinContent( h2d->FindBin(1.4,0.5)) );  bin_content_error_vector.push_back(h2d->GetBinError( h2d->FindBin(1.4,0.5))); cut_value.push_back(cut[l]);     bin_content_vector_sim.push_back(h2d1->GetBinContent( h2d1->FindBin(1.4,0.5)) );  bin_content_vector_sim_error.push_back(h2d1->GetBinError( h2d1->FindBin(1.4,0.5)));    bin_content_vector_sigma_method.push_back(h2d_sigma_method->GetBinContent( h2d_sigma_method->FindBin(1.4,0.5)) );  bin_content_error_vector_sigma_method.push_back(h2d_sigma_method->GetBinError( h2d_sigma_method->FindBin(1.4,0.5)));
      bin_content_vector_mis_match_removed.push_back(h2d_mismatch_removed->GetBinContent( h2d_mismatch_removed->FindBin(1.4,0.5)) );  bin_content_error_vector_mis_match_removed.push_back(h2d_mismatch_removed->GetBinError( h2d_mismatch_removed->FindBin(1.4,0.5)));
    }
    cout<<"BDT is "<< cut[l]<<endl;    
    h2d->SetTitle(Form("corrected yield/simulated yield at BDT%g", cut[l]));

    
    canvas->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf)",cut[l]));
   
   delete h2d; delete h2d1;
      } 
      


     //correction of histograms  
    
    Double_t yield_opt, yield_opt_error,yield_sigma_method_opt, yield_opt_sigma_method_error,yield_sim, yield_sim_error, yield_mismatch_removed, yield_mismatch_removed_err;
    TFile *file_data1 = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/2M/simulated_files/urqmd_prim_signal_m_200_400_bins_10.root");
    TH2F *Mc_sim_data1 = (TH2F*)file_data1->Get("SimParticles_McLambda200_400/SimParticles_rapidity_SimParticles_pT_McLambda200_400");
    Double_t true_yield = Mc_sim_data1->GetBinContent( Mc_sim_data1->FindBin(1.4,0.5));
    Double_t true_yield_error = Mc_sim_data1->GetBinError( Mc_sim_data1->FindBin(1.4,0.5));
    
    auto hist_corr_yield = corrected_yield_hist(bin_content_vector,(*min_element(bin_content_vector.begin(), bin_content_vector.end())) -18000, (*max_element(bin_content_vector.begin(), bin_content_vector.end()))+5000);
    hist_corr_yield->GetXaxis()->SetTitle("Corrected yield");
    auto hist_corr_yield_sim = corrected_yield_hist(bin_content_vector_sim, (*min_element(bin_content_vector_sim.begin(), bin_content_vector_sim.end())) , (*max_element(bin_content_vector_sim.begin(), bin_content_vector_sim.end()))); hist_corr_yield_sim->SetLineColor(kGreen);
    auto hist_corr_yield_sigma_method = corrected_yield_hist(bin_content_vector_sigma_method, (*min_element(bin_content_vector_sigma_method.begin(), bin_content_vector_sigma_method.end()))-25000 , (*max_element(bin_content_vector_sigma_method.begin(), bin_content_vector_sigma_method.end()))+5000);
    hist_corr_yield_sigma_method->GetXaxis()->SetTitle("Corrected yield");hist_corr_yield_sigma_method->SetLineStyle(2);
    
    auto hist_corr_yield_mis_match_removed = corrected_yield_hist(bin_content_vector_mis_match_removed,(*min_element(bin_content_vector_mis_match_removed.begin(), bin_content_vector_mis_match_removed.end()))-1000 , (*max_element(bin_content_vector_mis_match_removed.begin(), bin_content_vector_mis_match_removed.end()))+1000);hist_corr_yield_mis_match_removed->SetLineColor(kYellow);hist_corr_yield_mis_match_removed->SetLineStyle(2);
    hist_corr_yield_mis_match_removed->GetXaxis()->SetTitle("Corrected yield ");
    
    auto hist_mu_fit = corrected_yield_hist(mu, *min_element(mu.begin(), mu.end()) , *max_element(mu.begin(), mu.end()));
    hist_mu_fit->GetXaxis()->SetTitle("mean");
    auto hist_sigma_fit = corrected_yield_hist(sigma1, *min_element(sigma1.begin(), sigma1.end()) , *max_element(sigma1.begin(), sigma1.end()));
    hist_sigma_fit->GetXaxis()->SetTitle("sigma");
    auto hist_a1_fit = corrected_yield_hist(a1, *min_element(a1.begin(), a1.end()) , *max_element(a1.begin(), a1.end()));    hist_a1_fit->GetXaxis()->SetTitle("a1");
    auto hist_n1_fit = corrected_yield_hist(n1, *min_element(n1.begin(), n1.end()) , *max_element(n1.begin(), n1.end()));    hist_n1_fit->GetXaxis()->SetTitle("n1");
    auto hist_a2_fit = corrected_yield_hist(a2, *min_element(a2.begin(), a2.end()) , *max_element(a2.begin(), a2.end()));    hist_a2_fit->GetXaxis()->SetTitle("a2");
    auto hist_n2_fit = corrected_yield_hist(n2, *min_element(n2.begin(), n2.end()) , *max_element(n2.begin(), n2.end()));    hist_n2_fit->GetXaxis()->SetTitle("n2");
    
    auto hist_back0_fit = corrected_yield_hist(back0, *min_element(back0.begin(), back0.end()) , *max_element(back0.begin(), back0.end()));    hist_back0_fit->GetXaxis()->SetTitle("back0");
    auto hist_back1_fit = corrected_yield_hist(back1, *min_element(back1.begin(), back1.end()) , *max_element(back1.begin(), back1.end()));    hist_back1_fit->GetXaxis()->SetTitle("back1");
    auto hist_back2_fit = corrected_yield_hist(back2, *min_element(back2.begin(), back2.end()) , *max_element(back2.begin(), back2.end()));    hist_back2_fit->GetXaxis()->SetTitle("back2");
    auto hist_back3_fit = corrected_yield_hist(back3, *min_element(back3.begin(), back3.end()) , *max_element(back3.begin(), back3.end()));    hist_back3_fit->GetXaxis()->SetTitle("back3");
    auto hist_back4_fit = corrected_yield_hist(back4, *min_element(back4.begin(), back4.end()) , *max_element(back4.begin(), back4.end()));    hist_back4_fit->GetXaxis()->SetTitle("back4");
    auto hist_back5_fit = corrected_yield_hist(back5, *min_element(back5.begin(), back5.end()) , *max_element(back5.begin(), back5.end()));    hist_back5_fit->GetXaxis()->SetTitle("back5");
    
    
    
    
    Double_t std_of_hist= hist_corr_yield->GetStdDev();
    auto gme = bdt_vs_corrected_yield(bin_content_vector, cut_value,bin_content_error_vector);
    auto gme_sigma_method = bdt_vs_corrected_yield(bin_content_vector_sigma_method, cut_value,bin_content_error_vector_sigma_method); gme_sigma_method->SetTitle("sigma method");
    auto gme_mc_method = bdt_vs_corrected_yield(bin_content_vector_sim, cut_value,bin_content_vector_sim_error);gme_mc_method->SetTitle("MC method");
    auto gme_significance = bdt_vs_corrected_yield(significance, cut_value,significance_error);gme_significance->SetTitle("significance");gme_significance->GetYaxis()->SetTitle("#sigma");
    auto gme_raw_yield_only = bdt_vs_corrected_yield(raw_yield_only, cut_value,raw_yield_only_error);gme_raw_yield_only->SetTitle("raw yield from fit");gme_raw_yield_only->GetYaxis()->SetTitle("Raw yield");
    auto gme_mu = bdt_vs_corrected_yield(mu, cut_value,mu_error);gme_mu->GetYaxis()->SetTitle("mean");
    auto gme_sigma1 = bdt_vs_corrected_yield(sigma1, cut_value,sigma1_error);gme_sigma1->GetYaxis()->SetTitle("sigma");
    auto gme_n1 = bdt_vs_corrected_yield(n1, cut_value,n1_error);gme_n1->GetYaxis()->SetTitle("n1");
    auto gme_a1 = bdt_vs_corrected_yield(a1, cut_value,a1_error);gme_a1->GetYaxis()->SetTitle("a1");
    auto gme_n2 = bdt_vs_corrected_yield(n2, cut_value,n2_error);gme_n2->GetYaxis()->SetTitle("n2");
    auto gme_a2 = bdt_vs_corrected_yield(a2, cut_value,a2_error);gme_a2->GetYaxis()->SetTitle("a2");
    
    auto gme_back0 = bdt_vs_corrected_yield(back0, cut_value,back0_error);gme_back0->GetYaxis()->SetTitle("back0");
    auto gme_back1 = bdt_vs_corrected_yield(back1, cut_value,back1_error);gme_back1->GetYaxis()->SetTitle("back1");
    auto gme_back2 = bdt_vs_corrected_yield(back2, cut_value,back2_error);gme_back2->GetYaxis()->SetTitle("back2");
    auto gme_back3 = bdt_vs_corrected_yield(back3, cut_value,back3_error);gme_back3->GetYaxis()->SetTitle("back3");
    auto gme_back4 = bdt_vs_corrected_yield(back4, cut_value,back4_error);gme_back4->GetYaxis()->SetTitle("back4");
    auto gme_back5 = bdt_vs_corrected_yield(back5, cut_value,back5_error);gme_back5->GetYaxis()->SetTitle("back5");
    

   for (int i=0;i<=bin_content_vector.size();i++)
   {     
     if( (cut_value[i]>0.9359) && (cut_value[i]<0.9361) ){     yield_opt=bin_content_vector[i]; yield_opt_error=bin_content_error_vector[i];  yield_sim=bin_content_vector_sim[i]; yield_sim_error=bin_content_vector_sim_error[i]; yield_sigma_method_opt = bin_content_vector_sigma_method[i];yield_opt_sigma_method_error=bin_content_error_vector_sigma_method[i]; yield_mismatch_removed=bin_content_vector_mis_match_removed[i];yield_mismatch_removed_err=bin_content_error_vector_mis_match_removed[i];}
   }
   
      TH1F* central_yield_sigma = new TH1F("central_yield_sigma","central_yield_sigma",1,yield_sigma_method_opt - yield_opt_sigma_method_error,yield_sigma_method_opt + yield_opt_sigma_method_error);     central_yield_sigma->SetMarkerColor(kOrange);     central_yield_sigma->SetLineColor(kOrange);
      TH1F* central_yield_sim = new TH1F("central_yield_sim","central_yield_sim",1,yield_sim - yield_sim_error,yield_sim + yield_sim_error);     central_yield_sim->SetMarkerColor(kPink);     central_yield_sim->SetLineColor(kPink);
      TH1F* true_yield_hist = new TH1F("true_yield_hist","true_yield_hist",1,true_yield - true_yield_error,true_yield + true_yield_error);     true_yield_hist->SetMarkerColor(kCyan);     true_yield_hist->SetLineColor(kCyan);
      TH1F* mis_match_removed_central_hist = new TH1F("mis_match_removed_central_hist","mis_match_removed_central_hist",1,yield_mismatch_removed - yield_mismatch_removed_err,yield_mismatch_removed + yield_mismatch_removed_err);     mis_match_removed_central_hist->SetMarkerColor(kPink);     mis_match_removed_central_hist->SetLineColor(kPink);
      for (int i=0;i<1;i++)     {     central_yield_sigma->Fill(yield_sigma_method_opt); true_yield_hist->Fill(true_yield); central_yield_sim->Fill(yield_sim); mis_match_removed_central_hist->Fill(yield_mismatch_removed);   }

   
/*when using RawYield method uncomment
  TLine* line = new TLine (yield_opt,0 ,yield_opt ,6);  line -> SetLineColor ( kRed );    line -> SetLineWidth(2);
  TLine* line1 = new TLine (yield_opt-yield_opt_error,0 , yield_opt-yield_opt_error,6);    line1 -> SetLineColor ( kRed );line1->SetLineStyle(3);line1 -> SetLineWidth(2);
  TLine* line2 = new TLine (yield_opt+yield_opt_error,0 , yield_opt+yield_opt_error,6);    line2 -> SetLineColor ( kRed );line2->SetLineStyle(3);    line2 -> SetLineWidth(2);
  */
  TLine* line = new TLine (yield_sigma_method_opt,0 ,yield_sigma_method_opt ,6);  line -> SetLineColor ( kRed );    line -> SetLineWidth(6);
  TLine* line1 = new TLine (yield_sigma_method_opt-yield_opt_sigma_method_error,0 , yield_sigma_method_opt-yield_opt_sigma_method_error,6);    line1 -> SetLineColor ( kRed );line1->SetLineStyle(3);line1 -> SetLineWidth(6);
  TLine* line2 = new TLine (yield_sigma_method_opt+yield_opt_sigma_method_error,0 , yield_sigma_method_opt+yield_opt_sigma_method_error,6);    line2 -> SetLineColor ( kRed );line2->SetLineStyle(3);    line2 -> SetLineWidth(6);

    
    
    TLine* line3 = new TLine (true_yield,0 ,true_yield ,6);  line3 -> SetLineColor ( kBlack );    line3 -> SetLineWidth(6);
    
    TLine* line4 = new TLine (true_yield+sqrt(true_yield),0 , true_yield+sqrt(true_yield),6);    line4 -> SetLineColor ( kBlack );line4->SetLineStyle(3);line4 -> SetLineWidth(6);
    TLine* line7 = new TLine (true_yield-sqrt(true_yield),0 , true_yield-sqrt(true_yield),6);    line7 -> SetLineColor ( kBlack );line7->SetLineStyle(3);line7 -> SetLineWidth(6);
    
      TLine* line5 = new TLine (yield_sim,0 , yield_sim,6);    line5 -> SetLineColor ( kMagenta );    line5 -> SetLineWidth(6);
  TLine* line6 = new TLine (yield_sim-yield_sim_error,0 , yield_sim-yield_sim_error,6);    line6 -> SetLineColor ( kMagenta );line6->SetLineStyle(3);    line6 -> SetLineWidth(6);
    
    TLegend* legend = new TLegend(0.55,0.7,0.9,0.9);
    //legend->AddEntry(hist_corr_yield_sim,"#Lambda hyperon corrected yield MC","pe");
  //legend->AddEntry(line5,"default corrected yield","l");
  //legend->AddEntry(line6,"default yield #pm #sigma_{default yield}","l");
  legend->AddEntry(hist_corr_yield_sigma_method,"#Lambda hyperon corrected yield","pe");
  //legend->AddEntry(central_yield_sigma,"default corrected yield","pe");
  legend->AddEntry(true_yield_hist,"True yield","pe");
  legend->AddEntry(hist_corr_yield_mis_match_removed,"mis-match removed","pe");
  legend->AddEntry(mis_match_removed_central_hist,"default","pe");
  //legend->AddEntry(line,"default corrected yield","l");
  //legend->AddEntry(line1,"default yield #pm #sigma_{default yield}","l");
  //legend->AddEntry(line3,"True yield","l");
  //legend->AddEntry(line4,"True yield + #sigma_{True yield}","l");
  legend->SetTextSize(0.023);

   auto c = new TCanvas();
    c->Draw();
//
    //hist_corr_yield->Draw("pe");
    hist_corr_yield_sigma_method->Draw("pe");
    //hist_corr_yield_sim->Draw("pe,same");
    //central_yield_sigma->Draw("pe,same");
    //central_yield_sim->Draw("pe,same");
    true_yield_hist->Draw("pe,same");
    hist_corr_yield_mis_match_removed->Draw("pe,same");
    mis_match_removed_central_hist->Draw("pe,same");
    /*line->Draw("SAME");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");*/
  legend->Draw("same");

   //gme->Draw("APS ; Z ; 5 s=0.5");
   c->Update();
        TFile *myfile = TFile::Open("corrected_yield_DSCB_alice.root","recreate");
    gme->Write();
    gme_mc_method->Write();
    gme_raw_yield_only->Write();
    gme_sigma_method->Write();
    gme_significance->Write();
    gme_mu->Write();    gme_a1->Write();    gme_sigma1->Write();    gme_n1->Write();  gme_a2->Write(); gme_n2->Write();
    gme_back0->Write(); gme_back1->Write();gme_back2->Write();gme_back3->Write();gme_back4->Write();gme_back5->Write();
    hist_corr_yield->Write();  central_yield_sigma->Write(); central_yield_sim->Write();true_yield_hist->Write();
    hist_corr_yield_sigma_method->Write();   hist_corr_yield_mis_match_removed->Write();mis_match_removed_central_hist->Write();
    hist_corr_yield_sim->Write();
    hist_mu_fit->Write();    hist_sigma_fit->Write();    hist_n1_fit->Write();    hist_a1_fit->Write();hist_a2_fit->Write();hist_n2_fit->Write();
    hist_back0_fit->Write();    hist_back1_fit->Write();    hist_back2_fit->Write();    hist_back3_fit->Write();    hist_back4_fit->Write();
    hist_back5_fit->Write(); 
    
    cout<<"central yield is"<<yield_sigma_method_opt<<" plus minus "<<yield_opt_sigma_method_error<<endl;
    
    auto c11 = new TCanvas();
    c11->Draw();
    hist_mu_fit->Draw("pe");
    
    
    std::string name1="canvas1";  TCanvas* canvas1 = CreateCanvas(name1);
    canvas1->Clear();    canvas1->Draw();
    hist_corr_yield_sigma_method->Draw("pe");
    hist_corr_yield_sim->Draw("pe,same");
    central_yield_sigma->Draw("pe,same");
    central_yield_sim->Draw("pe,same");
    true_yield_hist->Draw("pe,same");
    canvas1->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/corrected_yield.pdf(");
    
    canvas1->Clear();    canvas1->Draw();
    hist_mu_fit->Draw("pe");canvas1->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/corrected_yield.pdf");
    
    canvas1->Clear();    canvas1->Draw();
    hist_sigma_fit->Draw("pe");canvas1->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/corrected_yield.pdf");
    
    canvas1->Clear();    canvas1->Draw();
    hist_n1_fit->Draw("pe");canvas1->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/corrected_yield.pdf");
    
    canvas1->Clear();    canvas1->Draw();
    hist_a1_fit->Draw("pe");canvas1->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/corrected_yield.pdf");
    
    canvas1->Clear();    canvas1->Draw();
    hist_n2_fit->Draw("pe");canvas1->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/corrected_yield.pdf");
    
    canvas1->Clear();    canvas1->Draw();
    hist_a2_fit->Draw("pe");canvas1->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/corrected_yield.pdf)");
    
  }
  
  int main(int argc,char** argv) {
    TString name = argv[1];
    TString name2 = argv[2];
    systematics_approach_2M_alice_dscb(name, name2);
    return 0;
  }
          

                

  
  
