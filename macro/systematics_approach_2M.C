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
/*0.08, 0.0825, 0.085, 0.0875, 0.09, 0.0925, 0.095, 0.0975, 0.1, 0.1025, 0.105, 0.1075, 0.11, 0.1125, 0.115, 0.1175, 0.12, 0.1225, 0.125, 0.1275, 0.13, 0.1325, 0.135, 0.1375, 0.1425, 0.145, 0.1475, 0.15, 0.1525, 0.155, 0.1575, 0.16, 0.1625, 0.165, 0.1675, 0.17, 0.1725, 0.175, 0.1775, 0.18, 0.1825, 0.185, 0.1875, 0.19, 0.1925, 0.195, 0.1975, 0.2, 0.2025, 0.205, 0.2075, 0.21, 0.2125, 0.215, 0.2175, 0.22, 0.2225, 0.225, 0.2275, 0.23, 0.2325, 0.235, 0.2375, 0.24, 0.2425, 0.245, 0.2475, 0.25, 0.2525, 0.255, 0.2575, 0.26, 0.2625, 0.265, 0.2675, 0.27, 0.2725, 0.275, 0.2775, 0.28, 0.2825, 0.285, 0.2875, 0.29, 0.2925, 0.295, 0.2975, 0.3, 0.3025, 0.305, 0.3075, 0.31, 0.3125, 0.315, 0.3175, 0.32, 0.3225, 0.325, 0.3275, 0.33, 0.3325, 0.335, 0.3375, 0.34, 0.3425, 0.345, 0.3475, 0.35, 0.3525, 0.355, 0.3575, 0.36, 0.3625, 0.365, 0.3675, 0.37, 0.3725, 0.375, 0.3775, 0.38, 0.3825, 0.385, 0.3875, 0.39, 0.3925, 0.395, 0.3975, 0.4025, 0.405, 0.4075, 0.41, 0.4125, 0.415, 0.4175, 0.42, 0.4225, 0.425, 0.4275, 0.43, 0.4325, 0.435, 0.4375, 0.44, 0.4425, 0.445, 0.4475, 0.45, 0.4525, 0.455, 0.4575, 0.46, 0.4625, 0.465, 0.4675, 0.47, 0.4725, 0.475, 0.4775, 0.48, 0.485, 0.4875, 0.49,  0.495, 0.4975, 0.5, 0.5025, 0.505, 0.5075, 0.51, 0.5125, 0.515, 0.5175, 0.52, 0.5225, 0.525, 0.5275, 0.53, 0.5325, 0.535, 0.5375, 0.54, 0.5425, 0.545, 0.5475, 0.55, 0.5525, 0.5575, 0.56, 0.5625, 0.565, 0.5675, 0.57, 0.5725, 0.575, 0.5775, 0.58, 0.5825, 0.585, 0.5875, 0.59, 0.5925, 0.595, 0.5975, 0.6, 0.6025, 0.605, 0.6075, 0.61, 0.6125, 0.615, 0.6175, 0.62, 0.6225, 0.625, 0.6275, 0.63, 0.6325, 0.635, 0.6375, 0.64, 0.6425, 0.645, 0.65, 0.6525, 0.655, 0.6575, 0.66, 0.6625, 0.665, 0.6675, 0.67, 0.6725, 0.675, 0.6775, 0.68, 0.6825, 0.685, 0.6875, 0.69, 0.6925, 0.695, 0.6975, 0.7, 0.7025, 0.705, 0.7075, 0.71, 0.7125, 0.715, 0.7175, 0.72, 0.7225, 0.725, 0.7275, 0.73, 0.7325, 0.735, 0.7375, 0.74, 0.7425, 0.745, 0.7475, 0.75, 0.7525, 0.755, 0.7575, 0.76, 0.7625, 0.765, 0.7675, 0.77, 0.7725,   0.78, 0.7825, 0.785, 0.7875, 0.79, 0.7925, 0.795, 0.7975, 0.8, 0.8025, 0.805, 0.81, 0.8125, 0.815, 0.8175, 0.82, 0.8225, 0.8275, 0.83, 0.835, 0.84, 0.8575, 0.86, 0.865, 0.8675, 0.87, 0.8725, 0.88, 0.8875, 0.8925, 0.895*/


enum { kGaus = 0,
  k2Gaus=1, kPol2=2, kNoBk=3, kPol3=6};
  

  Int_t typeb = kPol3;
  Int_t types = k2Gaus;
  
  using std::cout;
  using std::endl;
  
 
  
  void systematics_approach_2M(TString data,TString sim, TString fNameOut = "") 
    {
    TGaxis::SetMaxDigits(4);
    
    TFile* f = TFile::Open(data.Data(), "read");    
    TTree* t = f->Get<TTree>("t_BDT"); 
    
    TFile* file_sim = TFile::Open(sim.Data(), "read"); 
    TTree* t1 = file_sim->Get<TTree>("t_BDT");
    
    
    Float_t mass, pT, rapidity, MCmass, MCpT, MCrapidity, xgb_preds, MCxgb_preds;Int_t issignal, MCissignal; vector<Int_t> bins={500};  vector<Float_t> mu,mu_error,sigma1,sigma1_error,sigma2,sigma2_error,frac2gaus,frac2gaus_error, bin_content_vector, bin_content_error_vector, bin_content_vector_sim,bin_content_vector_sim_error,bin_content_vector_sigma_method,bin_content_error_vector_sigma_method, significance,significance_error, raw_yield_only,raw_yield_only_error,cut_value;
    vector<Float_t> cut={  0.02 ,0.022 ,0.024 ,0.026 ,0.028 ,0.03 ,0.032 ,0.034 ,0.036 ,0.038 ,0.04 ,0.042 ,0.044 ,0.046 ,0.048 ,0.05 ,0.052 ,0.054 ,0.056 ,0.058 ,0.06 ,0.062 ,0.064 ,0.066 ,0.068 ,0.07 ,0.072 ,0.074 ,0.076 ,0.078 ,0.08 ,0.082 ,0.084 ,0.086 ,0.088 ,0.09 ,0.092 ,0.094 ,0.096 ,0.098 ,0.1 ,0.102 ,0.104 ,0.106 ,0.108 ,0.11 ,0.112 ,0.114 ,0.116 ,0.118 ,0.12 ,0.122 ,0.124 ,0.126 ,0.128 ,0.13 ,0.132 ,0.134 ,0.136 ,0.138 ,0.14 ,0.142 ,0.144 ,0.146 ,0.148 ,0.15 ,0.152 ,0.154 ,0.156 ,0.158 ,0.16 ,0.162 ,0.164 ,0.166 ,0.168 ,0.17 ,0.172 ,0.174 ,0.176 ,0.178 ,0.18 ,0.182 ,0.184 ,0.186 ,0.188 ,0.19 ,0.192 ,0.194 ,0.196 ,0.198 ,0.2 ,0.202 ,0.204 ,0.206 ,0.208 ,0.21 ,0.212 ,0.214 ,0.216 ,0.218 ,0.22 ,0.222 ,0.224 ,0.226 ,0.228 ,0.23 ,0.232 ,0.234 ,0.236 ,0.238 ,0.24 ,0.242 ,0.244 ,0.246 ,0.248 ,0.25 ,0.252 ,0.254 ,0.256 ,0.258 ,0.26 ,0.262 ,0.264 ,0.266 ,0.268 ,0.27 ,0.272 ,0.274 ,0.276 ,0.278 ,0.28 ,0.282 ,0.284 ,0.286 ,0.288 ,0.29 ,0.292 ,0.294 ,0.296 ,0.298 ,0.3 ,0.302 ,0.304 ,0.306 ,0.308 ,0.31 ,0.312 ,0.314 ,0.316 ,0.318 ,0.32 ,0.322 ,0.324 ,0.326 ,0.328 ,0.33 ,0.332 ,0.334 ,0.336 ,0.338 ,0.34 ,0.342 ,0.344 ,0.346 ,0.348 ,0.35 ,0.352 ,0.354 ,0.356 ,0.358 ,0.36 ,0.362 ,0.364 ,0.366 ,0.368 ,0.37 ,0.372 ,0.374 ,0.376 ,0.378 ,0.38 ,0.382 ,0.384 ,0.386 ,0.388 ,0.39 ,0.392 ,0.394 ,0.396 ,0.398 ,0.4 ,0.402 ,0.404 ,0.406 ,0.408 ,0.41 ,0.412 ,0.414 ,0.416 ,0.418 ,0.42 ,0.422 ,0.424 ,0.426 ,0.428 ,0.43 ,0.432 ,0.434 ,0.436 ,0.438 ,0.44 ,0.442 ,0.444 ,0.446 ,0.448 ,0.45 ,0.452 ,0.454 ,0.456 ,0.458 ,0.46 ,0.462 ,0.464 ,0.466 ,0.468 ,0.47 ,0.472 ,0.474 ,0.476 ,0.478 ,0.48 ,0.482 ,0.484 ,0.486 ,0.488 ,0.49 ,0.492 ,0.494 ,0.496 ,0.498 ,0.5 ,0.502 ,0.504 ,0.506 ,0.508 ,0.51 ,0.512 ,0.514 ,0.516 ,0.518 ,0.52 ,0.522 ,0.524 ,0.526 ,0.528 ,0.532 ,0.534 ,0.536 ,0.538 ,0.54 ,0.542 ,0.544 ,0.546 ,0.548 ,0.55 ,0.552 ,0.554 ,0.556 ,0.558 ,0.56 ,0.562 ,0.564 ,0.566 ,0.568 ,0.57 ,0.572 ,0.574 ,0.576 ,0.578 ,0.58 ,0.582 ,0.584 ,0.586 ,0.588 ,0.59 ,0.592 ,0.594 ,0.596 ,0.598 ,0.6 ,0.602 ,0.604 ,0.606 ,0.608 ,0.61 ,0.612 ,0.614 ,0.616 ,0.618 ,0.62 ,0.622 ,0.624 ,0.626 ,0.628 ,0.63 ,0.632 ,0.634 ,0.636 ,0.638 ,0.64 ,0.642 ,0.644 ,0.646 ,0.648 ,0.65 ,0.652 ,0.654 ,0.656 ,0.658 ,0.66 ,0.662 ,0.664 ,0.666 ,0.668 ,0.67 ,0.672 ,0.674 ,0.676 ,0.678 ,0.68 ,0.682 ,0.684 ,0.686 ,0.688 ,0.69 ,0.692 ,0.694 ,0.696 ,0.698 ,0.7 ,0.702 ,0.704 ,0.706 ,0.708 ,0.71 ,0.712 ,0.714 ,0.716 ,0.718 ,0.72 ,0.722 ,0.724 ,0.726 ,0.728 ,0.73 ,0.732 ,0.734 ,0.736 ,0.738 ,0.74 ,0.742 ,0.744 ,0.746 ,0.748 ,0.75 ,0.752 ,0.754 ,0.756 ,0.758 ,0.76 ,0.762 ,0.764 ,0.766 ,0.768 ,0.77 ,0.772 ,0.774 ,0.776 ,0.778 ,0.78 ,0.782 ,0.784 ,0.786 ,0.788 ,0.79 ,0.792 ,0.794 ,0.796 ,0.798 ,0.8 ,0.802 ,0.804 ,0.806 ,0.808 ,0.81 ,0.812 ,0.814 , 0.816 ,0.818 ,0.82 ,0.822 ,0.824 ,0.826 ,0.828 ,0.83 ,0.832 ,0.834 ,0.836 ,0.838 ,0.84  };
    vector<Float_t> initial_mass={8}; vector<Float_t> final_mass={20};Double_t sigma_range = 3.85;
    
    t->SetBranchAddress("mass", &mass);t->SetBranchAddress("pT",&pT);t->SetBranchAddress("rapidity",&rapidity);t->SetBranchAddress("BDT_score",&xgb_preds);t->SetBranchAddress("issignal",&issignal);
    t1->SetBranchAddress("pT",&MCpT);t1->SetBranchAddress("rapidity",&MCrapidity);t1->SetBranchAddress("mass", &MCmass);t1->SetBranchAddress("BDT_score", &MCxgb_preds);t1->SetBranchAddress("issignal", &MCissignal);
    
    std::string name="canvas";  TCanvas* canvas = CreateCanvas(name);
    std::string name0="canvas0";  TCanvas* canvas0 = CreateCanvas(name0);
    
    for (int l=0;l<cut.size();l++){
      
      vector<Double_t> y_bin_low1, pt_bin_low1, yield, yield_error,  mc_yield, sim_mc_yield, mc_yield_for_bin_vec, pt_y_yield_bdt_0, yield_sigma_method, yield_error_sigma_method;
    
    Double_t bin_size = 0.3;
    float y_bin_low=0.9, y_bin_up =1.2;
    for (int i = 0; i<1; i++)
      {
      y_bin_low = y_bin_low+bin_size;      y_bin_up = y_bin_up+bin_size;
      float pt_bin_low =0.,         pt_bin_up =0.3;
      
      for (int i = 0; i<1; i++)
        {
        pt_bin_low = pt_bin_low+bin_size;  pt_bin_up = pt_bin_up+bin_size;
      
        vector<Double_t> yield_for_sigma, yield_for_sigma_error;
        Double_t mc_yield_for_bin, middle_corrected_yield;
        for (int k=0;k<bins.size();k++){ 
        TH1F* h1 = new TH1F("h1", Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
        h1 -> SetStats (0);
        TH1F* h0 = new TH1F("h0",Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
        Float_t sum = 0; Float_t sum_mc = 0; Float_t sum_bdt_0=0;
        for (int i =0; i < t->GetEntries(); i++)
          {
          t->GetEntry(i);
        //MC original cut less than 2 is applied remove this
          if ( (pT<pt_bin_up) && (pT>pt_bin_low) && (rapidity >y_bin_low) && (rapidity <y_bin_up) && (xgb_preds>cut[l]) &&(issignal<2) ){
            h1->Fill(mass);
            if(issignal==1){
            sum_mc+=1;              
            }
            }
          }
          
          for (int i =0; i < t1->GetEntries(); i++)
          {
          t1->GetEntry(i);
          //MC original cut less than 2 is applied remove this
          if ( (MCpT<pt_bin_up) && (MCpT>pt_bin_low) && (MCrapidity >y_bin_low) && (MCrapidity <y_bin_up) && (MCissignal==1) && (MCxgb_preds>cut[l]))
          {
            //if ((MCxgb_preds>cut[l]))
            //{
            h0->Fill(MCmass);
            sum +=1;
            //}
            /* if (MCxgb_preds>=0)    //lines for ML efficiency calculation, if want these lines remove also the (MCxgb_preds>cut[l]) from the first if condition
            //{
            sum_bdt_0 +=1;
            }*/  
          }
            
           }
          
              int j=0;
            
              if (sum>1){
                j++;
                
                Double_t minMassForFit0 = h0->GetMean()-sigma_range*(h0->GetRMS());         Double_t maxMassForFit0 = h0->GetMean()+sigma_range*(h0->GetRMS());              

                AliHFInvMassFitter* fitter0 = new AliHFInvMassFitter(h0, minMassForFit0, maxMassForFit0,   AliHFInvMassFitter::ETypeOfBkg::kNoBk, types);
                fitter0->SetInitialGaussianMean(1.1156);                                fitter0->SetInitialGaussianSigma(0.05);
                fitter0->SetInitialFrac2Gaus(0.2);                                      fitter0->SetInitialSecondGaussianSigma(0.0009);
                fitter0->SetParticlePdgMass(1.115683);
                //if (pt_bin_up<0.3){fitter0->SetInitialGaussianSigma(0.01);fitter0->SetInitialSecondGaussianSigma(0.00001);fitter0->SetInitialGaussianMean(1.1156);}
                //fitter->SetPolDegreeForBackgroundFit(0);
    
                Bool_t ok = fitter0->MassFitter(kFALSE);
    
                if (ok) {
                 
               TF1 *sig_func0 = fitter0->GetSignalFunc(); sig_func0 = PlotSignal(sig_func0);
                 h0 = PlotHistos(h0);
               
                canvas0->Draw();
                canvas0->Clear();
                TPad* pad0 = CreatePad1();
                pad0 -> Clear();
                pad0 -> SetLogy();
                h0->Draw("pe,X0");
                sig_func0->Draw("same");
                auto latex0 = new TLatex ();
                latex0 -> SetNDC ();
                latex0 -> SetTextSize (0.03);
                latex0 -> DrawLatex (0.5 ,0.55, Form(  "yield = %0.1f#pm %0.1f ; MC %0.2f", fitter0->GetRawYield(), fitter0->GetRawYieldError(), sum   ));
                latex0 -> DrawLatex (0.5 ,0.52, Form(  "Yield/MC = %0.3f#pm%0.3f", (fitter0->GetRawYield())/(sum) ,(  (fitter0->GetRawYield())/ (sum)   )*(TMath::Sqrt( (fitter0->GetRawYieldError()/fitter0->GetRawYield() )*(fitter0->GetRawYieldError()/fitter0->GetRawYield() ) + (TMath::Sqrt(sum)/sum)*(TMath::Sqrt(sum)/sum) ))    ));
                latex0 -> DrawLatex (0.5 ,0.45, Form("#sigma_{1} = %0.4f; #sigma_{2} =  %0.4f ",sig_func0->GetParameter(2),sig_func0->GetParameter(4)));
                latex0 -> DrawLatex (0.5 ,0.4, Form("#mu = %0.4f ", sig_func0->GetParameter(1)));
                latex0 -> DrawLatex (0.5 ,0.35, Form("#chi^{2}_{red} = %0.3f ", fitter0->GetReducedChiSquare()));
                latex0->Draw("same");
                canvas0->cd();
                TPad* pad02 = CreatePad2();
                pad02 -> Clear();
                pad02->cd();
                ratio_plot_hists(h0, sig_func0, bins[k], minMassForFit0, maxMassForFit0);
                
                  

                Double_t minMassForFit = sig_func0->GetParameter(1)- initial_mass[0] * ( sig_func0->GetParameter(2)+ sig_func0->GetParameter(4) );
                Double_t maxMassForFit = sig_func0->GetParameter(1)+ final_mass[0]   * ( sig_func0->GetParameter(2)+ sig_func0->GetParameter(4) );
                
                

                AliHFInvMassFitter* fitter = new AliHFInvMassFitter(h1, minMassForFit, maxMassForFit,  typeb, types);
                fitter->SetInitialGaussianMean(sig_func0->GetParameter(1));             fitter->SetInitialGaussianSigma(sig_func0->GetParameter(2));
                fitter->SetInitialFrac2Gaus(sig_func0->GetParameter(3));                fitter->SetInitialSecondGaussianSigma(sig_func0->GetParameter(4));
//                 fitter->SetFixGaussianMean(sig_func0->GetParameter(1));             fitter->SetFixGaussianSigma(sig_func0->GetParameter(2));
//                 fitter->SetFixFrac2Gaus(sig_func0->GetParameter(3));                fitter->SetFixSecondGaussianSigma(sig_func0->GetParameter(4));
                fitter->SetPolDegreeForBackgroundFit(4); fitter->SetParticlePdgMass(1.115683);
                Int_t sigma_sides=6;
                fitter->SetNSigma4SideBands(sigma_sides);
                
                Bool_t ok1 = fitter->MassFitter(kFALSE);//kTRUE to draw
                if (ok1) {
                                                      
                  TF1 *sig_func = fitter->GetSignalFunc(); sig_func = PlotSignal(sig_func);
                  TF1 *bac_func= fitter->GetBackgroundRecalcFunc(); bac_func = PlotBackground(bac_func);
                  TF1 *M_func = fitter->GetMassFunc();M_func = PlotMass(M_func);
                  TF1 *initial_bac = fitter->GetInitialBackgroundFit();  
                  TF1 *initial_bac_extend = fitter->GetBackgroundFullRangeFunc(); initial_bac_extend->SetLineStyle(3); initial_bac_extend->SetLineColor(kMagenta);
                  
                  //Double_t cov_mat = M_func->GetCovarianceMatrix();
                  h1 = PlotHistos(h1);
                  
                  canvas ->Draw();
                  canvas -> Clear();
                  
                  //fitter->DrawHere(canvas, 3, 2);
                  TPad* pad1 = CreatePad1();
                  pad1 -> Clear();
                  pad1->SetLogy();
                  pad1->cd();
                  h1->Draw("pe,X0");
                  bac_func->Draw("SAME");
                  M_func->Draw("SAME");//ESAME
                  sig_func->Draw("SAME");
                  
                  initial_bac->Draw("same");  initial_bac_extend->Draw("same");
                  
                  TLine* line = new TLine (sig_func0->GetParameter(1)-sigma_sides* fitter0->GetSigma(),0 ,sig_func0->GetParameter(1)-sigma_sides* fitter0->GetSigma() ,60000);  line -> SetLineColor ( kRed );    line -> SetLineWidth(2);
                  TLine* line1 = new TLine (sig_func0->GetParameter(1)+sigma_sides* fitter0->GetSigma(),0 ,sig_func0->GetParameter(1)+sigma_sides* fitter0->GetSigma() ,60000);  line1 -> SetLineColor ( kRed );    line1 -> SetLineWidth(2);
                  line->Draw("same");
                  line1->Draw("same");
                  
                  
                  Double_t sigma_for_yield = 1.5;
                  Double_t sig, errsig, s, errs, b, errb;
                  fitter->Signal(sigma_for_yield,s,errs);
                  fitter->Background(sigma_for_yield,b,errb);
                  fitter->Significance(sigma_for_yield,sig,errsig);
                  

                 TLegend* legend = new TLegend(0.5,0.7,0.8,0.85);
                  legend->SetTextSize(0.023);
                  legend->AddEntry(h1,"#Lambda hyperon","pe,X0");
                  legend->AddEntry(M_func,"Double gaussian+pol4","l");
                  legend->AddEntry(sig_func,"Double gaussian","l");
                  legend->AddEntry(bac_func,"pol4","l");
                  legend->AddEntry(initial_bac,"Initial background fit","l");
                  legend->AddEntry(initial_bac_extend,"Initial background extended","l");
                  legend->Draw("same");   
                  
                  auto latex = new TLatex ();
                  latex -> SetNDC ();
                  latex -> SetTextSize (0.03);
                  latex -> DrawLatex (0.5 ,0.65, Form("#mu = %0.4f; #sigma_{1} = %0.4f; #sigma_{2} = %0.4f", sig_func->GetParameter(1), sig_func->GetParameter(2), sig_func->GetParameter(4) ));
                  latex -> DrawLatex (0.5 ,0.25, Form("yield = %0.1f#pm %0.1f ; MC %0.2f", fitter->GetRawYield(), fitter->GetRawYieldError(), sum_mc ));
                  latex -> DrawLatex (0.5 ,0.2, Form("#color[2]{Yield/MC = %0.3f#pm%0.3f}", (fitter->GetRawYield())/(sum_mc),(  (fitter->GetRawYield())/ (sum_mc)   )*(TMath::Sqrt( (fitter->GetRawYieldError()/fitter->GetRawYield() )*(fitter->GetRawYieldError()/fitter->GetRawYield() ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc) ))  ));
                  latex -> DrawLatex (0.5 ,0.15, Form("yield from signal in %0.1f #sigma = %0.1f#pm %0.1f",sigma_for_yield, s, errs ));
                  latex -> DrawLatex (0.5 ,0.1, Form("#color[2]{Yield/MC = %0.3f#pm%0.3f}", (s)/(sum_mc),(  (s)/(sum_mc)   )*(TMath::Sqrt( (errs/s )*(errs/s ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc) ))  ));
                  latex -> DrawLatex (0.5 ,0.08, Form("#chi^{2}_{red} = %0.3f ", fitter->GetReducedChiSquare()));
                  latex->Draw("same");
                  
                  
                  canvas->cd();
                  TPad* pad2 = CreatePad2();
                  pad2 -> Clear();
                  pad2->cd();

                  ratio_plot_hists(h1, M_func, bins[k], minMassForFit, maxMassForFit);
                  
                  
                  if (j==1){
                    canvas0->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf(",cut[l]));
                    canvas->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                            }
                  else{
                    canvas0->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                      canvas->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                      }
                      
                  
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
                  yield.push_back(fitter->GetRawYield());
                  yield_error.push_back(fitter->GetRawYieldError());
                  significance.push_back(sig);
                  significance_error.push_back(errsig);
                  mc_yield.push_back(sum_mc);
                  pt_y_yield_bdt_0.push_back(sum_bdt_0);
                  raw_yield_only.push_back(fitter->GetRawYield());
                  raw_yield_only_error.push_back(fitter->GetRawYieldError());
                  
                  mu.push_back(sig_func->GetParameter(1));
                  mu_error.push_back(sig_func->GetParError(1));
                  sigma1.push_back(sig_func->GetParameter(2));
                  sigma1_error.push_back(sig_func->GetParError(2));
                  sigma2.push_back(sig_func->GetParameter(4));
                  sigma2_error.push_back(sig_func->GetParError(4));
                  frac2gaus.push_back(sig_func->GetParameter(3));
                  frac2gaus_error.push_back(sig_func->GetParError(3));
                  
                  
                        }
                        }
                else if (!ok) {
                  cout << ".........I am sorry" << endl;
                  fitter0->GetHistoClone()->Draw();
                  y_bin_low1.push_back(y_bin_low);
                  pt_bin_low1.push_back(pt_bin_low);
                  yield.push_back(0.0);
                  yield_error.push_back(0.0);
                  mc_yield.push_back(0.0);
                  sim_mc_yield.push_back(0.0);
                  pt_y_yield_bdt_0.push_back(0);
                              }
                delete fitter0;
      
                         }           
        
        delete h1; 
          delete h0;
        }
        
        }
                    
      }
      
      
      

    TFile *file_sim = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/2M/jul19/dcm_prim_signal_m_200_400_bins_10.root");
    TH2F *Mc_sim = (TH2F*)file_sim->Get("SimParticles_McLambda200_400/SimParticles_rapidity_SimParticles_pT_McLambda200_400");
    
    TFile *file_data = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/2M/jul19/urqmd_prim_signal_m_200_400_bins_10.root");
    TH2F *Mc_sim_data = (TH2F*)file_data->Get("SimParticles_McLambda200_400/SimParticles_rapidity_SimParticles_pT_McLambda200_400");
    
    TH2F* h2d = eff_creator_new(y_bin_low1, pt_bin_low1, yield, yield_error, sim_mc_yield, Mc_sim, cut[l], Mc_sim_data,pt_y_yield_bdt_0 ,0, mc_yield);
    TH2F* h2d1 = eff_creator_sim(y_bin_low1, pt_bin_low1, yield, yield_error, sim_mc_yield, Mc_sim, cut[l], Mc_sim_data,pt_y_yield_bdt_0 ,0, mc_yield);
    TH2F* h2d_sigma_method = eff_creator_new(y_bin_low1, pt_bin_low1, yield_sigma_method, yield_error_sigma_method, sim_mc_yield, Mc_sim, cut[l], Mc_sim_data,pt_y_yield_bdt_0 ,0, mc_yield);
    

    if (h2d->GetBinContent( h2d->FindBin(1.4,0.5))>0){bin_content_vector.push_back(h2d->GetBinContent( h2d->FindBin(1.4,0.5)) );  bin_content_error_vector.push_back(h2d->GetBinError( h2d->FindBin(1.4,0.5))); cut_value.push_back(cut[l]);     bin_content_vector_sim.push_back(h2d1->GetBinContent( h2d1->FindBin(1.4,0.5)) );  bin_content_vector_sim_error.push_back(h2d1->GetBinError( h2d1->FindBin(1.4,0.5)));    bin_content_vector_sigma_method.push_back(h2d_sigma_method->GetBinContent( h2d_sigma_method->FindBin(1.4,0.5)) );  bin_content_error_vector_sigma_method.push_back(h2d_sigma_method->GetBinError( h2d_sigma_method->FindBin(1.4,0.5)));}
    cout<<"BDT is "<< cut[l]<<endl;
    
    if (h2d1->GetBinContent( h2d1->FindBin(1.4,0.5))>0){ }
    
    h2d->SetTitle(Form("corrected yield/simulated yield at BDT%g", cut[l]));
    //h2d->SetName(Form("corrected BDT%g", cut[l])); h2d->GetYaxis()->SetRangeUser (0 ,1.8);  h2d->GetXaxis()->SetRangeUser (0 ,1.8); h2d->GetZaxis()->SetRangeUser (0.6 ,1);
    
    //auto c = new TCanvas(); c->Draw(); h2d->Draw("colz"); //to draw the corrected yield/sim yield
    //TFile *myfile = TFile::Open("hist.root","recreate");
    //h2d->Write();
    //Mc_sim_data->Write();
    
  /*  if (l==0.2)
    {
    TFile *bdt = TFile::Open("BDT.root","recreate");
    h2d->Write();
    }*/
    //if (l!=0.2){
    //TFile *bdt = TFile::Open("BDT.root","update");
    //h2d->Write();}
    
    canvas->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/pT_rapidity_BDT_%0.5f.pdf)",cut[l]));
   
   delete h2d; delete h2d1;
      } 
      


      
    
    Double_t yield_opt, yield_opt_error,yield_sigma_method_opt, yield_opt_sigma_method_error,yield_096, yield_096_error;
    auto hist_corr_yield = corrected_yield_hist(bin_content_vector,(*min_element(bin_content_vector.begin(), bin_content_vector.end())) -18000, (*max_element(bin_content_vector.begin(), bin_content_vector.end()))+5000);
    hist_corr_yield->GetXaxis()->SetTitle("Corrected yield");
    auto hist_corr_yield_sim = corrected_yield_hist(bin_content_vector_sim, (*min_element(bin_content_vector_sim.begin(), bin_content_vector_sim.end())) , (*max_element(bin_content_vector_sim.begin(), bin_content_vector_sim.end()))); hist_corr_yield_sim->SetLineColor(kGreen);
    auto hist_corr_yield_sigma_method = corrected_yield_hist(bin_content_vector_sigma_method, (*min_element(bin_content_vector_sigma_method.begin(), bin_content_vector_sigma_method.end()))-18000 , (*max_element(bin_content_vector_sigma_method.begin(), bin_content_vector_sigma_method.end()))+5000);
    hist_corr_yield_sigma_method->GetXaxis()->SetTitle("Corrected yield");
    
    auto hist_mu_fit = corrected_yield_hist(mu, *min_element(mu.begin(), mu.end()) , *max_element(mu.begin(), mu.end()));
    hist_mu_fit->GetXaxis()->SetTitle("mean");
    auto hist_sigma_fit = corrected_yield_hist(sigma1, *min_element(sigma1.begin(), sigma1.end()) , *max_element(sigma1.begin(), sigma1.end()));
    hist_sigma_fit->GetXaxis()->SetTitle("sigma1");
    auto hist_sigma2_fit = corrected_yield_hist(sigma2, *min_element(sigma2.begin(), sigma2.end()) , *max_element(sigma2.begin(), sigma2.end()));
    hist_sigma2_fit->GetXaxis()->SetTitle("sigma2");
    auto hist_frac2gaus_fit = corrected_yield_hist(frac2gaus, *min_element(frac2gaus.begin(), frac2gaus.end()) , *max_element(frac2gaus.begin(), frac2gaus.end()));    hist_frac2gaus_fit->GetXaxis()->SetTitle("frac2gaus");
    
    
    
    Double_t std_of_hist= hist_corr_yield->GetStdDev();
    auto gme = bdt_vs_corrected_yield(bin_content_vector, cut_value,bin_content_error_vector);
    auto gme_sigma_method = bdt_vs_corrected_yield(bin_content_vector_sigma_method, cut_value,bin_content_error_vector_sigma_method); gme_sigma_method->SetTitle("sigma method");
    auto gme_mc_method = bdt_vs_corrected_yield(bin_content_vector_sim, cut_value,bin_content_vector_sim_error);gme_mc_method->SetTitle("MC method");
    auto gme_significance = bdt_vs_corrected_yield(significance, cut_value,significance_error);gme_significance->SetTitle("significance");
    auto gme_raw_yield_only = bdt_vs_corrected_yield(raw_yield_only, cut_value,raw_yield_only_error);gme_raw_yield_only->SetTitle("raw yield from fit");
    auto gme_mu = bdt_vs_corrected_yield(mu, cut_value,mu_error);gme_mu->SetTitle("mean from fit");
    auto gme_sigma1 = bdt_vs_corrected_yield(sigma1, cut_value,sigma1_error);gme_sigma1->SetTitle("sigma1 from fit");
    auto gme_sigma2 = bdt_vs_corrected_yield(sigma2, cut_value,sigma2_error);gme_sigma2->SetTitle("sigma2 from fit");
    auto gme_frac2gaus = bdt_vs_corrected_yield(frac2gaus, cut_value,frac2gaus_error);gme_frac2gaus->SetTitle("frac2gaus from fit");
    
   for (int i=0;i<=bin_content_vector.size();i++)
   {     
     if( (cut_value[i]>0.5199) && (cut_value[i]<0.5201) ){     yield_opt=bin_content_vector[i]; yield_opt_error=bin_content_error_vector[i];  yield_096=bin_content_vector_sim[i]; yield_096_error=bin_content_vector_sim_error[i]; yield_sigma_method_opt = bin_content_vector_sigma_method[i];yield_opt_sigma_method_error=bin_content_error_vector_sigma_method[i];}
   }

   
/*when using RawYield method uncomment
  TLine* line = new TLine (yield_opt,0 ,yield_opt ,6);  line -> SetLineColor ( kRed );    line -> SetLineWidth(2);
  TLine* line1 = new TLine (yield_opt-yield_opt_error,0 , yield_opt-yield_opt_error,6);    line1 -> SetLineColor ( kRed );line1->SetLineStyle(3);line1 -> SetLineWidth(2);
  TLine* line2 = new TLine (yield_opt+yield_opt_error,0 , yield_opt+yield_opt_error,6);    line2 -> SetLineColor ( kRed );line2->SetLineStyle(3);    line2 -> SetLineWidth(2);
  */
  TLine* line = new TLine (yield_sigma_method_opt,0 ,yield_sigma_method_opt ,6);  line -> SetLineColor ( kRed );    line -> SetLineWidth(6);
  TLine* line1 = new TLine (yield_sigma_method_opt-yield_opt_sigma_method_error,0 , yield_sigma_method_opt-yield_opt_sigma_method_error,6);    line1 -> SetLineColor ( kRed );line1->SetLineStyle(3);line1 -> SetLineWidth(6);
  TLine* line2 = new TLine (yield_sigma_method_opt+yield_opt_sigma_method_error,0 , yield_sigma_method_opt+yield_opt_sigma_method_error,6);    line2 -> SetLineColor ( kRed );line2->SetLineStyle(3);    line2 -> SetLineWidth(6);

    
    TFile *file_data1 = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/2M/jul19/urqmd_prim_signal_m_200_400_bins_10.root");
    TH2F *Mc_sim_data1 = (TH2F*)file_data1->Get("SimParticles_McLambda200_400/SimParticles_rapidity_SimParticles_pT_McLambda200_400");
    Double_t true_yield = Mc_sim_data1->GetBinContent( Mc_sim_data1->FindBin(1.4,0.5));
    TLine* line3 = new TLine (true_yield,0 ,true_yield ,6);  line3 -> SetLineColor ( kBlack );    line3 -> SetLineWidth(6);
    
    TLine* line4 = new TLine (true_yield+sqrt(true_yield),0 , true_yield+sqrt(true_yield),6);    line4 -> SetLineColor ( kBlack );line4->SetLineStyle(3);line4 -> SetLineWidth(6);
    TLine* line7 = new TLine (true_yield-sqrt(true_yield),0 , true_yield-sqrt(true_yield),6);    line7 -> SetLineColor ( kBlack );line7->SetLineStyle(3);line7 -> SetLineWidth(6);
    
      TLine* line5 = new TLine (yield_096,0 , yield_096,6);    line5 -> SetLineColor ( kMagenta );    line5 -> SetLineWidth(6);
  TLine* line6 = new TLine (yield_096-yield_096_error,0 , yield_096-yield_096_error,6);    line6 -> SetLineColor ( kMagenta );line6->SetLineStyle(3);    line6 -> SetLineWidth(6);
    
    TLegend* legend = new TLegend(0.55,0.7,0.9,0.9);
    legend->AddEntry(hist_corr_yield_sim,"#Lambda hyperon corrected yield MC","pe");
  legend->AddEntry(line5,"default corrected yield","l");
  legend->AddEntry(line6,"default yield #pm #sigma_{default yield}","l");
  legend->AddEntry(hist_corr_yield,"#Lambda hyperon corrected yield","pe");
  legend->AddEntry(line,"default corrected yield","l");
  legend->AddEntry(line1,"default yield #pm #sigma_{default yield}","l");
  legend->AddEntry(line3,"True yield","l");
  legend->AddEntry(line4,"True yield + #sigma_{True yield}","l");
  legend->SetTextSize(0.023);

   auto c = new TCanvas();
    c->Draw();
//
    //hist_corr_yield->Draw("pe");
    hist_corr_yield_sigma_method->Draw("pe");
    hist_corr_yield_sim->Draw("pe,same");
    line->Draw("SAME");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  legend->Draw("same");

   //gme->Draw("APS ; Z ; 5 s=0.5");
   c->Update();
        TFile *myfile = TFile::Open("hist.root","recreate");
    gme->Write();
    gme_mc_method->Write();
    gme_raw_yield_only->Write();
    gme_sigma_method->Write();
    gme_significance->Write();
    gme_mu->Write();
    gme_frac2gaus->Write();
    gme_sigma1->Write();
    gme_sigma2->Write();
    hist_corr_yield->Write();
    hist_corr_yield_sigma_method->Write();
    hist_mu_fit->Write();
    hist_sigma_fit->Write();
    hist_sigma2_fit->Write();
    hist_frac2gaus_fit->Write();
    
    cout<<"central yield is"<<yield_sigma_method_opt<<" plus minus "<<yield_opt_sigma_method_error<<endl;
    
    auto c11 = new TCanvas();
    c11->Draw();
    hist_mu_fit->Draw("pe");
    
  }
  
  int main(int argc,char** argv) {
    TString name = argv[1];
    TString name2 = argv[2];
    systematics_approach_2M(name, name2);
    return 0;
  }
          

                

  
  
