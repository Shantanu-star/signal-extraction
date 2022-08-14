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

 
Double_t DSCB(Double_t *x, Double_t *par)
{
  Double_t mu_shift = 1.1157;
  Double_t factor = par[0];
  Double_t mu = par[1];// mu and sig are the parameters of the gaussians
  Double_t sigma = par[2];
  Double_t a1 = par[3];// a1, n1 the parameters of the left power law tail
  Double_t n1 = par[4];
  Double_t a2 = par[5];// a2, n2 the parameters of the right power law tail
  Double_t n2 = par[6];// a1, a2 >= 0; n1, n2 >= 1
  
  Double_t u = (x[0]-mu_shift-mu)/sigma;
  Double_t A1 = TMath::Power(n1/a1, n1)*TMath::Exp(-a1*a1/2);
  Double_t A2 = TMath::Power(n2/a2, n2)*TMath::Exp(-a2*a2/2);
  Double_t B1 = n1/a1 - a1;
  Double_t B2 = n2/a2 - a2;
  
  if(u<-a1)
    return factor*A1*TMath::Power((B1-u), -n1);
  else if(u>=-a1 && u<a2)
    return factor*TMath::Exp(-u*u/2);
  else if(u>=a2)
    return factor*A2*TMath::Power((B2+u), -n2);
  else
    return -1.;
}

//pol2 function
Double_t pol2(Double_t *x, Double_t *par)
{
   Double_t f2 = par[0]+par[1]*(x[0])+par[2]*(x[0]*x[0]);
  //Double_t f2 = par[0]/(fMaxMass-fMinMass)+par[1]*(x[0]-0.5*(fMaxMass+fMinMass))+par[2]*(x[0]*x[0]-1/3.*(fMaxMass*fMaxMass*fMaxMass-fMinMass*fMinMass*fMinMass)/(fMaxMass-fMinMass)); 
  return f2;
}

//pol4 function
Double_t pol4(Double_t *x, Double_t *par)
{
   Double_t f4 = par[0]+par[1]*(x[0])+par[2]*(x[0]*x[0])+par[3]*(x[0]*x[0]*x[0])+par[4]*(x[0]*x[0]*x[0]*x[0]);
  //Double_t f2 = par[0]/(fMaxMass-fMinMass)+par[1]*(x[0]-0.5*(fMaxMass+fMinMass))+par[2]*(x[0]*x[0]-1/3.*(fMaxMass*fMaxMass*fMaxMass-fMinMass*fMinMass*fMinMass)/(fMaxMass-fMinMass)); 
  return f4;
}
//pol3 function
Double_t pol3(Double_t *x, Double_t *par)
{
   Double_t f3 = par[0]+par[1]*(x[0])+par[2]*(x[0]*x[0])+par[3]*(x[0]*x[0]*x[0]);
  //Double_t f2 = par[0]/(fMaxMass-fMinMass)+par[1]*(x[0]-0.5*(fMaxMass+fMinMass))+par[2]*(x[0]*x[0]-1/3.*(fMaxMass*fMaxMass*fMaxMass-fMinMass*fMinMass*fMinMass)/(fMaxMass-fMinMass)); 
  return f3;
}

//total fit function
Double_t fit_func(Double_t *x, Double_t *par)
{
  return pol4(x,par) + DSCB(x,&par[5]); 
}  
  
  
 
  
  void systematics_approach_2M_crystall_ball(TString data,TString sim, TString fNameOut = "") 
{
    
    
    TFile* f = TFile::Open(data.Data(), "read");    
    TTree* t = f->Get<TTree>("t_BDT"); 
    
    TFile* file_sim = TFile::Open(sim.Data(), "read"); 
    TTree* t1 = file_sim->Get<TTree>("t_BDT");
    
    
    Float_t mass, pT, rapidity, MCmass, MCpT, MCrapidity, xgb_preds, MCxgb_preds;Int_t issignal, MCissignal; vector<Int_t> bins={500};  vector<Float_t> mu,mu_error,sigma1,sigma1_error,sigma2,sigma2_error,frac2gaus,frac2gaus_error, bin_content_vector, bin_content_error_vector, bin_content_vector_sim,bin_content_vector_sim_error,bin_content_vector_sigma_method,bin_content_error_vector_sigma_method, significance,significance_error, raw_yield_only,raw_yield_only_error,cut_value; vector<Float_t> cut={0.02 ,0.022 ,0.024 ,0.026 ,0.028 ,0.03 ,0.032 ,0.034 ,0.036 ,0.038 ,0.04 ,0.042 ,0.044 ,0.046 ,0.048 ,0.05 ,0.052 ,0.054 ,0.056 ,0.058 ,0.06 ,0.062 ,0.064 ,0.066 ,0.068 ,0.07 ,0.072 ,0.074 ,0.076 ,0.078 ,0.08 ,0.082 ,0.084 ,0.086 ,0.088 ,0.09 ,0.092 ,0.094 ,0.096 ,0.098 ,0.1 ,0.102 ,0.104 ,0.106 ,0.108 ,0.11 ,0.112 ,0.114 ,0.116 ,0.118 ,0.12 ,0.122 ,0.124 ,0.126 ,0.128 ,0.13 ,0.132 ,0.134 ,0.136 ,0.138 ,0.14 ,0.142 ,0.144 ,0.146 ,0.148 ,0.15 ,0.152 ,0.154 ,0.156 ,0.158 ,0.16 ,0.162 ,0.164 ,0.166 ,0.168 ,0.17 ,0.172 ,0.174 ,0.176 ,0.178 ,0.18 ,0.182 ,0.184 ,0.186 ,0.188 ,0.19 ,0.192 ,0.194 ,0.196 ,0.198 ,0.2 ,0.202 ,0.204 ,0.206 ,0.208 ,0.21 ,0.212 ,0.214 ,0.216 ,0.218 ,0.22 ,0.222 ,0.224 ,0.226 ,0.228 ,0.23 ,0.232 ,0.234 ,0.236 ,0.238 ,0.24 ,0.242 ,0.244 ,0.246 ,0.248 ,0.25 ,0.252 ,0.254 ,0.256 ,0.258 ,0.26 ,0.262 ,0.264 ,0.266 ,0.268 ,0.27 ,0.272 ,0.274 ,0.276 ,0.278 ,0.28 ,0.282 ,0.284 ,0.286 ,0.288 ,0.29 ,0.292 ,0.294 ,0.296 ,0.298 ,0.3 ,0.302 ,0.304 ,0.306 ,0.308 ,0.31 ,0.312 ,0.314 ,0.316 ,0.318 ,0.32 ,0.322 ,0.324 ,0.326 ,0.328 ,0.33 ,0.332 ,0.334 ,0.336 ,0.338 ,0.34 ,0.342 ,0.344 ,0.346 ,0.348 ,0.35 ,0.352 ,0.354 ,0.356 ,0.358 ,0.36 ,0.362 ,0.364 ,0.366 ,0.368 ,0.37 ,0.372 ,0.374 ,0.376 ,0.378 ,0.38 ,0.382 ,0.384 ,0.386 ,0.388 ,0.39 ,0.392 ,0.394 ,0.396 ,0.398 ,0.4 ,0.402 ,0.404 ,0.406 ,0.408 ,0.41 ,0.412 ,0.414 ,0.416 ,0.418 ,0.42 ,0.422 ,0.424 ,0.426 ,0.428 ,0.43 ,0.432 ,0.434 ,0.436 ,0.438 ,0.44 ,0.442 ,0.444 ,0.446 ,0.448 ,0.45 ,0.452 ,0.454 ,0.456 ,0.458 ,0.46 ,0.462 ,0.464 ,0.466 ,0.468 ,0.47 ,0.472 ,0.474 ,0.476 ,0.478 ,0.48 ,0.482 ,0.484 ,0.486 ,0.488 ,0.49 ,0.492 ,0.494 ,0.496 ,0.498 ,0.5 ,0.502 ,0.504 ,0.506 ,0.508 ,0.51 ,0.512 ,0.514 ,0.516 ,0.518 ,0.52 ,0.522 ,0.524 ,0.526 ,0.528 ,0.532 ,0.534 ,0.536 ,0.538 ,0.54 ,0.542 ,0.544 ,0.546 ,0.548 ,0.55 ,0.552 ,0.554 ,0.556 ,0.558 ,0.56 ,0.562 ,0.564 ,0.566 ,0.568 ,0.57 ,0.572 ,0.574 ,0.576 ,0.578 ,0.58 ,0.582 ,0.584 ,0.586 ,0.588 ,0.59 ,0.592 ,0.594 ,0.596 ,0.598 ,0.6 ,0.602 ,0.604 ,0.606 ,0.608 ,0.61 ,0.612 ,0.614 ,0.616 ,0.618 ,0.62 ,0.622 ,0.624 ,0.626 ,0.628 ,0.63 ,0.632 ,0.634 ,0.636 ,0.638 ,0.64 ,0.642 ,0.644 ,0.646 ,0.648 ,0.65 , 0.652 ,0.654 ,0.656 ,0.658 ,0.66 ,0.662 ,0.664 ,0.666 ,0.668 ,0.67 ,0.672 ,0.674 ,0.676 ,0.678 ,0.68 ,0.682 ,0.684 ,0.686 ,0.688 ,0.69 ,0.692 ,0.694 ,0.696 ,0.698 ,0.7 ,0.702 ,0.704 ,0.706 ,0.708 ,0.71 ,0.712 ,0.714 ,0.716 ,0.718 ,0.72 ,0.722 ,0.724 ,0.726 ,0.728 ,0.73 ,0.732 ,0.734 ,0.736 ,0.738 ,0.74 ,0.742 ,0.744 ,0.746 ,0.748 ,0.75 ,0.752 ,0.754 ,0.756 ,0.758 ,0.76 ,0.762 ,0.764 ,0.766 ,0.768 ,0.77 ,0.772 ,0.774 ,0.776 ,0.778 ,0.78 ,0.782 ,0.784 ,0.786 ,0.788 ,0.79 ,0.792 ,0.794 ,0.796 ,0.798 ,0.8 ,0.802 ,0.804 ,0.806 ,0.808 ,0.81 ,0.812 ,0.814 , 0.816 ,0.818 ,0.82 ,0.822 ,0.824 ,0.826 ,0.828 ,0.83 ,0.832 ,0.834 ,0.836 ,0.838 ,0.84  };
    vector<Float_t> initial_mass={21}; vector<Float_t> final_mass={72};Double_t sigma_for_yield = 2.5;Double_t sigma_range = 3.85;
    
    t->SetBranchAddress("mass", &mass);t->SetBranchAddress("pT",&pT);t->SetBranchAddress("rapidity",&rapidity);t->SetBranchAddress("BDT_score",&xgb_preds);t->SetBranchAddress("issignal",&issignal);
    t1->SetBranchAddress("pT",&MCpT);t1->SetBranchAddress("rapidity",&MCrapidity);t1->SetBranchAddress("mass", &MCmass);t1->SetBranchAddress("BDT_score", &MCxgb_preds);t1->SetBranchAddress("issignal", &MCissignal);
    
    std::string name="canvas";  TCanvas* canvas = CreateCanvas(name);
    std::string name0="canvas0";  TCanvas* canvas0 = CreateCanvas(name0);
    std::string name1="canvas01";  TCanvas* canvas1 = CreateCanvas(name1);
    
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
        //hist for background only outisde the lambda peak from URQMD model
        TH1F* hback = new TH1F("hback", Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
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
             if((mass<1.108)| (mass>1.13)){hback->Fill(mass);}//Fill background part
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

                //Fit DCM signal
                TF1* f0 = new TF1("DSCB", DSCB, minMassForFit0,maxMassForFit0, 7);
                f0->SetParameters(h0->GetBinContent(h0->FindBin(1.1157)),0,0.001,1,1,1,1);
                //f0->SetParameter(0, h0->GetBinContent(h0->FindBin(1.1157)));
                f0->SetParLimits(0, 0, 10*h0->GetBinContent(h0->FindBin(1.1157)));
                f0->SetParLimits(3, 0, 10);
                f0->SetParLimits(4, 1, 100);
                f0->SetParLimits(5, 0, 10);
                f0->SetParLimits(6, 1, 100);
                
                auto fitResult0 = h0->Fit(f0,"R,L,M,E,+,0", "",minMassForFit0,maxMassForFit0);
                //auto covMatrix0 = fitResult0->GetCovarianceMatrix();
                 auto par0 = f0 -> GetParameters();
                 auto par0_errors = f0 -> GetParErrors();
//                 
//                  Double_t min_for_sig0 = (1.1157+par0[1]) - 3*(f0->GetParameter(2));
//                  Double_t max_for_sig0 = (1.1157+par0[1]) + 3*(f0->GetParameter(2));
//  
//                 Double_t fRawYield0= (f0->Integral(min_for_sig0,max_for_sig0)) / h0->GetBinWidth(1);
//                 Double_t fRawYieldErr0 = (f0->IntegralError(min_for_sig0,max_for_sig0, fitResult0->GetParams() , covMatrix0.GetMatrixArray(), 0.0000000001) )/h0->GetBinWidth(1);
//                  Double_t signal0 = f0->Integral(min_for_sig0,max_for_sig0)/h1->GetBinWidth(1);
//                  Double_t errsignal0=(fRawYieldErr0/fRawYield0)*signal0;
                
                
                
                auto sig_func0 = PlotSignal(f0);
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
                //latex0 -> DrawLatex (0.5 ,0.55, Form(  "yield = %0.1f#pm %0.1f ; MC %0.2f", signal0, errsignal0, sum   ));
                //latex0 -> DrawLatex (0.5 ,0.52, Form(  "Yield/MC = %0.3f#pm%0.3f", (signal0)/(sum) ,(  (signal0)/ (sum)   )*(TMath::Sqrt( (errsignal0/signal0 )*(errsignal0/signal0 ) + (TMath::Sqrt(sum)/sum)*(TMath::Sqrt(sum)/sum) ))    ));
                //latex0 -> DrawLatex (0.5 ,0.45, Form("#mu-3*sigma = %0.3f + %0.3f = %0.3f ",sig_func0->GetParameter(1),sig_func0->GetParameter(2), sig_func0->GetParameter(1)-3* sig_func0->GetParameter(2)));
                //latex0 -> DrawLatex (0.5 ,0.4, Form("#mu+3*sigma = %0.3f ", sig_func0->GetParameter(1)+3* sig_func0->GetParameter(2)));
                latex0->Draw("same");
                canvas0->cd();
                TPad* pad02 = CreatePad2();
                pad02 -> Clear();
                pad02->cd();
                ratio_plot_hists(h0, sig_func0, bins[k], minMassForFit0, maxMassForFit0);
                

                
                  

                //back fit only
                Double_t minMassForFit = (1.1157+par0[1])-initial_mass[0]*par0[2];
                Double_t maxMassForFit = (1.1157+par0[1])+final_mass[0]*par0[2];
                Double_t integralHisto=hback->Integral(hback->FindBin(minMassForFit),hback->FindBin(maxMassForFit),"width");
                //Fit URQMD background
                TF1* f11 = new TF1("pol3", pol3, minMassForFit,maxMassForFit,4);
                f11->SetParameters(integralHisto,-10,5,1);
                hback->Fit(f11,"", "RLEMN",minMassForFit,minMassForFit);
                auto par11 = f11 -> GetParameters();
                auto par11_error = f11 -> GetParErrors();
                
                TF1* f1 = new TF1("pol4", pol4, minMassForFit,maxMassForFit,5);
                f1->SetParameters(par11[0],par11[1],par11[2],par11[3],1);
                //f1->SetParErrors(par11_error[0],par11_error[1],par11_error[2],par11_error[3],0.1);
                hback->Fit(f1,"", "RLEMN",minMassForFit,minMassForFit);
                auto par1 = f1 -> GetParameters();
                auto par1_error = f1 -> GetParErrors();
                canvas1->Draw();
                canvas1->Clear();
                TPad* pad11 = CreatePad1();
                pad11 -> Clear();
                pad11 -> SetLogy();
                auto back_pol4_side_bands = PlotBackground(f1);
                hback = PlotHistos(hback);
                hback->Draw("pe,X0");
                f1->Draw("same");
                
                
                canvas ->Draw();
                canvas -> Clear();

                TF1* f2 = new TF1("DSCB+pol2", fit_func, minMassForFit,maxMassForFit, 12);
                //f2->SetParameters(par1[0],par1[1],par1[2],par1[3],h1->GetBinContent(h1->FindBin(1.1157)),par0[1],par0[2],par0[3],par0[4],par0[5],par0[6]);
                
                f2->SetParameter(0,par1[0]);    f2->SetParError(0, par1_error[0]);
                f2->SetParameter(1,par1[1]);    f2->SetParError(1, par1_error[1]);
                f2->SetParameter(2,par1[2]);    f2->SetParError(2, par1_error[2]);
                f2->SetParameter(3,par1[3]);    f2->SetParError(3, par1_error[3]);
                f2->SetParameter(4,par1[4]);    f2->SetParError(4, par1_error[4]);
                f2->SetParameter(5,h1->GetBinContent(h1->FindBin(1.1157)));
                f2->SetParameter(6,par0[1]);    f2->SetParError(6, par0_errors[1]); 
                //f2->SetParLimits(6, par0[1]-par0_errors[1], par0[1]+par0_errors[1]);
                f2->SetParameter(7,par0[2]);    f2->SetParError(7, par0_errors[2]);
                //f2->SetParLimits(7,0.00000001,0.9);
                f2->SetParameter(8,par0[3]);    f2->SetParError(8, par0_errors[3]);
                f2->SetParameter(9,par0[4]);    f2->SetParError(9, par0_errors[4]);
                f2->SetParameter(10,par0[5]);   f2->SetParError(10, par0_errors[5]);
                f2->SetParameter(11,par0[6]);   f2->SetParError(11, par0_errors[6]);
//                 f2->FixParameter(6,par0[1]);
//                 f2->FixParameter(7,par0[2]);
//                 f2->FixParameter(8,par0[3]);
//                 f2->FixParameter(9,par0[4]);
//                 f2->FixParameter(10,par0[5]);
//                 f2->FixParameter(11,par0[6]);
                
                
//                 f2->SetStepSize(0.)

                //if (pt_bin_low<0.6){f2->SetParLimits(3, 0, 10*h1->GetBinContent(h1->FindBin(1.1157)));f2->SetParLimits(6, 0.1, 100);f2->SetParLimits(8, 0.1, 100);}
                
                //f2->SetParLimits(4, 0, 10*h1->GetBinContent(h1->FindBin(1.1157)));
                //f2->SetParLimits(7, 0, 10);
                //f2->SetParLimits(8, 1, 100);
                //f2->SetParLimits(9, 0, 10);
                //f2->SetParLimits(10, 1, 100);
                
                
                auto fitResult = h1->Fit(f2,"R,L,M,E,+,0,S", "",minMassForFit,maxMassForFit);
                //fitResult->Print();
                auto covMatrix = fitResult->GetCovarianceMatrix();
                auto par2 = f2 -> GetParameters();
                auto par2_error = f2 -> GetParErrors();
                
                
                TF1* sig_func = new TF1("fs", DSCB, minMassForFit,maxMassForFit, 7); 
                sig_func->SetParameters(par2[5],par2[6],par2[7],par2[8],par2[9],par2[10],par2[11]);
                TF1* bac_func_total = new TF1("background", pol4, minMassForFit,maxMassForFit, 5); 
                bac_func_total->SetParameters(par2[0],par2[1],par2[2],par2[3],par2[4]);
//                TF1* bac_func_total = new TF1("background", pol3, minMassForFit,maxMassForFit, 3); 
//                bac_func_total->SetParameters(par2[0],par2[1],par2[2],par2[3]);
                
                
                
                 h1 = PlotHistos(h1);
                 f2 = PlotMass(f2);
                 sig_func = PlotSignal(sig_func);
                 bac_func_total = PlotBackground(bac_func_total);


                  TPad* pad1 = CreatePad1();
                  pad1 -> Clear();
                  pad1->cd();
                  //h1->GetYaxis()->SetRangeUser(-100000,10000);
                  pad1->SetLogy();
                  h1->Draw("pe,X0");
                  bac_func->Draw("SAME");
                  M_func->Draw("SAME");//ESAME
                  sig_func->Draw("SAME");
                  
                  
                  Double_t min_for_sig = (1.1157+par2[6]) - sigma_for_yield*(f2->GetParameter(7));
                  cout<<"minimum = "<<min_for_sig<<endl;
                  Double_t max_for_sig = (1.1157+par2[6]) + sigma_for_yield*(f2->GetParameter(7));
                  
                  //float fRawYield = f2->GetParameter(2)/h1->GetBinWidth(1);
                  //float fRawYieldErr = f2->GetParError(2)/h1->GetBinWidth(1);
                  
                  Double_t fRawYield= (f2->Integral(min_for_sig,max_for_sig)) / h1->GetBinWidth(1);
                  double fRawYieldErr = (f2->IntegralError(min_for_sig,max_for_sig, fitResult->GetParams() , covMatrix.GetMatrixArray(), 0.0000000001) )/h1->GetBinWidth(1);
                  Double_t signal = sig_func->Integral(min_for_sig,max_for_sig)/h1->GetBinWidth(1);
                  Double_t errsignal=(fRawYieldErr/fRawYield)*signal;/*assume relative error is the same as for total integral*/
                  //Double_t errsignal= sig_func->IntegralError(minMassForFit,maxMassForFit)/h1->GetBinWidth(1);
                  
                  auto latex = new TLatex ();
                  latex -> SetNDC ();
                  latex -> SetTextSize (0.03);
                  latex -> DrawLatex (0.48 ,0.85, Form("s = %4.5f #pm %4.5f; mc = %4.3f", signal, errsignal, sum_mc));
                   latex -> DrawLatex (0.48 ,0.8, Form("#color[2]{Yield/MC = %0.3f#pm%0.3f}", (signal)/(sum_mc),(  (signal)/(sum_mc)   )*(TMath::Sqrt( (errsignal/signal )*(errsignal/signal ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc) ))  ));
                  latex -> DrawLatex (0.2 ,0.75, Form("#chi_{red}^{2} = %4.3f",f2->GetChisquare()/f2->GetNDF()));
                  latex -> DrawLatex (0.48 ,0.7, Form("m_{0}: = %4.5f #pm %4.3f", 1.1157+par2[6], f2->GetParError(6) ));
                  latex -> DrawLatex (0.48 ,0.65, Form("#sigma: = %4.5f #pm %4.5f", f2->GetParameter(7), f2->GetParError(7) ));
                  
                  //latex -> DrawLatex (0.48 ,0.5, Form("raw yield = %4.5f #pm %4.5f", fRawYield, fRawYieldErr));
                  latex -> DrawLatex (0.48 ,0.55, "CBM Performance");
                  latex -> DrawLatex (0.48 ,0.5, "URQMD, AuAu @ 12 #it{A}GeV/#it{c}");
                  latex -> Draw("SAME");
                  
                  TLine* line = new TLine (min_for_sig,0 ,min_for_sig ,60000);  line -> SetLineColor ( kRed );    line -> SetLineWidth(2);
                  TLine* line1 = new TLine (max_for_sig,0 ,max_for_sig ,60000);  line1 -> SetLineColor ( kRed );    line1 -> SetLineWidth(2);
                  line->Draw("same");
                  line1->Draw("same");
                  
                  
                  
                  
                  
                  canvas->cd();
                  TPad* pad2 = CreatePad2();
                  pad2 -> Clear();
                  pad2->cd();
                ratio_plot_hists(h1, f2, bins[k], minMassForFit, maxMassForFit);
                
                
                
                
                
                
                if (j==1){
                    canvas0->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/dscb/pT_rapidity_BDT_%0.5f.pdf(",cut[l]));
                    canvas->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/dscb/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                            }
                  else{
                    canvas0->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/dscb/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                      canvas->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/dscb/pT_rapidity_BDT_%0.5f.pdf",cut[l]));
                      }
                
                  y_bin_low1.push_back(y_bin_low);
                  pt_bin_low1.push_back(pt_bin_low);
                  sim_mc_yield.push_back(sum);
                  pt_y_yield_bdt_0.push_back(sum_bdt_0);
                  

                  yield_sigma_method.push_back(signal);
                  yield_error_sigma_method.push_back(errsignal);
                  mc_yield.push_back(sum_mc);
                  raw_yield_only.push_back(signal);
                  raw_yield_only_error.push_back(errsignal);
                  
                  mu.push_back((1.1157+par2[6]));
                  mu_error.push_back(par2_error[6]);
                  sigma1.push_back(par2[7]);
                  sigma1_error.push_back(par2_error[7]);
                
                

  }
        }
        }
      }
    
      
    TFile *file_sim = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/2M/jul19/dcm_prim_signal_m_200_400_bins_10.root");
    TH2F *Mc_sim = (TH2F*)file_sim->Get("SimParticles_McLambda200_400/SimParticles_rapidity_SimParticles_pT_McLambda200_400");
    
    TFile *file_data = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/2M/jul19/urqmd_prim_signal_m_200_400_bins_10.root");
    TH2F *Mc_sim_data = (TH2F*)file_data->Get("SimParticles_McLambda200_400/SimParticles_rapidity_SimParticles_pT_McLambda200_400");
    
    TH2F* h2d1 = eff_creator_sim(y_bin_low1, pt_bin_low1, yield_sigma_method, yield_error_sigma_method, sim_mc_yield, Mc_sim, cut[l], Mc_sim_data,pt_y_yield_bdt_0 ,0, mc_yield);
    TH2F* h2d_sigma_method = eff_creator_new(y_bin_low1, pt_bin_low1, yield_sigma_method, yield_error_sigma_method, sim_mc_yield, Mc_sim, cut[l], Mc_sim_data,pt_y_yield_bdt_0 ,0, mc_yield);
    

    if (h2d_sigma_method->GetBinContent( h2d_sigma_method->FindBin(1.4,0.5))>0){ cut_value.push_back(cut[l]);     bin_content_vector_sim.push_back(h2d1->GetBinContent( h2d1->FindBin(1.4,0.5)) );  bin_content_vector_sim_error.push_back(h2d1->GetBinError( h2d1->FindBin(1.4,0.5)));    bin_content_vector_sigma_method.push_back(h2d_sigma_method->GetBinContent( h2d_sigma_method->FindBin(1.4,0.5)) );  bin_content_error_vector_sigma_method.push_back(h2d_sigma_method->GetBinError( h2d_sigma_method->FindBin(1.4,0.5)));}
    cout<<"BDT is "<< cut[l]<<endl;
    
    if (h2d1->GetBinContent( h2d1->FindBin(1.4,0.5))>0){ }
    
    h2d_sigma_method->SetTitle(Form("corrected yield/simulated yield at BDT%g", cut[l]));
    

    canvas->Print (Form("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/fit_plots/dscb/pT_rapidity_BDT_%0.5f.pdf)",cut[l]));
    }
    
    
    Double_t  yield_sigma_method_opt, yield_opt_sigma_method_error,yield_096, yield_096_error;
    
    auto hist_corr_yield_sim = corrected_yield_hist(bin_content_vector_sim, (*min_element(bin_content_vector_sim.begin(), bin_content_vector_sim.end())) , (*max_element(bin_content_vector_sim.begin(), bin_content_vector_sim.end()))); hist_corr_yield_sim->SetLineColor(kGreen);
    auto hist_corr_yield_sigma_method = corrected_yield_hist(bin_content_vector_sigma_method, (*min_element(bin_content_vector_sigma_method.begin(), bin_content_vector_sigma_method.end()))-18000 , (*max_element(bin_content_vector_sigma_method.begin(), bin_content_vector_sigma_method.end()))+5000);
    hist_corr_yield_sigma_method->GetXaxis()->SetTitle("Corrected yield");
    
    auto hist_mu_fit = corrected_yield_hist(mu, *min_element(mu.begin(), mu.end()) , *max_element(mu.begin(), mu.end()));
    hist_mu_fit->GetXaxis()->SetTitle("mean");
    auto hist_sigma_fit = corrected_yield_hist(sigma1, *min_element(sigma1.begin(), sigma1.end()) , *max_element(sigma1.begin(), sigma1.end()));
    hist_sigma_fit->GetXaxis()->SetTitle("sigma1");
    
    
    
    Double_t std_of_hist= hist_corr_yield_sigma_method->GetStdDev();
    auto gme_sigma_method = bdt_vs_corrected_yield(bin_content_vector_sigma_method, cut_value,bin_content_error_vector_sigma_method); gme_sigma_method->SetTitle("sigma method");
    auto gme_mc_method = bdt_vs_corrected_yield(bin_content_vector_sim, cut_value,bin_content_vector_sim_error);gme_mc_method->SetTitle("MC method");
    auto gme_raw_yield_only = bdt_vs_corrected_yield(raw_yield_only, cut_value,raw_yield_only_error);gme_raw_yield_only->SetTitle("raw yield from fit");
    auto gme_mu = bdt_vs_corrected_yield(mu, cut_value,mu_error);gme_mu->SetTitle("mean from fit");
    auto gme_sigma1 = bdt_vs_corrected_yield(sigma1, cut_value,sigma1_error);gme_sigma1->SetTitle("sigma1 from fit");
    
   for (int i=0;i<=bin_content_vector_sigma_method.size();i++)
   {     
     if( (cut_value[i]>0.5199) && (cut_value[i]<0.5201) ){ yield_096=bin_content_vector_sim[i]; yield_096_error=bin_content_vector_sim_error[i]; yield_sigma_method_opt = bin_content_vector_sigma_method[i];yield_opt_sigma_method_error=bin_content_error_vector_sigma_method[i];}
   }

   

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
  legend->AddEntry(hist_corr_yield_sigma_method,"#Lambda hyperon corrected yield","pe");
  legend->AddEntry(line,"default corrected yield","l");
  legend->AddEntry(line1,"default yield #pm #sigma_{default yield}","l");
  legend->AddEntry(line3,"True yield","l");
  legend->AddEntry(line4,"True yield + #sigma_{True yield}","l");
  legend->SetTextSize(0.023);

   auto c = new TCanvas();
    c->Draw();

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
        TFile *myfile = TFile::Open("corrected_yield_dscb.root","recreate");
    gme_mc_method->Write();
    gme_raw_yield_only->Write();
    gme_sigma_method->Write();
    gme_mu->Write();
    gme_sigma1->Write();

    hist_corr_yield_sigma_method->Write();
    hist_mu_fit->Write();
    hist_sigma_fit->Write();

    
    cout<<"central yield is"<<yield_sigma_method_opt<<" plus minus "<<yield_opt_sigma_method_error<<endl;
    
    auto c11 = new TCanvas();
    c11->Draw();
    hist_mu_fit->Draw("pe");  
    
    
}
  int main(int argc,char** argv) {
    TString name = argv[1];
    TString name2 = argv[2];
    systematics_approach_2M_crystall_ball(name, name2);
    return 0;
  }
          

                

  
  
