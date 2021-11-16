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
//using namespace std;

// #include "AliHFInvMassFitter.h"
// #include "AliHFMassFitter.h"
#endif

enum { kGaus = 0,
  k2Gaus=1, kPol2=2, kNoBk=3, kPol3=6};
  

  Int_t typeb = kNoBk;
  Int_t types = k2Gaus;
  
  using std::cout;
  using std::endl;
  
  vector<Double_t> FitXicZerotoXiPiInvMass2(TString fname, TString fNameOut = "") 
    {
    
    //open files
    TFile* f = TFile::Open(fname.Data(), "read");    
    TTree* t = f->Get<TTree>("t1"); TTree* t1 = f->Get<TTree>("t2");
    TFile *f1 = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/lambda_qa_urqmd.root");
    TDirectory* dir = (TDirectory*)f1->Get("SimParticles_McLambda");
    TH2D *mc_spectra = (TH2D*)dir->Get("SimParticles_rapidity_SimParticles_pT_McLambda");
    TFile* f2 = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/new_dcm_100_efficiency_pt_y_yield_bdt_cut_0.8.root");
    TH2D *efficiency = (TH2D*)f2->Get("Efficiency");
    TFile *f3 = TFile::Open("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/lambda_qa_dcm.root");
    TDirectory* dir1 = (TDirectory*)f3->Get("SimParticles_McLambda");
    TH2D *Mc = (TH2D*)dir1->Get("SimParticles_rapidity_SimParticles_pT_McLambda");
    Mc->SetName("Mc"); Mc->SetTitle("Mc");
    
    
    vector<Double_t> red_chi2, y_bin_low1, pt_bin_low1, sigma_1, sigma_1_error, sigma_2, sigma_2_error, sigma_1_2_error, yield, yield_error, mc_yield, sim_mc_yield;float mass, pT, rapidity, MCmass, MCpT, MCrapidity;Double_t xgb_preds, MCxgb_preds, issignal, MCissignal;
    
    t->SetBranchAddress("mass", &mass);t->SetBranchAddress("pT",&pT);t->SetBranchAddress("rapidity",&rapidity);t->SetBranchAddress("xgb_preds",&xgb_preds);t->SetBranchAddress("issignal",&issignal);
    t1->SetBranchAddress("MCpT",&MCpT);t1->SetBranchAddress("MCrapidity",&MCrapidity);t1->SetBranchAddress("MCmass", &MCmass);t1->SetBranchAddress("MCxgb_preds", &MCxgb_preds);t1->SetBranchAddress("MCissignal", &MCissignal);
    TCanvas* canvas = CreateCanvas();
    float y_bin_low=-0.2, y_bin_up =0;
    for (int i = 0; i<15; i++)
      {
      y_bin_low = y_bin_low+0.2;
      y_bin_up = y_bin_up+0.2;
      float pt_bin_low =-0.2, pt_bin_up =0;
      for (int i = 0; i<15; i++)
        {
        pt_bin_low = pt_bin_low+0.2;
        pt_bin_up = pt_bin_up+0.2;
        
        TH1F* h1 = new TH1F("h1", Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), 200, 1.08, 1.2);
        h1 -> SetStats (0);
        TH1F* h0 = new TH1F("h0",Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), 200, 1.08, 1.2);
        Float_t sum = 0; Float_t sum_mc = 0; Float_t cut = 0.9;
        for (int i =0; i < t->GetEntries(); i++)
          {
          t->GetEntry(i);
          if ( (pT<pt_bin_up) && (pT>pt_bin_low) && (rapidity >y_bin_low) && (rapidity <y_bin_up) && (xgb_preds>cut) ){
            //fill h1 with lambda candidates
            h1->Fill(mass);
            if(issignal>0){
            sum_mc+=1;              
            }
            }
          }
          
          for (int i =0; i < t1->GetEntries(); i++)
          {
          t1->GetEntry(i);
          if ( (MCpT<pt_bin_up) && (MCpT>pt_bin_low) && (MCrapidity >y_bin_low) && (MCrapidity <y_bin_up) && (MCxgb_preds>cut) & (MCissignal>0) ){
            //fill h0 with signal only candidates
            h0->Fill(MCmass);
            sum +=1;
            }
          }
          
              int j=0;
            
              if (sum>1){
                j++;
                Double_t minMassForFit0 = h0->GetMean()-4*(h0->GetRMS());
                Double_t maxMassForFit0 = h0->GetMean()+4*(h0->GetRMS());              
//signal only fitter creation
                AliHFInvMassFitter* fitter0 = new AliHFInvMassFitter(h0, minMassForFit0, maxMassForFit0,   AliHFInvMassFitter::ETypeOfBkg::kNoBk, types);
                fitter0->SetInitialGaussianMean(1.1156);
                fitter0->SetInitialGaussianSigma(0.001);
                fitter0->SetInitialFrac2Gaus(0.2);
                fitter0->SetInitialSecondGaussianSigma(0.0009);
                //fitter->SetPolDegreeForBackgroundFit(0);
    
                Bool_t ok = fitter0->MassFitter(kFALSE);
    
                if (ok) {
                Double_t minMassForFit = 1.09;
                Double_t maxMassForFit = 1.2; 
                
                TF1 *sig_func0 = fitter0->GetSignalFunc();
//lambda candidates fitter creation
                AliHFInvMassFitter* fitter = new AliHFInvMassFitter(h1, minMassForFit, maxMassForFit,   AliHFInvMassFitter::ETypeOfBkg::kPol2, types);
                fitter->SetInitialGaussianMean(sig_func0->GetParameter(1));
                fitter->SetInitialGaussianSigma(sig_func0->GetParameter(2));
                fitter->SetInitialFrac2Gaus(sig_func0->GetParameter(3));
                fitter->SetInitialSecondGaussianSigma(sig_func0->GetParameter(4));
                //fitter->SetPolDegreeForBackgroundFit(3);
                
                Bool_t ok1 = fitter->MassFitter(kFALSE);//KTRUE to draw
                if (ok1) {
                  
                  
                  TF1 *sig_func = fitter->GetSignalFunc(); sig_func = PlotSignal(sig_func);
                  TF1 *bac_func= fitter->GetBackgroundRecalcFunc(); bac_func = PlotBackground(bac_func);
                  TF1 *M_func = fitter->GetMassFunc();M_func = PlotMass(M_func);
                  
                  //plotting
                  h1 = PlotHistos(h1);
                  
                  canvas ->Draw();
                  canvas -> Clear();
                  
                  //fitter->DrawHere(canvas, 3, 2);
                  TPad* pad1 = CreatePad1();
                  pad1 -> Clear();
                  pad1->cd();
                  h1->Draw("pe,X0");
                  bac_func->Draw("SAME");
                  M_func->Draw("SAME");//ESAME
                  sig_func->Draw("SAME");
                  

                  TLegend* legend = legend_plot(h1, M_func, sig_func, bac_func);
                  
                  Float_t sigma1, sigma1e, sigma2,sigma2e;
                  if (sig_func->GetParameter(2) < sig_func->GetParameter(4)){sigma_1.push_back(sig_func->GetParameter(2)); sigma_1_error.push_back(sig_func->GetParError(2)); sigma1= sig_func->GetParameter(2); sigma1e=sig_func->GetParError(2); }
                  else if (sig_func->GetParameter(4) < sig_func->GetParameter(2)){sigma_1.push_back(sig_func->GetParameter(4));sigma_1_error.push_back(sig_func->GetParError(4)); sigma1= sig_func->GetParameter(4); sigma1e=sig_func->GetParError(4);}
                  if(sig_func->GetParameter(4) > sig_func->GetParameter(2)){sigma_2.push_back(sig_func->GetParameter(4));sigma_2_error.push_back(sig_func->GetParError(4)); sigma2= sig_func->GetParameter(4); sigma2e=sig_func->GetParError(4);}
                  if (sig_func->GetParameter(2) > sig_func->GetParameter(4)){sigma_2.push_back(sig_func->GetParameter(2));sigma_2_error.push_back(sig_func->GetParError(2)); sigma2= sig_func->GetParameter(4); sigma2e=sig_func->GetParError(4);}
                  
                  /*TF1 *gaus1 = new TF1("gaus1",gaus1_func,minMassForFit,maxMassForFit,4);
                  gaus1->SetParameters(sig_func->GetParameter(0), sig_func->GetParameter(1), sig_func->GetParameter(2),sig_func->GetParameter(3));
                  
                  TF1 *gaus2 = new TF1("gaus2",gaus2_func,minMassForFit,maxMassForFit,4);
                  gaus1->SetParameters(sig_func->GetParameter(0), sig_func->GetParameter(1), sig_func->GetParameter(3), sig_func->GetParameter(4));
                  
                  TF1 *double_gaus1 = new TF1("double_gaus",double_gaus,minMassForFit,maxMassForFit,5);
                  double_gaus1->SetParameters(sig_func->GetParameter(0), sig_func->GetParameter(1), sig_func->GetParameter(2),sig_func->GetParameter(3),sig_func->GetParameter(4));
                  gaus2->Draw();
                  */
                  auto latex = new TLatex ();
                  latex -> SetNDC ();
                  latex -> SetTextSize (0.03);
                  //latex -> DrawLatex (0.5 ,0.3, Form("#sigma_{1} = %4.5f, sigma_{2} = %4.5f, mean = %4.5f",f2->GetParameter(1), f2->GetParameter(4),f2->GetParameter(2)));
                 // latex -> DrawLatex (0.2 ,0.8, Form("#chi_{red}^{2} = %4.5f",fitter0->GetReducedChiSquare()));
                  latex -> DrawLatex (0.48 ,0.8, Form("m_{0}: = %4.5f #pm %4.5f", sig_func->GetParameter(1), sig_func->GetParError(1) ));
                  latex -> DrawLatex (0.48 ,0.75, Form("#sigma_{1}: = %4.5f #pm %4.5f; #sigma_{2}: = %4.5f #pm %4.5f", sigma1, sigma1e,sigma2,sigma2e));
                  //latex -> DrawLatex (0.18 ,0.6, Form("#sigma_{2}: = %4.5f #pm %4.5f", sig_func->GetParameter(4), sig_func->GetParError(4)));
                  latex -> DrawLatex (0.48 ,0.85, Form("s = %4.5f #pm %4.5f; mc = %4.5f", fitter->GetRawYield(), fitter->GetRawYieldError(), sum_mc));
                  latex -> Draw();

                  

                  canvas->cd();
                  TPad* pad2 = CreatePad2();
                  pad2 -> Clear();
                  pad2->cd();
                  TH1F *func_hist = new TH1F("func_hist", "", 200, 1.08, 1.2);
                  TH1F* diff_hist = new TH1F("diff_hist", "",200, 1.08, 1.2);
                  for (int i =h1->FindBin(minMassForFit); i<h1->FindBin(maxMassForFit); i++){
                    Double_t f_value= M_func->Eval(h1->GetBinCenter(i));
                    Double_t t_value = h1->GetBinContent(i);
                    func_hist->SetBinContent(i,f_value);
                    if (h1->GetBinError(i) > 0){
                      diff_hist->SetBinContent(i,(t_value-f_value)/h1->GetBinError(i));}
                  }
                  diff_hist = dif_hist(diff_hist);
                  diff_hist->Draw();
                  TLine* line = new TLine (1.08,0 ,1.2 ,0);
                  line -> SetLineColor ( kRed );
                  line->Draw("SAME");

                  //save some of the data
                  y_bin_low1.push_back(y_bin_low);
                  pt_bin_low1.push_back(pt_bin_low);
                  red_chi2.push_back(fitter->GetReducedChiSquare());
                  sim_mc_yield.push_back(sum);

                  
                  
                  yield.push_back(fitter->GetRawYield());
                  yield_error.push_back(fitter->GetRawYieldError());
                  mc_yield.push_back(sum_mc);
                  
                  //sigma_1_2_error.push_back((sig_func->GetParError(2)*sig_func->GetParError(2))+(sig_func->GetParError(4)*sig_func->GetParError(4)+2*));
  
                  if (j==1){
                    canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf(","pdf");
                            }
                  else{
                      canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf","pdf");
                      }
                        }
                        }
                else if (!ok) {
                  cout << ".........I am sorry" << endl;
                  fitter0->GetHistoClone()->Draw();
                  red_chi2.push_back(0.0);
                  y_bin_low1.push_back(y_bin_low);
                  pt_bin_low1.push_back(pt_bin_low);
                  sigma_1.push_back(0.0);
                  sigma_1_error.push_back(0.0);
                  sigma_2.push_back(0.0);
                  sigma_2_error.push_back(0.0);
                  yield.push_back(0.0);
                  yield_error.push_back(0.0);
                  mc_yield.push_back(0.0);
                  sim_mc_yield.push_back(0.0);
                              }
                delete fitter0;
      
                         }           
        
        delete h1;
          delete h0;
        }
      }
      
//this part is just double checking everything
    canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf)","pdf");
    TH2F* h4 = hist2d_pt_y();
    h4->SetTitle("#sigma_{1}+#sigma_{2}");
    TH2F* h6 = hist2d_pt_y();
    h6->SetTitle("#sigma_{2}");
    TH2F* h7 = hist2d_pt_y();
    h7->SetTitle("#frac{#Lambda}{#epsilon_{DCM}} x #frac{1}{#Lambda_{sim}}");
    h7->GetYaxis()->SetRangeUser (0 ,3);
    h7->GetXaxis()->SetRangeUser (0 ,3);
    //h4->SetTitle("#chi^{2}_{red}");
    //h4->SetTitle("yield - MC/ #sigma");
    //h4->GetZaxis()->SetRangeUser (-2 ,2);
    TH1F* h5 = new TH1F("h1","",300,0,350);
    TH2F* recons = new TH2F("recons", "recons", 15,0,3,15,0,3);
    for (int i=0;i<225;i++){
    Double_t y= y_bin_low1[i];
    Double_t pT=pt_bin_low1[i];
    Double_t y_bin = int((y+0.1)/0.2 + 1);
    Double_t pT_bin = int((pT+0.1)/0.2 + 1);
    //h4->SetBinContent(y_bin, pT_bin, red_chi2[i]);
    if (sim_mc_yield[i]>0){recons->SetBinContent(y_bin,pT_bin,sim_mc_yield[i]);}
    
    if (TMath::Abs(sigma_1_error[i])>0){h4->SetBinContent(y_bin, pT_bin, sigma_1[i]+sigma_2[i]);}
    if (TMath::Abs(sigma_2_error[i])>0){h6->SetBinContent(y_bin, pT_bin, sigma_2[i]);}
    h7->SetBinContent(y_bin, pT_bin, yield[i]);
    h7->SetBinError(y_bin, pT_bin,yield_error[i]);
    if(((TMath::Abs(yield_error[i]))>0)  && (mc_yield[i]>0) && (yield[i]>0)){
      //h4->SetBinContent(y_bin, pT_bin, (yield[i]-mc_yield[i])/TMath::Abs(yield_error[i]));
      h5->SetBinContent(int((i+0.1)/0.2 + 1),yield[i]/mc_yield[i]);
      h5->SetBinError(int((i+0.1)/0.2 + 1),yield_error[i]/yield[i]);
      //cout<<"i is "<<i<<endl;
      //cout<<"bin is "<<int((i+0.1)/0.2 + 1)<<endl;
    }
    
    }
    auto c0 = new TCanvas(" c0 ","", 500,500);
    c0->Draw();
    h5->SetStats(0);
    h5->GetYaxis()->SetRangeUser (0.8 ,1.1);
    h5->Draw("pe");
    TLine* line = new TLine (0,1 ,350 ,1);
    line -> SetLineColor ( kRed );
    line->Draw("SAME");
    auto c = new TCanvas(" c ","", 500,500);
    c->Draw();
    h4->Draw("colz1");
    gStyle->SetPaintTextFormat("4.4f");
    h4->Draw("TEXT SAME");
    auto c6 = new TCanvas(" c6 ","", 500,500);
    c6->Draw();
    //h6->Draw("colz1");
    //h6->Draw("TEXT SAME");
    
    auto c7 = new TCanvas(" c7 ","", 500,500);
    gStyle->SetPaintTextFormat("4.2f");
    c7->Draw();
    TH1* hc = (TH1*)h7->Clone();
    hc->SetName("efficiency_dcm");
    auto divide = hc->Divide(efficiency);
        for (int i=0;i<225;i++){Double_t y= y_bin_low1[i]; Double_t pT=pt_bin_low1[i]; Double_t y_bin = int((y+0.1)/0.2 + 1);  Double_t pT_bin = int((pT+0.1)/0.2 + 1);
    h7->SetBinError(y_bin, pT_bin,yield_error[i]/yield[i]);}
    auto divide_eff = hc->Divide(mc_spectra);
    auto eff_dcm_ratio = h7->Divide(hc);
    auto eff_corr_sim = h7->Divide(mc_spectra);
    //h7->GetZaxis()->SetRangeUser (0 ,1.2);
    //h7->Draw("colz1");
    //h7->Draw("TEXT SAME");
    
    //this function is for creating efficiency plot
    eff_creator(y_bin_low1, pt_bin_low1, yield, yield_error, sim_mc_yield);
    
    auto c8 = new TCanvas("c8","",500,500);
    c8->Draw();
    auto ratio_sim_mc = recons->Divide(Mc);
    hc->Draw("colz1");
    hc->Draw("TEXT SAME");

    
    
    c->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal_hist.png","png");
    c6->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal_hist1.png","png");
    c7->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal_hist2.png","png");
    return red_chi2; 
  }
  
  int main(int argc,char** argv) {
    TString name = argv[1];
    vector<Double_t> red_chi2 = FitXicZerotoXiPiInvMass2(name);
    return 0;
  }
          

                

  
  
