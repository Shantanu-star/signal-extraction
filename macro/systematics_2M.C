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
  Int_t types = k2Gaus;
  
  using std::cout;
  using std::endl;
  
  void systematics_2M(TString data,TString sim, TString fNameOut = "") 
    {
    
    
    TFile* f = TFile::Open(data.Data(), "read");    
    TTree* t = f->Get<TTree>("t_BDT"); 
    
    TFile* file_sim = TFile::Open(sim.Data(), "read"); 
    TTree* t1 = file_sim->Get<TTree>("t_BDT");
    
    
    vector<Double_t> red_chi2, y_bin_low1, pt_bin_low1, sigma_1, sigma_1_error, sigma_2, sigma_2_error, sigma_1_2_error, yield, yield_error, mc_yield, sim_mc_yield, mc_yield_for_bin_vec;Float_t mass, pT, rapidity, MCmass, MCpT, MCrapidity, xgb_preds, MCxgb_preds;Int_t issignal, MCissignal; vector<Int_t> bins={250,500,1000};  vector<Float_t> cut={0.78};vector<Float_t> initial_mass={9,8,7}; vector<Float_t> final_mass={19,20,21}; vector<Int_t>back_function{3,4,5}; vector<Float_t> sigma_for_yield_calculation={1,1.5,2};vector<Float_t> sigma_sides={5,6,7};Double_t sigma_range = 3.85;
    
    t->SetBranchAddress("mass", &mass);t->SetBranchAddress("pT",&pT);t->SetBranchAddress("rapidity",&rapidity);t->SetBranchAddress("BDT_score",&xgb_preds);t->SetBranchAddress("issignal",&issignal);
    t1->SetBranchAddress("pT",&MCpT);t1->SetBranchAddress("rapidity",&MCrapidity);t1->SetBranchAddress("mass", &MCmass);t1->SetBranchAddress("BDT_score", &MCxgb_preds);t1->SetBranchAddress("issignal", &MCissignal);
    std::string name="canvas";  TCanvas* canvas = CreateCanvas(name);
    std::string name0="canvas0";  TCanvas* canvas0 = CreateCanvas(name0);
    
    Double_t bin_size = 0.3;
    float y_bin_low=0.9, y_bin_up =1.2;
    for (int i = 0; i<1; i++)
      {
        
      y_bin_low = y_bin_low+bin_size;      y_bin_up = y_bin_up+bin_size;
      float pt_bin_low =0.,         pt_bin_up =0.3;
      
      for (int pt_i = 0; pt_i<1; pt_i++)
        {
          
        pt_bin_low = pt_bin_low+bin_size;  pt_bin_up = pt_bin_up+bin_size;
      
        vector<Double_t> yield_for_sigma, yield_for_sigma_error;
        Double_t mc_yield_for_bin, middle_yield; 
        for (int k=0;k<bins.size();k++){ for (int l=0;l<cut.size();l++){ for (int m=0; m<initial_mass.size();m++){ for (int n=0;n<final_mass.size();n++){for (int nnn=0;nnn<=2;nnn++){for (int mmm=0;mmm<=2;mmm++){ for (int kkk=0;kkk<=2;kkk++){ 
        TH1F* h1 = new TH1F("h1", Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
        h1 -> SetStats (0);
        TH1F* h0 = new TH1F("h0",Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
        Float_t sum = 0; Float_t sum_mc = 0;
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
          if ( (MCpT<pt_bin_up) && (MCpT>pt_bin_low) && (MCrapidity >y_bin_low) && (MCrapidity <y_bin_up) && (MCxgb_preds>cut[l]) &&(MCissignal==1)){
            h0->Fill(MCmass);
            sum +=1;
            }
          }
          
              int j=0;
            
              if (sum>1){
                j++;
                Double_t minMassForFit0 = h0->GetMean()-sigma_range*(h0->GetRMS());         Double_t maxMassForFit0 = h0->GetMean()+sigma_range*(h0->GetRMS());             

                AliHFInvMassFitter* fitter0 = new AliHFInvMassFitter(h0, minMassForFit0, maxMassForFit0,   AliHFInvMassFitter::ETypeOfBkg::kNoBk, types);
                fitter0->SetInitialGaussianMean(1.1156);                                fitter0->SetInitialGaussianSigma(0.001);
                fitter0->SetInitialFrac2Gaus(0.2);                                      fitter0->SetInitialSecondGaussianSigma(0.0009);
                if (pt_bin_up<0.3){fitter0->SetInitialGaussianSigma(0.01);fitter0->SetInitialSecondGaussianSigma(0.00001);fitter0->SetInitialGaussianMean(1.1156);}
                //fitter->SetPolDegreeForBackgroundFit(0);
    
                Bool_t ok = fitter0->MassFitter(kFALSE);
    
                if (ok) {
                TF1 *sig_func0 = fitter0->GetSignalFunc();
                
                Double_t minMassForFit = sig_func0->GetParameter(1)- initial_mass[m] * ( sig_func0->GetParameter(2)+ sig_func0->GetParameter(4) );
                Double_t maxMassForFit = sig_func0->GetParameter(1)+ final_mass[n]   * ( sig_func0->GetParameter(2)+ sig_func0->GetParameter(4) );
                
/*                cout<<"min mass for fit is "<<minMassForFit<<endl;
                cout<<"max mass for fit is "<<maxMassForFit<<endl;
                cout<<"mean is  "<<sig_func0->GetParameter(1)<<endl;
                cout<<"sum of sigma is "<<sig_func0->GetParameter(2)+ sig_func0->GetParameter(4)<<endl;
                cout<<"back pol "<<back_function[nnn]<<endl;
*/

                AliHFInvMassFitter* fitter = new AliHFInvMassFitter(h1, minMassForFit, maxMassForFit,  typeb, types);
                fitter->SetInitialGaussianMean(sig_func0->GetParameter(1));             fitter->SetInitialGaussianSigma(sig_func0->GetParameter(2));
                fitter->SetInitialFrac2Gaus(sig_func0->GetParameter(3));                fitter->SetInitialSecondGaussianSigma(sig_func0->GetParameter(4));
                
                fitter->SetPolDegreeForBackgroundFit(back_function[nnn]); fitter->SetParticlePdgMass(1.115683);
                fitter->SetNSigma4SideBands(sigma_sides[kkk]);
                
                Bool_t ok1 = fitter->MassFitter(kFALSE);//KTRUE to draw
                if (ok1) {
                  
                  TF1 *sig_func = fitter->GetSignalFunc(); sig_func = PlotSignal(sig_func);
                  TF1 *bac_func= fitter->GetBackgroundRecalcFunc(); bac_func = PlotBackground(bac_func);
                  TF1 *M_func = fitter->GetMassFunc();M_func = PlotMass(M_func); 
                  
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
                  
                  Double_t sig, errsig, s, errs, b, errb;
                  fitter->Signal(sigma_for_yield_calculation[mmm],s,errs);
                  fitter->Background(sigma_for_yield_calculation[mmm],b,errb);
                  fitter->Significance(sigma_for_yield_calculation[mmm],sig,errsig);
                  
                  
                  auto latex = new TLatex ();
                  latex -> SetNDC ();
                  latex -> SetTextSize (0.03);
                  latex -> DrawLatex (0.5 ,0.5, Form("m_{min} = #mu - %0.2f#sigma_{1}+#sigma_{2}; m_{max} = #mu+ %0.2f#sigma_{1}+#sigma_{2}", initial_mass[m], final_mass[k] ));
                  latex -> DrawLatex (0.5 ,0.45, Form("bins = %d;  pol %0d", bins[k], back_function[nnn] ));                   
                  latex -> DrawLatex (0.5 ,0.25, Form("yield = %0.1f#pm %0.1f ; MC %0.2f", fitter->GetRawYield(), fitter->GetRawYieldError(), sum_mc ));
                  latex -> DrawLatex (0.5 ,0.2, Form("#color[2]{Yield/MC = %0.3f#pm%0.3f}", (fitter->GetRawYield())/(sum_mc),(  (fitter->GetRawYield())/ (sum_mc)   )*(TMath::Sqrt( (fitter->GetRawYieldError()/fitter->GetRawYield() )*(fitter->GetRawYieldError()/fitter->GetRawYield() ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc) ))  ));
                  latex -> DrawLatex (0.5 ,0.15, Form("yield from signal in 3#sigma = %0.1f#pm %0.1f", s, errs ));
                  latex -> DrawLatex (0.5 ,0.1, Form("#color[2]{Yield/MC = %0.3f#pm%0.3f}", (s)/(sum_mc),(  (s)/(sum_mc)   )*(TMath::Sqrt( (errs/s )*(errs/s ) + (TMath::Sqrt(sum_mc)/sum_mc)*(TMath::Sqrt(sum_mc)/sum_mc) ))  ));
                  latex -> DrawLatex (0.5 ,0.08, Form("#chi^{2}_{red} = %0.3f ", fitter->GetReducedChiSquare()));
                  latex->Draw("same");
                  
                  
                  TLegend* legend = new TLegend(0.5,0.6,0.8,0.8);
                  legend->SetTextSize(0.023);
                  legend->AddEntry(h1,"#Lambda hyperon","pe,X0");
                  legend->AddEntry(M_func,"Double gaussian+pol2","l");
                  legend->AddEntry(sig_func,"Double gaussian","l");
                  legend->AddEntry(bac_func,Form("pol %d",back_function[nnn]),"l");
                  legend->Draw("same");
                  
                  
                  
                  canvas->cd();
                  TPad* pad2 = CreatePad2();
                  pad2 -> Clear();
                  pad2->cd();

                  ratio_plot_hists(h1, M_func, bins[k], minMassForFit, maxMassForFit);
                  
                  
                  
                                    if (j==1){
                    canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf(","pdf");
                            }
                  else{
                      canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf","pdf");
                      }
                    
                  
                  sim_mc_yield.push_back(sum);
                  
                  yield_for_sigma.push_back(s);
                  yield_for_sigma_error.push_back(errs);
                  mc_yield.push_back(sum_mc);
                  if ( (bins[k]==500) && (initial_mass[m]==7) && (final_mass[n]==17) && (back_function[nnn]==4) && (sigma_for_yield_calculation[mmm]==1.5) & (sigma_sides[kkk]==6)){
                  middle_yield = s;
                  
                  }

                        }
                        }
                else if (!ok) {
                  cout << ".........I am sorry" << endl;
                  fitter0->GetHistoClone()->Draw();
                              }
                delete fitter0;
      
                         }           
        
        delete h1; 
          delete h0;
          mc_yield_for_bin=sum;
          
        }
        }
        }
        }
        }
        }


        
        }
        double max = *max_element(yield_for_sigma.begin(), yield_for_sigma.end());
        double min = *min_element(yield_for_sigma.begin(), yield_for_sigma.end());
        TH1F* hist_yield = new TH1F("hist_yield","hist_yield",yield_for_sigma.size(),min,max);
        //hist_yield->GetXaxis()->SetTitle("variation in ");
      for (int i=0;i<yield_for_sigma.size();i++)
      {
        //hist_yield->SetBinContent(i+1,yield_for_sigma[i]);
        //hist_yield->SetBinError(i+1,yield_for_sigma_error[i]);
        hist_yield->Fill(yield_for_sigma[i]);
      }
     TH1F* central_yield = new TH1F("central_yield","central_yield",1,middle_yield - sqrt(middle_yield),middle_yield + sqrt(middle_yield));
     central_yield->SetMarkerColor(kRed);
     central_yield->SetLineColor(kRed);
     for (int i=0;i<1;i++)
     {
     central_yield->Fill(middle_yield);
     }
     cout<<" middle yield is "<<middle_yield<<endl;
     
        TFile *myfile = TFile::Open("raw_yield_systematics.root","recreate");
        hist_yield->Write();
        central_yield->Write();
        /*hist_yield->Rebin();
        TF1* yiel_bdt_fit = new TF1("gaus","gaus");
        yiel_bdt_fit->SetParameters( 5,hist_yield->GetMean(),hist_yield->GetRMS()/2);
        hist_yield->Fit(yiel_bdt_fit,"ELM");
        hist_yield->GetXaxis()->SetTitle("Yield");
        hist_yield->SetYTitle("Counts");*/

        auto c1 = new TCanvas(" c1 ","", 500,500);
        c1->Draw();
        hist_yield->Draw("pe");
        //yiel_bdt_fit->Draw("same");
        central_yield->Draw("pe,same");
      
        
        y_bin_low1.push_back(y_bin_low);
        pt_bin_low1.push_back(pt_bin_low);
        mc_yield_for_bin_vec.push_back(mc_yield_for_bin);
        yield.push_back(hist_yield->GetMean());
        yield_error.push_back(hist_yield->GetRMS());
        
        
      }
      
      
      }
      canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf)","pdf");
    
    
    
  }
  
  int main(int argc,char** argv) {
    TString name = argv[1];
    TString name2 = argv[2];
    systematics_2M(name, name2);
    return 0;
  }
          

                

  
  
