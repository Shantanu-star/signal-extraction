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
#endif

// Double sided crystal ball function
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

//total fit function
Double_t fit_func(Double_t *x, Double_t *par)
{
  return pol2(x,par) + DSCB(x,&par[3]); 
}

  //Fitting starts here
  void crystal_ball(TString fname, TString fNameOut = "") 
    {
    
    //retrieve the tree and branches
    TFile* f = TFile::Open(fname.Data(), "read");    
    TTree* t = f->Get<TTree>("t1"); TTree* t1 = f->Get<TTree>("t2");
    
    vector<Double_t> red_chi2, y_bin_low1, pt_bin_low1, sigma_1, sigma_1_error, sigma_2, sigma_2_error, sigma_1_2_error, yield, yield_error, mc_yield, sim_mc_yield;
    float mass, pT, rapidity, MCmass, MCpT, MCrapidity;Double_t xgb_preds, MCxgb_preds, issignal, MCissignal, original_issignal, MCoriginal_issignal;vector<Int_t> bins={220}; vector<Float_t> cut={0.95};
    
    t->SetBranchAddress("mass", &mass);t->SetBranchAddress("pT",&pT);t->SetBranchAddress("rapidity",&rapidity);t->SetBranchAddress("xgb_preds",&xgb_preds);t->SetBranchAddress("issignal",&issignal);t->SetBranchAddress("original_issignal",&original_issignal);
    t1->SetBranchAddress("MCpT",&MCpT);t1->SetBranchAddress("MCrapidity",&MCrapidity);t1->SetBranchAddress("MCmass", &MCmass);t1->SetBranchAddress("MCxgb_preds", &MCxgb_preds);t1->SetBranchAddress("MCissignal", &MCissignal);t1->SetBranchAddress("MCoriginal_issignal", &MCoriginal_issignal);
    TCanvas* canvas = CreateCanvas();
    
    //start multi-differential fitting
    float y_bin_low=-0.3, y_bin_up =0;
    for (int i = 0; i<10; i++)
      {
      y_bin_low = y_bin_low+0.3;
      y_bin_up = y_bin_up+0.3;
      float pt_bin_low =-0.3, pt_bin_up =0;
      for (int i = 0; i<10; i++)
        {
        pt_bin_low = pt_bin_low+0.3;
        pt_bin_up = pt_bin_up+0.3;
        for (int k=0;k<bins.size();k++){
          for (int l=0;l<cut.size();l++){
          //hist for URQMD (data)
        TH1F* h1 = new TH1F("h1", Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
      
        //hist for DCM (simulated data)
        TH1F* h0 = new TH1F("h0",Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
      
        //hist for background only outisde the lambda peak from URQMD model
        TH1F* hback = new TH1F("hback", Form("rapidity=[%g,%g] & p_{T}=[%g,%g]",y_bin_low,y_bin_up,pt_bin_low,pt_bin_up), bins[k], 1.08, 1.2);
        Float_t sum = 0; Float_t sum_mc = 0; 
        for (int i =0; i < t->GetEntries(); i++)
          {
          t->GetEntry(i);
          if ( (pT<pt_bin_up) && (pT>pt_bin_low) && (rapidity >y_bin_low) && (rapidity <y_bin_up) && (xgb_preds>cut[l]) && (original_issignal<2)){
            h1->Fill(mass);//Fill total mass URQMD
            if(issignal>0){
            sum_mc+=1;              
            }
            if((mass<1.108)| (mass>1.13)){hback->Fill(mass);}//Fill background part
            }
          }
          
          for (int i =0; i < t1->GetEntries(); i++)
          {
          t1->GetEntry(i);
          if ( (MCpT<pt_bin_up) && (MCpT>pt_bin_low) && (MCrapidity >y_bin_low) && (MCrapidity <y_bin_up) && (MCxgb_preds>cut[l]) && (MCissignal>0) &&(MCoriginal_issignal==1)){
            h0->Fill(MCmass);//Fill signal only from DCM
            sum +=1;
            }
          }
          
              int j=0;
            
              if (sum>100){
                j++;
                Double_t minMassForFit0 = h0->GetMean()-4*(h0->GetRMS());
                Double_t maxMassForFit0 = h0->GetMean()+4*(h0->GetRMS());              
                //Fit DCM signal
                TF1* f0 = new TF1("DSCB", DSCB, minMassForFit0,maxMassForFit0, 7);
                f0->SetParameters(h0->GetBinContent(h0->FindBin(1.1157)),0,0.001,1,1,1,1);
                //f0->SetParameter(0, h0->GetBinContent(h0->FindBin(1.1157)));
                f0->SetParLimits(0, 0, 10*h0->GetBinContent(h0->FindBin(1.1157)));
                f0->SetParLimits(3, 0, 10);
                f0->SetParLimits(4, 1, 100);
                f0->SetParLimits(5, 0, 10);
                f0->SetParLimits(6, 1, 100);
                
                h0->Fit(f0,"R,M,E,+,0", "",minMassForFit0,maxMassForFit0);
                auto par0 = f0 -> GetParameters();
    
                Double_t minMassForFit = 1.08;
                Double_t maxMassForFit = 1.2; 
                //Fit URQMD background
                TF1* f1 = new TF1("pol2", pol2, minMassForFit,maxMassForFit,3);
                hback->Fit(f1,"", "REMN",minMassForFit,minMassForFit);
                auto par1 = f1 -> GetParameters();

                //Fit URQMD total
                TF1* f2 = new TF1("DSCB+pol2", fit_func, minMassForFit,maxMassForFit, 10);
                f2->SetParameters(par1[0],par1[1],par1[2],h1->GetBinContent(h1->FindBin(1.1157)),par0[1],par0[2],par0[3],par0[4],par0[5],par0[6]);
                /*f2->SetParError(0, 0.01);
                f2->SetParError(1, 0.01);
                f2->SetParError(2, 0.01);
                f2->SetParError(3, 0.01);
                f2->SetParError(4, 0.01);
                f2->SetParError(5, 0.01);
                f2->SetParError(6, 0.01);*/
                //f2->SetParError(7, 0.001);
                //f2->SetParError(8, 0.01);
                //f2->SetParError(9, 0.01);
                //f2->SetStepSize(0.)

                //if (pt_bin_low<0.6){f2->SetParLimits(3, 0, 10*h1->GetBinContent(h1->FindBin(1.1157)));f2->SetParLimits(6, 0.1, 100);f2->SetParLimits(8, 0.1, 100);}
                f2->SetParLimits(3, 0, 10*h1->GetBinContent(h1->FindBin(1.1157)));
                f2->SetParLimits(6, 0, 10);
                f2->SetParLimits(7, 1, 100);
                f2->SetParLimits(8, 0, 10);
                f2->SetParLimits(9, 1, 100);
                
                auto fitResult = h1->Fit(f2,"R,M,E,+,0,S", "",minMassForFit,maxMassForFit);
                //fitResult->Print();
                auto covMatrix = fitResult->GetCovarianceMatrix();
                auto par2 = f2 -> GetParameters();
                
                
                TF1* sig_func = new TF1("fs", DSCB, minMassForFit,maxMassForFit, 7); 
                sig_func->SetParameters(par2[3],par2[4],par2[5],par2[6],par2[7],par2[8],par2[9]);
                TF1* bac_func = new TF1("background", pol2, minMassForFit,maxMassForFit, 3); 
                bac_func->SetParameters(par2[0],par2[1],par2[2]);
                
                
                
                 h1 = PlotHistos(h1);
                 f2 = PlotMass(f2);
                 sig_func = PlotSignal(sig_func);
                 bac_func = PlotBackground(bac_func);
                  

                  

                  TH1F *func_hist = new TH1F("func_hist", "", bins[k], minMassForFit, maxMassForFit);
                  TH1F* diff_hist = new TH1F("diff_hist", "",bins[k], minMassForFit, maxMassForFit);
                  
                  for (int i =h1->FindBin(minMassForFit); i<h1->FindBin(maxMassForFit); i++)
                  {
                    Double_t f_value= f2->Eval(h1->GetBinCenter(i));
                    Double_t t_value = h1->GetBinContent(i);
                    func_hist->SetBinContent(i,f_value);
                    if (h1->GetBinError(i) > 0){
                      diff_hist->SetBinContent(i,(t_value-f_value)/h1->GetBinError(i));}
                  }
                  diff_hist = dif_hist(diff_hist);
                  
                  TLine* line = new TLine (minMassForFit,0 ,maxMassForFit ,0);
                  line -> SetLineColor ( kRed );
                
                  
                  canvas ->Draw();
                  canvas -> Clear();
                  TPad* pad1 = CreatePad1();
                  pad1 -> Clear();
                  pad1->cd();
                  pad1->SetLogy();
                  h1->Draw("pe,X0");
                  f2->Draw("SAME");//ESAME
                  sig_func->Draw("SAME");
                  bac_func->Draw("SAME");
                  TLegend* legend = legend_plot(h1, f2, sig_func, bac_func);
                  
                  
                  

                  Double_t min_for_sig = (1.1157+par2[4]) - 5*(f2->GetParameter(5));
                  cout<<"minimum = "<<min_for_sig<<endl;
                  Double_t max_for_sig = (1.1157+par2[4]) + 5*(f2->GetParameter(5));
                  Double_t fRawYield= (f2->Integral(min_for_sig,max_for_sig)) / h1->GetBinWidth(1);
                  double fRawYieldErr = (f2->IntegralError(min_for_sig,max_for_sig, fitResult->GetParams() , covMatrix.GetMatrixArray(), 0.0000000001) )/h1->GetBinWidth(1);
                  Double_t signal = sig_func->Integral(min_for_sig,max_for_sig)/h1->GetBinWidth(1);
                  //Double_t errsignal= sig_func->IntegralError(minMassForFit,maxMassForFit)/h1->GetBinWidth(1);
                  Double_t errsignal=(fRawYieldErr/fRawYield)*signal;/*assume relative error is the same as for total integral*/
                  
                  auto latex = new TLatex ();
                  latex -> SetNDC ();
                  latex -> SetTextSize (0.03);
                  latex -> DrawLatex (0.2 ,0.8, Form("#chi_{red}^{2} = %4.3f",f2->GetChisquare()/f2->GetNDF()));
                  latex -> DrawLatex (0.48 ,0.8, Form("m_{0}: = %4.5f #pm %4.3f", 1.1157+par2[4], f2->GetParError(4) ));
                  latex -> DrawLatex (0.48 ,0.75, Form("#sigma: = %4.5f #pm %4.5f", f2->GetParameter(5), f2->GetParError(5) ));
                  latex -> DrawLatex (0.48 ,0.85, Form("s = %4.5f #pm %4.5f; mc = %4.3f", signal, errsignal, sum_mc));
                  latex -> DrawLatex (0.48 ,0.5, Form("raw yield = %4.5f #pm %4.5f", fRawYield, fRawYieldErr));
                  latex -> DrawLatex (0.48 ,0.4, "CBM Performance");
                  latex -> DrawLatex (0.48 ,0.35, "URQMD, AuAu @ 12 #it{A}GeV/#it{c}");
                  latex -> Draw("SAME");
                  
  
                  
                  canvas->cd();
                  TPad* pad2 = CreatePad2();
                  pad2 -> Clear();
                  pad2->cd();   
                  diff_hist->Draw();
                  line->Draw("SAME");
                  
                  red_chi2.push_back(f2->GetChisquare()/f2->GetNDF());
                  y_bin_low1.push_back(y_bin_low);
                  pt_bin_low1.push_back(pt_bin_low);
                  yield.push_back(signal);
                  if (errsignal>0){yield_error.push_back(errsignal);}
                  if (errsignal==0){yield_error.push_back(sqrt(signal));}
                  mc_yield.push_back(sum_mc);
                  sim_mc_yield.push_back(sum);
                  
                  if (j==1){
                    canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf(","pdf");
                            }
                  else{
                      canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf","pdf");
                      }
                       
                delete func_hist, diff_hist;
              }
                        
                else {
                  red_chi2.push_back(0.0);
                  y_bin_low1.push_back(y_bin_low);
                  pt_bin_low1.push_back(pt_bin_low);
                  yield.push_back(0.0);
                  yield_error.push_back(0.0);
                  mc_yield.push_back(0.0);
                  sim_mc_yield.push_back(0.0);
                  
                  //cout << ".........low bin counts" << endl;

                              }
      
            delete h0;delete h1; delete hback;            } } }          
        
        }
      
      

    canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf)","pdf");
    
    TH2F* h7 = hist2d_pt_y();  h7->SetTitle("#chi^{2}_{red}");h7->GetYaxis()->SetRangeUser (0 ,3);h7->GetXaxis()->SetRangeUser (0 ,3);
    for(int i=0;i<225;i++)
    {Double_t y= y_bin_low1[i];    Double_t pT=pt_bin_low1[i];  Double_t y_bin = int((y+0.1)/0.2 + 1); Double_t pT_bin = int((pT+0.1)/0.2 + 1);
    h7->SetBinContent(y_bin,pT_bin,red_chi2[i]);    
    }
    auto c0 = new TCanvas(" c0 ","", 500,500);
    c0->Draw();
    //gStyle->SetPaintTextFormat("4.2f");
    //h7->Draw("colz1");
    //h7->Draw("TEXT SAME");    
    //c0->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/red_chi2.png","png");
    
    const char data_path [] = "/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/urqmd/m_200_400_prim_z_minus_10_80_urqmd.root";
    const char MC_path   [] = "/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/dcm/m_200_400_prim_z_minus_10_80_dcm.root";
    eff_creator(y_bin_low1, pt_bin_low1, yield, yield_error, sim_mc_yield, data_path,MC_path);
    
    TH1F* h5 = new TH1F("h1","",yield.size(),0,yield.size()); 
    for (int i=0;i<yield.size();i++)
    {
      if(((TMath::Abs(yield_error[i]))>0)  && (mc_yield[i]>0) && (yield[i]>0))
      {
      h5->SetBinContent(i+1,yield[i]/mc_yield[i]);
      h5->SetBinError(i+1,(yield[i]/mc_yield[i]) * (TMath::Sqrt( (yield_error[i]/yield[i])*((yield_error[i]/yield[i])) + ((TMath::Sqrt(mc_yield[i])/mc_yield[i]))*((TMath::Sqrt(mc_yield[i])/mc_yield[i])))));
      }
    }
    h5->SetStats(0);
    //h5->GetYaxis()->SetRangeUser (0.9 ,1.1); h5->GetYaxis()->SetTitle("yield/MC");
    h5->Draw("pe");
    TLine* line = new TLine (0,1 ,yield.size() ,1);
    line -> SetLineColor ( kRed );
    line->Draw("SAME");
    
  }
  
  int main(int argc,char** argv) {
    TString name = argv[1];
    crystal_ball(name);
    return 0;
  }
          

                

  
  
