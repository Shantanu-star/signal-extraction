#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include "plot_tools.C"

// #include "AliHFInvMassFitter.h"
// #include "AliHFMassFitter.h"
#endif

enum { kGaus = 0,
  k2Gaus=1, kPol2=2, kNoBk=3, kPol3=6};
  

  

  Int_t typeb = kNoBk;
  Int_t types = k2Gaus;
  
  using std::cout;
  using std::endl;
  
  vector<Double_t> FitXicZerotoXiPiInvMass1(TString fname, TString fNameOut = "") 
    {
    
    vector<Double_t> red_chi2, y_bin_low1, pt_bin_low1;
    TFile* f = TFile::Open(fname.Data(), "read");    
    TTree* t = f->Get<TTree>("t1");
    
    float mass, pT, rapidity, MCmass, MCpT, MCrapidity;Double_t issignal;
    t->SetBranchAddress("mass", &mass);t->SetBranchAddress("pT",&pT);t->SetBranchAddress("rapidity",&rapidity);t->SetBranchAddress("issignal",&issignal);
    auto* canvas = new TCanvas(" canvas ","", 500,500);
    canvas->SetBottomMargin(0.15);
    canvas->SetLeftMargin (0.15);
    //canvas->SetLogy();
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
        
        for (int i =0; i < t->GetEntries(); i++)
          {
          t->GetEntry(i);
          if ( (pT<pt_bin_up) && (pT>pt_bin_low) && (rapidity >y_bin_low) && (rapidity <y_bin_up) ){
            h1->Fill(mass);
            if(issignal>0){
              h0->Fill(mass);
              }
            }
          }
            Float_t sum = 0;
            for (Int_t i=1;i<=200;i++) {
              sum += h0->GetBinContent(i);
                }
              int j=0;
            
              if (sum>1000){
                j++;
                Double_t minMassForFit0 = h0->GetMean()-4*(h0->GetRMS());
                Double_t maxMassForFit0 = h0->GetMean()+4*(h0->GetRMS());              

                AliHFInvMassFitter* fitter0 = new AliHFInvMassFitter(h0, minMassForFit0, maxMassForFit0,   typeb, types);
                fitter0->SetInitialGaussianMean(1.1156);
                fitter0->SetInitialGaussianSigma(0.0009);
                //fitter0->SetInitialFrac2Gaus(0.2);
                fitter0->SetInitialSecondGaussianSigma(0.0009);
                //fitter->SetPolDegreeForBackgroundFit(0);
                TF1 *f2 = new TF1("double_gaus",double_gaus,minMassForFit0,maxMassForFit0,5);
                //TF1* f2 = new TF1("full","[0]*((1-[3])*exp(-0.5*((x-[2])/[1])^2)+[3]*exp(-0.5*((x-[2])/[4])^2)) ",minMassForFit0,maxMassForFit0);
                f2->SetParameters(0,1.1156,0.001,0.2,0.0009);
                h0->Fit(f2,"MNRI","L",minMassForFit0,maxMassForFit0);
                
                //TF1* pointer = fitter0->GetBackgroundFullRangeFunc();
                //pointer->FixParameter(0,0);
                
                Bool_t ok = fitter0->MassFitter(kFALSE);
    
                if (ok) {
                Double_t minMassForFit = h1->GetMean()-1.5*(h1->GetRMS());
                Double_t maxMassForFit = h1->GetMean()+3*(h1->GetRMS()); 
                
                TF1 *sig_func0 = fitter0->GetSignalFunc();

              
                  
                  
                  
                  //fitter->DrawHere(canvas, 3, 1);
                  canvas -> SetLogy();
                  canvas ->Draw();
                  h0 = PlotHistos(h0);
                  h0->GetXaxis()->SetRangeUser(1.105, 1.13);
                  h0->Draw("pe,X0");

                  //h1 = PlotHistos(h1);
                  //h1->Draw("pe,X0");

                  TF1 *sig_func = fitter0->GetSignalFunc();
                  sig_func = PlotSignal(sig_func);

                  //cout<<"This is the mean"<<sig_func->GetParameter(1)<<endl;
                  
                  
                  TF1 *bac_func= fitter0->GetBackgroundRecalcFunc();
                  bac_func->SetNpx(100);
                  //bac_func->SetLineStyle(3);
                  bac_func->SetLineColor(kYellow);
                  bac_func->SetLineWidth(4);
                  //bac_func->Draw("SAME");
                  
                  TF1 *M_func = fitter0->GetMassFunc();
                  M_func->SetNpx(1000);
                  //M_func->SetLineStyle(3);
                  M_func->SetLineWidth(4);
                  M_func->SetLineColor(kRed);
                  M_func->Draw("SAME");//ESAME
                  f2->SetLineColor(kBlue);
                  f2->Draw("SAME");
                  sig_func->Draw("SAME");
                  bac_func->Draw("SAME");

                  
                  
                  auto legend = new TLegend(0.5,0.6,0.8,0.9);
                  legend->AddEntry(h0,"#Lambda hyperon","pe,X0");
                  legend->AddEntry(M_func,"A[#frac{(1-B)}{#sqrt{2#pi}/#sigma_{1}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{(#sigma_{1})^{2}}}+#frac{(B)}{#sqrt{2#pi}/#sigma_{2}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{#sigma_{2}^{2}}}]+C","l");
                  
                  legend->AddEntry(sig_func,"A[#frac{(1-B)}{#sqrt{2#pi}/#sigma_{1}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{(#sigma_{1})^{2}}}+#frac{(B)}{#sqrt{2#pi}/#sigma_{2}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{#sigma_{2}^{2}}}]","l");
                  
                  legend->AddEntry(f2,"A[#frac{(1-B)}{#sqrt{2#pi}/#sigma_{1}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{(#sigma_{1})^{2}}}+#frac{(B)}{#sqrt{2#pi}/#sigma_{2}}e^{#frac{-1}{2} #frac{(m-m_{0})^{2}}{#sigma_{2}^{2}}}]","l");
                  
                  //legend->AddEntry(bac_func,"C+Dx+Ex^{2}","l");
                  legend->AddEntry(bac_func,"C","l");
                  //legend->AddEntry(bac_func,"#frac{C}{m_{f}-m_{i}}+D(x-0.5(m_{f}-m_{i}))","l");
                  legend -> SetLineWidth (0);
                  legend->SetTextSize(0.02);
                  legend->Draw();
                  
                  auto latex = new TLatex ();
                  latex -> SetNDC ();
                  latex -> SetTextSize (0.02);
                  //latex -> DrawLatex (0.5 ,0.3, Form("#sigma_{1} = %4.5f, sigma_{2} = %4.5f, mean = %4.5f",f2->GetParameter(1), f2->GetParameter(4),f2->GetParameter(2)));
                  latex -> DrawLatex (0.18 ,0.85, Form("C = %4.5f #pm %4.5f", bac_func->GetParameter(0), bac_func->GetParError(0)));
                  latex -> DrawLatex (0.18 ,0.45, "#color[4]{Root basic}               vs                            #color[3]{Alice (signal function from total fit)}");
                  latex -> DrawLatex (0.18 ,0.4, Form("#color[4]{m_{0}: = %4.5f #pm %4.5f;}                         #color[3]{= %4.5f #pm %4.5f}",f2->GetParameter(1),f2->GetParError(1), sig_func->GetParameter(1), sig_func->GetParError(1) ));
                  latex -> DrawLatex (0.18 ,0.35, Form("#color[4]{#sigma_{1}: = %4.5f #pm %4.5f;}                         #color[3]{= %4.5f #pm %4.5f}",f2->GetParameter(2), f2->GetParError(2), sig_func->GetParameter(2), sig_func->GetParError(2)));
                  latex -> DrawLatex (0.18 ,0.3, Form("#color[4]{#sigma_{2}: = %4.5f #pm %4.5f;}                           #color[3]{= %4.5f #pm %4.5f}",f2->GetParameter(4), f2->GetParError(4), sig_func->GetParameter(4), sig_func->GetParError(4)));
                  latex -> DrawLatex (0.18 ,0.25, Form("#color[4]{#chi_{red}^{2} = %4.5f;}                                                #color[3]{=  %4.5f}",f2->GetChisquare() / f2->GetNDF(),fitter0->GetReducedChiSquare()));
                  latex -> Draw();
                  
                  y_bin_low1.push_back(y_bin_low);
                  pt_bin_low1.push_back(pt_bin_low);
                  red_chi2.push_back(fitter0->GetReducedChiSquare());
                  //red_chi2.push_back(f2->GetChisquare() / f2->GetNDF());
  
                  if (j==1){
                    canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf(","pdf");
                            }
                  else{
                      canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf","pdf");
                      }
                        
                        }
                else if (!ok) {
                  cout << ".........I am sorry" << endl;
                  fitter0->GetHistoClone()->Draw();
                              }
                             
                            //Double_t p1 = fitter->GetSignalFunc()->GetParameter(0);
                              //cout<<"This is the 1st parameter"<<par<<endl;
                delete fitter0;
      
                         } 
            else         {
                red_chi2.push_back(0);
                y_bin_low1.push_back(y_bin_low);
                pt_bin_low1.push_back(pt_bin_low);
                          }
            
        
        delete h1;
          delete h0;
        }
      }
      

    canvas->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal.pdf)","pdf");
    TH2F* h4 = new TH2F("reduced chi2", "reduced chi2", 15,0,3,15,0,3);
    for (int i=0;i<225;i++){
    Double_t y= y_bin_low1[i];
    Double_t pT=pt_bin_low1[i];
    Double_t y_bin = int((y+0.1)/0.2 + 1);
    Double_t pT_bin = int((pT+0.1)/0.2 + 1);
    h4->SetBinContent(y_bin, pT_bin, red_chi2[i]);
    }
    auto c = new TCanvas(" c ","", 500,500);
    c->Draw();
    h4->SetStats(0);
    h4->Draw("colz");
    gStyle->SetPaintTextFormat("4.1f");
    h4->Draw("TEXT SAME");
    //h4->GetYaxis()->SetRangeUser (0 ,2);
    //h4->GetXaxis()->SetRangeUser (0.4 ,2.6);
    c->Print ("/home/shahid/cbmsoft/Cut_optimization/uncut_data/Project/pT_rapidity_distribution_XGB_extracted_signal_hist.pdf","pdf");
    return red_chi2; 
  }
  
  int main(int argc,char** argv) {
    TString name = argv[1];
    vector<Double_t> red_chi2 = FitXicZerotoXiPiInvMass1(name);
    return 0;
  }
          

                

  
  
