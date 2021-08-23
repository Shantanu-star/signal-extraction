#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include <TTree.h>

#include "AliHFInvMassFitter.h"
#include "AliHFMassFitter.h"
#endif

enum {kGaus=0, kDoubleGaus};

Double_t minMassForFit=1.08;
Double_t maxMassForFit=2.0;

//const Int_t nPtBins=5;
//Double_t ptlims[nPtBins+1]={1.,2.,3.,5.,8.,12.};
//Int_t rebin[nPtBins]={2,2,1,1,2};
const Int_t nPtBins=15;
Double_t ptlims[nPtBins+1]={0.,3.};
Int_t rebin[nPtBins]={.2};
Int_t typeb=2;
Int_t types=kGaus;

void FitXicZerotoXiPiInvMass(TString fname, TString fNameOut=""){
  TFile *f = TFile::Open(fname.Data());
  if (!f) {
    printf("ERROR: file %s does not exist",fname.Data());
    return;
  } 

  
    TCanvas *cInvMass[nPtBins];
    for(int i=0; i<nPtBins; i++){
        cInvMass[i]=new TCanvas(Form("c%d",i), Form("c%d",i));
        cInvMass[i]->Divide(5,2);
    }
    
  float         mass;
  float         pT;
  float         mva;

  TBranch        *mass;
  TBranch        *pT;
  TBranch        *b_mva;

  TTree *t = (TTree*)f->Get("t1");
    
  t->SetBranchAddress("mass", &mass, &mass);
  t->SetBranchAddress("pT", &pT, &pT);
  t->SetBranchAddress("MVA2", &mva, &b_mva);

    TH2F *hMassVsPt[20];
   
    for(int j =0; j<10; j++){
        hMassVsPt[j] = new TH2F(Form("hMassVsPt%d",j),Form("hMassVsPt%d",j),80,2.468-0.2,2.468+0.2, 20, 0, 20);
        for (Int_t i=0; i<t->GetEntries(); i++){
            t->GetEntry(i);
            if (mva>j/10.) hMassVsPt[j]->Fill(mass,pT);
        }
    }
  TH1F** hMass=new TH1F*[nPtBins];
  
  for(Int_t i=0;i<nPtBins;i++) {
    hMass[i]=0x0;
  }
  AliHFInvMassFitter.h** fitter=new AliHFInvMassFitter.h*[10];

    cout<< "1"<<endl;
    
 Double_t sig,errsig,s,errs,b,errb;
 
 Double_t mxicPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    
for(int j =0; j<10; j++){
    
    for (Int_t i=0; i<nPtBins; i++){
    Printf("bin %d ", i);
    Int_t ptminbin = hMassVsPt[j]->GetYaxis()->FindBin(ptlims[i]+0.0001);
    Int_t ptmaxbin = hMassVsPt[j]->GetYaxis()->FindBin(ptlims[i+1]-0.0001);
    Printf("====== pt limits %f %f , bins %d %d", ptlims[i],ptlims[i+1], ptminbin, ptmaxbin );
    hMass[i]=(TH1F*)hMassVsPt[j]->ProjectionX(Form("projection_%d",i), ptminbin, ptmaxbin) ;
    cInvMass[i]->cd(j+1);
        hMass[i]->SetTitle(Form(" %0.0f < p_{T} < %0.0f GeV/c, MVA cut %f", ptlims[i],ptlims[i+1],j/10.));
    hMass[i]->Rebin(rebin[i]);
    hMass[i]->Sumw2();
    //hMass[i]->Draw();
   
    Int_t nMassBins=hMass[i]->GetNbinsX();
    Double_t hmin=TMath::Max(minMassForFit,hMass[i]->GetBinLowEdge(2));
    Double_t hmax=TMath::Min(maxMassForFit,hMass[i]->GetBinLowEdge(nMassBins-2)+hMass[i]->GetBinWidth(nMassBins-2));

    cout<<"LIMITS FOR FITS = " << hmin << " " << hmax << endl;

    fitter[i]=new AliHFMassFitter(hMass[i],hmin, hmax,1,typeb,types);
    //
    //    fitter[i]->SetInitialGaussianMean(2.467);
    fitter[i]->SetInitialGaussianMean(2.470);
//    fitter[i]->SetFixGaussianSigma(0.010);
  fitter[i]->SetInitialGaussianSigma(0.010);

    Printf("mean = %f\n", fitter[i][j]->GetMean());
    //Printf("mean = %f, sigma = %f \n", fitter[i]->GetMean(), fitter[i]->GetSigma());

    //fitter[i]->SetInitialGaussianSigma(0.012);
    Bool_t ok=fitter[i]->MassFitter(kFALSE);
    if(ok){
      fitter[i]->DrawHere(cInvMass[i]->cd(j+1),3,1);
    }
    else if(!ok) {
      cout << ".........I am sorry" <<endl;
      fitter[i]->GetHistoClone()->Draw();
    }
    // Double_t mass=fitter[i]->GetMean();
    //Double_t massErr= fitter[i]->GetMeanUncertainty();
    // Double_t sigma=fitter[i]->GetSigma();
    // Double_t sigmaErr=fitter[i]->GetSigmaUncertainty();
    //fitter[i]->DrawHere(gPad);    
    //fitter[i]->Signal(3,s,errs);
    //fitter[i]->Background(3,b,errb);
    //fitter[i]->Significance(3,sig,errsig);

    
  }//loop on pt
}
  if(fNameOut!=""){
      //cInvMass->SaveAs(Form("%s.png",fNameOut.Data()));
    }
}

