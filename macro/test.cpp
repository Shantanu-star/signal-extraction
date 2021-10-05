#if !defined(__CLING__)

#include "AliHFInvMassFitter.h"
#include <TH1F.h>
#include <TF1.h>

#endif

TH1F* CreateTestHisto(){
  TH1F* h = new TH1F("h_test", "", 100, 1.05, 1.2);
  auto f = new TF1("f","gaus(0) + gaus(3) + pol1(5)",1.05, 1.2);
  f->SetParameters(2,1.115,0.002,1,1.115,0.005, 0.1, 0.1);

  h->FillRandom("f", 20000);
  return h;
}

void test(){
  auto h = CreateTestHisto();
//  h->Draw();

  auto fitter = new AliHFInvMassFitter(h, 1.05, 1.2, AliHFInvMassFitter::ETypeOfBkg::kLin, AliHFInvMassFitter::ETypeOfSgn::k2Gaus);
  fitter->SetInitialFrac2Gaus(0.5);
  fitter->SetInitialGaussianMean(1.115);
  fitter->SetInitialGaussianSigma(0.002);
  fitter->SetInitialSecondGaussianSigma(0.005);
  Bool_t ok = fitter->MassFitter(kTRUE);
  if(!ok){
    std::cout << "Not OK!" << std::endl;
  }
}

int main(){
  test();
  return 0;
}