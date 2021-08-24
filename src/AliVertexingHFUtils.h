#ifndef SIGNALEXTRACTION_SRC_ALIVERTEXINGHFUTILS_H_
#define SIGNALEXTRACTION_SRC_ALIVERTEXINGHFUTILS_H_

namespace AliVertexingHFUtils{

//______________________________________________________________________
inline void ComputeSignificance(Double_t signal,
                                              Double_t errsignal,
                                              Double_t background,
                                              Double_t errbackground,
                                              Double_t& significance,
                                              Double_t& errsignificance) {
  /// calculate significance from S, B and errors
  Double_t errSigSq=errsignal*errsignal;
  Double_t errBkgSq=errbackground*errbackground;
  Double_t sigPlusBkg=signal+background;
  if(sigPlusBkg>0. && signal>0.){
    significance =  signal/TMath::Sqrt(signal+background);
    errsignificance = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
  }else{
    significance=0.;
    errsignificance=0.;
  }
  return;
}

}

#endif //SIGNALEXTRACTION_SRC_ALIVERTEXINGHFUTILS_H_
