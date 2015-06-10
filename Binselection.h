#ifndef BINSELECTION_H
#define BINSELECTION_H

class Binselection {
 public:
  // Returns result of search region bin
  static unsigned int searchBin(double met, double mt2, unsigned int nBjets, unsigned int nTop);

};

unsigned int Binselection::searchBin(double met, double mt2, unsigned int nBjets, unsigned int nTop) {
  unsigned int bin = 0;//this will contain total events
  if(nBjets ==1 && nTop == 1){
    if(met== && mt2==)bin == ;

}
  if(nBjets ==2 && nTop == 1){
    if(met== && mt2==)bin == ;

  }
  if(nBjets >=3 && nTop == 1){
    if(met== && mt2==)bin == ;

  }
  if(nBjets ==1 && nTop == 2){
    if(met== && mt2==)bin == ;

  }
  if(nBjets ==2 && nTop == 2){
    if(met== && mt2==)bin == ;

  }
  if(nBjets >=3 && nTop == 2){
    if(met== && mt2==)bin == ;

  }
  if(nBjets ==1 && nTop == 3){
    if(met== && mt2==)bin == ;

  }
  if(nBjets ==2 && nTop == 3){
    if(met== && mt2==)bin == ;

  }
  if(nBjets >=3 && nTop == 3){
    if(met== && mt2==)bin == ;

  }
  return bin;
}

#endif
