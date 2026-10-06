#include "RampdemReader.h"
#include "TCanvas.h"

int main(){
  TProfile2D* myMap = RampdemReader::getMap(RampdemReader::surface, 4);

  TCanvas c1("","", 1920, 1080);
  myMap->Draw("COLZ");
  c1.SaveAs("foo.pdf");



  return 0;
}
