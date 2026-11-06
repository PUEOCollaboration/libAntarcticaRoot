#include "RampdemReader.h"
#include "TCanvas.h"

void rampdemReader_test(){
  TProfile2D* myMap = RampdemReader::getMap(RampdemReader::surface, 4);

  TCanvas c1("","", 1920, 1080);
  myMap->Draw("COLZ");
  c1.SaveAs("foo.pdf");


}

int main(){
  rampdemReader_test();
}
