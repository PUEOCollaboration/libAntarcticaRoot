#include "TH2DAntarctica.h"
#include "RampdemReader.h"
#include "AntarcticaBackground.h"
#include "TRandom3.h"
#include "TVirtualPad.h"
#include "TPaletteAxis.h"
#include "TList.h"

ClassImp(TProfile2DAntarctica)
ClassImp(TH2DAntarctica)

const double palX1 = AntarcticaBackgroundDefaults::zAxisRightMargin;
const double palX2 = AntarcticaBackgroundDefaults::zAxisRightMargin + AntarcticaBackgroundDefaults::zAxisWidth;
const double palY1 = 1 - AntarcticaBackgroundDefaults::zAxisTopBottomMargin - AntarcticaBackgroundDefaults::zAxisHeight;
const double palY2 = 1 - AntarcticaBackgroundDefaults::zAxisTopBottomMargin;

TProfile2DAntarctica::TProfile2DAntarctica(Int_t nx, Int_t ny)
    : TProfile2D(), fAntarcticaBackground(NULL){

  AntarcticaBackground* b = getBackground();

  nx = nx <= 0 ? b->GetNbinsX() : nx;
  ny = ny <= 0 ? b->GetNbinsY() : ny;
  SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
          ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));
}


TProfile2DAntarctica::TProfile2DAntarctica(const char* name, const char* title, Int_t nx, Int_t ny)
    : TProfile2D(), fAntarcticaBackground(NULL){

  SetNameTitle(name, title);
  AntarcticaBackground* b = getBackground();

  nx = nx <= 0 ? b->GetNbinsX() : nx;
  ny = ny <= 0 ? b->GetNbinsY() : ny;
  SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
          ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));
}



void TProfile2DAntarctica::Draw(Option_t* opt){

  AntarcticaBackground* b = getBackground();
  b->Draw();
  b->SetBit(kCanDelete, false);
  TString sameOpt = TString::Format("%s same", opt);
  TProfile2D::Draw(sameOpt);
  ResetColorAxis();
  RescaleBackground();
  gPad->Modified();
  gPad->Update();
}


Int_t TProfile2DAntarctica::Fill(Double_t lon, Double_t lat, Double_t val){

  Double_t easting, northing;
  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  return TProfile2D::Fill(easting, northing, val);
}


/** 
 * This is just to try out drawing one 2D histogram on top of another background
 * 
 * @param nTimes 
 */
void TProfile2DAntarctica::FillRandomly(Int_t nTimes){

  TRandom3 randy;

  for(int i=0; i < nTimes; i++){

    double lon = randy.Uniform(0, 360);
    double lat = randy.Uniform(-90, -60);

    Fill(lon, lat);
  }
}



/**
 * Helper function which prettifies the z-axis
 */
void TProfile2DAntarctica::ResetColorAxis(){
  gPad->Modified();
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*) GetListOfFunctions()->FindObject("palette");
  if(palette){
    palette->SetX1NDC(palX1);
    palette->SetX2NDC(palX2);
    palette->SetY1NDC(palY1);
    palette->SetY2NDC(palY2);
    
    TAxis* zAxis = GetZaxis();
    // zAxis->SetTitle();
    zAxis->SetTitleSize(AntarcticaBackgroundDefaults::zAxisTextSize);
    zAxis->SetLabelSize(AntarcticaBackgroundDefaults::zAxisTextSize);
    // std::cout << zAxis->GetTitleOffset() << std::endl;
    // zAxis->SetTitleOffset(0.1);
    gPad->Modified();
    gPad->Update();
  }
}


void TProfile2DAntarctica::ExecuteEvent(int event, int px, int py){
  if(event==kButton1Double){
    this->UnZoom();
  }
  getBackground()->ExecuteEvent(event, px, py);
  TProfile2D::ExecuteEvent(event, px, py);
}
























TH2DAntarctica::TH2DAntarctica(Int_t nx, Int_t ny)
    : TH2D(), fAntarcticaBackground(NULL){

  AntarcticaBackground* b = getBackground();

  nx = nx <= 0 ? b->GetNbinsX() : nx;
  ny = ny <= 0 ? b->GetNbinsY() : ny;
  SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
          ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));
}


TH2DAntarctica::TH2DAntarctica(const char* name, const char* title, Int_t nx, Int_t ny)
    : TH2D(), fAntarcticaBackground(NULL){

  SetNameTitle(name, title);
  AntarcticaBackground* b = getBackground();

  nx = nx <= 0 ? b->GetNbinsX() : nx;
  ny = ny <= 0 ? b->GetNbinsY() : ny;
  SetBins(nx, b->GetXaxis()->GetBinLowEdge(1), b->GetXaxis()->GetBinUpEdge(b->GetNbinsX()),
          ny, b->GetYaxis()->GetBinLowEdge(1), b->GetYaxis()->GetBinUpEdge(b->GetNbinsY()));
}



void TH2DAntarctica::Draw(Option_t* opt){

  AntarcticaBackground* b = getBackground();
  b->Draw();
  b->SetBit(kCanDelete, false);  
  TString sameOpt = TString::Format("%s same", opt);
  TH2D::Draw(sameOpt);
  ResetColorAxis();
  RescaleBackground();
  gPad->Modified();
  gPad->Update();
}


Int_t TH2DAntarctica::Fill(Double_t lon, Double_t lat, Double_t val){

  Double_t easting, northing;
  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
  return TH2D::Fill(easting, northing, val);  
}


/** 
 * This is just to try out drawing one 2D histogram on top of another background
 * 
 * @param nTimes 
 */
void TH2DAntarctica::FillRandomly(Int_t nTimes){

  TRandom3 randy;

  for(int i=0; i < nTimes; i++){

    double lon = randy.Uniform(0, 360);
    double lat = randy.Uniform(-90, -60);

    Fill(lon, lat);
  }
}


/**
 * Helper function which prettifies the z-axis
 */
void TH2DAntarctica::ResetColorAxis(){
  gPad->Modified();
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*) GetListOfFunctions()->FindObject("palette");
  if(palette){
    palette->SetX1NDC(palX1);
    palette->SetX2NDC(palX2);
    palette->SetY1NDC(palY1);
    palette->SetY2NDC(palY2);
    
    TAxis* zAxis = GetZaxis();
    // zAxis->SetTitle();
    zAxis->SetTitleSize(AntarcticaBackgroundDefaults::zAxisTextSize);
    zAxis->SetLabelSize(AntarcticaBackgroundDefaults::zAxisTextSize);
    // std::cout << zAxis->GetTitleOffset() << std::endl;
    // zAxis->SetTitleOffset(0.1);
    gPad->Modified();
    gPad->Update();
  }
}


void TH2DAntarctica::ExecuteEvent(int event, int px, int py){
  if(event==kButton1Double){
    this->UnZoom();
  }
  getBackground()->ExecuteEvent(event, px, py);
  TH2D::ExecuteEvent(event, px, py);
}
