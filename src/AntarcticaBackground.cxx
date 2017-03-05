#include "AntarcticaBackground.h"
#include "TGraphAntarctica.h"
#include "TVirtualPad.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"


ClassImp(AntarcticaBackground)



// simple global variable to increment integer in name so we can have as manny of these background as we like.
static int numAntarcticaBackgrounds = 0;


AntarcticaBackground::AntarcticaBackground(RampdemReader::dataSet dataSet, Int_t coarseness)
  : TProfile2D() {
  init(dataSet, coarseness);
}


void AntarcticaBackground::init(RampdemReader::dataSet dataSet, Int_t coarseness){
  SetDirectory(0);
  fName = TString::Format("fAntarctica%d", numAntarcticaBackgrounds);
  numAntarcticaBackgrounds++;
  fDataSet = dataSet;
  fCoarseness = coarseness;
  needRemakeHist = true;  // force updateHist() to read in data by setting impossible coarseness
  fAlreadyDrawn = false; // set by Draw(), needed by updateHist()
  updateHist();

  fGridPoints = 1000;
  fDeltaLon = 45; // degrees
  fDeltaLat = 5; // degrees
  fGrid = false;
  needRemakeGrid = true; // force updateGrid() to make grid TGraphs on first call
  updateGrid();

  fUseToolTip = true;
  fToolTip = NULL;
}



/**
 * Update the map of Antarctica as different options are selected
 */
void AntarcticaBackground::updateHist(){


  if(needRemakeHist){

    // Here I save the current drawn axis range so I can set it to be the same after updating contents
    Int_t firstX = fXaxis.GetFirst();
    Int_t lastX = fXaxis.GetLast();
    Int_t firstY = fYaxis.GetFirst();
    Int_t lastY = fYaxis.GetLast();

    Double_t lowX = fXaxis.GetBinLowEdge(firstX);
    Double_t highX = fXaxis.GetBinUpEdge(lastX);
    Double_t lowY = fXaxis.GetBinLowEdge(firstY);
    Double_t highY = fXaxis.GetBinUpEdge(lastY);

    // std::cout << firstX << "\t" << lastX << "\t" << lowX << "\t" << highX << std::endl;


    // now I get the new histogram binning
    Int_t nx, ny;
    RampdemReader::getNumXY(nx, ny, fDataSet);
    Double_t xMin, xMax, yMin, yMax;
    RampdemReader::getMapCoordinates(xMin, yMin, xMax, yMax, fDataSet);

    // accounting for coarseness
    nx /= fCoarseness;
    ny /= fCoarseness;


    // Get rid of old the AntarcticaBackground content
    fBinEntries.Reset();
    for(int by=0; by <= GetNbinsY() + 1; by++){
      for(int bx=0; bx <= GetNbinsX() + 1; bx++){
	SetBinContent(bx, by, 0);
      }
    }
    // change the histogram dimensions
    SetBins(nx, xMin, xMax, ny, yMin, yMax);

    // insert new data
    RampdemReader::fillThisHist(this, fDataSet);

    if(fAlreadyDrawn){

      // now set the viewing range to the same as before if we have a gPad instance
      fXaxis.SetRangeUser(lowX, highX);
      fYaxis.SetRangeUser(lowY, highY);

      // and update the z-axis title if needed
      prettifyPalette();
    }

  }

  // prettification
  GetXaxis()->SetNdivisions(0, kFALSE);
  GetYaxis()->SetNdivisions(0, kFALSE);

  needRemakeGrid = false;


}






Int_t AntarcticaBackground::GetCoarseness(){
  return fCoarseness;
}



void AntarcticaBackground::ToolTip(Bool_t toolTip){
  fUseToolTip = toolTip;
  if(!fUseToolTip){
    fToolTip->Hide();
  }
}

Bool_t AntarcticaBackground::GetToolTip(){
  return fUseToolTip;
}


void AntarcticaBackground::SetCoarseness(Int_t coarseness){

  // sanity check
  if(coarseness < 1){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", coarsenesss must be >= 1. Setting coareness = "
	      << defaultCoarseness << std::endl;
    coarseness = defaultCoarseness;
  }

  needRemakeHist = fCoarseness == coarseness ? false : true;
  fCoarseness = coarseness;
  updateHist();
}


RampdemReader::dataSet AntarcticaBackground::GetDataSet(){
  return fDataSet;
}



void AntarcticaBackground::SetDataSet(RampdemReader::dataSet dataSet){

  needRemakeHist = dataSet == fDataSet ? false : true;
  fDataSet = dataSet;
  updateHist();
}





void AntarcticaBackground::Grid(Bool_t grid){
  fGrid = grid;

  updateGrid();

  if(gPad){
    TList* prims = gPad->GetListOfPrimitives();

    if(fGrid){
      // Manually add the grid TGraphs into the list of pad primitives...
      for(UInt_t grInd=0; grInd < grGrids.size(); grInd++){
	TGraph* gr = grGrids.at(grInd);
	prims->AddAfter(this, gr);
      }
    }
    else{
      // ...or manually add the grid TGraphs into the list of pad primitives
      for(UInt_t grInd=0; grInd < grGrids.size(); grInd++){
	TGraph* gr = grGrids.at(grInd);
	prims->RecursiveRemove(gr);
      }
    }
  }

}

Bool_t AntarcticaBackground::GetGrid(){
  return fGrid;
}



void AntarcticaBackground::updateGrid(){

  if(needRemakeGrid){

    const int minLat = -90 + fDeltaLat;
    const int maxLat = -50;

    deleteGrid();

    std::cout << "here" << std::endl;

    // make circles of constant latitude
    for(Int_t lat = minLat; lat<= maxLat; lat += fDeltaLat){
      std::cout << "lat\t" << lat << "\t" << fDeltaLat << std::endl;
      TGraph* gr = new TGraph();
      gr->SetLineColor(kGray);
      const Double_t deltaLon = 360./fGridPoints;
      for(int i=0; i < fGridPoints; i++){
	Double_t theLat = lat;
	Double_t theLon = i*deltaLon;
	Double_t easting, northing;
	RampdemReader::LonLatToEastingNorthing(theLon, theLat, easting, northing);
	gr->SetPoint(gr->GetN(), easting, northing);
	// std::cout << gr << "\t" << gr->GetN() << "\t" << easting << "\t" << northing << std::endl;
      }
      gr->SetEditable(false);
      gr->SetName(Form("latitude=%d", lat)); // descriptive name
      grGrids.push_back(gr);
    }

    // make lines of constant longitude
    for(Int_t lon = 0; lon < 360; lon+= fDeltaLon){
      std::cout << "lon\t" << lon << "\t" << fDeltaLat << std::endl;
      TGraph* gr = new TGraph();
      gr->SetLineColor(kGray);
      const Double_t deltaLat = double(maxLat - -90)/fGridPoints;
      for(int i=0; i < fGridPoints; i++){
	Double_t theLat = -90 + deltaLat*i;
	Double_t theLon = lon;
	Double_t easting, northing;
	RampdemReader::LonLatToEastingNorthing(theLon, theLat, easting, northing);
	gr->SetPoint(gr->GetN(), easting, northing);
      }
      gr->SetEditable(false);
      gr->SetName(Form("longitude=%d", lon)); // descriptive name
      grGrids.push_back(gr);
    }

    needRemakeGrid = false;
  }

}




void AntarcticaBackground::deleteGrid(){

  // Grid(false); // avoid segfault?
  while(grGrids.size() > 0){
    TGraph* gr = grGrids.back();
    delete gr;
    grGrids.pop_back();
  }

}






void AntarcticaBackground::SetGridDivisions(Int_t deltaLon, Int_t deltaLat){

  needRemakeGrid = fDeltaLon == deltaLon && fDeltaLat == deltaLat ? false : true;
  fDeltaLon = deltaLon;
  fDeltaLat = deltaLat;
  updateGrid();
}



/**
 * Draw function for Antarctica Background which also prettifies the pad/canvas.
 *
 * @param opt is the draw option (default is colz)
 */
void AntarcticaBackground::Draw(Option_t* opt){
  TProfile2D::Draw(opt);
  setPadMargins();
  prettifyPalette();
  fXaxis.SetAxisColor(kWhite);
  fYaxis.SetAxisColor(kWhite);

  // gPad->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
  // 		"AntarcticaBackground::Interactive(Int_t,Int_t,Int_t,TObject*)");
  gPad->Update();


  fAlreadyDrawn = true;
  fToolTip = new TGToolTip();

}


/**
 * Helper function which prettifies the pad
 */
void AntarcticaBackground::setPadMargins(){
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.02);
  gPad->SetLeftMargin(0.02);
  gPad->SetRightMargin(0.1);
  // gPad->SetTopMargin(0.05);
  // gPad->SetBottomMargin(0.05);
  // gPad->SetLeftMargin(0.05);
  // gPad->SetRightMargin(0.1);
  // gPad->SetRightMargin(0.05);
  gPad->SetFrameLineColor(0);
  gPad->SetFrameLineWidth(0);
  gPad->SetFrameBorderSize(0);

}



/**
 * Helper function which prettifies the z-axis
 */
void AntarcticaBackground::prettifyPalette(){
  gPad->Modified();
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*) GetListOfFunctions()->FindObject("palette");
  if(palette){
    palette->SetX1NDC(0.91);
    palette->SetX2NDC(0.95);
    palette->SetY1NDC(0.55);
    palette->SetY2NDC(0.95);
    // palette->SetX1NDC(0.03);
    // palette->SetX2NDC(0.06);
    // palette->SetY1NDC(0.03);
    // palette->SetY2NDC(0.16);
    // palette->SetTitleSize(0.001);
    // palette->SetTitleOffset(0.1);

    TAxis* zAxis = GetZaxis();
    zAxis->SetTitle(RampdemReader::dataSetToAxisTitle(fDataSet));
    zAxis->SetTitleSize(0.001);
    std::cout << zAxis->GetTitleOffset() << std::endl;
    zAxis->SetTitleOffset(25);
    gPad->Modified();
    gPad->Update();
  }
}




// Interactive functions...

void AntarcticaBackground::Rampdem(bool useRampdem){
  SetDataSet(RampdemReader::rampdem);
}

Bool_t AntarcticaBackground::GetRampdem(){
  return fDataSet == RampdemReader::rampdem;
}

void AntarcticaBackground::Bed(bool useBed){
  SetDataSet(RampdemReader::bed);
};

Bool_t AntarcticaBackground::GetBed(){
  return fDataSet == RampdemReader::bed;
}

void AntarcticaBackground::Icemask(bool useIcemask){
  SetDataSet(RampdemReader::icemask_grounded_and_shelves);
}

Bool_t AntarcticaBackground::GetIcemask(){
  return fDataSet == RampdemReader::icemask_grounded_and_shelves;
}

void AntarcticaBackground::Surface(bool useSurface){
  SetDataSet(RampdemReader::surface);
}

Bool_t AntarcticaBackground::GetSurface(){
  return fDataSet == RampdemReader::surface;
}

void AntarcticaBackground::Thickness(bool useThickness){
  SetDataSet(RampdemReader::thickness);
}

Bool_t AntarcticaBackground::GetThickness(){
  return fDataSet == RampdemReader::thickness;
}



void AntarcticaBackground::ExecuteEvent(Int_t event, Int_t x, Int_t y)
{
//   static int keyWasPressed=0;
//    switch (event) {
//    case kKeyPress:
//      //     std::cout << "kKeyPress" << std::endl;
//      keyWasPressed=1;
//      break;
//    case kButtonPress:
//      //     cout << "kButtonPress" << endl;
//      break;

//    case kButton1Double:
//      //     std::cout << "kButtonDoubleClick" << std::endl;
//      //     new TCanvas();
//      break;

//    case kButton1Down:
//      //     std::cout << "kButton1Down" << std::endl;
//      if(!keyWasPressed) {
//        if(!fNewCanvas) drawInNewCanvas();
//        else this->TGraph::ExecuteEvent(event,px,py);
//      }
//      else {
//        //       std::cout << "ctrl + click\n";
//        CorrelationFactory::Instance()->addWaveformToCorrelation(this);
//        keyWasPressed=0;
//      }

//      break;

//    default:
//        this->TGraph::ExecuteEvent(event,px,py);
//        break;
//    }
// }

// void AntarcticaBackground::Interactive(Int_t event, Int_t x, Int_t y, TObject* selected){

  // std::cout << event << "\t" << x << "\t" << y << std::endl;

  if(fUseToolTip){
    Double_t easting = gPad->AbsPixeltoX(x);
    Double_t northing = gPad->AbsPixeltoY(y);
    Double_t val = GetBinContent(FindBin(easting, northing));
    // gPad->AbsPixeltoX(y);
    Double_t lon, lat;
    RampdemReader::EastingNorthingToLonLat(easting, northing, lon, lat);
    fToolTip->SetText(Form("Lon %4.2lf \nLat %4.2lf\n%4.2f", lon, lat, val));
    fToolTip->Show(x, y);
  }


  TProfile2D::ExecuteEvent(event, x, y);
}
