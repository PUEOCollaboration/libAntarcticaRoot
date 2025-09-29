#include "Geoid.h"

#include "TMath.h"
#include "TClass.h"

#include <iostream>

ClassImp(Geoid::Position)

/*-------------------- Free Functions (some of them used by NiceMC or pueoAnalysisFramework) --------------------*/

Double_t Geoid::getGeoidRadiusAtCosTheta(Double_t cosTheta) {
  /**
   * I discovered an approximately ~0.3 meter discrepancy at the poles between
   * methods setting lon/lat/alt=0 and getGeoidRadiusAtCosTheta.
   * Call this function with higherOrderCorrection = false to restore previous behaviour.
   */
  return GEOID_MIN*GEOID_MAX/TMath::Sqrt(GEOID_MIN*GEOID_MIN-(GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*cosTheta*cosTheta);
}

Double_t Geoid::getGeoidRadiusAtLatitude(Double_t latitude) {
  Position v;
  v.SetLonLatAlt(0, latitude, 0);
  return getGeoidRadiusAtCosTheta(v.CosTheta());
}

Double_t Geoid::getGeoidRadiusAtTheta(Double_t theta) {
  return getGeoidRadiusAtCosTheta(TMath::Cos(theta));
}

void Geoid::getCartesianCoords(Double_t lat, Double_t lon, Double_t alt, Double_t p[3]){

  // see page 71 onwards of https://web.archive.org/web/20120118224152/http://mercator.myzen.co.uk/mercator.pdf  
  lat *= TMath::DegToRad();
  lon *= TMath::DegToRad();
  //calculate x,y,z coordinates

  // Double_t C2 = pow(TMath::Cos(lat)*TMath::Cos(lat)+(1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR)*TMath::Sin(lat)*TMath::Sin(lat),-0.5);
  Double_t C2 = 1./TMath::Sqrt(TMath::Cos(lat)*TMath::Cos(lat)+(1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR)*TMath::Sin(lat)*TMath::Sin(lat));
  Double_t Q2 = (1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR)*C2;

  // Swapping x/y and inverting z
  p[1]=(R_EARTH*C2+alt)*TMath::Cos(lat)*TMath::Cos(lon);
  p[0]=(R_EARTH*C2+alt)*TMath::Cos(lat)*TMath::Sin(lon);
  p[2]=-(R_EARTH*Q2+alt)*TMath::Sin(lat);

}

void Geoid::getLatLonAltFromCartesian(const Double_t p[3], Double_t &lat, Double_t &lon, Double_t &alt){


// swapping x,y and inverting z
  Double_t x=p[1]; 
  Double_t y=p[0];
  Double_t z=-p[2];

  // see page 71 onwards of https://web.archive.org/web/20120118224152/http://mercator.myzen.co.uk/mercator.pdf  

  constexpr Double_t cosaeSq=(1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR);
  
  const Double_t lonVal   = TMath::ATan2(y,x);
  const Double_t xySq     = TMath::Sqrt(x*x+y*y);
  const Double_t tanPsit  = z/xySq;
  Double_t latGuess = TMath::ATan(tanPsit/cosaeSq);
  Double_t nextLat  = latGuess;
  Double_t geomBot  = R_EARTH*R_EARTH*xySq;

  const double deltaLatCloseEnough = 1e-6; // this corresponds to < 1m at the equator
  do {
    latGuess=nextLat;
    Double_t N      = R_EARTH/TMath::Sqrt(cos(latGuess)*cos(latGuess)+cosaeSq*sin(latGuess)*sin(latGuess));
    Double_t top    = (R_EARTH*R_EARTH*z + (1-cosaeSq)*cosaeSq*TMath::Power(N*TMath::Sin(latGuess),3));
    Double_t bottom = geomBot-(1-cosaeSq)*TMath::Power(N*TMath::Cos(latGuess),3);
    nextLat = TMath::ATan(top/bottom);
    // std::cout << latGuess << "\t" << nextLat << "\n";
  } while(TMath::Abs(nextLat-latGuess) > deltaLatCloseEnough);
  latGuess=nextLat;

  Double_t N = R_EARTH/TMath::Sqrt(cos(latGuess)*cos(latGuess)+cosaeSq*sin(latGuess)*sin(latGuess));
  Double_t height=(xySq/TMath::Cos(nextLat))-N;
  
  lat = latGuess*TMath::RadToDeg();
  lon = lonVal*TMath::RadToDeg();
  alt = height;
}

Double_t Geoid::getDistanceToCentreOfEarth(Double_t lat)
{
  Position v;
  v.SetLonLatAlt(0, lat, 0);
  return v.Mag();
}

  /**
   * Convert longitude and latitude to easting and northing using the geoid model
   *
   * @param lon is the longitude in degrees
   * @param lat is the latitude in degrees
   * @param easting in meters
   * @param northing in meters
   */
  void Geoid::LonLatToEastingNorthing(Double_t lon,Double_t lat,Double_t &easting,Double_t &northing){
    Double_t lon_rad = lon * TMath::DegToRad(); //convert to radians
    Double_t lat_rad = -lat * TMath::DegToRad();
    const double R_factor = scale_factor*c_0 * pow(( (1 + eccentricity*sin(lat_rad)) / (1 - eccentricity*sin(lat_rad)) ),eccentricity/2) * tan((TMath::Pi()/4) - lat_rad/2);
    easting = R_factor * sin(lon_rad);///(x_max-x_min);
    northing = R_factor * cos(lon_rad);///(y_max-y_min);
  }

  /**
   * Convert from easting/northing to longitude and latitude
   *
   * @param easting in meters
   * @param northing in meters
   * @param lon is the longitude
   * @param lat is the latitude
   */
  void Geoid::EastingNorthingToLonLat(Double_t easting,Double_t northing,Double_t &lon,Double_t &lat){

    double lon_rad = atan2(easting,northing);
    lon = lon_rad * TMath::RadToDeg();
    double R_factor = sqrt(easting*easting+northing*northing);
    double isometric_lat = (TMath::Pi()/2) - 2*atan(R_factor/(scale_factor*c_0));
    lat = isometric_lat + a_bar*sin(2*isometric_lat) + b_bar*sin(4*isometric_lat) + c_bar*sin(6*isometric_lat) + d_bar*sin(8*isometric_lat);
    lat =  -lat*TMath::RadToDeg(); //convert to degrees, with -90 degrees at the south pole
    return;
  }


std::ostream& operator<<(std::ostream& os, const TVector3& v){  
  os << "(" << v.X() << "," << v.Y() << "," << v.Z() << ")";
  return os;
}

/*-------------------- End of Free Functions  --------------------*/

/** 
 * @brief Custom streamer for ROOT to handle the laziness of coordinate conversions.
 * 
 * Since phi/theta, lon/lat/alt and easting/northing are calculated lazily they might
 * not be correct at the moment an object is written to disk if we used the default streamer.
 * Here we implement our own streamer which injects the proper coordinate conversions before
 * writing and initialises "AtLast" members (which aren't stored) to values that imply the
 * stored values are correct.
 * 
 * @param R__b the ROOT buffer
 */
void Geoid::Position::Streamer(TBuffer &R__b){
  if (R__b.IsReading()) {
    Position::Class()->ReadBuffer(R__b, this);

    fCartAtLastAngleCalc[0] = X();
    fCartAtLastAngleCalc[1] = Y();
    fCartAtLastAngleCalc[2] = Z();

    fCartAtLastGeoidCalc[0] = X();
    fCartAtLastGeoidCalc[1] = Y();
    fCartAtLastGeoidCalc[2] = Z();

    fLonLatAtLastEastNorthCalc[0] = longitude;
    fLonLatAtLastEastNorthCalc[1] = latitude;
  }
  else {
    // force all cached values to be correct *before* write
    updateAnglesFromCartesian();
    updateGeoidFromCartesian();
    updateEastingNorthingFromLonLat();

    Position::Class()->WriteBuffer(R__b, this);
  }
}

/*-------------------- PUBLIC Class Method Defintions --------------------*/

Double_t Geoid::Position::Longitude() const {
  updateGeoidFromCartesian();
  return longitude;
}

Double_t Geoid::Position::Latitude() const {
  updateGeoidFromCartesian();
  return latitude;
}

Double_t Geoid::Position::Altitude() const {
  updateGeoidFromCartesian();
  return altitude;
}

Double_t Geoid::Position::Easting() const {
  updateEastingNorthingFromLonLat();
  return easting;
}

Double_t Geoid::Position::Northing() const {
  updateEastingNorthingFromLonLat();
  return northing;
}

Double_t Geoid::Position::Theta() const {
  updateAnglesFromCartesian();
  return theta;
}

Double_t Geoid::Position::Phi() const {
  updateAnglesFromCartesian();
  return phi;
}

Double_t Geoid::Position::EllipsoidSurface() const {
  return getGeoidRadiusAtCosTheta(CosTheta());
}

void Geoid::Position::SetLongitude(double lon) {
  if(lon > 180){
    lon -= 360;
  }
  longitude = lon;
  updateCartesianFromGeoid();
}

void Geoid::Position::SetLatitude(double lat){
  latitude = lat;
  updateCartesianFromGeoid();
}

void Geoid::Position::SetAltitude(double alt) {
  altitude = alt;
  updateCartesianFromGeoid();
}

void Geoid::Position::SetLonLatAlt(double lon, double lat, double alt) {
  latitude = lat;
  altitude = alt;
  SetLongitude(lon);
}

void Geoid::Position::SetEasting(double east) {
  northing = Northing();
  easting = east;
  updateLonLatFromEastingNorthing(false);
}

void Geoid::Position::SetNorthing(double north) {    
  easting = Easting();
  northing = north;
  updateLonLatFromEastingNorthing(false);
}

void Geoid::Position::SetEastingNorthing(double east, double north) {
  easting = east;
  northing = north;
  updateLonLatFromEastingNorthing(true);
}

void Geoid::Position::SetEastingNorthingAlt(double east, double north, double alt) {
  easting = east;
  northing = north;
  altitude = alt;
  updateLonLatFromEastingNorthing(false);
}

Geoid::Pole Geoid::Position::nearerPole() const {
  return __getPole(Z());
}

double Geoid::Position::surfaceZ(Pole pole){
  // ellipse defined by: p^{2}/(geoid_max^{2}) + z^{2}/(geoid_min^{2}) = 1
  const double pSq = X()*X() + Y()*Y(); //lateral width of the geoid
  const double zSq = GEOID_MIN*GEOID_MIN*(1 - pSq/(GEOID_MAX*GEOID_MAX));
  if(zSq < 0){
    std::cerr <<  "Error in" << __PRETTY_FUNCTION__ << " can't find z if outside Geoid in x/y plane! "
  << " x = " << X() << ", y = " << Y() << ", GEOID_MAX = " << GEOID_MAX << std::endl;
    return TMath::QuietNaN();
  }
  else {
    return __signOfZ(pole)*TMath::Sqrt(zSq);
  }
}

double Geoid::Position::surfaceZ(){
  return surfaceZ(nearerPole()); 
}

void Geoid::Position::moveToGeoidZ(Pole pole){
  SetZ(surfaceZ(pole));
}

void Geoid::Position::moveToGeoidZ(){
  moveToGeoidZ(nearerPole());
}

Double_t Geoid::Position::Distance(const Position& p2) const{
  return (*this - p2).Mag();
}

/*-------------------- PRIVATE Class Method Defintions --------------------*/

void Geoid::Position::copyState(const Position& other){
  SetXYZ(other.X(), other.Y(), other.Z());// maybe redundant, but oh well
  longitude = other.longitude;
  latitude = other.latitude;

  altitude = other.altitude;
  theta = other.theta;
  phi = other.phi;
  easting = other.easting;
  northing = other.northing;
  for(size_t i=0; i < fCartAtLastGeoidCalc.size(); i++){
    fCartAtLastGeoidCalc[i] = other.fCartAtLastGeoidCalc[i];
  }
  for(size_t i=0; i < fCartAtLastAngleCalc.size(); i++){
    fCartAtLastAngleCalc[i] = other.fCartAtLastAngleCalc[i];
  }
  for(size_t i=0; i < fLonLatAtLastEastNorthCalc.size(); i++){
    fLonLatAtLastEastNorthCalc[i] = other.fLonLatAtLastEastNorthCalc[i];
  }    
}

void Geoid::Position::updateCartesianFromGeoid() {
  // always called after any lon/lat/alt has been updated
  Geoid::getCartesianCoords(latitude, longitude, altitude, fCartAtLastGeoidCalc.data());
  SetXYZ(fCartAtLastGeoidCalc[0], fCartAtLastGeoidCalc[1], fCartAtLastGeoidCalc[2]);
}

void Geoid::Position::updateGeoidFromCartesian() const {
  // called when Longitude(), Latitude(), Altitude() is requested
  if(X() != fCartAtLastGeoidCalc[0] ||
     Y() != fCartAtLastGeoidCalc[1] ||
     Z() != fCartAtLastGeoidCalc[2]){

    // Then we've moved, so must recalculate lon, lat alt;
    GetXYZ(fCartAtLastGeoidCalc.data());
    Geoid::getLatLonAltFromCartesian(fCartAtLastGeoidCalc.data(), latitude, longitude, altitude);
  }  
}

void Geoid::Position::updateAnglesFromCartesian() const {

  bool xDirty = X() != fCartAtLastAngleCalc[0];
  bool yDirty = Y() != fCartAtLastAngleCalc[1];
  
  if(xDirty || yDirty){
    phi = TVector3::Phi();
    fCartAtLastAngleCalc[0] = X();
    fCartAtLastAngleCalc[1] = Y();
  }

  if(xDirty || yDirty || Z() != fCartAtLastAngleCalc[2]){
    theta = TVector3::Theta();
    // if x or y was dirty, already stored them in fCartAtLastAngleCalc
    fCartAtLastAngleCalc[2] = Z();
  }
}

void Geoid::Position::updateEastingNorthingFromLonLat() const {

  if(longitude != fLonLatAtLastEastNorthCalc[0] ||
     latitude != fLonLatAtLastEastNorthCalc[1]){

    fLonLatAtLastEastNorthCalc[0] = Longitude();
    fLonLatAtLastEastNorthCalc[1] = Latitude();
    
    LonLatToEastingNorthing(fLonLatAtLastEastNorthCalc[0],
			    fLonLatAtLastEastNorthCalc[1],
			    easting, northing);
    
  }
}

void Geoid::Position::updateLonLatFromEastingNorthing(bool mustRecalcuateAltitudeFirst) {
  // Since we are going to eventually update cartesian from lon/lat/alt
  // we must make sure altitude is up to date...
  //
  double lon, lat;
  double alt = mustRecalcuateAltitudeFirst ? Altitude() : altitude;
  
  EastingNorthingToLonLat(easting, northing, lon, lat);
  SetLonLatAlt(lon, lat, alt);
}

int Geoid::Position::__signOfZ(const Geoid::Pole& pole) const {
  switch(pole){
    case Geoid::Pole::North: return -1;
  default:
    case Geoid::Pole::South: return 1;
  }
}

Geoid::Pole Geoid::Position::__getPole(double z) const{
  Geoid::Pole p = z >= 0 ? Geoid::Pole::South : Geoid::Pole::North;
  return p;
}
