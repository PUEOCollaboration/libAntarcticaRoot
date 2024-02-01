#ifndef GEOID_MODEL_H
#define GEOID_MODEL_H

#include "TVector3.h"
#include "TMath.h"
#include <iostream>
#include "TObject.h"
#include "TBuffer.h"

/**
 * @namespace Geoid
 * @brief Get positions, radii, latitudes, longitudes, and other goodies when modelling the Earth
 * 
 * A note on Cartesian coordinates: We don't use the WGS84 convention!
 * 
 * The cartesian coordinate system in WGS84 has:
 * The origin at the center of mass of the Earth.
 * The +ve z-axis running through the north pole.
 * The +ve x-axis running through the prime merian (0 longitude) at the equator
 * The +ve y-axis is picked such that the coordinate system is right handed.
 *
 * ANITA and icemc, for historical reasons, do things a little differently.
 * The origin is in the same place but:
 * The +ve z-axis runs through the SOUTH pole.
 * Then the x-axis and y-axis are swapped relative to WGS84.
 * That x/y swap maintains a right handed coordinate system.
 * 
 * i.e. the +ve y-axis comes out of the Earth at the equator at 0 longitude.
 * The ANITA/icemc coordinate system is still right handed, and has the property
 * that +ve x-axis aligns with easting, +ve y-axis aligns with northing,
 * and quantities like elevation are more +ve in z for higher altitude in Antarctica.
 * I suppose it does make Cartesian plots of Antarctica a bit easier to look at.
 */

class TBuffer;

namespace Geoid {

  /**
   * Ellipsoid Constants
   */
 // parameters of geoid model
  static constexpr double FLATTENING_FACTOR = (1./298.257223563);
  static constexpr double GEOID_MAX = 6.378137E6;
  static constexpr double R_EARTH = GEOID_MAX;
  static constexpr double GEOID_MIN = GEOID_MAX*(1 - FLATTENING_FACTOR); // parameters of geoid model



  enum class Pole {North,South}; // for choosing solutions for Geoid z as a function of x,y
  inline int signOfZ(Pole pole){ // See namespace comments on the coordinate system for context
    switch(pole){
    case Pole::North: return -1;
    default:
    case Pole::South: return 1;
    }
  }
  inline Pole getPole(double z){ // See namespace comments on the coordinate system for context
    Pole p = z >= 0 ? Pole::South : Pole::North;
    return p;
  }


  

  inline Double_t getGeoidRadiusAtCosTheta(Double_t cosTheta);
  Double_t getGeoidRadiusAtLatitude(Double_t lat);
  inline Double_t getGeoidRadiusAtTheta(Double_t theta);
  void getCartesianCoords(Double_t lat, Double_t lon, Double_t alt, Double_t p[3]);
  void getLatLonAltFromCartesian(const Double_t p[3], Double_t &lat, Double_t &lon, Double_t &alt);
  Double_t getDistanceToCentreOfEarth(Double_t lat);

  
  /**
   * Variables for conversion between polar stereographic coordinates and lat/lon.
   * i.e. Easting/Northing from Longitude/Latitude
   * Conversion equations from ftp://164.214.2.65/pub/gig/tm8358.2/TM8358_2.pdf  
   */
  static constexpr double scale_factor=0.97276901289;
  static constexpr double ellipsoid_inv_f = 1./FLATTENING_FACTOR;
  static constexpr double ellipsoid_b = R_EARTH*(1-(1/ellipsoid_inv_f));
  static const double eccentricity = sqrt((1/ellipsoid_inv_f)*(2-(1/ellipsoid_inv_f)));
  static const double a_bar = pow(eccentricity,2)/2 + 5*pow(eccentricity,4)/24 + pow(eccentricity,6)/12 + 13*pow(eccentricity,8)/360;
  static const double b_bar = 7*pow(eccentricity,4)/48 + 29*pow(eccentricity,6)/240 + 811*pow(eccentricity,8)/11520;
  static const double c_bar = 7*pow(eccentricity,6)/120 + 81*pow(eccentricity,8)/1120;
  static const double d_bar = 4279*pow(eccentricity,8)/161280;
  static const double c_0 = (2*R_EARTH / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);

  inline void LonLatToEastingNorthing(Double_t lon,Double_t lat,Double_t &easting,Double_t &northing);
  inline void EastingNorthingToLonLat(Double_t easting,Double_t northing,Double_t &lon,Double_t &lat);

  /**
   * @class Position
   * 
   * @brief The ultimate way to represent a position on the Earth.
   * 
   * Very much in the spirit of https://xkcd.com/927/
   * 
   * This class is supposed to combine the best features of a TVector, an icemc::Position, 
   * and the proper Geoid transformations.
   * 
   * Features:
   * Lazy, behind the scenes conversion of x,y,z to lon, lat, alt, easting, northing
   * 
   */
  class Position : public TVector3 {

  public:

    Position(Double_t x=0, Double_t y=0) : TVector3(x, y, 0) {
      moveToGeoidZ();
    }
    Position(Pole pole) : TVector3(0, 0, 0) {
      moveToGeoidZ(pole);
    }

    Position(Double_t x, Double_t y, Double_t z) : TVector3(x, y, z) {};
    Position(const TVector3& v) : TVector3(v) {};
    Position(const Position& p)
      : TVector3(p) 
    {
      copyState(p);
    }
    virtual ~Position(){;}

    template <class T> Position(const T& t);
    template <class T> Position(const T* t);

    /**
     * Longitude/Latitude/Altitude Getter functions
     */
    inline Double_t Latitude() const;
    inline Double_t Longitude() const;
    inline Double_t Altitude() const;
    inline Double_t Easting() const;
    inline Double_t Northing() const;

    Double_t EllipsoidSurface() const;
    inline void GetLonLatAlt(Double_t& lon, Double_t& lat, Double_t& alt) const;
    template <class T> void GetLonLatAlt(T& t) const;
    template <class T> void GetLonLatAlt(T* t) const;

    inline void SetLongitude(double longitude);
    inline void SetLatitude(double latitude);
    inline void SetAltitude(double altitude);
    inline void SetLonLatAlt(double lon, double lat, double alt);

    inline void SetEasting(double easting);
    inline void SetNorthing(double northing);
    inline void SetEastingNorthing(double easting, double northing);
    inline void SetEastingNorthingAlt(double easting, double northing, double alt);

    inline Double_t Theta() const;
    inline Double_t Phi() const;

    inline Pole nearerPole() const {
      return getPole(Z());
    }


    /** 
     * Find the value of z on Geoid surface given the values X(), Y()
     * Note: There are two solutions here, pole chooses which.
     * 
     * @param pole solution of z to find (Pole::North or Pole::South), default is Pole::South (near Antarctica).
     * 
     * @return The value of Z
     */
    inline double surfaceZ(Pole pole);

    /** 
     * Find the value of z on Geoid surface given the values X(), Y()
     * Note: There are two solutions here, the nearer solution (based  on the value of Z()) is chosen.
     * @see surfaceZ(Pole pole)
     * 
     * @return The value of Z
     */
    inline double surfaceZ(){
      return surfaceZ(nearerPole());
    }

    /** 
     * Change z so that we are on the surface at this X(),  Y()
     * @see surfaceZ(Pole pole)
     * 
     * @param signZ pole chose the sign of Z()
     */
    inline void moveToGeoidZ(Pole pole);    

    /** 
     * Change z so that we are on the surface at this X(),  Y()
     * @see surfaceZ(Pole pole)
     * 
     */
    inline void moveToGeoidZ(){
      moveToGeoidZ(nearerPole());
    }



    inline Double_t Distance(const Position& p2) const; ///@todo make this better

  private:


    /**
     * How the class actually works:
     *
     * The fX, fY, fZ in the base TVector class are the definitive 
     * representation of the position. In other words, any non-const 
     * functions ALWAYS update the fX, fY, fZ immediately
     * (fX, fY, fZ are the variable names in TVector3.)
     * However, the longitude, latitude, altitude, and theta, phi,
     * and easting, northing are not updated, and are recalculated on
     * demand. When calculated they are stored and if fX, fY, fZ didn't change
     * The last answer will be returned again.
     *
     * The model for the object constness is that only changing the vector position
     * is non-const. Therefore the derived quantities (which are calculated and
     * stored on demand) are mutable.
     *
     * The derived coordiante caches are non-ROOT persistent,
     * (marked by comments after variable beginning with //!).
     * i.e. they won't be stored in a TTree.
     */
    void copyState(const Position& other);

    void updateGeoidFromCartesian() const;
    void updateAnglesFromCartesian() const;
    void updateEastingNorthingFromLonLat() const;

    // can't make const as the cartesian x,y,z are represented in TVector3
    // and accessors aren't virtual and so can't be overloaded. 
    void updateLonLatFromEastingNorthing(bool mustRecalcuateAltitudeFirst);
    void updateCartesianFromGeoid();

    mutable Double_t longitude = 0;
    mutable Double_t latitude  = 0;
    mutable Double_t altitude  = 0;
    mutable Double_t theta     = 0;
    mutable Double_t phi       = 0;
    mutable Double_t easting   = 0;
    mutable Double_t northing  = 0;

    ClassDef(Position, 1);
    
    mutable std::array<Double_t, 3> fCartAtLastGeoidCalc       = {-1, -1, -1};    //! fX, fY, fZ when the geoid was last updated. Is not stored!
    mutable std::array<Double_t, 3> fCartAtLastAngleCalc       = {0, 0, 0};       //! fX, fY, fZ when the angles were last updated. Is not stored!
    mutable std::array<Double_t, 3> fLonLatAtLastEastNorthCalc = {-9999, -9999};  //! Longitude(), Latitude() when the Easting/Northing were last calculated. Is not stored!

  };






















  template <class T> Position::Position(const T& t){
    SetLonLatAlt(t.longitude,  t.latitude, t.altitude);
    updateCartesianFromGeoid();
  }

  template <class T> Position::Position(const T* t) : Position(*t){;}


  inline void Position::SetLongitude(double lon) {
    if(lon > 180){
      lon -= 360;
    }
    longitude = lon;
    updateCartesianFromGeoid();
  }
  inline void Position::SetLatitude(double lat){
    latitude = lat;
    updateCartesianFromGeoid();
  }
  inline void Position::SetAltitude(double alt) {
    altitude = alt;
    updateCartesianFromGeoid();
  }
  inline void Position::SetLonLatAlt(double lon, double lat, double alt) {
    latitude = lat;
    altitude = alt;
    SetLongitude(lon);
  }


  inline void Position::SetEasting(double easting) {
    northing = Northing();
    easting = easting;
    updateLonLatFromEastingNorthing(false);
  }
  inline void Position::SetNorthing(double northing) {    
    easting = Easting();
    northing = northing;
    updateLonLatFromEastingNorthing(false);
  }
  inline void Position::SetEastingNorthing(double easting, double northing) {
    easting = easting;
    northing = northing;
    updateLonLatFromEastingNorthing(true);
  }
  inline void Position::SetEastingNorthingAlt(double easting, double northing, double alt) {
    easting = easting;
    northing = northing;
    altitude = alt;
    updateLonLatFromEastingNorthing(false);
  }



  inline Double_t Position::Latitude() const {
    updateGeoidFromCartesian();
    return latitude;
  }
  inline Double_t Position::Longitude() const {
    updateGeoidFromCartesian();
    return longitude;
  }
  inline Double_t Position::Altitude() const {
    updateGeoidFromCartesian();
    return altitude;
  }

  inline Double_t Position::Theta() const {
    updateAnglesFromCartesian();
    return theta;
  }
  inline Double_t Position::Phi() const {
    updateAnglesFromCartesian();
    return phi;
  }


  inline Double_t Position::Easting() const {
    updateEastingNorthingFromLonLat();
    return easting;
  }

  inline Double_t Position::Northing() const {
    updateEastingNorthingFromLonLat();
    return northing;
  }

  inline Double_t Position::EllipsoidSurface() const {
    return getGeoidRadiusAtCosTheta(CosTheta());
  }



  void Position::GetLonLatAlt(Double_t& lon, Double_t& lat, Double_t& alt) const {
    // call functions rather than direct access so cache is updated!
    lon = Longitude();
    lat = Latitude();
    alt = Altitude();
  }

  template <class T>
  void Position::GetLonLatAlt(T& t) const {
    GetLonLatAlt(t.longitude,  t.latitude,  t.altitude);
  }

  template <class T>
  void Position::GetLonLatAlt(T* t) const {
    GetLonLatAlt(*t);
  }




  inline double Position::surfaceZ(Pole pole){
    // ellipse defined by: p^{2}/(geoid_max^{2}) + z^{2}/(geoid_min^{2}) = 1
    const double pSq = X()*X() + Y()*Y(); //lateral width of the geoid
    const double zSq = GEOID_MIN*GEOID_MIN*(1 - pSq/(GEOID_MAX*GEOID_MAX));
    if(zSq < 0){
      std::cerr <<  "Error in" << __PRETTY_FUNCTION__ << " can't find z if outside Geoid in x/y plane! "
		<< " x = " << X() << ", y = " << Y() << ", GEOID_MAX = " << GEOID_MAX << std::endl;
      return TMath::QuietNaN();
    }
    else {
      return signOfZ(pole)*TMath::Sqrt(zSq);
    }
  }

  inline void Position::moveToGeoidZ(Pole pole){
    SetZ(surfaceZ(pole));
  }

  inline Double_t Position::Distance(const Position& p2) const{
    return (*this - p2).Mag();
  }



  inline void Position::copyState(const Position& other){
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



  inline Double_t getGeoidRadiusAtCosTheta(Double_t cosTheta) {
    /**
     * I discovered an approximately ~0.3 meter discrepancy at the poles between
     * methods setting lon/lat/alt=0 and getGeoidRadiusAtCosTheta.
     * Call this function with higherOrderCorrection = false to restore previous behaviour.
     */
    return GEOID_MIN*GEOID_MAX/TMath::Sqrt(GEOID_MIN*GEOID_MIN-(GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*cosTheta*cosTheta);
  }
  inline Double_t getGeoidRadiusAtTheta(Double_t theta) {
    return getGeoidRadiusAtCosTheta(TMath::Cos(theta));
  }



  

  
  /**
   * Convert longitude and latitude to easting and northing using the geoid model
   *
   * @param lon is the longitude in degrees
   * @param lat is the latitude in degrees
   * @param easting in meters
   * @param northing in meters
   */
  void LonLatToEastingNorthing(Double_t lon,Double_t lat,Double_t &easting,Double_t &northing){
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
  void EastingNorthingToLonLat(Double_t easting,Double_t northing,Double_t &lon,Double_t &lat){

    double lon_rad = atan2(easting,northing);
    lon = lon_rad * TMath::RadToDeg();
    double R_factor = sqrt(easting*easting+northing*northing);
    double isometric_lat = (TMath::Pi()/2) - 2*atan(R_factor/(scale_factor*c_0));
    lat = isometric_lat + a_bar*sin(2*isometric_lat) + b_bar*sin(4*isometric_lat) + c_bar*sin(6*isometric_lat) + d_bar*sin(8*isometric_lat);
    lat =  -lat*TMath::RadToDeg(); //convert to degrees, with -90 degrees at the south pole
    return;
  }
  

  
}

/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param v is the TVector3
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const TVector3& v);



#endif
