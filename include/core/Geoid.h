#ifndef GEOID_MODEL_H
#define GEOID_MODEL_H

#include "TVector3.h"
#include "TMath.h"
#include "TObject.h"
#include "TBuffer.h"

#include <iostream>

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

  enum class Pole {North,South}; // for choosing solutions for Geoid z as a function of x,y
  /**
   * Ellipsoid Constants
   */
 // parameters of geoid model
  constexpr double FLATTENING_FACTOR = (1./298.257223563);
  constexpr double GEOID_MAX = 6.378137E6;
  constexpr double R_EARTH = GEOID_MAX;
  constexpr double GEOID_MIN = GEOID_MAX*(1 - FLATTENING_FACTOR); // parameters of geoid model
  /**
   * Variables for conversion between polar stereographic coordinates and lat/lon.
   * i.e. Easting/Northing from Longitude/Latitude
   * Conversion equations from ftp://164.214.2.65/pub/gig/tm8358.2/TM8358_2.pdf  
   */
  constexpr double scale_factor=0.97276901289;
  constexpr double ellipsoid_inv_f = 1./FLATTENING_FACTOR;
  constexpr double ellipsoid_b = R_EARTH*(1-(1/ellipsoid_inv_f));
  const double eccentricity = sqrt((1/ellipsoid_inv_f)*(2-(1/ellipsoid_inv_f)));
  const double a_bar = pow(eccentricity,2)/2 + 5*pow(eccentricity,4)/24 + pow(eccentricity,6)/12 + 13*pow(eccentricity,8)/360;
  const double b_bar = 7*pow(eccentricity,4)/48 + 29*pow(eccentricity,6)/240 + 811*pow(eccentricity,8)/11520;
  const double c_bar = 7*pow(eccentricity,6)/120 + 81*pow(eccentricity,8)/1120;
  const double d_bar = 4279*pow(eccentricity,8)/161280;
  const double c_0 = (2*R_EARTH / sqrt(1-pow(eccentricity,2))) * pow(( (1-eccentricity) / (1+eccentricity) ),eccentricity/2);

  // namespace scoped free functions; implementation in Geoid.cxx
  Double_t getGeoidRadiusAtCosTheta(Double_t cosTheta);
  Double_t getGeoidRadiusAtLatitude(Double_t lat);
  Double_t getGeoidRadiusAtTheta(Double_t theta);
  void getCartesianCoords(Double_t lat, Double_t lon, Double_t alt, Double_t p[3]);
  void getLatLonAltFromCartesian(const Double_t p[3], Double_t &lat, Double_t &lon, Double_t &alt);
  Double_t getDistanceToCentreOfEarth(Double_t lat);
  
  void LonLatToEastingNorthing(Double_t lon,Double_t lat,Double_t &easting,Double_t &northing);
  void EastingNorthingToLonLat(Double_t easting,Double_t northing,Double_t &lon,Double_t &lat);

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
    Position & operator=(const Position & p) 
    {
      copyState(p); 
      return *this; 
    }
    virtual ~Position(){;}

    template <class T> Position(const T& t);
    template <class T> Position(const T* t);

    /**
     * Longitude/Latitude/Altitude Getter functions
     */
    Double_t Longitude() const;
    Double_t Latitude() const;
    Double_t Altitude() const;
    Double_t Easting() const;
    Double_t Northing() const;
    Double_t Theta() const;
    Double_t Phi() const;
    Double_t EllipsoidSurface() const;

    inline void GetLonLatAlt(Double_t& lon, Double_t& lat, Double_t& alt) const;
    template <class T> void GetLonLatAlt(T& t) const;
    template <class T> void GetLonLatAlt(T* t) const;

    void SetLongitude(double longitude);
    void SetLatitude(double latitude);
    void SetAltitude(double altitude);
    void SetLonLatAlt(double lon, double lat, double alt);

    void SetEasting(double easting);
    void SetNorthing(double northing);
    void SetEastingNorthing(double easting, double northing);
    void SetEastingNorthingAlt(double easting, double northing, double alt);

    Pole nearerPole() const;

    /** 
     * Find the value of z on Geoid surface given the values X(), Y()
     * Note: There are two solutions here, pole chooses which.
     * 
     * @param pole solution of z to find (Pole::North or Pole::South), default is Pole::South (near Antarctica).
     * 
     * @return The value of Z
     */
    double surfaceZ(Pole pole);

    /** 
     * Find the value of z on Geoid surface given the values X(), Y()
     * Note: There are two solutions here, the nearer solution (based  on the value of Z()) is chosen.
     * @see surfaceZ(Pole pole)
     * 
     * @return The value of Z
     */
    double surfaceZ();

    /** 
     * Change z so that we are on the surface at this X(),  Y()
     * @see surfaceZ(Pole pole)
     * 
     * @param signZ pole chose the sign of Z()
     */
    void moveToGeoidZ(Pole pole);    

    /** 
     * Change z so that we are on the surface at this X(),  Y()
     * @see surfaceZ(Pole pole)
     * 
     */
    void moveToGeoidZ();

    Double_t Distance(const Position& p2) const; ///@todo make this better

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

    // can't make const as the cartesian x,y,z are represented in TVector3
    // and accessors aren't virtual and so can't be overloaded. 
    void updateCartesianFromGeoid();
    void updateGeoidFromCartesian() const;
    void updateAnglesFromCartesian() const;
    void updateEastingNorthingFromLonLat() const;
    void updateLonLatFromEastingNorthing(bool mustRecalcuateAltitudeFirst);

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

    // See namespace comments on the coordinate system for context
    int __signOfZ(const Pole& pole) const;

    // See namespace comments on the coordinate system for context
    Pole __getPole(double z) const;

}; // end class def


  template <class T> Position::Position(const T& t){
    SetLonLatAlt(t.longitude,  t.latitude, t.altitude);
    updateCartesianFromGeoid();
  }

  template <class T> Position::Position(const T* t) : Position(*t){;}

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
