#ifndef _ANTARCTIC_ATMOSPHERE_H 
#define _ANTARCTIC_ATMOSPHERE_H 

/* Various routines for atmospheric modeling
 *
 * Some things depend on GeographicLib 
 *
 **/


class TGraph; 
class Adu5Pat; 
namespace AntarcticAtmosphere
{

  /** Geoid model for converting from WGS84 to MSL and vice versa 
   * These require GeographicLib and the corresponding data files. See GeographicLib documentation for more details. 
   * */
  enum Geoid
  {
    EGM96_5, 
    EGM2008_1
  };
    
  /** Convert from mean sea level to WGS84 altitude. Requires GeographicLib */ 
  double MSLtoWGS84( double h, double lat, double lon, Geoid g = EGM96_5)  ; 

  /** Convert from WGS84 altitude to height above mean sea level. Requires GeographicLib*/ 
  double WGS84toMSL( const Adu5Pat * pat,  Geoid g = EGM96_5) ; 


  /** Atmospheric Parameters. 
   * */ 
  struct Pars
  {
    double rho;  //kg / m^3 
    double P;    //mb 
    double T;    //Kelvin 
    double N;    //refractivity in units of (n-1) * 1e6 
  }; 

  enum Par
  {
    DENSITY,
    PRESSURE,
    TEMPERATURE,
    REFRACTIVITY
  };



  class AtmosphericModel 
  {
    public: 
      virtual int computeAtmosphere(double h, Pars * p) const= 0; 
      virtual double get(double h, Par p) const; 
      TGraph * makeGraph(double hmin, double hmax, int nh, Par P) const; 
      virtual ~AtmosphericModel() { ; }
  }; 


  class StandardUS : public AtmosphericModel
  {
    public: 
      StandardUS(double sea_level_T_kelvin = 265, double sea_level_P_mbar = 970)
      {
        sea_level_T = sea_level_T_kelvin; 
        sea_level_P = sea_level_P_mbar; 
      }
      virtual ~StandardUS() { ;} 

      virtual int computeAtmosphere(double h, Pars * p) const; 
    private: 
      double sea_level_T; 
      double sea_level_P; 
  };


}



#endif