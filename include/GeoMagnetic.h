#ifndef GEOMAGNETIC_H
#define GEOMAGNETIC_H

#include <vector>
#include "TArrow.h"
#include "Geoid.h"
#include "TVector3.h"

/** 
 * @namespace GeoMagentic Functions to calculate the Earth's geo-magnetic field in our Geoid coordinate system
 *
 * Currently uses the IGRF model
 * 
 */

class TCanvas;

namespace GeoMagnetic{

  
  // I want these in the GeoMagnetic namespace so I can just call them
  // I read xMax for a 1e19 eV proton off a plot in an Auger paper
  const double xMaxP19 = 0.8e4; // kg / m^{2}
  // The minimum for ANITA, an 1e18 Fe, would be more like:
  const double xMaxFe18 = 0.65e4; // kg / m^{2}


  /** 
   * Get the GeoMagnetic field at a certain time and position
   * @see Geoid.h
   * 
   * @param unixTime the time
   * @param position the place, in the Geoid cartesian coordinate system
   * 
   * @return the GeoMagnetic field, represented in the Geoid cartesian coordinate system
   */
  TVector3 getField(UInt_t unixTime, const Geoid::Position& position);
  

  
  /** 
   * @class FieldPoint Drawable pair of position/field values.
   */
  class FieldPoint : public TArrow {
  public:
    FieldPoint (UInt_t unixTime, double lon, double lat, double alt);
    FieldPoint (UInt_t unixTime, const Geoid::Position& position);
    virtual ~FieldPoint(){}
    virtual void Draw(Option_t* opt = "");

    inline double posX() const {return fPosition.X();}
    inline double posY() const {return fPosition.Y();}
    inline double posZ() const {return fPosition.Z();}
    inline double posR() const {return fPosition.Mag();}
    inline double posTheta() const {return fPosition.Theta();}
    inline double posPhi() const {return fPosition.Phi();}
    inline double componentX() const {return fField.X();}
    inline double componentY() const {return fField.Y();}
    inline double componentZ() const {return fField.Z();}

    const TVector3& field(){return fField;}
    const Geoid::Position& position(){return fPosition;}  
    UInt_t getUnixTime(){return fUnixTime;}
  
 private:
    double fDrawScaleFactor; ///< Conversion factor (metres/nT) to draw on an AntarcticaBackground
    TVector3 fField; ///< Cartesian components of the geomagnetic field in the Geoid coordinate system
    Geoid::Position fPosition; ///< Location of the magnetic field in the Geoid coordinate system
    UInt_t fUnixTime; ///< Time at which to calculate the field
    inline void calculateFieldAtPosition(){
      fField = getField(fUnixTime, fPosition);
      SetLineColor(fField.Z() > 0 ? kRed : kBlue);
    }
  };

  const double n_air = 1;
  const double n_ice = 1.31; // might need to check this

  TVector3 specularReflection(const TVector3& reflectionPointToSource, const TVector3& surfaceNormal);
  TVector3 fresnelReflection(const TVector3& sourceToReflection, const TVector3& surfaceNormal, TVector3& electricFieldVec, double n1=n_air, double n2=n_ice);
  TCanvas* plotFresnelReflection();

  TCanvas* plotFieldAtAltitude(UInt_t unixTime, double altitude);
  TCanvas* plotAtmosphere();

  double getAtmosphericDensity(double altitude);
  TVector3 getXMaxPosition(const TVector3& initialPosition, const TVector3& cosmicRayDirection, double xMax);
  TVector3 getInitialPosition(const TVector3& destination, const TVector3& destinationToSource);
  void setDebug(bool db);

}
#endif
