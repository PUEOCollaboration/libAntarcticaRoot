#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace BaseList;
#pragma link C++ class BaseList::base+;
#pragma link C++ class BaseList::path+;
#pragma link C++ class BaseList::abstract_base+;
#pragma link C++ class AntarcticaBackground+;
#pragma link C++ class TGraphAntarctica+;
#pragma link C++ class TArrowAntarctica+;
#pragma link C++ class TH2DAntarctica+;
#pragma link C++ class TProfile2DAntarctica+;
#pragma link C++ class SkyMap+;

#pragma link C++ class AntarcticSegmentationScheme;
#pragma link C++ class AntarcticCoord;
#pragma link C++ class StereographicGrid;
#pragma link C++ class PayloadParameters;
#pragma link C++ class CartesianSurfaceMap;
#pragma link C++ namespace AntarcticAtmosphere; 
#pragma link C++ class AntarcticAtmosphere::AtmosphericModel; 
#pragma link C++ class AntarcticAtmosphere::StandardUS; 
#pragma link C++ class AntarcticAtmosphere::ExponentialRefractivity; 
#pragma link C++ class AntarcticAtmosphere::ArtificialInversion;

#pragma link C++ namespace Geoid;
#pragma link C++ class Geoid::Position-;

#pragma link C++ namespace Refraction;
#pragma link C++ class Refraction::Model; 
#pragma link C++ class Refraction::PositionIndependentModel; 
#pragma link C++ class Refraction::PGFit; 
#pragma link C++ class Refraction::RaytracerSpherical; 
#pragma link C++ class Refraction::SphRay; 

#pragma link C++ namespace pueo;
#pragma link C++ class pueo::UsefulAttitude+;
#endif
