target_sources(${PROJECT_NAME} PRIVATE core/Geoid.cxx core/RampdemReader.cxx)

if(NOT CORE_ONLY)
  target_sources(${PROJECT_NAME} PRIVATE
    misc/AntarcticaAtmosphere.cxx
    misc/AntarcticaBackground.cxx
    misc/AntarcticaGeometry.cxx
    misc/BaseList.cxx
    misc/GeoMagnetic.cxx
    misc/RefractionModel.cxx
    misc/SkyMap.cxx
    misc/TArrowAntarctica.cxx
    misc/TGraphAntarctica.cxx
    misc/TH2DAntarctica.cxx
  )
endif()
