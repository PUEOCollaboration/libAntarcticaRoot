list(APPEND HEADER_FILES core/Geoid.h core/RampdemReader.h)

if(NOT CORE_ONLY)
  list(APPEND HEADER_FILES
    misc/AntarcticAtmosphere.h
    misc/AntarcticaBackground.h
    misc/AntarcticaGeometry.h
    misc/BaseList.h
    misc/GeoMagnetic.h
    misc/RefractionModel.h
    misc/SkyMap.h
    misc/TArrowAntarctica.h
    misc/TGraphAntarctica.h
    misc/TH2DAntarctica.h
  )

  target_compile_definitions(${PROJECT_NAME} PRIVATE INCLUDE_MISC)
endif()

target_sources(${PROJECT_NAME} PUBLIC
  FILE_SET HEADERS
  BASE_DIRS core misc   # <-- Everything in BASE_DIRS is available during the build
  FILES ${HEADER_FILES} # <-- But only FILES will be installed
)

set(DICTNAME G__${PROJECT_NAME})
root_generate_dictionary(${DICTNAME} ${HEADER_FILES} MODULE ${PROJECT_NAME} LINKDEF LinkDef.h)

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
