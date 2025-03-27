cmake_minimum_required(VERSION 3.23.0)
project(AntarcticaRoot VERSION 1.0.0)
add_library(${PROJECT_NAME} SHARED)

#================================================================================================
#                                       HOUSEKEEPING
#================================================================================================

if(NOT CMAKE_BUILD_TYPE) 
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE) 
endif()

find_package(RootFftwWrapper CONFIG REQUIRED)
find_package(ROOT CONFIG REQUIRED COMPONENTS HistPainter Physics Ged)
find_package(ZLIB REQUIRED) # TODO: figure out if zlib is required

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF) # turns -std=gnu++17 -> -std=c++17
#================================================================================================
#                                       OPTIONS
#================================================================================================

option(USE_GEOGRAPHIC_LIB "Use GeographicLib for geodesic calculations" OFF)
# note: For Fedora Linux users, run `dnf install GeographicLib-devel` first
if(USE_GEOGRAPHIC_LIB)
  find_package(GeographicLib REQUIRED)
  target_compile_definitions(${PROJECT_NAME} PRIVATE USE_GEOGRAPHIC_LIB)
  target_link_libraries(${PROJECT_NAME} PRIVATE ${GeographicLib_LIBRARIES} )
endif()

#================================================================================================
#                                       BUILDING
#================================================================================================

set(HEADER_FILES 
  include/AntarcticAtmosphere.h
  include/AntarcticaBackground.h
  include/AntarcticaGeometry.h
  include/BaseList.h
  include/GeoMagnetic.h
  include/Geoid.h
  include/RampdemReader.h
  include/RefractionModel.h
  include/SkyMap.h
  include/TArrowAntarctica.h
  include/TGraphAntarctica.h
  include/TH2DAntarctica.h
)
target_sources(${PROJECT_NAME} PRIVATE
  src/AntarcticaAtmosphere.cxx
  src/AntarcticaBackground.cxx
  src/AntarcticaGeometry.cxx
  src/BaseList.cxx
  src/GeoMagnetic.cxx
  src/Geoid.cxx
  src/RampdemReader.cxx
  src/RefractionModel.cxx
  src/SkyMap.cxx
  src/TArrowAntarctica.cxx
  src/TGraphAntarctica.cxx
  src/TH2DAntarctica.cxx
)
target_sources(${PROJECT_NAME} PUBLIC
  FILE_SET HEADERS
  BASE_DIRS include     # <-- Everything in BASE_DIRS is available during the build
  FILES ${HEADER_FILES} # <-- But only FILES will be installed, although these two are currently the same
)

# use warning flags in RelWithDebInfo mode, which is set to be the default mode
target_compile_options(${PROJECT_NAME} PRIVATE $<$<CONFIG:RelWithDebInfo>:-Wall -Wextra>)

target_link_libraries(${PROJECT_NAME} PUBLIC  RootFftwWrapper::RootFftwWrapper
                                      PRIVATE ZLIB::ZLIB
                                      PRIVATE ROOT::HistPainter ROOT::Ged ROOT::Physics)

set(DICTNAME G__${PROJECT_NAME})
root_generate_dictionary(${DICTNAME} ${HEADER_FILES} MODULE ${PROJECT_NAME} LINKDEF LinkDef.h)

#================================================================================================
#                                       INSTALLING
#================================================================================================

# Generate AntarcticaRootConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Install config files
include(GNUInstallDirs)
install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)

# Install this library and its headers as well as exporting as a CMake Target
install(
  TARGETS ${PROJECT_NAME} 
  EXPORT ${PROJECT_NAME}Targets 
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

# installing CERN ROOT's _rdict.pcm file
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Generate a Targets.cmake file and directly install it
install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE AntarcticaRoot::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)

# installing datafiles (ie. stuff in data/)
install(DIRECTORY data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/anitaCalib)
    # note: trailing "/" is important here
