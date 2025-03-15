cmake_minimum_required(VERSION 3.19.0)
project(AntarcticaRoot VERSION 1.0.0)

# note: do not glob files with wildcards
set(HEADER_FILES 
  ${CMAKE_CURRENT_SOURCE_DIR}/include/AntarcticAtmosphere.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/AntarcticaBackground.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/AntarcticaGeometry.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/BaseList.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/GeoMagnetic.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/Geoid.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/RampdemReader.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/RefractionModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/SkyMap.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/TArrowAntarctica.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/TGraphAntarctica.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/TH2DAntarctica.h
)
add_library(${PROJECT_NAME} SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/src/AntarcticaAtmosphere.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/AntarcticaBackground.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/AntarcticaGeometry.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BaseList.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/GeoMagnetic.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Geoid.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RampdemReader.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RefractionModel.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/SkyMap.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TArrowAntarctica.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TGraphAntarctica.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TH2DAntarctica.cxx
)
                                  
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

target_include_directories(${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
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

install(FILES ${HEADER_FILES} DESTINATION include/${PROJECT_NAME})

# Config files and such
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Generating AntarcticaRootConfigVersion.cmake
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
# Installing config files
install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "lib/cmake/${PROJECT_NAME}"
)

# installing the libraries
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets LIBRARY DESTINATION lib)

# installing CERN ROOT's _rdict.pcm file
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm"
  DESTINATION DESTINATION lib
)

# installing the target file
install(EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE AntarcticaRoot::
    DESTINATION lib/cmake/${PROJECT_NAME}
)

# installing datafiles (ie. stuff in data/)
install(DIRECTORY data/ DESTINATION share/anitaCalib)
    # note: trailing "/" is important here
