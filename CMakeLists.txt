cmake_minimum_required(VERSION 3.19.0)
project(AntarcticaRoot VERSION 1.0.0)

file(GLOB SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx")
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
                                  
find_package(RootFftwWrapper CONFIG REQUIRED)
find_package(ROOT CONFIG REQUIRED COMPONENTS HistPainter Physics Ged)
find_package(ZLIB REQUIRED) # TODO: figure out if zlib is required

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

option(USE_HEALPIX
  "Enable Healpix Segmentation, requires healpix to be installed on your system in a place that can be found easily."
  OFF)

if(USE_HEALPIX)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  find_package(HealPix REQUIRED)
  # TODO: this preprocessor definition doesn't seem to be implemented in the source code?
  target_compile_definitions(${PROJECT_NAME} PRIVATE USE_HEALPIX)
  # Are these really necessary?!?? 
  target_link_libraries (${PROJECT_NAME} PRIVATE healpix_cxx cxxsupport sharp fftpack c_utils cfitsio)
endif()


#================================================================================================
#                                       BUILDING
#================================================================================================

target_include_directories(${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

# turn on warnings only in debug mode
target_compile_options(${PROJECT_NAME} PRIVATE
  $<$<CONFIG:Debug>:-Wall -Wextra>
)

# TODO: probably harmless, but do downstream projects need this flag?
target_compile_definitions(${PROJECT_NAME} PUBLIC USE_FFT_TOOLS)

target_link_libraries(${PROJECT_NAME} PUBLIC  RootFftwWrapper::RootFftwWrapper
                                      PRIVATE ZLIB::ZLIB
                                      PRIVATE ROOT::HistPainter ROOT::Ged ROOT::Physics)


set(DICTNAME G__${PROJECT_NAME})
file(GLOB HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
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
# TODO: not sure if pcm files should be installed, or how (https://root-forum.cern.ch/t/cmake-not-installing-pcm-and-rootmap/50997)

# installing the target file
install(EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE AntarcticaRoot::
    DESTINATION lib/cmake/${PROJECT_NAME}
)

# installing datafiles (ie. stuff in data/)
install(DIRECTORY data/ DESTINATION share/anitaCalib)
    # note: trailing "/" is important here
