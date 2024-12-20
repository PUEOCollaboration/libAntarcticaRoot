cmake_minimum_required(VERSION 3.19.0)
project(AntarcticaRoot VERSION 1.0.0)

file(GLOB SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx")
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
                                  
find_package(ROOT CONFIG REQUIRED)# COMPONENTS MathMore MathCore RIO Hist Tree Net Minuit Gui HistPainter)
# find_package(FFTW REQUIRED)     # TODO: compiled just fine without these
# find_package(ZLIB REQUIRED)

#================================================================================================
#                                       OPTIONS
#================================================================================================

option(FORCE_OLD_GPP_ABI
  "Force old g++ ABI; this might be necessary if using new g++ with ROOT compiled with older g++ or other similar situations"
  OFF)

if(FORCE_OLD_GPP_ABI)
  target_compile_options(${PROJECT_NAME} PRIVATE -D_GLIBCXX_USE_CXX11_ABI=0)
endif()

option(USE_GEOGRAPHIC_LIB "Use GeographicLib for geodesic calculations" OFF)
# note: For Fedora Linux users, run `dnf install GeographicLib-devel` first
if(USE_GEOGRAPHIC_LIB)
  find_package(GeographicLib REQUIRED)
  target_compile_options(${PROJECT_NAME} PRIVATE -DUSE_GEOGRAPHIC_LIB)
  target_link_libraries (${PROJECT_NAME} PRIVATE ${GeographicLib_LIBRARIES} )
endif()

option(USE_HEALPIX
  "Enable Healpix Segmentation, requires healpix to be installed on your system in a place that can be found easily."
  OFF)
# TODO: this is broken af, the module can't seem to find healpix, even after dnf install healpix-devel and/or chealpix-devel
if(USE_HEALPIX)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/module")
  find_package(HealPix REQUIRED)

  target_compile_options(${PROJECT_NAME} PRIVATE -DUSE_HEALPIX)
  # Are these really necessary?!?? 
  target_link_libraries (${PROJECT_NAME} PRIVATE healpix_cxx cxxsupport sharp fftpack c_utils cfitsio)
endif()


#================================================================================================
#                                       BUILDING
#================================================================================================

target_include_directories(${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>
)

# Transform ROOT_CXX_FLAGS into a list of individual flags, and add them to the compilation command
separate_arguments(
  ROOT_CXX_FLAGS_LIST 
  UNIX_COMMAND "${ROOT_CXX_FLAGS} -Wall -Wextra -march=native -DUSE_FFT_TOOLS"
)
target_compile_options(${PROJECT_NAME} PRIVATE ${ROOT_CXX_FLAGS_LIST})

# TODO: compiled just fine without linking
# target_link_libraries(${PROJECT_NAME} PRIVATE ${ZLIB_LIBRARIES} RootFftwWrapper ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})

set(DICTNAME G__${PROJECT_NAME})
file(GLOB HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
root_generate_dictionary(${DICTNAME} ${HEADER_FILES} MODULE ${PROJECT_NAME} LINKDEF LinkDef.h)

#================================================================================================
#                                       INSTALLING
#================================================================================================

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)

# Config files and such
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Generating RootFftwWrapperConfig.cmake
configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
)
# Generating RootFftwWrapperConfigVersion.cmake
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
# Installing these two config files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "lib/cmake/${PROJECT_NAME}"
)

# installing the libraries
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets LIBRARY DESTINATION lib)
# TODO: not sure if pcm files should be installed, or how (https://root-forum.cern.ch/t/cmake-not-installing-pcm-and-rootmap/50997)

# installing the target file
install(EXPORT ${PROJECT_NAME}Targets
    NAMESPACE AntarcticaRoot::
    FILE ${PROJECT_NAME}Targets.cmake
    DESTINATION lib/cmake/${PROJECT_NAME}
)

# installing datafiles (ie. stuff in data/)
install(DIRECTORY data/ DESTINATION share/anitaCalib)
