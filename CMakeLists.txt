cmake_minimum_required(VERSION 3.19.0)
project(AntarcticaRoot VERSION 1.0.0)

    # set(UTIL_LIB_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/lib)
    # set(UTIL_INC_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/include)
    # set(UTIL_BIN_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/bin)
    # set(UTIL_SHARE_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/share)
    # set(LD_UTIL $ENV{PUEO_UTIL_INSTALL_DIR}/lib)
    # set(INC_UTIL $ENV{PUEO_UTIL_INSTALL_DIR}/include)


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net
  Minuit Gui HistPainter)
find_package(FFTW REQUIRED)
find_package(ZLIB REQUIRED)


if(TARGET RootFftwWrapper)
  message("RootFftwWrapper is a target")
else()
 link_directories(${UTIL_LIB_DIR})
endif()




#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)

message("ROOT_VERSION is set to ${ROOT_VERSION}")
if( ${ROOT_VERSION} VERSION_GREATER "5.99/99")
    message("Using ROOT_VERSION 6")
    include(${ROOT_USE_FILE})
else()
#    include(RootNewMacros)
     add_definitions(${ROOT_DEFINITIONS})
endif()

message("ROOT_INCLUDE_DIRS is set to ${ROOT_INCLUDE_DIRS}")

# Transform ROOT_CXX_FLAGS into a list of individual flags
separate_arguments(ROOT_CXX_FLAGS_LIST UNIX_COMMAND ${ROOT_CXX_FLAGS})
target_compile_options(${PROJECT_NAME} PRIVATE ${ROOT_CXX_FLAGS_LIST})

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include ${FFTW_INCLUDES} ${ZLIB_INCLUDES} ${UTIL_INC_DIR})

file(GLOB HEADER_FILES ${PROJECT_SOURCE_DIR}/include/*.h ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB_RECURSE SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cxx ${PROJECT_SOURCE_DIR}/src/*.cc)

set(DICTNAME G__${PROJECT_NAME})

ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
#target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES} RootFftwWrapper ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES} ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})

if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})
endif()


message("UTIL_INC_DIR is set to ${UTIL_INC_DIR}")
message("DICTNAME is set to ${DICTNAME}" )
message("PROJECT_BINARY_DIR is set to ${PROJECT_BINARY_DIR}")



install (FILES ${HEADER_FILES} DESTINATION ${UTIL_INC_DIR})
install (TARGETS ${PROJECT_NAME}
	ARCHIVE DESTINATION ${UTIL_LIB_DIR}
	LIBRARY DESTINATION ${UTIL_LIB_DIR}
	RUNTIME DESTINATION ${UTIL_BIN_DIR})
      #Only needed for ROOT6
if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
  install (FILES ${PROJECT_BINARY_DIR}/${DICTNAME}_rdict.pcm DESTINATION ${UTIL_LIB_DIR})
#  install (FILES ${PROJECT_BINARY_DIR}/lib${PROJECT_NAME}.rootmap DESTINATION ${UTIL_LIB_DIR})
endif()

install (DIRECTORY data/ DESTINATION ${UTIL_SHARE_DIR}/anitaCalib)

#################################################################################3
### The default build configuration is INSANE. No optimization? What is this 1971?

message (" Adding new build type")

set(CMAKE_CXX_FLAGS_DEFAULT
  "-Os -g -Wall -Wextra -march=native"
  CACHE STRING "c++ Flags used during default eventReaderRoot builds"
  FORCE )

set(CMAKE_C_FLAGS_DEFAULT
  "-Os -g -Wall -march=native"
  CACHE STRING "c Flags used during default eventReaderRoot builds"
  FORCE )

set(CMAKE_EXE_LINKER_FLAGS_DEFAULT
  "-g"
  CACHE STRING "ld Flags used during default eventReaderRoot builds"
  FORCE )

set(CMAKE_SHARED_LINKER_FLAGS_DEFAULT
  "-g"
  CACHE STRING "ld Flags used during default eventReaderRoot builds"
  FORCE )


mark_as_advanced ( CMAKE_CXX_FLAGS_DEFAULT CMAKE_C_FLAGS_DEFAULT CMAKE_EXE_LINKER_FLAGS_DEFAULT CMAKE_SHARED_LINKER_FLAGS_DEFAULT)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Default
    CACHE STRING "Choose tye type of build: None Debug Release RelWithDebInfo MinSizeRel Default"
    FORCE )
endif()


add_definitions("-DUSE_FFT_TOOLS")

option ( FORCE_OLD_GPP_ABI " Force old g++ ABI; this might be necessary if using new g++ with ROOT compiled with older g++ or other similar situations" OFF )
if (FORCE_OLD_GPP_ABI)
  add_definitions( -D_GLIBCXX_USE_CXX11_ABI=0 )
endif()

option ( USE_GEOGRAPHIC_LIB "Use GeographicLib for geodesic calculations" ) 
if (USE_GEOGRAPHIC_LIB)
  add_definitions  (-DUSE_GEOGRAPHIC_LIB )
  find_package(GeographicLib REQUIRED) 
  target_link_libraries (${PROJECT_NAME} ${GeographicLib_LIBRARIES} )
endif()


option (USE_HEALPIX "Enable Healpix Segmentation, requires healpix to be
installed on your system in a place that can be found easily." OFF)

if (USE_HEALPIX)
  add_definitions  (-DUSE_HEALPIX )
  # Are these really necessary?!?? 
  target_link_libraries (${PROJECT_NAME} healpix_cxx cxxsupport sharp fftpack c_utils cfitsio)
endif()



