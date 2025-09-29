cmake_minimum_required(VERSION 3.23.0)
project(AntarcticaRoot VERSION 1.0.0)
add_library(${PROJECT_NAME} SHARED)

#================================================================================================
#                                       HOUSEKEEPING
#================================================================================================
include(GNUInstallDirs)
if(NOT CMAKE_BUILD_TYPE) 
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE) 
endif()

find_package(RootFftwWrapper CONFIG REQUIRED)
find_package(ROOT CONFIG REQUIRED COMPONENTS HistPainter Physics Ged)

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

option(CORE_ONLY "Compile only Geoid and RampdemReader (the only ones used by NiceMC)" OFF)

#================================================================================================
#                                       BUILDING
#================================================================================================
add_subdirectory(include)
add_subdirectory(src)

# use warning flags in RelWithDebInfo mode, which is set to be the default mode
target_compile_options(${PROJECT_NAME} PRIVATE $<$<CONFIG:RelWithDebInfo>:-Wall -Wextra>)

target_link_libraries(${PROJECT_NAME} PUBLIC  RootFftwWrapper::RootFftwWrapper
                                              ROOT::HistPainter ROOT::Ged ROOT::Physics)

#================================================================================================
#                                       INSTALLING
#================================================================================================

# Generate AntarcticaRootConfigVersion.cmake in the build/ directory
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Install the config file generated above and the Config file in the cmake/ directory
install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)

# Install this library and its headers, exporting as a CMake Target that downstream projects can import
install(
  TARGETS ${PROJECT_NAME} 
  EXPORT ${PROJECT_NAME}Targets 
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Generate a Targets.cmake file and directly install it
# this is so that downstream projects can find the target created by this project.
install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE AntarcticaRoot::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)

# installing datafiles (ie. stuff in data/)
install(DIRECTORY data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/anitaCalib)
    # note: trailing "/" is important here
