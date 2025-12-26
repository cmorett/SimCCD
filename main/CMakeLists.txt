# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# When vcpkg is available configure its toolchain upfront so dependencies such
# as assimp are picked up automatically without additional manual flags.
if(NOT SIMCCD_VCPKG_ROOT AND DEFINED ENV{VCPKG_ROOT})
  set(SIMCCD_VCPKG_ROOT "$ENV{VCPKG_ROOT}" CACHE PATH "Detected vcpkg root for this build")
endif()
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE OR NOT CMAKE_TOOLCHAIN_FILE)
  set(_VCPKG_CANDIDATES)
  if(DEFINED ENV{VCPKG_ROOT})
    list(APPEND _VCPKG_CANDIDATES "$ENV{VCPKG_ROOT}")
  endif()
  list(APPEND _VCPKG_CANDIDATES "C:/vcpkg")
  foreach(_vcpkg_root IN LISTS _VCPKG_CANDIDATES)
    if(_vcpkg_root AND EXISTS "${_vcpkg_root}/scripts/buildsystems/vcpkg.cmake")
      set(CMAKE_TOOLCHAIN_FILE "${_vcpkg_root}/scripts/buildsystems/vcpkg.cmake"
          CACHE FILEPATH "Toolchain file for vcpkg integration")
      if(NOT DEFINED VCPKG_TARGET_TRIPLET AND EXISTS "${_vcpkg_root}/installed/x64-windows")
        set(VCPKG_TARGET_TRIPLET "x64-windows"
            CACHE STRING "Default vcpkg triplet used for dependency discovery")
      endif()
      set(SIMCCD_VCPKG_ROOT "${_vcpkg_root}" CACHE PATH "Detected vcpkg root for this build")
      break()
    endif()
  endforeach()
  unset(_vcpkg_root)
  unset(_VCPKG_CANDIDATES)
endif()

project(B02)

if(SIMCCD_VCPKG_ROOT AND DEFINED VCPKG_TARGET_TRIPLET)
  set(_VCPKG_INSTALLED_DIR "${SIMCCD_VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}")
  if(EXISTS "${_VCPKG_INSTALLED_DIR}")
    list(APPEND CMAKE_PREFIX_PATH "${_VCPKG_INSTALLED_DIR}")
  endif()
  unset(_VCPKG_INSTALLED_DIR)
endif()

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(_SIMCCD_USING_MINGW_WITH_MSVC_LIBS OFF)
  foreach(_g4lib ${Geant4_LIBRARIES})
    if(_g4lib MATCHES "\\.lib$")
      set(_SIMCCD_USING_MINGW_WITH_MSVC_LIBS ON)
      break()
    endif()
  endforeach()
  if(_SIMCCD_USING_MINGW_WITH_MSVC_LIBS)
    message(FATAL_ERROR
      "Geant4 at ${Geant4_DIR} was built with MSVC (.lib files). "
      "Please configure this project with the Visual Studio toolchain "
      "(e.g. run 'cmake --preset vs2022-release' or reuse the VS Developer "
      "Command Prompt) instead of MinGW.")
  endif()
  unset(_SIMCCD_USING_MINGW_WITH_MSVC_LIBS)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
set(SIMCCD_ASSETS_DIR "${PROJECT_SOURCE_DIR}/assets")
set(SIMCCD_BUILD_ASSETS_DIR "${PROJECT_BINARY_DIR}/assets")
configure_file(${PROJECT_SOURCE_DIR}/include/SimCCDConfig.hh.in
  ${PROJECT_BINARY_DIR}/generated/SimCCDConfig.hh @ONLY)
include_directories(${PROJECT_BINARY_DIR}/generated)

find_package(assimp CONFIG REQUIRED)


add_definitions(-DG4ANALYSIS_USE)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
                                                               
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
# Remove analysis.cc from the sources list
list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/src/B02AnalysisManager.cc)
file(COPY ${PROJECT_SOURCE_DIR}/assets DESTINATION ${PROJECT_BINARY_DIR})



#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(b02_executable B02.cc ${sources} ${headers})
target_include_directories(b02_executable PRIVATE ${Geant4_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/generated)
target_link_libraries(b02_executable PRIVATE assimp::assimp ${Geant4_LIBRARIES})


add_custom_command(TARGET b02_executable POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${PROJECT_SOURCE_DIR}/assets"
    "$<TARGET_FILE_DIR:b02_executable>/assets"
)


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(B02_SCRIPTS
  b02.in
  b02.out
  init_vis.mac
  macros/geom_check.mac
  macros/run_tierB_cosmic_muons.mac
  macros/run_tierB_cosmic_muons_20cm.mac
  macros/run_tierB_paper_default_1M.mac
  macros/run_with_cad.mac
  macros/run_with_cad_quiet.mac
  macros/run_without_cad.mac
  run1.mac
  run2.mac
  vis.mac
  )

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/macros)

foreach(_script ${B02_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(B02 DEPENDS b02_executable)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS b02_executable DESTINATION bin)
