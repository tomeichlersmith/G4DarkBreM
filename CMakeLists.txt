cmake_minimum_required(VERSION 3.12)

project(G4DarkBreM
  VERSION 1.0.0
  DESCRIPTION "Dark brem simulation integrated within Geant4"
  LANGUAGES CXX)

# Search for Geant4 and load its settings
find_package(Geant4 REQUIRED gdml ui_all vis_all)

# Create an imported Geant4 target if it hasn't been done yet.
if(NOT TARGET Geant4::Interface)
  # Geant4_DEFINITIONS already include -D, this leads to the error
  # <command-line>:0:1: error: macro names must be identifiers if not removed.
  set(Geant4_COMPILE_DEFINITIONS "")
  foreach(def ${Geant4_DEFINITIONS})
    string(REPLACE "-D" "" def ${def})
    list(APPEND Geant4_COMPILE_DEFINITIONS ${def})
  endforeach()

  # Create the Geant4 target
  add_library(Geant4::Interface INTERFACE IMPORTED GLOBAL)

  # Set the target properties
  set_target_properties(
    Geant4::Interface
    PROPERTIES INTERFACE_LINK_LIBRARIES "${Geant4_LIBRARIES}"
               INTERFACE_COMPILE_OPTIONS "${Geant4_Flags}"
               INTERFACE_COMPILE_DEFINITIONS "${Geant4_COMPILE_DEFINITIONS}"
               INTERFACE_INCLUDE_DIRECTORIES "${Geant4_INCLUDE_DIRS}")

  message(STATUS "Found Geant4 version ${Geant4_VERSION}")
endif()

add_library(G4DarkBreM SHARED
  src/G4DarkBreM/APrimePhysics.cxx
  src/G4DarkBreM/G4DarkBremsstrahlung.cxx
  src/G4DarkBreM/G4DarkBremsstrahlung.cxx)
target_link_libraries(G4DarkBreM PUBLIC Geant4::Interface)
target_include_directories(include)

find_package(Boost COMPONENTS program_options)

add_executable(test-xsec-calc app/test_xsec_calc.cxx)
target_link_libraries(test-xsec-calc PRIVATE G4DarkBreM Boost::program_options)
install(TARGETS test-xsec-calc DESTINATION bin)

add_executable(test-scample app/test_scample.cxx)
target_link_libraries(test-scample PRIVATE G4DarkBreM Boost::program_options)
install(TARGETS test-scample DESTINATION bin)

#add_executable(simulate-dark-brem app/simulate_dark_brem.cxx)
#target_link_libraries(simulate-dark-brem PRIVATE G4DarkBreM Boost::program_options)
#install(TARGETS simulate-dark-brem DESTINATION bin)

