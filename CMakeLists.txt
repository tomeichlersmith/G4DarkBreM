cmake_minimum_required(VERSION 3.12)

project(G4DarkBreM
  VERSION 1.1.0
  DESCRIPTION "Dark brem simulation integrated within Geant4"
  LANGUAGES CXX)

# Search for Geant4 and load its settings
find_package(Geant4 10.2.3 REQUIRED)
find_package(Boost 1.68 REQUIRED)

include(${Geant4_USE_FILE})

add_library(G4DarkBreM SHARED
  src/G4DarkBreM/ElementXsecCache.cxx
  src/G4DarkBreM/G4APrime.cxx
  src/G4DarkBreM/G4DarkBreMModel.cxx
  src/G4DarkBreM/G4DarkBremsstrahlung.cxx)
target_link_libraries(G4DarkBreM PUBLIC ${Geant4_LIBRARIES} Boost::boost)
target_include_directories(G4DarkBreM PUBLIC include)
install(TARGETS G4DarkBreM DESTINATION lib)

add_executable(g4db-xsec-calc app/xsec_calc.cxx)
target_link_libraries(g4db-xsec-calc PRIVATE G4DarkBreM)
install(TARGETS g4db-xsec-calc DESTINATION bin)

add_executable(g4db-scample app/scample.cxx)
target_link_libraries(g4db-scample PRIVATE G4DarkBreM)
install(TARGETS g4db-scample DESTINATION bin)

add_executable(g4db-simulate app/simulate.cxx)
target_link_libraries(g4db-simulate PRIVATE G4DarkBreM ${Geant4_LIBRARIES})
install(TARGETS g4db-simulate DESTINATION bin)

