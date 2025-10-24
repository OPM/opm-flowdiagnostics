# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

cmake_minimum_required(VERSION 3.27)

project(opm-flowdiagnostics
  VERSION   2024.04
  LANGUAGES C CXX
)

enable_testing()

set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

set(CMAKE_C_STANDARD          99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS        OFF)

set(CMAKE_DEBUG_POSTFIX "-d")

add_library(OpmFDEngine)

set_target_properties (OpmFDEngine PROPERTIES
  SOVERSION 1
  VERSION   2024.04
)

target_sources(OpmFDEngine
  PRIVATE
    opm/flowdiagnostics/CellSet.cpp
    opm/flowdiagnostics/ConnectionValues.cpp
    opm/flowdiagnostics/ConnectivityGraph.cpp
    opm/flowdiagnostics/DerivedQuantities.cpp
    opm/flowdiagnostics/Solution.cpp
    opm/flowdiagnostics/Toolbox.cpp
    opm/flowdiagnostics/TracerTofSolver.cpp
    opm/utility/graph/AssembledConnections.cpp
    opm/utility/numeric/RandomVector.cpp
    opm/utility/graph/tarjan.c

  PUBLIC FILE_SET engine_public_headers
    TYPE HEADERS
    FILES
      opm/flowdiagnostics/CellSet.hpp
      opm/flowdiagnostics/CellSetValues.hpp
      opm/flowdiagnostics/ConnectionValues.hpp
      opm/flowdiagnostics/ConnectivityGraph.hpp
      opm/flowdiagnostics/DerivedQuantities.hpp
      opm/flowdiagnostics/Solution.hpp
      opm/flowdiagnostics/Toolbox.hpp
      opm/flowdiagnostics/TracerTofSolver.hpp
      opm/utility/graph/AssembledConnectionsIteration.hpp
      opm/utility/graph/AssembledConnections.hpp
      opm/utility/numeric/RandomVector.hpp

  PRIVATE FILE_SET engine_internal_headers
    TYPE HEADERS
    FILES
      opm/utility/graph/tarjan.h
)

target_include_directories(OpmFDEngine PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

write_basic_package_version_file(OpmFDEngineConfigVersion.cmake
  VERSION ${OpmFDEngine_VERSION}
  COMPATIBILITY ExactVersion
)

install(TARGETS OpmFDEngine EXPORT OpmFDEngineTargets
  RUNTIME           # Following options apply to runtime artifacts.
    COMPONENT OpmFDEngine_Runtime
  LIBRARY           # Following options apply to library artifacts.
    COMPONENT OpmFDEngine_Runtime
    NAMELINK_COMPONENT OpmFDEngine_Development
  ARCHIVE           # Following options apply to archive artifacts.
    COMPONENT OpmFDEngine_Development
  FILE_SET engine_public_headers # Following options apply to file set HEADERS.
    COMPONENT OpmFDEngine_Development
)

install(EXPORT OpmFDEngineTargets
  FILE OpmFDEngineTargets.cmake
  NAMESPACE OpmFDEngine::
  DESTINATION share/cmake/OpmFDEngine
)

install(FILES
    "OpmFDEngineConfig.cmake"
    "${CMAKE_BINARY_DIR}/OpmFDEngineConfigVersion.cmake"
  DESTINATION
    share/cmake/OpmFDEngine
)

find_package(Boost
  1.44
  REQUIRED
  COMPONENTS "unit_test_framework"
)

add_executable(test_assembledconnections tests/test_assembledconnections.cpp)
target_link_libraries(test_assembledconnections
  PRIVATE OpmFDEngine Boost::unit_test_framework
  )

add_executable(test_cellset tests/test_cellset.cpp)
target_link_libraries(test_cellset
  PRIVATE OpmFDEngine Boost::unit_test_framework
  )

add_executable(test_connectionvalues tests/test_connectionvalues.cpp)
target_link_libraries(test_connectionvalues
  PRIVATE OpmFDEngine Boost::unit_test_framework
  )

add_executable(test_connectivitygraph tests/test_connectivitygraph.cpp)
target_link_libraries(test_connectivitygraph
  PRIVATE OpmFDEngine Boost::unit_test_framework
  )

add_executable(test_derivedquantities tests/test_derivedquantities.cpp)
target_link_libraries(test_derivedquantities
  PRIVATE OpmFDEngine Boost::unit_test_framework
  )

add_executable(test_flowdiagnosticstool tests/test_flowdiagnosticstool.cpp)
target_link_libraries(test_flowdiagnosticstool
  PRIVATE OpmFDEngine Boost::unit_test_framework
  )

add_executable(test_tarjan tests/test_tarjan.cpp)
target_link_libraries(test_tarjan
  PRIVATE OpmFDEngine Boost::unit_test_framework
  )

add_test(AssembledConnections test_assembledconnections)
add_test(CellSet test_cellset)
add_test(ConnectionValues test_connectionvalues)
add_test(ConnectivityGraph test_connectivitygraph)
add_test(DerivedQuantities test_derivedquantities)
add_test(FlowDiagnsticsTool test_flowdiagnosticstool)
add_test(Tarjan test_tarjan)
