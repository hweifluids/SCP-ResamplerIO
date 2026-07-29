cmake_minimum_required(VERSION 3.24)

project(SCPResamplerIO VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(SCP_RESAMPLER_IO_BUILD_TOOL "Build the Resample command-line tool." ${PROJECT_IS_TOP_LEVEL})
option(SCP_RESAMPLER_IO_BUILD_TESTS "Build SCP-ResamplerIO tests." ${PROJECT_IS_TOP_LEVEL})
option(SCP_RESAMPLER_IO_BUILD_FULL_LIBRARY
  "Build the Qt-facing inspector and resampling library in addition to the lightweight native readers."
  ON)

if(SCP_RESAMPLER_IO_BUILD_TESTS)
  include(CTest)
endif()

set(SCP_RESAMPLER_IO_VTK_COMPONENTS
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersGeneral
  FiltersPoints
  FiltersVerdict
  IOEnSight
  IOGeometry
  IOHDF
  IOLegacy
  IOParallel
  IOXML
)

if(NOT TARGET SCP::SolverSDK)
  set(scp_solver_sdk_sibling "${CMAKE_CURRENT_SOURCE_DIR}/../SCP-SolverSDK")
  set(scp_solver_sdk_local "${CMAKE_CURRENT_SOURCE_DIR}/.deps/SCP-SolverSDK")
  if(EXISTS "${scp_solver_sdk_sibling}/CMakeLists.txt")
    set(scp_solver_sdk_source "${scp_solver_sdk_sibling}")
  elseif(EXISTS "${scp_solver_sdk_local}/CMakeLists.txt")
    set(scp_solver_sdk_source "${scp_solver_sdk_local}")
  endif()

  if(scp_solver_sdk_source)
    file(STRINGS "${scp_solver_sdk_source}/VERSION" scp_solver_sdk_version LIMIT_COUNT 1)
    if(NOT scp_solver_sdk_version STREQUAL "0.2.0")
      message(FATAL_ERROR
        "SCP-ResamplerIO requires SCP-SolverSDK 0.2.0; found '${scp_solver_sdk_version}' "
        "at ${scp_solver_sdk_source}.")
    endif()
    add_subdirectory(
      "${scp_solver_sdk_source}"
      "${CMAKE_CURRENT_BINARY_DIR}/_deps/SCP-SolverSDK"
      EXCLUDE_FROM_ALL
    )
  else()
    file(GLOB scp_solver_sdk_installed_configs
      "${scp_solver_sdk_sibling}/release/*/lib/cmake/SCPSolverSDK/SCPSolverSDKConfig.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/.deps/SCP-SolverSDK/release/*/lib/cmake/SCPSolverSDK/SCPSolverSDKConfig.cmake"
    )
    list(SORT scp_solver_sdk_installed_configs)
    if(NOT scp_solver_sdk_installed_configs)
      message(FATAL_ERROR
        "SCP-SolverSDK was not found at the fixed sibling or module-local .deps locations.")
    endif()
    list(GET scp_solver_sdk_installed_configs 0 scp_solver_sdk_config)
    get_filename_component(SCPSolverSDK_DIR "${scp_solver_sdk_config}" DIRECTORY)
    find_package(SCPSolverSDK 0.2.0 EXACT CONFIG REQUIRED NO_DEFAULT_PATH)
  endif()
endif()

find_package(VTK CONFIG REQUIRED COMPONENTS ${SCP_RESAMPLER_IO_VTK_COMPONENTS})

add_library(SCP_ResamplerReaders STATIC
  src/DataSetReader.cpp
)
add_library(SCP::ResamplerReaders ALIAS SCP_ResamplerReaders)
set_target_properties(SCP_ResamplerReaders PROPERTIES EXPORT_NAME ResamplerReaders)
target_compile_features(SCP_ResamplerReaders PUBLIC cxx_std_17)
target_include_directories(SCP_ResamplerReaders PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(SCP_ResamplerReaders PUBLIC
  SCP::SolverSDK
  ${VTK_LIBRARIES}
)
if(WIN32)
  target_compile_definitions(SCP_ResamplerReaders PUBLIC NOMINMAX WIN32_LEAN_AND_MEAN)
endif()
if(MSVC)
  target_compile_options(SCP_ResamplerReaders PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
endif()
vtk_module_autoinit(
  TARGETS SCP_ResamplerReaders
  MODULES ${VTK_LIBRARIES}
)

if(NOT SCP_RESAMPLER_IO_BUILD_FULL_LIBRARY)
  if(SCP_RESAMPLER_IO_BUILD_TOOL OR SCP_RESAMPLER_IO_BUILD_TESTS)
    message(FATAL_ERROR
      "SCP_RESAMPLER_IO_BUILD_TOOL/TESTS require SCP_RESAMPLER_IO_BUILD_FULL_LIBRARY=ON.")
  endif()
  install(TARGETS SCP_ResamplerReaders
    EXPORT SCPResamplerIOTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  )
  install(DIRECTORY include/streamcenterplus/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/streamcenterplus"
  )
  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SCPResamplerIOConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/SCPResamplerIOConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SCPResamplerIO"
  )
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/SCPResamplerIOConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
  )
  install(EXPORT SCPResamplerIOTargets
    FILE SCPResamplerIOTargets.cmake
    NAMESPACE SCP::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SCPResamplerIO"
  )
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SCPResamplerIOConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SCPResamplerIOConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SCPResamplerIO"
  )
  return()
endif()

find_package(Qt6 6.5 REQUIRED COMPONENTS Core)

add_library(SCP_ResamplerIO STATIC
  src/FlowDataInspector.cpp
  resample/src/Parameters.cpp
  resample/src/Solver.cpp
  resample/src/VtkIO.cpp
)
add_library(SCP::ResamplerIO ALIAS SCP_ResamplerIO)
set_target_properties(SCP_ResamplerIO PROPERTIES EXPORT_NAME ResamplerIO)

target_compile_features(SCP_ResamplerIO PUBLIC cxx_std_17)
target_include_directories(SCP_ResamplerIO PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/resample/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(SCP_ResamplerIO PUBLIC
  SCP::ResamplerReaders
  SCP::SolverSDK
  Qt6::Core
  ${VTK_LIBRARIES}
)

if(WIN32)
  target_compile_definitions(SCP_ResamplerIO PUBLIC NOMINMAX WIN32_LEAN_AND_MEAN)
endif()
if(MSVC)
  target_compile_options(SCP_ResamplerIO PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
endif()

vtk_module_autoinit(
  TARGETS SCP_ResamplerIO
  MODULES ${VTK_LIBRARIES}
)

if(SCP_RESAMPLER_IO_BUILD_TOOL)
  include("${SCP_SOLVER_SDK_CMAKE_DIR}/StreamcenterPlusSolverInfo.cmake")
  include("${SCP_SOLVER_SDK_CMAKE_DIR}/StreamcenterPlusSolverBundle.cmake")

  set(resampler_runtime_dirs)
  foreach(runtime_target IN ITEMS Qt6::Core ${VTK_LIBRARIES})
    if(NOT TARGET "${runtime_target}")
      continue()
    endif()
    foreach(location_property IN ITEMS
        IMPORTED_LOCATION_RELEASE
        IMPORTED_LOCATION_RELWITHDEBINFO
        IMPORTED_LOCATION)
      get_target_property(runtime_location "${runtime_target}" "${location_property}")
      if(runtime_location AND NOT runtime_location MATCHES "-NOTFOUND$")
        get_filename_component(runtime_dir "${runtime_location}" DIRECTORY)
        if(EXISTS "${runtime_dir}")
          list(APPEND resampler_runtime_dirs "${runtime_dir}")
        endif()
      endif()
    endforeach()
  endforeach()
  list(APPEND STREAMCENTERPLUS_PACKAGE_RUNTIME_DIRS ${resampler_runtime_dirs})
  list(REMOVE_DUPLICATES STREAMCENTERPLUS_PACKAGE_RUNTIME_DIRS)

  add_executable(streamcenterplus_resample_solver
    resample/src/main.cpp
    resample/src/Cli.cpp
  )
  add_executable(scp_nek5000_smoke
    tools/src/Nek5000Smoke.cpp
  )
  target_link_libraries(streamcenterplus_resample_solver PRIVATE SCP::ResamplerIO)
  target_link_libraries(scp_nek5000_smoke PRIVATE SCP::ResamplerIO)
  streamcenterplus_configure_solver_info(
    streamcenterplus_resample_solver "Resample" ON OFF "tools"
    "single_static,two_zone_static,two_zone_nested,two_zone_multigrid,two_zone_sliding"
  )
  vtk_module_autoinit(
    TARGETS streamcenterplus_resample_solver scp_nek5000_smoke
    MODULES ${VTK_LIBRARIES}
  )
  streamcenterplus_install_solver_bundle(
    streamcenterplus_resample_solver DESTINATION "tools"
  )
  streamcenterplus_install_solver_bundle(
    scp_nek5000_smoke DESTINATION "tools"
  )
endif()

if(SCP_RESAMPLER_IO_BUILD_TESTS AND BUILD_TESTING)
  add_executable(scp_resampler_io_native_readers
    tests/native-readers/main.cpp
  )
  target_link_libraries(scp_resampler_io_native_readers PRIVATE SCP::ResamplerReaders)
  target_compile_definitions(scp_resampler_io_native_readers PRIVATE
    SCP_REPOSITORY_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/../.."
  )
  vtk_module_autoinit(
    TARGETS scp_resampler_io_native_readers
    MODULES ${VTK_LIBRARIES}
  )
  add_test(NAME SCPResamplerIO.NativeReaders
    COMMAND scp_resampler_io_native_readers
  )

  add_executable(scp_resampler_io_nek_resample
    tests/nek-resample/main.cpp
  )
  target_link_libraries(scp_resampler_io_nek_resample PRIVATE SCP::ResamplerIO)
  vtk_module_autoinit(
    TARGETS scp_resampler_io_nek_resample
    MODULES ${VTK_LIBRARIES}
  )
  add_test(NAME SCPResamplerIO.NekSpectralResample
    COMMAND scp_resampler_io_nek_resample
  )
  set_tests_properties(SCPResamplerIO.NativeReaders PROPERTIES TIMEOUT 120)

  add_executable(scp_resampler_io_multizone
    tests/multizone/main.cpp
  )
  target_link_libraries(scp_resampler_io_multizone PRIVATE
    SCP::ResamplerIO
    VTK::hdf5
  )
  target_compile_definitions(scp_resampler_io_multizone PRIVATE
    SCP_REPOSITORY_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/../.."
  )
  vtk_module_autoinit(
    TARGETS scp_resampler_io_multizone
    MODULES ${VTK_LIBRARIES}
  )
  add_test(NAME SCPResamplerIO.MultiZoneManifest
    COMMAND scp_resampler_io_multizone
  )
  set_tests_properties(SCPResamplerIO.MultiZoneManifest PROPERTIES TIMEOUT 120)

  add_executable(scp_resampler_io_nek5000_integration
    tests/nek5000-integration/main.cpp
  )
  target_link_libraries(scp_resampler_io_nek5000_integration PRIVATE SCP::ResamplerIO)
  target_compile_definitions(scp_resampler_io_nek5000_integration PRIVATE
    SCP_RESAMPLER_IO_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
  )
  vtk_module_autoinit(
    TARGETS scp_resampler_io_nek5000_integration
    MODULES ${VTK_LIBRARIES}
  )
  add_test(NAME SCPResamplerIO.Nek5000Integration
    COMMAND scp_resampler_io_nek5000_integration
  )
  set_tests_properties(SCPResamplerIO.Nek5000Integration PROPERTIES
    SKIP_RETURN_CODE 77
    TIMEOUT 600
  )
endif()

install(TARGETS SCP_ResamplerReaders SCP_ResamplerIO
  EXPORT SCPResamplerIOTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
install(DIRECTORY include/
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(DIRECTORY resample/include/
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SCPResamplerIOConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/SCPResamplerIOConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SCPResamplerIO"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SCPResamplerIOConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)

install(EXPORT SCPResamplerIOTargets
  FILE SCPResamplerIOTargets.cmake
  NAMESPACE SCP::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SCPResamplerIO"
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/SCPResamplerIOConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/SCPResamplerIOConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SCPResamplerIO"
)
