if(NOT DEFINED RESAMPLE_EXECUTABLE OR NOT EXISTS "${RESAMPLE_EXECUTABLE}")
  message(FATAL_ERROR "RESAMPLE_EXECUTABLE must name the built Release Resample tool.")
endif()
if(NOT DEFINED FIXTURE_ROOT OR NOT FIXTURE_ROOT MATCHES "[/.]tests[/\\]")
  message(FATAL_ERROR "The CLI regression fixture must be below the workspace .tests directory.")
endif()
file(MAKE_DIRECTORY "${FIXTURE_ROOT}/input")
file(WRITE "${FIXTURE_ROOT}/input/frame.vti" [=[<?xml version="1.0"?>
<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian">
  <ImageData WholeExtent="0 1 0 1 0 1" Origin="0 0 0" Spacing="1 1 1">
    <Piece Extent="0 1 0 1 0 1">
      <PointData Scalars="Signal"><DataArray type="Float64" Name="Signal" format="ascii">0 1 2 3 4 5 6 7</DataArray></PointData>
      <CellData/>
    </Piece>
  </ImageData>
</VTKFile>
]=])
file(WRITE "${FIXTURE_ROOT}/input/series.pvd" [=[<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
    <DataSet timestep="4.257" group="" part="0" file="frame.vti"/>
    <DataSet timestep="4.258" group="" part="0" file="frame.vti"/>
  </Collection>
</VTKFile>
]=])
file(WRITE "${FIXTURE_ROOT}/input.resample.scpo"
  "input_path = ${FIXTURE_ROOT}/input/series.pvd\n"
  "output_dir = ${FIXTURE_ROOT}/output\n"
  "output_basename = cli\nsolver_status = error\ndimensions = [3, 3, 3]\n"
  "method = linear\nbounds_mode = global\nphysical_dt = 0.001\n")
execute_process(
  COMMAND "${RESAMPLE_EXECUTABLE}" "${FIXTURE_ROOT}/input.resample.scpo"
  WORKING_DIRECTORY "${FIXTURE_ROOT}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error
  TIMEOUT 90
)
file(WRITE "${FIXTURE_ROOT}/stdout.log" "${output}")
file(WRITE "${FIXTURE_ROOT}/stderr.log" "${error}")
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Resample CLI failed (${result}):\n${output}\n${error}")
endif()
foreach(name IN ITEMS cli.pvd cli_f000000.vti cli_f000001.vti)
  if(NOT EXISTS "${FIXTURE_ROOT}/output/${name}")
    message(FATAL_ERROR "Resample CLI did not publish ${name}.")
  endif()
endforeach()
file(READ "${FIXTURE_ROOT}/output/cli.pvd" collection)
string(REGEX MATCHALL "timestep=\"[^\"]+\"" times "${collection}")
list(LENGTH times time_count)
if(NOT time_count EQUAL 2)
  message(FATAL_ERROR "Resample CLI did not publish both source times.")
endif()
list(GET times 0 first_time)
list(GET times 1 second_time)
string(REGEX REPLACE "^timestep=\"([^\"]+)\"$" "\\1" first_time "${first_time}")
string(REGEX REPLACE "^timestep=\"([^\"]+)\"$" "\\1" second_time "${second_time}")
if(NOT first_time EQUAL 4.257 OR NOT second_time EQUAL 4.258)
  message(FATAL_ERROR "Resample CLI changed the explicit source PVD timestamps.")
endif()
message(STATUS "Resample CLI completed both PVD frames using the available VTK workers.")
