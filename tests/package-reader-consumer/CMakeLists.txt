cmake_minimum_required(VERSION 3.24)

project(SCPResamplerIOReaderConsumer LANGUAGES CXX)

find_package(SCPResamplerIO CONFIG REQUIRED COMPONENTS Readers)

add_executable(scp_resampler_io_reader_consumer main.cpp)
target_link_libraries(scp_resampler_io_reader_consumer
  PRIVATE SCP::ResamplerReaders
)
