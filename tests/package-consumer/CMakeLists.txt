cmake_minimum_required(VERSION 3.24)

project(SCPResamplerIOConsumer LANGUAGES CXX)

find_package(SCPResamplerIO CONFIG REQUIRED)

add_executable(scp_resampler_io_package_consumer main.cpp)
target_link_libraries(scp_resampler_io_package_consumer PRIVATE SCP::ResamplerIO)
