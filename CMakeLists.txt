cmake_minimum_required(VERSION 3.16)
project(B0Trackers)

find_package(EICrecon REQUIRED)
find_package(spdlog REQUIRED)

set(CMAKE_CXX_STANDARD ${EICrecon_CXX_STANDARD})

# Explicit source list
set(B0Trackers_PLUGIN_SOURCES
    B0Trackers.cc
    B0Trackers.h
    variables.h
)

# Build the shared plugin
add_library(B0Trackers_plugin SHARED ${B0Trackers_PLUGIN_SOURCES})
target_link_libraries(B0Trackers_plugin EICrecon::rootfile_plugin)
target_link_libraries(B0Trackers_plugin ${JANA_LIB} ${ROOT_LIBRARIES} spdlog::spdlog)

# Set output name and format
set_target_properties(B0Trackers_plugin
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "B0Trackers"
        SUFFIX ".so"
)

# Install location
install(TARGETS B0Trackers_plugin
        DESTINATION ${USER_PLUGIN_OUTPUT_DIRECTORY})
