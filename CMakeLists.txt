cmake_minimum_required(VERSION 3.20)
project(B0Trackers)

option(B0Trackers_BUILD_PLUGIN "Build the JANA/EICrecon B0Trackers plugin" ON)

include(CTest)

if(B0Trackers_BUILD_PLUGIN)
    find_package(EICrecon REQUIRED)
    find_package(spdlog REQUIRED)

    set(CMAKE_CXX_STANDARD ${EICrecon_CXX_STANDARD})

    # Explicit source list
    set(B0Trackers_PLUGIN_SOURCES
        B0Trackers.cc
        B0Trackers.h
        B0TrackersHelpers.h
    )

    # Build the shared plugin
    add_library(B0Trackers_plugin SHARED ${B0Trackers_PLUGIN_SOURCES})
    target_link_libraries(B0Trackers_plugin EICrecon::rootfile_plugin)
    target_link_libraries(B0Trackers_plugin ${JANA_LIB} ${ROOT_LIBRARIES} spdlog::spdlog)

    # Set output name and format
    # JANA loads plugins by bare name, so no "lib" prefix; the suffix follows the
    # platform rather than being pinned to .so.
    set_target_properties(B0Trackers_plugin
        PROPERTIES
            PREFIX ""
            OUTPUT_NAME "B0Trackers"
            SUFFIX "${CMAKE_SHARED_LIBRARY_SUFFIX}"
    )

    # Install location
    install(TARGETS B0Trackers_plugin
            DESTINATION ${USER_PLUGIN_OUTPUT_DIRECTORY})
endif()

if(BUILD_TESTING)
    add_executable(B0TrackersHelpers_test tests/test_helpers.cc)
    target_include_directories(B0TrackersHelpers_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_features(B0TrackersHelpers_test PRIVATE cxx_std_17)
    # The checks report failures explicitly rather than through assert(), so they
    # stay live under -DNDEBUG.
    add_test(NAME B0TrackersHelpers COMMAND B0TrackersHelpers_test)
endif()
