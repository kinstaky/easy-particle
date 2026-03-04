cmake_minimum_required(VERSION 3.10)

project(
	easy-particle
	VERSION 0.1.0
	DESCRIPTION "A project for easy use physics nuclear particle."
	LANGUAGES CXX
)

include(GNUInstallDirs)
include(FetchContent)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
	message(FATAL_ERROR "You cannot build in a source directory (or any directory with"
						"CmakeLists.txt file). Please make a build subdirectory. Feel free"
						"to remove CMakeCache.txt and CMakeFiles.")
endif()

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Fetching catima...")
FetchContent_Declare(
	catima
	GIT_REPOSITORY https://github.com/kinstaky/catima.git
	GIT_TAG v1.90
	GIT_SHALLOW TRUE
	GIT_SHALLOW_DEPTH 1
)
FetchContent_MakeAvailable(catima)

message(STATUS "Fetching cxxopts...")
FetchContent_Declare(
	cxxopts
	GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
	GIT_TAG v3.3.1
	GIT_SHALLOW TRUE
	GIT_SHALLOW_DEPTH 1
)
FetchContent_MakeAvailable(cxxopts)

message(STATUS "Fetch finished.")

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	set(CMAKE_CXX_EXTENSIONS OFF)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra")
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)
	include(CTest)
endif()

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist Graf GenVector MathCore)

option(BUILD_TESTING "Build test" OFF)
option(USE_HARDCODE_DATA "Use hardcode AMDC2020 atomic mass data" OFF)

function(show_config_options)
	message("")
	message("Build Properties")
	message("  c++ compiler.................: " ${CMAKE_CXX_COMPILER})
	message("  c++ stardard.................: " ${CMAKE_CXX_STANDARD})
	message("  cmake........................: " ${CMAKE_VERSION})
	message("  cmake build type.............: " ${CMAKE_BUILD_TYPE})
	message("")
	message("Build Options")
	message("  build testing................: " ${BUILD_TESTING})
	message("User-defined Flags")
	message("  USE_HARDCODE_DATA............: " ${USE_HARDCODE_DATA})
endfunction()

show_config_options()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_compile_definitions(DEBUG)
endif()

set(EASY_PARTICLE_DATA_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(EASY_PARTICLE_DATA_FILE_SOURCE "${EASY_PARTICLE_DATA_SOURCE_DIR}/amdc_ion_2020.txt")
set(EASY_PARTICLE_DATA_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/share/easy-particle")
set(EASY_PARTICLE_DATA_INSTALL_DIR "${CMAKE_INSTALL_FULL_DATADIR}/easy-particle")

file(MAKE_DIRECTORY "${EASY_PARTICLE_DATA_BUILD_DIR}")
file(COPY "${EASY_PARTICLE_DATA_FILE_SOURCE}" DESTINATION "${EASY_PARTICLE_DATA_BUILD_DIR}")

install(
	FILES "${EASY_PARTICLE_DATA_FILE_SOURCE}"
	DESTINATION "${CMAKE_INSTALL_DATADIR}/easy-particle"
)

add_subdirectory(src)
add_subdirectory(bin)

if (BUILD_TESTING)
	add_subdirectory(tests)
endif()
