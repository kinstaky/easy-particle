add_library(
	particle
	particle.cpp
	nuclear_data.cpp
	path_manager.cpp
)

target_compile_features(particle PUBLIC cxx_std_17)

target_include_directories(
	particle
	PUBLIC
		$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
		$<INSTALL_INTERFACE:include>
)

if (USE_HARDCODE_DATA)
	target_compile_definitions(particle PUBLIC USE_HARDCODE_DATA)
endif()

target_compile_definitions(
	particle
	PUBLIC
		EASY_PARTICLE_SOURCE_DATA_DIR="${EASY_PARTICLE_DATA_SOURCE_DIR}"
		EASY_PARTICLE_DATA_DIR="${EASY_PARTICLE_DATA_INSTALL_DIR}"
)

target_link_libraries(
	particle
	PUBLIC
		catima
		ROOT::Core
		ROOT::GenVector
		ROOT::MathCore
)

install(TARGETS particle)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)
