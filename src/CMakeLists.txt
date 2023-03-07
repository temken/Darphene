# OPTION A: Build an executable and a static library
# Executable
add_executable(DarPhene
	main.cpp)

target_compile_options(DarPhene PUBLIC -Wall -pedantic)

target_link_libraries(DarPhene
	PUBLIC
	coverage_config
	DarPhene_library)

target_include_directories(DarPhene
	PRIVATE
	${GENERATED_DIR})

install(TARGETS DarPhene DESTINATION ${BIN_DIR})

# Static library
file(GLOB FILES "*.cpp")
list(REMOVE_ITEM FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_library(DarPhene_library STATIC
	${FILES})

set_target_properties(DarPhene_library PROPERTIES PREFIX "")

target_compile_options(DarPhene_library PUBLIC -Wall -pedantic)

target_link_libraries(DarPhene_library
	PRIVATE
	${Boost_LIBRARIES}
	"${ARB}"
	PUBLIC
	${MPI_CXX_LIBRARIES}
	Eigen3::Eigen
	coverage_config
	libobscura
	libphysica)

target_include_directories(DarPhene_library
	PRIVATE
	${GENERATED_DIR}
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
	${INCLUDE_DIR})

install(TARGETS DarPhene_library DESTINATION ${LIB_DIR})
