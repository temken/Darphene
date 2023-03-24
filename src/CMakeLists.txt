# OPTION A: Build an executable and a static library
# Executable
add_executable(Darphene
	main.cpp)

target_compile_options(Darphene PUBLIC -Wall -pedantic)

target_link_libraries(Darphene
	PUBLIC
	coverage_config
	Darphene_library)

target_include_directories(Darphene
	PRIVATE
	${GENERATED_DIR})

install(TARGETS Darphene DESTINATION ${BIN_DIR})

# Static library
file(GLOB FILES "*.cpp")
list(REMOVE_ITEM FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_library(Darphene_library STATIC
	${FILES})

set_target_properties(Darphene_library PROPERTIES PREFIX "")

target_compile_options(Darphene_library PUBLIC -Wall -pedantic)

target_link_libraries(Darphene_library
	PRIVATE
	${Boost_LIBRARIES}
	"${ARB}"
	PUBLIC
	${MPI_CXX_LIBRARIES}
	Eigen3::Eigen
	coverage_config
	libobscura
	libphysica)

target_include_directories(Darphene_library
	PRIVATE
	${GENERATED_DIR}
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
	${INCLUDE_DIR})

install(TARGETS Darphene_library DESTINATION ${LIB_DIR})
