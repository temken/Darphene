# OPTION A: Build an executable and a static library
	# Executable
	add_executable(graphene
		main.cpp)

	target_compile_options(graphene PUBLIC -Wall -pedantic)

	target_link_libraries(graphene
		PUBLIC
			coverage_config
			graphene_library)

	target_include_directories(graphene
		PRIVATE
			${GENERATED_DIR} )

	install(TARGETS graphene DESTINATION ${BIN_DIR})

	# Static library
	add_library(graphene_library STATIC
	    Graphene.cpp
	    Hydrogenic_Wavefunctions.cpp
		Multidimensional_Integration.cpp)

	set_target_properties(graphene_library PROPERTIES PREFIX "")

	target_compile_options(graphene_library PUBLIC -Wall -pedantic)

	target_link_libraries(graphene_library
		PUBLIC
			${Boost_LIBRARIES}
			coverage_config
			Eigen3::Eigen
			libobscura
			libphysica)

	target_include_directories(graphene_library
		PRIVATE
			${GENERATED_DIR}
		PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}
			${INCLUDE_DIR})

	install(TARGETS graphene_library DESTINATION ${LIB_DIR})

#####################################################################
# # OPTION B: Build only an executable.
# 	# Executable
# 	add_executable(executable
# 		main.cpp
# 		module.cpp)

# 	target_compile_options(executable PUBLIC -Wall -pedantic)

# 	target_link_libraries(executable
# 		PUBLIC
# 			coverage_config
# 			libobscura)

# 	target_include_directories(executable
# 		PRIVATE
# 			${GENERATED_DIR} )

# 	install(TARGETS executable DESTINATION ${BIN_DIR})
