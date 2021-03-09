# 1. Unit tests: Graphene.xpp
add_executable(test_Graphene
	test_Graphene.cpp)
target_link_libraries(test_Graphene 
	PRIVATE
		graphene_library
		gtest_main ) #contains the main function
target_include_directories(test_Graphene
	PRIVATE
		${GENERATED_DIR} )
target_compile_options(test_Graphene PUBLIC -Wall -pedantic)
# install(TARGETS tests DESTINATION ${TESTS_DIR})
add_test(NAME test_Graphene COMMAND test_Graphene
	WORKING_DIRECTORY ${TESTS_DIR} )

# 2. Unit tests: Hydrogenic_Wavefunctions.xpp
add_executable(test_Hydrogenic_Wavefunctions
	test_Hydrogenic_Wavefunctions.cpp)
target_link_libraries(test_Hydrogenic_Wavefunctions 
	PRIVATE
		graphene_library
		gtest_main ) #contains the main function
target_include_directories(test_Hydrogenic_Wavefunctions
	PRIVATE
		${GENERATED_DIR} )
target_compile_options(test_Hydrogenic_Wavefunctions PUBLIC -Wall -pedantic)
# install(TARGETS tests DESTINATION ${TESTS_DIR})
add_test(NAME test_Hydrogenic_Wavefunctions COMMAND test_Hydrogenic_Wavefunctions
	WORKING_DIRECTORY ${TESTS_DIR} )

# 3. Unit tests: Multidimensional_Integration.xpp
add_executable(test_Multidimensional_Integration
test_Multidimensional_Integration.cpp)
target_link_libraries(test_Multidimensional_Integration 
PRIVATE
	graphene_library
	gtest_main ) #contains the main function
target_include_directories(test_Multidimensional_Integration
PRIVATE
	${GENERATED_DIR} )
target_compile_options(test_Multidimensional_Integration PUBLIC -Wall -pedantic)
# install(TARGETS tests DESTINATION ${TESTS_DIR})
add_test(NAME test_Multidimensional_Integration COMMAND test_Multidimensional_Integration
WORKING_DIRECTORY ${TESTS_DIR} )