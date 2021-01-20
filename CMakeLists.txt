
cmake_minimum_required(VERSION 2.8)
project(MERIDA_ILP)

set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_CXX_COMPILER g++)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -Wall -pedantic")

find_package(Boost COMPONENTS filesystem regex) 

if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS}) 
else()
	message(SEND_ERROR "Boost was not found on your system")
endif()

set(GT2_INCLUDE libraries)
set(GT2_CORE  libraries/genetrail2/core)
set(GT2_applications  applications)
set(MERIDA_BASE MERIDA_ILP)
set(LOGICALMODEL_BASE LogicalModel)
set(CV_BASE CV)
set(PATHWAY_BASE Pathway)

set(SOURCES
        ${GT2_CORE}/AbstractMatrix.cpp
        ${GT2_CORE}/DenseMatrix.cpp
        ${GT2_CORE}/DenseMatrixReader.cpp
        ${GT2_CORE}/Exception.cpp
        ${GT2_CORE}/DenseMatrixWriter.cpp
        ${GT2_CORE}/MatrixWriter.cpp
        ${GT2_CORE}/DenseMatrixSubset.cpp
        ${GT2_applications}/matrixTools.cpp
	${MERIDA_BASE}/MERIDA_ILP.cpp
	${LOGICALMODEL_BASE}/LogicalModel.cpp
	${PATHWAY_BASE}/Pathway.cpp
        )


set(CV_SOURCES
	${CV_BASE}/CrossValidator.cpp
	${CV_BASE}/LogicalModelCrossValidator.cpp
        )



include_directories(
	${CPLEX_PATH}/cplex/include
	${CPLEX_PATH}/concert/include
        ${GT2_CORE}
        ${GT2_INCLUDE}
        ${GT2_applications}
        /usr/include/boost
	${EIGEN_PATH}
)

link_directories(
	${CPLEX_PATH}/concert/lib/x86-64_linux/static_pic
	${CPLEX_PATH}/cplex/lib/x86-64_linux/static_pic
)


set(LINK_FLAGS
        "-pthread"
)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")

# Create library for cross validation
add_library(CrossValidation_lib SHARED ${CV_SOURCES})
set_target_properties(CrossValidation_lib PROPERTIES CXX_STANDARD 14)

include_directories(
	${CV_BASE}
	${LOGICALMODEL_BASE}
	${PATHWAY_BASE}
)

# Create library (for testing)
add_library(MERIDA_ILP_lib SHARED ${SOURCES})
set_target_properties(MERIDA_ILP_lib PROPERTIES CXX_STANDARD 14)
target_link_libraries(MERIDA_ILP_lib
        CrossValidation_lib
	ilocplex
	concert
	cplex
        -pthread
	${Boost_LIBRARIES}
)
#set(COMPILE_FLAGS
#	${CXX_FLAGS}
#)

# Create executable
add_executable(MERIDA_ILP ${SOURCES} ${MERIDA_BASE}/main.cpp)
set_target_properties(MERIDA_ILP PROPERTIES CXX_STANDARD 14)
target_link_libraries(MERIDA_ILP
        CrossValidation_lib
	ilocplex
	concert
	cplex
        -pthread
	${Boost_LIBRARIES}
)
# Aktiviere die Unittestfunktionalität von CMake.
#enable_testing()

# Füge die Google Test Bibliothek zum Projekt hinzu.
#add_subdirectory(${GTEST_SRC_DIR} "gtest")

# Füge die Unittests zum Project hinzu.
#add_subdirectory(test)
