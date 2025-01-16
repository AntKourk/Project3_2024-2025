cmake_minimum_required(VERSION 3.1...3.23)

project(opt_triangulation)

# Set the build type
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# CGAL and its components with Qt5
find_package(CGAL REQUIRED COMPONENTS Qt5)

# Boost and its components
find_package(Boost REQUIRED COMPONENTS system filesystem)

# Qt5
find_package(Qt5 REQUIRED COMPONENTS Widgets)

# OpenGL
find_package(OpenGL REQUIRED)

# GLEW
find_package(GLEW REQUIRED)

# Check if CGAL is found
if (NOT CGAL_FOUND)
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
endif()

# Check if Boost is found
if (NOT Boost_FOUND)
    message(STATUS "This project requires the Boost library, and will not be compiled.")
    return()
endif()

# Create the executable with both source files
add_executable(opt_triangulation center.cpp projection.cpp centroid.cpp circumcenter.cpp output.cpp flipEdges.cpp inputs.cpp main.cpp inside_convex_polygon_centroid.cpp utils.cpp local_search.cpp simulated_annealing.cpp ant_colony.cpp)

# Link the executable to CGAL, Boost, Qt5, OpenGL, and GLEW libraries
target_link_libraries(opt_triangulation 
    PRIVATE 
        CGAL::CGAL 
        Boost::system 
        Boost::filesystem 
        Qt5::Widgets 
        OpenGL::GL 
        GLEW::GLEW
)

# If CGAL Qt5 is found, define CGAL_USE_BASIC_VIEWER and link to CGAL_Qt5
if (CGAL_Qt5_FOUND)
    add_definitions(-DCGAL_USE_BASIC_VIEWER)
    target_link_libraries(opt_triangulation PRIVATE CGAL::CGAL_Qt5)
endif()
