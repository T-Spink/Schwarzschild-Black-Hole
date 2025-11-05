cmake_minimum_required(VERSION 3.15)
project(BlackHole)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

# add include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include/
    ${PROJECT_SOURCE_DIR}/integrator/
)

# Add source files
set(SOURCES
    ${PROJECT_SOURCE_DIR}/main.cpp
    ${PROJECT_SOURCE_DIR}/src/geodesics.cpp
    ${PROJECT_SOURCE_DIR}/src/constants.cpp
)

add_executable(BlackHole ${SOURCES})

target_link_libraries(BlackHole OpenGL::GL glfw glm Eigen3::Eigen)


