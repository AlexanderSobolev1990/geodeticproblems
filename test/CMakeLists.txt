cmake_minimum_required(VERSION 3.7)
project(test_geodetic LANGUAGES CXX)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}")

enable_testing(true)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt5Test REQUIRED)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

add_executable(${PROJECT_NAME} test_geodetic.cpp)
add_test(NAME test_geodetic COMMAND test_geodetic)

target_link_libraries(${PROJECT_NAME}
    geodetic
    Qt5::Test
    )
