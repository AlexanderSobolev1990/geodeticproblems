cmake_minimum_required(VERSION 3.7)
project(geodetic LANGUAGES CXX)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HEADERS
    include/geodesy.h
    include/consts.h
    include/convert.h
    include/compare.h
    include/units.h
)

set(SOURCES
    src/geodesy.cpp
    src/convert.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADERS} ${SOURCES}) # Динамическая библиотека

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}
PUBLIC
    ${ARMADILLO_LIBRARIES}
)

target_compile_options(${PROJECT_NAME} PUBLIC -fPIC) # Position Independent Code

# Включить замеры времени сборки
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CMAKE_COMMAND} -E time")
