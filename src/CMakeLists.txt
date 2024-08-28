cmake_minimum_required(VERSION 3.5)

project(dedup VERSION 0.17 LANGUAGES CXX)

find_package(PkgConfig REQUIRED)
pkg_check_modules(avisynth REQUIRED avisynth)

add_library(dedup SHARED dedup.cpp)

target_include_directories(dedup PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${AVISYNTH_INCLUDE_DIRS})

include(GNUInstallDirs)

INSTALL(TARGETS dedup
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/avisynth")

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
