cmake_minimum_required(VERSION 3.12)
project(PAIPAI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(mc_paipai
    ${SRC_DIR}/mcpaipai.cpp
    ${SRC_DIR}/constants.cpp
    ${SRC_DIR}/element.cpp
    ${SRC_DIR}/structure.cpp

)

target_include_directories(mc_paipai
    PRIVATE
        ${SRC_DIR}
)

find_package(Threads REQUIRED)
target_link_libraries(mc_paipai PRIVATE Threads::Threads)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
        target_link_libraries(mc_paipai PRIVATE stdc++fs)
    endif()
endif()


include(GNUInstallDirs) 
find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(PAIPAI_SCRIPTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/scripts)

configure_file(
    ${PAIPAI_SCRIPTS_DIR}/paipai.in
    ${CMAKE_CURRENT_BINARY_DIR}/paipai
    @ONLY
)

install(TARGETS mc_paipai
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(PROGRAMS
    ${CMAKE_CURRENT_BINARY_DIR}/paipai
    ${PAIPAI_SCRIPTS_DIR}/fast_worker.py
    ${PAIPAI_SCRIPTS_DIR}/slow_worker.py
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)