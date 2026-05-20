cmake_minimum_required(VERSION 3.15)

project(PAIPAI VERSION 2.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_executable(mc_paipai
    src/constants.cpp
    src/mcpaipai.cpp
    src/element.cpp
    src/structure.cpp
)

add_executable(findinter
    src/constants.cpp
    src/element.cpp
    src/findinter.cpp
    src/findinter_main.cpp
)

target_include_directories(mc_paipai PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(findinter PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_options(mc_paipai PRIVATE -O2 -pthread)
target_compile_options(findinter PRIVATE -O2)
target_link_libraries(mc_paipai PRIVATE stdc++fs pthread)

install(TARGETS mc_paipai
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(TARGETS findinter
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(PROGRAMS
    scripts/fast_worker.py
    scripts/slow_worker.py
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/paipai.in
    ${CMAKE_CURRENT_BINARY_DIR}/paipai
    @ONLY
)

install(PROGRAMS
    ${CMAKE_CURRENT_BINARY_DIR}/paipai
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)
