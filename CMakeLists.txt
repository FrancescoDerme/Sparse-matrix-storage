cmake_minimum_required(VERSION 3.8)
project(SparseMatrix)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Change path from /src if needed, or add more directories
file(GLOB_RECURSE sources
        "${CMAKE_SOURCE_DIR}/src/*.cpp"
        "${CMAKE_SOURCE_DIR}/tests/*.cpp"
        )
# Add precompiler definitions like that:
# add_definitions(-DSOME_DEFINITION)

add_executable(executable ${sources})

# Add more include directories if needed
target_include_directories(executable PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_include_directories(executable PUBLIC "${CMAKE_SOURCE_DIR}/tests")
target_include_directories(executable PUBLIC "${CMAKE_SOURCE_DIR}/include/Compressed")
target_include_directories(executable PUBLIC "${CMAKE_SOURCE_DIR}/include/Dynamic")
target_include_directories(executable PUBLIC "${CMAKE_SOURCE_DIR}/include/Implementation")
target_include_directories(executable PUBLIC "${CMAKE_SOURCE_DIR}/include/Utils")

# If you have precompiled headers you can add them like this
# target_precompiled_headers(my_app PRIVATE "${CMAKE_SOURCE_DIR}/src/pch.h")
