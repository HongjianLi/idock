# Requires cmake v3.20 for Boost 1.75.0
cmake_minimum_required(VERSION 3.20)
project(jdock LANGUAGES CXX)

# Enforce C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Create the target and add source files
add_executable(${PROJECT_NAME}
  src/array.cpp
  src/io_service_pool.cpp
  src/main.cpp
  src/random_forest.cpp
  src/random_forest_y.cpp
  src/residue.cpp
  src/safe_counter.cpp
  src/stopwatch.cpp
  src/atom.cpp
  src/ligand.cpp
  src/pka.cpp
  src/random_forest_x.cpp
  src/receptor.cpp
  src/result.cpp
  src/scoring_function.cpp
)

# https://cmake.org/cmake/help/latest/module/FindThreads.html
# Use posix thread lib if the system doesn't provide the thread functions
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

# https://cmake.org/cmake/help/latest/module/FindBoost.html
# Allow running without Boost .dll or .so
set(Boost_USE_STATIC_LIBS TRUE)

# Since Boost binaries compiled with static C++ runtime are only downloadable on Windows,
#   here we have to limit to MSVC or other build tools will throw the following error:
#   "No suitable build variant has been found."
# We can get rid of this error by building our own Boost binaries but that will take a lot
#   more of time for the GitHub workflow.
# For GCC we will enforce static linking everything later in this file.
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # Allow running without C++ Runtime
  set(Boost_USE_STATIC_RUNTIME TRUE)
endif()

find_package(Boost REQUIRED COMPONENTS
  program_options
)

# Set include path for the target only
target_include_directories(${PROJECT_NAME} PRIVATE
  ${Boost_INCLUDE_DIRS}
)

# Set lib path for the target only
target_link_libraries(${PROJECT_NAME}
  Threads::Threads
  Boost::program_options
)

# Setup static linking C++ runtime for GCC and MSVC runtime
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # using GCC
  # Static linking of libgcc and libstdc++ is not enough because some Linux distro like Alpine
  #   ships with an older version of libc, which cannot opt-in for static linking individually.
  target_link_options(${PROJECT_NAME} PRIVATE
    -static
  )
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # using Visual Studio C++
  set_property(TARGET ${PROJECT_NAME} PROPERTY
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"
  )
endif()

# Enable cmake --install to copy the binary to system dir
install(
  TARGETS ${PROJECT_NAME}
)