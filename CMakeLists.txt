cmake_minimum_required(VERSION 3.12)
project(Biosaur2TOPP)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMS
find_package(OpenMS REQUIRED)

# Include OpenMS directories
include_directories(${OpenMS_INCLUDE_DIRS})

# Add the executable
add_executable(Biosaur2 src/topp/Biosaur2.cpp)

# Link against OpenMS
target_link_libraries(Biosaur2 ${OpenMS_LIBRARIES})

# Install target
install(TARGETS Biosaur2
        RUNTIME DESTINATION bin)

# Set compiler flags
if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# Add optimization flags for release builds
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
