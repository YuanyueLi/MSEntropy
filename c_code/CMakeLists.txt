# CMakeList.txt : CMake project for SpectralEntropy_C, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("SpectralEntropy")

# Add source to this project's executable.
ADD_LIBRARY(SpectralEntropy STATIC "SpectralEntropy.c" "SpectralEntropy.h" "CleanSpectrum.c" "CleanSpectrum.h")
add_executable (Example "Example.c")
TARGET_LINK_LIBRARIES(Example PRIVATE SpectralEntropy m)

set_target_properties(Example PROPERTIES
    LINK_SEARCH_START_STATIC 1
    LINK_SEARCH_END_STATIC 1
)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")


if (CMAKE_VERSION VERSION_GREATER 3.12)
    set_property(TARGET Example PROPERTY C_STANDARD 99)
endif()
