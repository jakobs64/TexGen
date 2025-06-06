# ---------------------------------------------------------------------------
#  TexGen ‒ Core/CMakeLists.txt
# ---------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15)
project(TexGenCore)

# ---------------------------------------------------------------------------
# 1.  Gather sources
# ---------------------------------------------------------------------------

# Every .cpp that sits in this directory (flat layout)
file(GLOB CORE_SOURCE_FILES
     *.cpp *.cxx *.cc
     */*.cpp */*.cxx */*.cc)
list(APPEND CORE_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/Mesh.cpp)

# Octree-specific sources are excluded unless the feature is requested
option(BUILD_OCTREE "Enable Octree refinement classes" OFF)
option(TRACE_CORE    "Emit list of Core source files during configure" OFF)

if(NOT BUILD_OCTREE)
    list(REMOVE_ITEM CORE_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/MeshOctreeClasses.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/OctreeVoxelMesh.cpp)
endif()

# C helper that lives one level up
list(APPEND CORE_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/../units/units.c)

# Optional trace so we can see exactly what ends up in the library
if(TRACE_CORE)
    message(STATUS "[DEBUG] Final Core sources:")
    foreach(src ${CORE_SOURCE_FILES})
        message(STATUS "  ${src}")
    endforeach()
endif()

# ---------------------------------------------------------------------------
# 2.  Define the library (STATIC unless the user enabled BUILD_SHARED_LIBS)
# ---------------------------------------------------------------------------

add_library(TexGenCore ${CORE_SOURCE_FILES})

# Public compile flags / definitions
find_package(VTK 9.0 REQUIRED COMPONENTS CommonCore IOXML)

target_compile_definitions(TexGenCore PUBLIC USE_VTK
	 BUILD_MESH USE_MESH_IO BUILD_MESH_IO)
if(WIN32)
    target_compile_definitions(TexGenCore PRIVATE EXPORT)
endif()

# ---------------------------------------------------------------------------
# 3.  Public include directories
# ---------------------------------------------------------------------------

target_include_directories(TexGenCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(BUILD_OCTREE)
    target_include_directories(TexGenCore PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/../OctreeRefinement/include)
endif()

# ---------------------------------------------------------------------------
# 4.  Link dependencies
# ---------------------------------------------------------------------------

target_link_libraries(TexGenCore PUBLIC
    ${VTK_LIBRARIES}
    triangle triangle-api TinyXML CSparse tetgenlib)

if(BUILD_OCTREE)
    target_link_libraries(TexGenCore PUBLIC Octree)
endif()

# ---------------------------------------------------------------------------
# 5.  Installation rules
# ---------------------------------------------------------------------------

include(../Python/FindPythonSitePackages.cmake)

if(WIN32)
    install(TARGETS TexGenCore RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
    install(TARGETS TexGenCore RUNTIME DESTINATION ${PYTHON_SITEPACKAGES_DIR}/TexGen)
    install(TARGETS TexGenCore RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../Installer)
else()
    if(BUILD_SHARED_LIBS)
        install(TARGETS TexGenCore LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
    else()
        install(TARGETS TexGenCore ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
    endif()
endif()