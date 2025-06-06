# ---------------------------------------------------------------------------
#  TexGen – Python/CMakeLists.txt
# ---------------------------------------------------------------------------

if (BUILD_PYTHON_INTERFACE)

    include(FindPythonSitePackages.cmake)

    # -----------------------------------------------------------------------
    # 1.  Collect headers that SWIG depends on
    # -----------------------------------------------------------------------
    file(GLOB CORE_HEADER_FILES
         ../Core/*.h)

    # -----------------------------------------------------------------------
    # 2.  Locate SWIG and the chosen Python
    # -----------------------------------------------------------------------
    find_package(SWIG REQUIRED)
    find_package(PythonLibs REQUIRED)

    include_directories(
        ${PYTHON_INCLUDE_PATH}
        ${CMAKE_CURRENT_SOURCE_DIR}/../OctreeRefinement/include
    )

    # -----------------------------------------------------------------------
    # 3.  Generate the SWIG wrapper
    # -----------------------------------------------------------------------
    add_custom_command(
        OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/Core_wrap.cxx"
        COMMAND "${SWIG_EXECUTABLE}" -c++ -python
                -outdir "${CMAKE_CURRENT_SOURCE_DIR}"
                "${CMAKE_CURRENT_SOURCE_DIR}/Core.i"
        MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/Core.i"
        DEPENDS ${CORE_HEADER_FILES}
        COMMENT "Generating SWIG wrapper code..."
    )

    # CMP0003: link interface libraries via plain targets
    if (COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
    endif ()

    # -----------------------------------------------------------------------
    # 4.  Build the Python module _Core
    # -----------------------------------------------------------------------
    add_library(_Core MODULE Core.i Core_wrap.cxx)
    set_source_files_properties(Core_wrap.cxx PROPERTIES GENERATED 1)

    # Force-load the static archive on macOS so every symbol is present
    #if (APPLE)
     #   set(TEXGEN_CORE_ARCHIVE "${CMAKE_BINARY_DIR}/libTexGenCore.a")
      #  target_link_options(_Core PRIVATE
      #      "LINKER:-force_load,${TEXGEN_CORE_ARCHIVE}"
      #  )
    #endif ()

    # Make sure we build TexGenCore first
    add_dependencies(_Core TexGenCore)

    set_target_properties(_Core PROPERTIES
        PREFIX ""
        DEBUG_POSTFIX "_d"
    )

    if (WIN32)
        set_target_properties(_Core PROPERTIES SUFFIX ".pyd")
    endif ()

    if (UNIX)
        set_target_properties(_Core PROPERTIES
            INSTALL_RPATH_USE_LINK_PATH 1
            BUILD_WITH_INSTALL_RPATH   1
            INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
        )
    endif ()

    target_link_libraries(_Core
        PRIVATE
            TexGenCore
            ${VTK_LIBRARIES}
        optimized ${PYTHON_LIBRARIES}
    )

    if (PYTHON_DEBUG_LIBRARIES)
        target_link_libraries(_Core debug ${PYTHON_DEBUG_LIBRARIES})
    endif ()

    # -----------------------------------------------------------------------
    # 5.  Package pure-Python helpers
    # -----------------------------------------------------------------------
    # Core package initialiser
    set(PYTHON_MODULES
        "'Abaqus', 'Ansys', 'Completer', 'FlowTex', 'GridFile', 'WiseTex'"
    )

    if (BUILD_RENDERER)
        set(PYTHON_MODULES "${PYTHON_MODULES}, 'Renderer'")
    endif ()

    if (BUILD_CASCADE_EXPORT)
        set(PYTHON_MODULES "${PYTHON_MODULES}, 'Export'")
    endif ()

    file(WRITE "__init__.py" "# __all__ = [${PYTHON_MODULES}]")

    install(
        FILES __init__.py Core.py
        DESTINATION ${PYTHON_SITEPACKAGES_DIR}/TexGen
    )

    file(GLOB PYTHON_MODULE_FILES Lib/*.py)
    install(
        FILES ${PYTHON_MODULE_FILES}
        DESTINATION ${PYTHON_SITEPACKAGES_DIR}/TexGen
    )

    # Sample scripts
    file(GLOB PYTHON_SAMPLE_SCRIPTS Scripts/*.py)

    install(
        TARGETS _Core
        LIBRARY DESTINATION ${PYTHON_SITEPACKAGES_DIR}/TexGen
    )

    if (WIN32)
        install(
            FILES ${PYTHON_SAMPLE_SCRIPTS}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/Scripts
        )
        install(
            TARGETS _Core
            LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../Installer
        )
    else ()
        install(
            FILES ${PYTHON_SAMPLE_SCRIPTS}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/share/TexGen/Scripts
        )
    endif ()

    # -----------------------------------------------------------------------
    # 6.  Optional modules (Renderer / Export)
    # -----------------------------------------------------------------------
    # ---------- Renderer ----------
    if (BUILD_RENDERER)

        file(GLOB RENDERER_HEADER_FILES ../Renderer/*.h)
        find_package(VTK 6.3 REQUIRED NO_MODULE)
        include_directories(${VTK_INCLUDE_DIRS})

        add_custom_command(
            OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/Renderer_wrap.cxx"
            COMMAND "${SWIG_EXECUTABLE}" -c++ -python
                    -outdir "${CMAKE_CURRENT_SOURCE_DIR}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/Renderer.i"
            MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/Renderer.i"
            DEPENDS ${RENDERER_HEADER_FILES}
            COMMENT "Generating SWIG wrapper code (Renderer)…"
        )

        add_library(_Renderer MODULE Renderer.i Renderer_wrap.cxx)
        set_source_files_properties(Renderer_wrap.cxx PROPERTIES GENERATED 1)

        set_target_properties(_Renderer PROPERTIES
            PREFIX ""
            DEBUG_POSTFIX "_d"
        )
        if (WIN32)
            set_target_properties(_Renderer PROPERTIES SUFFIX ".pyd")
        endif ()
        if (UNIX)
            set_target_properties(_Renderer PROPERTIES
                COMPILE_FLAGS "-Wno-deprecated"
                INSTALL_RPATH_USE_LINK_PATH 1
                BUILD_WITH_INSTALL_RPATH   1
                INSTALL_RPATH "${VTK_LIBRARY_DIRS};${CMAKE_INSTALL_PREFIX}/lib"
            )
        endif ()

        target_link_libraries(_Renderer
            PRIVATE
                TexGenRenderer
                ${VTK_LIBRARIES}
            optimized ${PYTHON_LIBRARIES}
        )
        if (PYTHON_DEBUG_LIBRARIES)
            target_link_libraries(_Renderer debug ${PYTHON_DEBUG_LIBRARIES})
        endif ()

        install(
            TARGETS _Renderer
            LIBRARY DESTINATION ${PYTHON_SITEPACKAGES_DIR}/TexGen
        )
        install(
            FILES Renderer.py
            DESTINATION ${PYTHON_SITEPACKAGES_DIR}/TexGen
        )
        if (WIN32)
            install(
                TARGETS _Renderer
                LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../Installer
            )
        endif ()
    endif ()

    # ---------- Exporter ----------
    if (BUILD_CASCADE_EXPORT)

        set(EXPORT_HEADER_FILES ../Export/Exporter.h)

        add_custom_command(
            OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/Export_wrap.cxx"
            COMMAND "${SWIG_EXECUTABLE}" -c++ -python
                    -outdir "${CMAKE_CURRENT_SOURCE_DIR}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/Export.i"
            MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/Export.i"
            DEPENDS ${EXPORT_HEADER_FILES}
            COMMENT "Generating SWIG wrapper code (Export)…"
        )

        add_library(_Export MODULE Export.i Export_wrap.cxx)
        set_source_files_properties(Export_wrap.cxx PROPERTIES GENERATED 1)

        set_target_properties(_Export PROPERTIES
            PREFIX ""
            DEBUG_POSTFIX "_d"
        )
        if (WIN32)
            set_target_properties(_Export PROPERTIES SUFFIX ".pyd")
        endif ()
        if (UNIX)
            set_target_properties(_Export PROPERTIES
                COMPILE_FLAGS "-Wno-deprecated"
                INSTALL_RPATH_USE_LINK_PATH 1
                BUILD_WITH_INSTALL_RPATH   1
                INSTALL_RPATH "${OPENCASCADE_LIBRARY_DIR};${CMAKE_INSTALL_PREFIX}/lib"
            )
        endif ()

        target_link_libraries(_Export
            PRIVATE TexGenExport
            optimized ${PYTHON_LIBRARIES}
        )
        if (PYTHON_DEBUG_LIBRARIES)
            target_link_libraries(_Export debug ${PYTHON_DEBUG_LIBRARIES})
        endif ()

        install(
            TARGETS _Export
            LIBRARY DESTINATION ${PYTHON_SITEPACKAGES_DIR}/TexGen
        )
        install(
            FILES Export.py
            DESTINATION ${PYTHON_SITEPACKAGES_DIR}/TexGen
        )
        if (WIN32)
            install(
                TARGETS _Export
                LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../Installer
            )
        endif ()
    endif ()

endif ()  # BUILD_PYTHON_INTERFACE