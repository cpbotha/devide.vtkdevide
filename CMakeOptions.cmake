#
# Configure output paths for libraries and executables.
#
SET(LIBRARY_OUTPUT_PATH ${VTKDEVIDE_BINARY_DIR}/bin CACHE PATH
    "Single output directory for building all libraries.")
SET(EXECUTABLE_OUTPUT_PATH ${VTKDEVIDE_BINARY_DIR}/bin CACHE PATH
    "Single output directory for building all executables.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)


#
# Try to find VTK and include its settings (otherwise complain)
#
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})


#
# Build shared libs ?
#
# Defaults to the same VTK setting.
#

# Standard CMake option for building libraries shared or static by default.
OPTION(BUILD_SHARED_LIBS
       "Build with shared libraries."
       ${VTK_BUILD_SHARED_LIBS})
# Copy the CMake option to a setting with VTKDEVIDE_ prefix for use in
# our project.  This name is used in vtkmyConfigure.h.in.
SET(VTKDEVIDE_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

# If this is a build tree, provide an option for putting
# this project's executables and libraries in with VTK's.
IF (EXISTS ${VTK_DIR}/bin)
  OPTION(USE_VTK_OUTPUT_PATHS
         "Use VTK's output directory for this project's executables and libraries."
         OFF)
  MARK_AS_ADVANCED (USE_VTK_OUTPUT_PATHS)
  IF (USE_VTK_OUTPUT_PATHS)
    SET (LIBRARY_OUTPUT_PATH ${VTK_DIR}/bin)
    SET (EXECUTABLE_OUTPUT_PATH ${VTK_DIR}/bin)
  ENDIF (USE_VTK_OUTPUT_PATHS)
ENDIF (EXISTS ${VTK_DIR}/bin)



#
# Wrap Tcl, Java, Python
#
# Rational: even if your VTK was wrapped, it does not mean that you want to
# wrap your own local classes.
# Default value is OFF as the VTK cache might have set them to ON but
# the wrappers might not be present (or yet not found).
#

#
# Tcl
#

IF (VTK_WRAP_TCL)

  OPTION(VTKDEVIDE_WRAP_TCL
         "Wrap classes into the TCL interpreted language."
         ON)

  IF(VTKDEVIDE_WRAP_TCL)
    SET(VTK_WRAP_TCL3_INIT_DIR "${VTKDEVIDE_SOURCE_DIR}/Wrapping")
    INCLUDE(${VTK_CMAKE_DIR}/vtkWrapTcl.cmake)
  ENDIF(VTKDEVIDE_WRAP_TCL)

ELSE (VTK_WRAP_TCL)

  IF (VTKDEVIDE_WRAP_TCL)
    MESSAGE("Warning. VTKDEVIDE_WRAP_TCL is ON but the VTK version you have "
            "chosen has not support for Tcl (VTK_WRAP_TCL is OFF).  "
            "Please set VTKDEVIDE_WRAP_TCL to OFF.")
    SET (VTKDEVIDE_WRAP_TCL OFF)
  ENDIF (VTKDEVIDE_WRAP_TCL)

ENDIF (VTK_WRAP_TCL)

#
# Python
#

IF (VTK_WRAP_PYTHON)

  OPTION(VTKDEVIDE_WRAP_PYTHON
         "Wrap classes into the Python interpreted language."
         ON)

  IF (VTKDEVIDE_WRAP_PYTHON)
    SET(VTK_WRAP_PYTHON3_INIT_DIR "${VTKDEVIDE_SOURCE_DIR}/Wrapping")
    INCLUDE(${VTK_CMAKE_DIR}/vtkWrapPython.cmake)
    IF (WIN32)
      IF (NOT BUILD_SHARED_LIBS)
        MESSAGE(FATAL_ERROR "Python support requires BUILD_SHARED_LIBS to be ON.")
        SET (VTKDEVIDE_CAN_BUILD 0)
      ENDIF (NOT BUILD_SHARED_LIBS)
    ENDIF (WIN32)
  ENDIF (VTKDEVIDE_WRAP_PYTHON)

ELSE (VTK_WRAP_PYTHON)

  IF (VTKDEVIDE_WRAP_PYTHON)
    MESSAGE("Warning. VTKDEVIDE_WRAP_PYTHON is ON but the VTK version you have "
            "chosen has not support for Python (VTK_WRAP_PYTHON is OFF).  "
            "Please set VTKDEVIDE_WRAP_PYTHON to OFF.")
    SET (VTKDEVIDE_WRAP_PYTHON OFF)
  ENDIF (VTKDEVIDE_WRAP_PYTHON)

ENDIF (VTK_WRAP_PYTHON)

#
# Java
#

IF (VTK_WRAP_JAVA)

  OPTION(VTKDEVIDE_WRAP_JAVA
         "Wrap classes into the Java interpreted language."
         ON)

  IF (VTKDEVIDE_WRAP_JAVA)
    SET(VTK_WRAP_JAVA3_INIT_DIR "${VTKDEVIDE_SOURCE_DIR}/Wrapping")
    INCLUDE(${VTK_CMAKE_DIR}/vtkWrapJava.cmake)
    IF (WIN32)
      IF (NOT BUILD_SHARED_LIBS)
        MESSAGE(FATAL_ERROR "Java support requires BUILD_SHARED_LIBS to be ON.")
        SET (VTKDEVIDE_CAN_BUILD 0)
      ENDIF (NOT BUILD_SHARED_LIBS)
    ENDIF (WIN32)

    # Tell the java wrappers where to go.
    SET(VTK_JAVA_HOME ${VTKDEVIDE_BINARY_DIR}/java/vtkmy)
    MAKE_DIRECTORY(${VTK_JAVA_HOME})
  ENDIF (VTKDEVIDE_WRAP_JAVA)

ELSE (VTK_WRAP_JAVA)

  IF (VTKDEVIDE_WRAP_JAVA)
    MESSAGE("Warning. VTKDEVIDE_WRAP_JAVA is ON but the VTK version you have "
            "chosen has not support for Java (VTK_WRAP_JAVA is OFF).  "
            "Please set VTKDEVIDE_WRAP_JAVA to OFF.")
    SET (VTKDEVIDE_WRAP_JAVA OFF)
  ENDIF (VTKDEVIDE_WRAP_JAVA)

ENDIF (VTK_WRAP_JAVA)



# let's also get DCMTK in here...

FIND_PATH(DCMTK_INCLUDE_PATH dcmdata/dctk.h
/usr/include
/usr/local/include
)

FIND_PATH(DCMTK_LIB_PATH libdcmdata.a
/usr/lib
/usr/local/lib
)

IF(${DCMTK_INCLUDE_PATH})
  IF(${DCMTK_LIB_PATH})
     SET (HAS_DCMTK 1 CACHE INTERNAL "DCMTK available.")
  ENDIF(${DCMTK_LIB_PATH})
ENDIF(${DCMTK_INCLUDE_PATH})
