#
# Try to find VTK and include its settings (otherwise complain)
#

INCLUDE (${CMAKE_ROOT}/Modules/FindVTK.cmake)

IF (USE_VTK_FILE)
  INCLUDE (${USE_VTK_FILE})
ELSE (USE_VTK_FILE)
  SET (VTKDEVIDE_CAN_BUILD 0)
ENDIF (USE_VTK_FILE)

#
# Try to find ITK
#INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
#IF (USE_ITK_FILE)
#  INCLUDE (${USE_ITK_FILE})
#ELSE (USE_ITK_FILE)
#  SET (VTKDEVIDE_CAN_BUILD 0)
#ENDIF (USE_ITK_FILE)

#
# Build shared libs ?
#
# Defaults to the same VTK setting.
#

IF (USE_VTK_FILE)

  OPTION(BUILD_SHARED_LIBS 
         "Build with shared libraries." 
         ${VTK_BUILD_SHARED_LIBS})
	 
  # This value has to be set so that it can be use in vtkmyConfigure.h.in
  # otherwise the BUILD_SHARED_LIB from VTK's vtkConfigure.h file is picked
  # first :(

  SET(VTKDEVIDE_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS} CACHE INTERNAL 
      "Is this VTKDEVIDE built with shared libraries.")
      
  #
  # Output path(s)
  #

  SET (LIBRARY_OUTPUT_PATH ${VTKDEVIDE_BINARY_DIR}/bin/ CACHE PATH
       "Single output directory for building all libraries.")

  SET (EXECUTABLE_OUTPUT_PATH ${VTKDEVIDE_BINARY_DIR}/bin/ CACHE PATH 
       "Single output directory for building all executables.")
       
  MARK_AS_ADVANCED (
    LIBRARY_OUTPUT_PATH 
    EXECUTABLE_OUTPUT_PATH
  )

ENDIF (USE_VTK_FILE)

#
# Wrap Tcl, Java, Python
#
# Rational: even if your VTK was wrapped, it does not mean that you want to 
# wrap your own local classes. 
# Default value is OFF as the VTK cache might have set them to ON but 
# the wrappers might not be present (or yet not found).
#

#
# Python
# 

IF (VTK_WRAP_PYTHON)

  OPTION(VTKDEVIDE_WRAP_PYTHON 
         "Wrap classes into the Python interpreted language." 
         ON)
         
  IF (VTKDEVIDE_WRAP_PYTHON)

    IF (NOT VTK_WRAP_PYTHON_EXE)

      MESSAGE("Error. Unable to find VTK_WRAP_PYTHON_EXE, please edit this value to specify the correct location of the VTK Python wrapper.")
      MARK_AS_ADVANCED(CLEAR VTK_WRAP_PYTHON_EXE)
      SET (VTKDEVIDE_CAN_BUILD 0)

    ELSE (NOT VTK_WRAP_PYTHON_EXE)

      FIND_FILE(VTK_WRAP_HINTS hints ${VTKDEVIDE_SOURCE_DIR}/Wrapping )
      MARK_AS_ADVANCED(VTK_WRAP_HINTS)

      IF (USE_INSTALLED_VTK)
        INCLUDE (${CMAKE_ROOT}/Modules/FindPythonLibs.cmake)
      ENDIF (USE_INSTALLED_VTK)

      IF (PYTHON_INCLUDE_PATH)
        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
      ENDIF (PYTHON_INCLUDE_PATH)

      IF (WIN32)
        IF (NOT BUILD_SHARED_LIBS)
          MESSAGE("Error. Python support requires BUILD_SHARED_LIBS to be ON.")
          SET (VTKDEVIDE_CAN_BUILD 0)
        ENDIF (NOT BUILD_SHARED_LIBS)  
      ENDIF (WIN32)

    ENDIF (NOT VTK_WRAP_PYTHON_EXE)
  ENDIF (VTKDEVIDE_WRAP_PYTHON)

ELSE (VTK_WRAP_PYTHON)

  IF (VTKDEVIDE_WRAP_PYTHON)
    MESSAGE("Warning. VTKDEVIDE_WRAP_PYTHON is ON but the VTK version you have chosen has not support for Python (VTK_WRAP_PYTHON is OFF). Please set VTKDEVIDE_WRAP_PYTHON to OFF.")
    SET (VTKDEVIDE_WRAP_PYTHON_OFF)
  ENDIF (VTKDEVIDE_WRAP_PYTHON)

ENDIF (VTK_WRAP_PYTHON)

# let's also get DCMTK in here...

FIND_PATH(DCMTK_INCLUDE_PATH dcmdata/dctk.h
/usr/include
/usr/local/include
)

FIND_PATH(DCMTK_LIB_PATH libdcmdata.a
/usr/lib
/usr/local/lib
)

IF(DCMTK_INCLUDE_PATH)
  IF(DCMTK_LIB_PATH)
     SET (HAS_DCMTK 1 CACHE INTERNAL "DCMTK available.")
  ENDIF(DCMTK_LIB_PATH)
ENDIF(DCMTK_INCLUDE_PATH)
