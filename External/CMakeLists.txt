#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SET (External_SRCS
     vtkPVExtractVOI.cxx)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdevideExternal C++ library
#

ADD_LIBRARY (vtkdevideExternal External_SRCS)

#
# Create the vtkdevideExternal Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL (vtkdevideExternalTCL ExternalTCL_SRCS External_SRCS)
  ADD_LIBRARY (vtkdevideExternalTCL SHARED ExternalTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideExternalTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideExternal Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdevideExternalPython ExternalPython_SRCS External_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdevideExternalPython MODULE ExternalPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdevideExternalPython SHARED ExternalPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdevideExternalPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideExternalPython 
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideExternalPython 
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideExternal Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideExternalJava ExternalJava_SRCS External_SRCS)
  ADD_LIBRARY (vtkdevideExternalJava SHARED ExternalJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideExternalJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdevideExternal
  vtkCommon
)
