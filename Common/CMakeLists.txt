#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SET (Common_SRCS 
     vtkEventOutputWindow.cxx)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdevideCommon C++ library
#

ADD_LIBRARY (vtkdevideCommon Common_SRCS)

#
# Create the vtkdevideCommon Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL (vtkdevideCommonTCL CommonTCL_SRCS Common_SRCS)
  ADD_LIBRARY (vtkdevideCommonTCL SHARED CommonTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideCommonTCL 
                         vtkCommonTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideCommon Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdevideCommonPython CommonPython_SRCS Common_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdevideCommonPython MODULE CommonPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdevideCommonPython SHARED CommonPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdevideCommonPython
			   vtkCommonPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideCommonPython 
	                   vtkCommonPython
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideCommonPython 
	                   vtkCommonPython
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideCommon Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideCommonJava CommonJava_SRCS Common_SRCS)
  ADD_LIBRARY (vtkdevideCommonJava SHARED CommonJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideCommonJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdevideCommon
  vtkCommon
)
