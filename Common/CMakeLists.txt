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
# Create the vtkdscasCommon C++ library
#

ADD_LIBRARY (vtkdscasCommon Common_SRCS)

#
# Create the vtkdscasCommon Tcl library
#

IF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

  VTK_WRAP_TCL (vtkdscasCommonTCL CommonTCL_SRCS Common_SRCS)
  ADD_LIBRARY (vtkdscasCommonTCL SHARED CommonTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasCommonTCL 
                         vtkCommonTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

#
# Create the vtkdscasCommon Python library
#

IF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdscasCommonPython CommonPython_SRCS Common_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdscasCommonPython MODULE CommonPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdscasCommonPython SHARED CommonPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdscasCommonPython
			   vtkCommonPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasCommonPython 
	                   vtkCommonPython
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasCommonPython 
	                   vtkCommonPython
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

#
# Create the vtkdscasCommon Java library
#

IF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdscasCommonJava CommonJava_SRCS Common_SRCS)
  ADD_LIBRARY (vtkdscasCommonJava SHARED CommonJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasCommonJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdscasCommon
  vtkCommon
)
