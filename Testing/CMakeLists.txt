#
# Source files
#
# Here is where you can add the name of your local Testing classes.
#

SOURCE_FILES (Testing_SRCS
	vtkImagePlaneWidget2
)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdscasTesting C++ library
#

ADD_LIBRARY (vtkdscasTesting Testing_SRCS)

#
# Create the vtkdscasTesting Tcl library
#

IF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

  VTK_WRAP_TCL (vtkdscasTestingTCL TestingTCL_SRCS Testing_SRCS)
  ADD_LIBRARY (vtkdscasTestingTCL SHARED TestingTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasTestingTCL 
                         vtkRenderingTCL
			 vtkFilteringTCL
			 vtkCommonTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

#
# Create the vtkdscasTesting Python library
#

IF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdscasTestingPython TestingPython_SRCS Testing_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdscasTestingPython MODULE TestingPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdscasTestingPython SHARED TestingPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdscasTestingPython
	                   vtkRenderingPython
			   vtkFilteringPython
			   vtkCommonPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasTestingPython 
	                   vtkRenderingPython
			   vtkFilteringPython
			   vtkCommonPython
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasTestingPython 
	                   vtkRenderingPython
			   vtkFilteringPython
			   vtkCommonPython
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

#
# Create the vtkdscasTesting Java library
#

IF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdscasTestingJava TestingJava_SRCS Testing_SRCS)
  ADD_LIBRARY (vtkdscasTestingJava SHARED TestingJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasTestingJava
                         vtkRenderingJava
			 vtkFilteringJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdscasTesting
  vtkRendering
  vtkFiltering
  vtkCommon
)
