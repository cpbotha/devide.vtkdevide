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
# Create the vtkdevideTesting C++ library
#

ADD_LIBRARY (vtkdevideTesting Testing_SRCS)

#
# Create the vtkdevideTesting Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL (vtkdevideTestingTCL TestingTCL_SRCS Testing_SRCS)
  ADD_LIBRARY (vtkdevideTestingTCL SHARED TestingTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideTestingTCL 
                         vtkRenderingTCL
			 vtkFilteringTCL
			 vtkCommonTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideTesting Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdevideTestingPython TestingPython_SRCS Testing_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdevideTestingPython MODULE TestingPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdevideTestingPython SHARED TestingPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdevideTestingPython
	                   vtkRenderingPython
			   vtkFilteringPython
			   vtkCommonPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideTestingPython 
	                   vtkRenderingPython
			   vtkFilteringPython
			   vtkCommonPython
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideTestingPython 
	                   vtkRenderingPython
			   vtkFilteringPython
			   vtkCommonPython
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideTesting Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideTestingJava TestingJava_SRCS Testing_SRCS)
  ADD_LIBRARY (vtkdevideTestingJava SHARED TestingJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideTestingJava
                         vtkRenderingJava
			 vtkFilteringJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdevideTesting
  vtkRendering
  vtkFiltering
  vtkCommon
)
