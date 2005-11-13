#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SET (Imaging_SRCS
     vtkHistogramLookupTable.cxx
     vtkImageBorderMask.cxx
     vtkImageGreyscaleReconstruct3D.cxx
     vtkImageHistogram2D.cxx
     vtkSelectConnectedComponents.cxx)


# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdevideImaging C++ library
#

ADD_LIBRARY (vtkdevideImaging Imaging_SRCS)

#
# Create the vtkdevideImaging Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL (vtkdevideImagingTCL ImagingTCL_SRCS Imaging_SRCS)
  ADD_LIBRARY (vtkdevideImagingTCL SHARED ImagingTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideImagingTCL 
                         vtkImagingTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideImaging Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdevideImagingPython ImagingPython_SRCS Imaging_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdevideImagingPython MODULE ImagingPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdevideImagingPython SHARED ImagingPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdevideImagingPython
			   vtkImagingPythonD)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideImagingPython 
	                   vtkImagingPythonD
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideImagingPython 
	                   vtkImagingPythonD
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideImaging Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideImagingJava ImagingJava_SRCS Imaging_SRCS)
  ADD_LIBRARY (vtkdevideImagingJava SHARED ImagingJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideImagingJava
			 vtkImagingJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdevideImaging
  vtkImaging
)
