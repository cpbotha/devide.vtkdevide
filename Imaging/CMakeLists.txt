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

ADD_LIBRARY (vtkdevideImaging ${Imaging_SRCS})

TARGET_LINK_LIBRARIES (
  vtkdevideImaging
  vtkImaging
)


#
# Create the vtkdevideImaging Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL3 (vtkdevideImagingTCL ImagingTCL_SRCS 
                 "${Imaging_SRCS}" "")
  ADD_LIBRARY (vtkdevideImagingTCL ${ImagingTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideImagingTCL 
                         vtkdevideImaging
                         vtkImagingTCL)

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideImaging Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON3 (vtkdevideImagingPython ImagingPython_SRCS 
                    "${Imaging_SRCS}")

  ADD_LIBRARY (vtkdevideImagingPythonD ${ImagingPython_SRCS})
  ADD_LIBRARY (vtkdevideImagingPython MODULE vtkdevideImagingPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtkdevideImagingPython
                         vtkdevideImagingPythonD)
  TARGET_LINK_LIBRARIES (vtkdevideImagingPythonD vtkImagingPythonD
                         vtkdevideImaging)
  # on Python 2.5 on Windows does not recognise DLL extensions
  # anymore, only PYD !  Previous versions support this.
  IF(WIN32 AND NOT CYGWIN)
      SET_TARGET_PROPERTIES(vtkdevideImagingPython 
          PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideImaging Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideImagingJava ImagingJava_SRCS "${Imaging_SRCS}")
  ADD_LIBRARY (vtkdevideImagingJava SHARED ${ImagingJava_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideImagingJava
			 vtkImagingJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)
