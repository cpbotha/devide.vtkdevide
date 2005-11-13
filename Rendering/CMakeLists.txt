#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SET (Rendering_SRCS 
     vtkShellExtractor.cxx
     vtkOpenGLVolumeShellSplatMapper.cxx
     vtkInteractorStyleTrackballActorConstrained.cxx)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdevideRendering C++ library
#

ADD_LIBRARY (vtkdevideRendering Rendering_SRCS)

#
# Create the vtkdevideRendering Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL (vtkdevideRenderingTCL RenderingTCL_SRCS Rendering_SRCS)
  ADD_LIBRARY (vtkdevideRenderingTCL SHARED RenderingTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideRenderingTCL 
                         vtkRenderingTCL
			 vtkVolumeRenderingTCL
			 vtkFilteringTCL
			 vtkCommonTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideRendering Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdevideRenderingPython RenderingPython_SRCS Rendering_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdevideRenderingPython MODULE RenderingPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdevideRenderingPython SHARED RenderingPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdevideRenderingPython
	                   vtkRenderingPythonD
			   vtkVolumeRenderingPythonD
			   vtkFilteringPythonD
			   vtkCommonPythonD)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideRenderingPython 
	                   vtkRenderingPythonD
			   vtkVolumeRenderingPythonD
			   vtkFilteringPythonD
			   vtkCommonPythonD
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideRenderingPython 
	                   vtkRenderingPythonD
			   vtkVolumeRenderingPythonD
			   vtkFilteringPythonD
			   vtkCommonPythonD
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideRendering Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideRenderingJava RenderingJava_SRCS Rendering_SRCS)
  ADD_LIBRARY (vtkdevideRenderingJava SHARED RenderingJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideRenderingJava
                         vtkRenderingJava
			 #vtkVolumeRenderingPython
			 vtkFilteringJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdevideRendering
  vtkRendering
  vtkVolumeRendering
  vtkFiltering
  vtkCommon
)
