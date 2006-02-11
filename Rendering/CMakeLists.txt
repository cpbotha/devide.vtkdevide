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

ADD_LIBRARY (vtkdevideRendering ${Rendering_SRCS})

TARGET_LINK_LIBRARIES (
  vtkdevideRendering
  vtkRendering
  vtkVolumeRendering
  vtkFiltering
  vtkCommon
)


#
# Create the vtkdevideRendering Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL3 (vtkdevideRenderingTCL RenderingTCL_SRCS 
                 "${Rendering_SRCS}" "")
  ADD_LIBRARY (vtkdevideRenderingTCL ${RenderingTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideRenderingTCL 
                         vtkRenderingTCL
			 vtkVolumeRenderingTCL
			 vtkFilteringTCL
			 vtkCommonTCL)

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideRendering Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON3 (vtkdevideRenderingPython RenderingPython_SRCS 
                    "${Rendering_SRCS}")

  ADD_LIBRARY (vtkdevideRenderingPythonD ${RenderingPython_SRCS})
  ADD_LIBRARY (vtkdevideRenderingPython MODULE 
               vtkdevideRenderingPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtkdevideRenderingPython
                         vtkdevideRenderingPythonD)
  TARGET_LINK_LIBRARIES (vtkdevideRenderingPythonD
	                   vtkRenderingPythonD
			   vtkVolumeRenderingPythonD
			   vtkFilteringPythonD
			   vtkCommonPythonD)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideRendering Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA3 (vtkdevideRenderingJava RenderingJava_SRCS ${Rendering_SRCS})
  ADD_LIBRARY (vtkdevideRenderingJava SHARED ${RenderingJava_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideRenderingJava
                         vtkRenderingJava
			 #vtkVolumeRenderingPython
			 vtkFilteringJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

