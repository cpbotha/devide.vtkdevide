#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SOURCE_FILES (Rendering_SRCS
   vtkShellExtractor
   vtkOpenGLVolumeShellSplatMapper
)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdscasRendering C++ library
#

ADD_LIBRARY (vtkdscasRendering Rendering_SRCS)

#
# Create the vtkdscasRendering Tcl library
#

IF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

  VTK_WRAP_TCL (vtkdscasRenderingTCL RenderingTCL_SRCS Rendering_SRCS)
  ADD_LIBRARY (vtkdscasRenderingTCL SHARED RenderingTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasRenderingTCL 
                         vtkRenderingTCL
			 vtkFilteringTCL
			 vtkCommonTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

#
# Create the vtkdscasRendering Python library
#

IF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdscasRenderingPython RenderingPython_SRCS Rendering_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdscasRenderingPython MODULE RenderingPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdscasRenderingPython SHARED RenderingPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdscasRenderingPython
	                   vtkRenderingPython
			   vtkFilteringPython
			   vtkCommonPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasRenderingPython 
	                   vtkRenderingPython
			   vtkFilteringPython
			   vtkCommonPython
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasRenderingPython 
	                   vtkRenderingPython
			   vtkFilteringPython
			   vtkCommonPython
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

#
# Create the vtkdscasRendering Java library
#

IF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdscasRenderingJava RenderingJava_SRCS Rendering_SRCS)
  ADD_LIBRARY (vtkdscasRenderingJava SHARED RenderingJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasRenderingJava
                         vtkRenderingJava
			 vtkFilteringJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdscasRendering
  vtkRendering
  vtkFiltering
  vtkCommon
)
