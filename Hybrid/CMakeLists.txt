#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SET (Hybrid_SRCS 
     vtkBoxWidgetConstrained.cxx
     vtkImageBacktracker.cxx
     vtkPolyLineWidget.cxx)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdevideRendering C++ library
#

ADD_LIBRARY (vtkdevideHybrid Hybrid_SRCS)

#
# Create the vtkdevideHybrid Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL (vtkdevideHybridTCL HybridTCL_SRCS Hybrid_SRCS)
  ADD_LIBRARY (vtkdevideHybridTCL SHARED HybridTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideHybridTCL 
                         vtkHybridTCL
			 vtkCommonTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideHybrid Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdevideHybridPython HybridPython_SRCS Hybrid_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdevideHybridPython MODULE HybridPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdevideHybridPython SHARED HybridPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdevideHybridPython
	                   vtkHybridPython
			   vtkCommonPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideHybridPython 
	                   vtkHybridPython
			   vtkCommonPython
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideHybridPython 
	                   vtkHybridPython
			   vtkCommonPython
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideHybrid Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideHybridJava HybridJava_SRCS Hybrid_SRCS)
  ADD_LIBRARY (vtkdevideHybridJava SHARED HybridJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideHybridJava
                         vtkHybridJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdevideHybrid
  vtkHybrid
  vtkCommon
)
