#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SET (Hybrid_SRCS 
     vtkBoxWidgetConstrained.cxx)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdscasRendering C++ library
#

ADD_LIBRARY (vtkdscasHybrid Hybrid_SRCS)

#
# Create the vtkdscasHybrid Tcl library
#

IF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

  VTK_WRAP_TCL (vtkdscasHybridTCL HybridTCL_SRCS Hybrid_SRCS)
  ADD_LIBRARY (vtkdscasHybridTCL SHARED HybridTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasHybridTCL 
                         vtkHybridTCL
			 vtkCommonTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

#
# Create the vtkdscasHybrid Python library
#

IF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdscasHybridPython HybridPython_SRCS Hybrid_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdscasHybridPython MODULE HybridPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdscasHybridPython SHARED HybridPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdscasHybridPython
	                   vtkHybridPython
			   vtkCommonPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasHybridPython 
	                   vtkHybridPython
			   vtkCommonPython
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasHybridPython 
	                   vtkHybridPython
			   vtkCommonPython
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

#
# Create the vtkdscasHybrid Java library
#

IF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdscasHybridJava HybridJava_SRCS Hybrid_SRCS)
  ADD_LIBRARY (vtkdscasHybridJava SHARED HybridJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasHybridJava
                         vtkHybridJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdscasHybrid
  vtkHybrid
  vtkCommon
)
