#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SET (Hybrid_SRCS 
     vtkBoxWidgetConstrained.cxx
     vtkPolyLineWidget.cxx)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdevideRendering C++ library
#

ADD_LIBRARY (vtkdevideHybrid ${Hybrid_SRCS})

TARGET_LINK_LIBRARIES (
  vtkdevideHybrid
  vtkHybrid
  vtkCommon
  vtkRendering
  vtkWidgets
)


#
# Create the vtkdevideHybrid Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL3 (vtkdevideHybridTCL HybridTCL_SRCS "${Hybrid_SRCS}" "")
  ADD_LIBRARY (vtkdevideHybridTCL ${HybridTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideHybridTCL
                         vtkdevideHybrid 
                         vtkHybridTCL
			 vtkCommonTCL
			 vtkWidgetsTCL)

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideHybrid Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON3 (vtkdevideHybridPython HybridPython_SRCS "${Hybrid_SRCS}")
  ADD_LIBRARY (vtkdevideHybridPythonD ${HybridPython_SRCS})
  ADD_LIBRARY (vtkdevideHybridPython MODULE vtkdevideHybridPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtkdevideHybridPython vtkdevideHybridPythonD)
  TARGET_LINK_LIBRARIES (vtkdevideHybridPythonD
	                 vtkWidgetsPythonD
			 vtkCommonPythonD
                         vtkdevideHybrid)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideHybrid Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA3 (vtkdevideHybridJava HybridJava_SRCS "${Hybrid_SRCS}")
  ADD_LIBRARY (vtkdevideHybridJava SHARED ${HybridJava_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideHybridJava
                         vtkHybridJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)
