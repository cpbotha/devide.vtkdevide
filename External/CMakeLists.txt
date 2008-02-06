#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SET (External_SRCS
     vtkPVExtractVOI.cxx
     vtkPVGlyphFilter.cxx)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdevideExternal C++ library
#

ADD_LIBRARY (vtkdevideExternal ${External_SRCS})

TARGET_LINK_LIBRARIES (
  vtkdevideExternal
  vtkFiltering
  vtkGraphics
  vtkImaging
  vtkCommon
)


#
# Create the vtkdevideExternal Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL3 (vtkdevideExternalTCL ExternalTCL_SRCS
                 "${External_SRCS}" "")
  ADD_LIBRARY (vtkdevideExternalTCL ${ExternalTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideExternalTCL vtkdevideExternal 
                         vtkdevideCommonTCL vtkFilteringTCL vtkGraphicsTCL)

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideExternal Python library
#

# and this is more like the new vtkMy (200410)

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON3 (vtkdevideExternalPython ExternalPython_SRCS 
                    "${External_SRCS}")
  ADD_LIBRARY (vtkdevideExternalPythonD ${ExternalPython_SRCS})
  ADD_LIBRARY (vtkdevideExternalPython MODULE
                    vtkdevideExternalPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtkdevideExternalPython vtkdevideExternalPythonD)
  TARGET_LINK_LIBRARIES (vtkdevideExternalPythonD 
                         vtkdevideCommonPythonD 
			 vtkFilteringPythonD
                         vtkGraphicsPythonD
			 vtkdevideExternal)

  # on Python 2.5 on Windows does not recognise DLL extensions
  # anymore, only PYD !  Previous versions support this.
  IF(WIN32 AND NOT CYGWIN)
      SET_TARGET_PROPERTIES(vtkdevideExternalPython 
          PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)


ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideExternal Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA3 (vtkdevideExternalJava ExternalJava_SRCS "${External_SRCS}")
  ADD_LIBRARY (vtkdevideExternalJava SHARED ${ExternalJava_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideExternalJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

