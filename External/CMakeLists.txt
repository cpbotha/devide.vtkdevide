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

ADD_LIBRARY (vtkdevideExternal External_SRCS)

#
# Create the vtkdevideExternal Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL (vtkdevideExternalTCL ExternalTCL_SRCS External_SRCS)
  ADD_LIBRARY (vtkdevideExternalTCL SHARED ExternalTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideExternalTCL
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideExternal Python library
#

# and this is more like the new vtkMy (200410)

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdevideExternalPython ExternalPython_SRCS ${External_SRCS})
  ADD_LIBRARY (vtkdevideExternalPython MODULE ${ExternalPython_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideExternalPython 
                         vtkdevideCommonPython 
			 vtkFilteringPython 
			 vtkdevideExternal)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideExternal Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideExternalJava ExternalJava_SRCS External_SRCS)
  ADD_LIBRARY (vtkdevideExternalJava SHARED ExternalJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideExternalJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdevideExternal
  vtkFiltering
  vtkGraphics
  vtkImaging
  vtkCommon
)
