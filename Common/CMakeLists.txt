#
# Source files
#
# Here is where you can add the name of your local Rendering classes.
#

SET (Common_SRCS vtkEventOutputWindow.cxx)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdevideCommon C++ library
#

ADD_LIBRARY (vtkdevideCommon ${Common_SRCS})

TARGET_LINK_LIBRARIES (
  vtkdevideCommon
  vtkCommon
)

#
# Create the vtkdevideCommon Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL3 (vtkdevideCommonTCL CommonTCL_SRCS "${Common_SRCS}" "")
  ADD_LIBRARY (vtkdevideCommonTCL ${CommonTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideCommonTCL vtkdevideCommon vtkCommonTCL)

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideCommon Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON3 (vtkdevideCommonPython CommonPython_SRCS
                    "${Common_SRCS}")
  # the D python lib contains the actual code
  ADD_LIBRARY (vtkdevideCommonPythonD ${CommonPython_SRCS})
  # and the non-D is just the Python import stub thingy
  ADD_LIBRARY (vtkdevideCommonPython MODULE
               vtkdevideCommonPythonInit.cxx)
  # stub thingy has to link with the true contents of course
  TARGET_LINK_LIBRARIES (vtkdevideCommonPython vtkdevideCommonPythonD)
  TARGET_LINK_LIBRARIES (vtkdevideCommonPythonD vtkCommonPythonD)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideCommon Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA3 (vtkdevideCommonJava CommonJava_SRCS "${Common_SRCS}")
  ADD_LIBRARY (vtkdevideCommonJava SHARED ${CommonJava_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideCommonJava
			 vtkCommonJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

