#
# Source files
#
# Here is where you can add the name of your local IO classes.
#

IF (${HAS_DCMTK})
   SET (IO_SRCS ${IO_SRCS} vtkDICOMVolumeReader.cxx)
ENDIF (${HAS_DCMTK})

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

IF (${HAS_DCMTK})
   # INCLUDE_ and LINK_ dirs **HAVE TO GO BEFORE THE TARGETS THAT THEY
   # APPLY TO**.
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH})
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/ofstd)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/dcmdata)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/dcmjpeg)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/dcmimgle)
   LINK_DIRECTORIES (${DCMTK_LIB_PATH})
ENDIF (${HAS_DCMTK})

#
# Create the vtkdevideIO C++ library
#

IF (${WIN32})
   ADD_LIBRARY (vtkdevideIO ${IO_SRCS})
ELSE (${WIN32})
   # under unix, we also want it to include a built-in dicom.dic
   # under windows, the Visual C++ builds already do this
   ADD_LIBRARY (vtkdevideIO ${IO_SRCS} dcdictzz.cxx)
ENDIF (${WIN32})

IF (${HAS_DCMTK})
   # from dcmimgle onwards is for JPEG support
   TARGET_LINK_LIBRARIES(vtkdevideIO dcmdata ofstd)
   # dcmimgle dcmjpeg ijg12 ijg16 ijg8)
ENDIF (${HAS_DCMTK})


#GET_DIRECTORY_PROPERTY(MY_LL, INCLUDE_DIRECTORIES)
#MESSAGE("${MY_LL}")

IF (${WIN32})
   TARGET_LINK_LIBRARIES(vtkdevideIO wsock32 netapi32)
ENDIF (${WIN32})

TARGET_LINK_LIBRARIES (
  vtkdevideIO
  vtkIO
  vtkCommon
)

#
# Create the vtkdevideIO Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL3 (vtkdevideIOTCL IOTCL_SRCS "${IO_SRCS}" "")
  ADD_LIBRARY (vtkdevideIOTCL SHARED ${IOTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideIOTCL
                         vtkdevideIO
                         vtkIOTCL )


ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideIO Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON3 (vtkdevideIOPython IOPython_SRCS "${IO_SRCS}")
  ADD_LIBRARY (vtkdevideIOPythonD ${IOPython_SRCS})
  ADD_LIBRARY (vtkdevideIOPython MODULE vtkdevideIOPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtkdevideIOPython vtkdevideIOPythonD)
  TARGET_LINK_LIBRARIES (vtkdevideIOPythonD vtkIOPythonD vtkdevideIO)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideIO Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideIOJava IOJava_SRCS "${IO_SRCS}")
  ADD_LIBRARY (vtkdevideIOJava SHARED ${IOJava_SRCS})
  TARGET_LINK_LIBRARIES (vtkdevideIOJava 
                         vtkIOJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)
