#
# Source files
#
# Here is where you can add the name of your local IO classes.
#

# if DCMTK is not found, IO won't be built (there are currently no other
# classes in the IO package)
IF (${HAS_DCMTK})
  SET (IO_SRCS ${IO_SRCS} vtkDICOMVolumeReader.cxx)
ENDIF (${HAS_DCMTK})

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# start setup of vtkdevideIO library
# IO_LIB_SRCS will be a list of objects that need to be linked together
SET(IO_LIB_SRCS ${IO_SRCS})


# it might be better to have the user selecting the directory
# containing all the dcmtk include files, and then to strip off the
# optional dcmtk...
IF (HAS_DCMTK)
   # INCLUDE_ and LINK_ dirs **HAVE TO GO BEFORE THE TARGETS THAT THEY
   # APPLY TO**.
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH})
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/dcmtk)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/dcmtk/config)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/dcmtk/ofstd)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/dcmtk/dcmdata)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/dcmtk/dcmjpeg)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/dcmtk/dcmimgle)
   LINK_DIRECTORIES (${DCMTK_LIB_PATH})

   IF (NOT WIN32)
      # under unix, we also want it to include a built-in dicom.dic
      # under windows, the Visual C++ builds already do this
      SET(IO_LIB_SRCS ${IO_LIB_SRCS} dcdictzz.cxx)
   ELSE (NOT WIN32)
      SET(ENV{HAVE_CONFIG_H} 1)
   ENDIF (NOT WIN32)

ENDIF (HAS_DCMTK)

#
# Create the vtkdevideIO C++ library
#


ADD_LIBRARY (vtkdevideIO ${IO_LIB_SRCS})

IF (HAS_DCMTK)
   # from dcmimgle onwards is for JPEG support
   TARGET_LINK_LIBRARIES(vtkdevideIO dcmdata ofstd)
   # dcmimgle dcmjpeg ijg12 ijg16 ijg8)
ENDIF (HAS_DCMTK)

IF (WIN32)
   TARGET_LINK_LIBRARIES(vtkdevideIO wsock32 netapi32)
ENDIF (WIN32)

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
  # on Python 2.5 on Windows does not recognise DLL extensions
  # anymore, only PYD !  Previous versions support this.
  IF(WIN32 AND NOT CYGWIN)
      SET_TARGET_PROPERTIES(vtkdevideIOPython 
          PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)

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
