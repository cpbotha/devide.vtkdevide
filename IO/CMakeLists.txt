#
# Source files
#
# Here is where you can add the name of your local IO classes.
#

IF (HAS_DCMTK)
   SET (IO_SRCS ${IO_SRCS} vtkDICOMVolumeReader.cxx)
ENDIF (HAS_DCMTK)

# include vtkIVCONDriver.cxx
#SET (IO_SRCS ${IO_SRCS} vtkIVCONDriver.cxx)

#INCLUDE_DIRECTORIES(${VTKDEVIDE_SOURCE_DIR}/IO/ivcon)
#FILE(GLOB IVCON_SRCS ivcon/*.cpp)


# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdevideIO C++ library
#

IF (WIN32)
   ADD_LIBRARY (vtkdevideIO IO_SRCS)
ELSE (WIN32)
   # under unix, we also want it to include a built-in dicom.dic
   # under windows, the Visual C++ builds already do this
   ADD_LIBRARY (vtkdevideIO IO_SRCS dcdictzz)
ENDIF (WIN32)


#
# Create the vtkdevideIO Tcl library
#

IF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

  VTK_WRAP_TCL (vtkdevideIOTCL IOTCL_SRCS IO_SRCS)
  ADD_LIBRARY (vtkdevideIOTCL SHARED IOTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideIOTCL
                         vtkIOTCL 
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDEVIDE_WRAP_TCL)

#
# Create the vtkdevideIO Python library
#

IF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdevideIOPython IOPython_SRCS IO_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdevideIOPython MODULE IOPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdevideIOPython SHARED IOPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdevideIOPython 
	                   vtkIOPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideIOPython 
	                   vtkIOPython
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdevideIOPython 
	                   vtkIOPython
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDEVIDE_WRAP_PYTHON)

#
# Create the vtkdevideIO Java library
#

IF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdevideIOJava IOJava_SRCS IO_SRCS)
  ADD_LIBRARY (vtkdevideIOJava SHARED IOJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdevideIOJava 
                         vtkIOJava)

ENDIF (VTK_WRAP_JAVA AND VTKDEVIDE_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdevideIO
  vtkIO
)

IF (HAS_DCMTK)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH})
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/ofstd)
   LINK_DIRECTORIES(${DCMTK_LIB_PATH})
   LINK_LIBRARIES(dcmdata ofstd)
ENDIF (HAS_DCMTK)

IF (WIN32)
   LINK_LIBRARIES(wsock32 netapi32)
ENDIF (WIN32)
