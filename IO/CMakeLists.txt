#
# Source files
#
# Here is where you can add the name of your local IO classes.
#

IF (HAS_HDF)
   SET (IO_SRCS vtkHDFVolumeReader.cxx)
ENDIF (HAS_HDF)

IF (HAS_DCMTK)
   SET (IO_SRCS ${IO_SRCS} vtkDICOMVolumeReader.cxx)
ENDIF (HAS_DCMTK)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

#
# Create the vtkdscasIO C++ library
#

IF (WIN32)
   ADD_LIBRARY (vtkdscasIO IO_SRCS)
ELSE (WIN32)
   # under unix, we also want it to include a built-in dicom.dic
   # under windows, the Visual C++ builds already do this
   ADD_LIBRARY (vtkdscasIO IO_SRCS dcdictzz)
ENDIF (WIN32)


#
# Create the vtkdscasIO Tcl library
#

IF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

  VTK_WRAP_TCL (vtkdscasIOTCL IOTCL_SRCS IO_SRCS)
  ADD_LIBRARY (vtkdscasIOTCL SHARED IOTCL_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasIOTCL 
                         vtkIOTCL 
                         ${TCL_LIBRARY})

ENDIF (VTK_WRAP_TCL AND VTKDSCAS_WRAP_TCL)

#
# Create the vtkdscasIO Python library
#

IF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

  VTK_WRAP_PYTHON (vtkdscasIOPython IOPython_SRCS IO_SRCS)

  IF (APPLE)
    ADD_LIBRARY (vtkdscasIOPython MODULE IOPython_SRCS)
  ELSE (APPLE)
    ADD_LIBRARY (vtkdscasIOPython SHARED IOPython_SRCS)
    TARGET_LINK_LIBRARIES (vtkdscasIOPython 
	                   vtkIOPython)
  ENDIF (APPLE)

  IF(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasIOPython 
	                   vtkIOPython
                           debug ${PYTHON_DEBUG_LIBRARY} 
                           optimized ${PYTHON_LIBRARY})
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES (vtkdscasIOPython 
	                   vtkIOPython
                           ${PYTHON_LIBRARY})
  ENDIF(WIN32)

ENDIF (VTK_WRAP_PYTHON AND VTKDSCAS_WRAP_PYTHON)

#
# Create the vtkdscasIO Java library
#

IF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

  VTK_WRAP_JAVA (vtkdscasIOJava IOJava_SRCS IO_SRCS)
  ADD_LIBRARY (vtkdscasIOJava SHARED IOJava_SRCS)
  TARGET_LINK_LIBRARIES (vtkdscasIOJava 
                         vtkIOJava)

ENDIF (VTK_WRAP_JAVA AND VTKDSCAS_WRAP_JAVA)

#
# Link to VTK and local
#

LINK_LIBRARIES (
  vtkdscasIO
  vtkIO
)

IF (HAS_HDF)
   INCLUDE_DIRECTORIES(${HDF_INCLUDE_PATH})
   LINK_DIRECTORIES(${HDF_LIB_PATH})
   LINK_LIBRARIES (mfhdf df)
ENDIF (HAS_HDF)

IF (HAS_DCMTK)
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH})
   INCLUDE_DIRECTORIES(${DCMTK_INCLUDE_PATH}/ofstd)
   LINK_DIRECTORIES(${DCMTK_LIB_PATH})
   LINK_LIBRARIES(dcmdata ofstd)
ENDIF (HAS_DCMTK)

IF (WIN32)
   LINK_LIBRARIES(wsock32 netapi32)
ENDIF (WIN32)
