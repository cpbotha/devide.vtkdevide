#
# Set the project name.
#

PROJECT (VTKDSCAS)

#
# Assume everything is set up correctly for build.
#
# If some requirements are not met, set it to 0.
# In that case, the build won't be done.
#

SET (VTKDSCAS_CAN_BUILD 1)

#
# Load CMake commands that you probably should not modify.
#

INCLUDE (${VTKDSCAS_SOURCE_DIR}/CMakeOptions.cmake)

#
# If everything is OK, then go into the sub directories and build.
#

IF (VTKDSCAS_CAN_BUILD)

  #
  # Here is where you can list the sub-directories holding your local 
  # classes. Sorting classes by 'package' type like VTK does (Common, 
  # Rendering, Filtering, Imaging, IO, etc.) is a good thing and prevents 
  # numerous dependencies problems.
  #

  SUBDIRS (
    Hybrid
    IO
    Rendering
  )

  #
  # You can put your include path(s) here
  #

  INCLUDE_DIRECTORIES(${VTKDSCAS_SOURCE_DIR}/IO)
  INCLUDE_DIRECTORIES(${VTKDSCAS_SOURCE_DIR}/Testing)  
  
  #
  # are we going to build the examples?
  #
  
  OPTION(BUILD_EXAMPLES "Build examples." ON)
  IF (BUILD_EXAMPLES)
    SUBDIRS(Examples)
  ENDIF (BUILD_EXAMPLES)
  

  #
  # Common configuration settings
  #
  # Do not worry about this one.
  #

  CONFIGURE_FILE(
    ${VTKDSCAS_SOURCE_DIR}/vtkdscasConfigure.h.in
    ${VTKDSCAS_BINARY_DIR}/vtkdscasConfigure.h
  )

  INCLUDE_DIRECTORIES(${VTKDSCAS_BINARY_DIR})

ENDIF (VTKDSCAS_CAN_BUILD)
