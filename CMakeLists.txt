#
# Set the project name.
#

PROJECT (VTKCPBOTHA)

#
# Assume everything is set up correctly for build.
#
# If some requirements are not met, set it to 0.
# In that case, the build won't be done.
#

SET (VTKCPBOTHA_CAN_BUILD 1)

#
# Load CMake commands that you probably should not modify.
#

INCLUDE (${VTKCPBOTHA_SOURCE_DIR}/CMakeOptions.cmake)

#
# If everything is OK, then go into the sub directories and build.
#

IF (VTKCPBOTHA_CAN_BUILD)

  #
  # Here is where you can list the sub-directories holding your local 
  # classes. Sorting classes by 'package' type like VTK does (Common, 
  # Rendering, Filtering, Imaging, IO, etc.) is a good thing and prevents 
  # numerous dependencies problems.
  #

  SUBDIRS (
    IO
    Rendering
#    Testing
  )

  #
  # You can put your include path(s) here
  #

  INCLUDE_DIRECTORIES(${VTKCPBOTHA_SOURCE_DIR}/IO)
  INCLUDE_DIRECTORIES(${VTKCPBOTHA_SOURCE_DIR}/Testing)  
  
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
    ${VTKCPBOTHA_SOURCE_DIR}/vtkcpbothaConfigure.h.in
    ${VTKCPBOTHA_BINARY_DIR}/vtkcpbothaConfigure.h
  )

  INCLUDE_DIRECTORIES(${VTKCPBOTHA_BINARY_DIR})

ENDIF (VTKCPBOTHA_CAN_BUILD)
