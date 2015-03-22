#ifndef __vtkdevideCommonWin32Header_h
#define __vtkdevideCommonWin32Header_h

#include <vtkdevideConfigure.h>

#if defined(WIN32) && !defined(VTK_DEVIDE_STATIC)
#if defined(vtkdevideCommon_EXPORTS)
#define VTK_DEVIDE_COMMON_EXPORT __declspec( dllexport ) 
#else
#define VTK_DEVIDE_COMMON_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DEVIDE_COMMON_EXPORT
#endif

#endif
