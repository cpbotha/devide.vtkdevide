#ifndef __vtkdscasCommonWin32Header_h
#define __vtkdscasCommonWin32Header_h

#include <vtkdscasConfigure.h>

#if defined(WIN32) && !defined(VTK_DSCAS_STATIC)
#if defined(vtkdscasCommon_EXPORTS)
#define VTK_DSCAS_COMMON_EXPORT __declspec( dllexport ) 
#else
#define VTK_DSCAS_COMMON_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DSCAS_COMMON_EXPORT
#endif

#endif
