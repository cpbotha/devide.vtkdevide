// vtkdscasIOWin32Header - manage Windows system differences
// $Id: vtkdscasIOWin32Header.h,v 1.1 2003/01/08 14:07:29 cpbotha Exp $
// The vtkdscasIOWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkdscasIOWin32Header_h
#define __vtkdscasIOWin32Header_h

#include <vtkdscasConfigure.h>

#if defined(WIN32) && !defined(VTK_DSCAS_STATIC)
#if defined(vtkdscasIO_EXPORTS)
#define VTK_DSCAS_IO_EXPORT __declspec( dllexport ) 
#else
#define VTK_DSCAS_IO_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DSCAS_IO_EXPORT
#endif

#endif
