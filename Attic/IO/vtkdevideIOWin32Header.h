// vtkdevideIOWin32Header - manage Windows system differences
// $Id: vtkdevideIOWin32Header.h,v 1.1 2004/01/15 11:00:55 cpbotha Exp $
// The vtkdevideIOWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkdevideIOWin32Header_h
#define __vtkdevideIOWin32Header_h

#include <vtkdevideConfigure.h>

#if defined(WIN32) && !defined(VTK_DEVIDE_STATIC)
#if defined(vtkdevideIO_EXPORTS)
#define VTK_DEVIDE_IO_EXPORT __declspec( dllexport ) 
#else
#define VTK_DEVIDE_IO_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DEVIDE_IO_EXPORT
#endif

#endif
