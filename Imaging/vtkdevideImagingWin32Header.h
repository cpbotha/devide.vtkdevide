#ifndef __vtkdevideImagingWin32Header_h
#define __vtkdevideImagingWin32Header_h

#include <vtkdevideConfigure.h>

#if defined(WIN32) && !defined(VTK_DEVIDE_STATIC)
#if defined(vtkdevideImaging_EXPORTS)
#define VTK_DEVIDE_IMAGING_EXPORT __declspec( dllexport ) 
#else
#define VTK_DEVIDE_IMAGING_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DEVIDE_IMAGING_EXPORT
#endif

#endif
