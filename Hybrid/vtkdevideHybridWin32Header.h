#ifndef __vtkdevideHybridWin32Header_h
#define __vtkdevideHybridWin32Header_h

#include <vtkdevideConfigure.h>

#if defined(WIN32) && !defined(VTK_DEVIDE_STATIC)
#if defined(vtkdevideHybrid_EXPORTS)
#define VTK_DEVIDE_HYBRID_EXPORT __declspec( dllexport ) 
#else
#define VTK_DEVIDE_HYBRID_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DEVIDE_HYBRID_EXPORT
#endif

#endif
