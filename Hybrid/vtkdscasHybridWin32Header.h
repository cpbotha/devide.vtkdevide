#ifndef __vtkdscasHybridWin32Header_h
#define __vtkdscasHybridWin32Header_h

#include <vtkdscasConfigure.h>

#if defined(WIN32) && !defined(VTK_DSCAS_STATIC)
#if defined(vtkdscasHybrid_EXPORTS)
#define VTK_DSCAS_HYBRID_EXPORT __declspec( dllexport ) 
#else
#define VTK_DSCAS_HYBRID_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DSCAS_HYBRID_EXPORT
#endif

#endif
