// vtkdevideTestingWin32Header - manage Windows system differences
// $Id: vtkdscasTestingWin32Header.h,v 1.2 2004/01/15 11:00:56 cpbotha Exp $
// The vtkdevideTestingWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkdevideTestingWin32Header_h
#define __vtkdevideTestingWin32Header_h

#include <vtkdevideConfigure.h>

#if defined(WIN32) && !defined(VTK_DEVIDE_STATIC)
#if defined(vtkdevideTesting_EXPORTS)
#define VTK_DEVIDE_TESTING_EXPORT __declspec( dllexport ) 
#else
#define VTK_DEVIDE_TESTING_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DEVIDE_TESTING_EXPORT
#endif

#endif
