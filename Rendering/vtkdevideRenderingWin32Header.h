/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkdevideRenderingWin32Header.h,v $
  Language:  C++
  Date:      $Date: 2004/01/15 11:00:55 $
  Version:   $Revision: 1.1 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkdevideIOWin32Header - manage Windows system differences
// .SECTION Description
// The vtkdevideIOWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkdevideRenderingWin32Header_h
#define __vtkdevideRenderingWin32Header_h

#include <vtkdevideConfigure.h>

#if defined(WIN32) && !defined(VTK_DEVIDE_STATIC)
#if defined(vtkdevideRendering_EXPORTS)
#define VTK_DEVIDE_RENDERING_EXPORT __declspec( dllexport ) 
#else
#define VTK_DEVIDE_RENDERING_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DEVIDE_RENDERING_EXPORT
#endif

#endif
