/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkdscasRenderingWin32Header.h,v $
  Language:  C++
  Date:      $Date: 2003/01/08 14:07:29 $
  Version:   $Revision: 1.1 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkdscasIOWin32Header - manage Windows system differences
// .SECTION Description
// The vtkdscasIOWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkdscasRenderingWin32Header_h
#define __vtkdscasRenderingWin32Header_h

#include <vtkdscasConfigure.h>

#if defined(WIN32) && !defined(VTK_DSCAS_STATIC)
#if defined(vtkdscasRendering_EXPORTS)
#define VTK_DSCAS_RENDERING_EXPORT __declspec( dllexport ) 
#else
#define VTK_DSCAS_RENDERING_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_DSCAS_RENDERING_EXPORT
#endif

#endif
