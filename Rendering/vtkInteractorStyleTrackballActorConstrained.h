// vtkInteractorStyleTrackballActorConstrained copyright (c) 2003 Charl P. Botha <cpbotha@ieee.org>
// $Id: vtkInteractorStyleTrackballActorConstrained.h,v 1.1 2003/06/23 16:20:47 cpbotha Exp $
// interactorstyle that can constrain object manipulation to planes, lines and objects

#ifndef __vtkInteractorStyleTrackballActorConstrained_h
#define __vtkInteractorStyleTrackballActorConstrained_h

#include <vtkInteractorStyleTrackballActor.h>
#include <vtkObjectFactory.h>
#include <vtkstd/vector>
#include "vtkdscasRenderingWin32Header.h"

class vtkProp3D;

class VTK_DSCAS_RENDERING_EXPORT vtkInteractorStyleTrackballActorConstrained : public vtkInteractorStyleTrackballActor
{
public:
   static vtkInteractorStyleTrackballActorConstrained *New();
   vtkTypeRevisionMacro(vtkInteractorStyleTrackballActorConstrained, vtkInteractorStyleTrackballActor);

   virtual void OnLeftButtonDown();
   virtual void OnMiddleButtonDown();
   virtual void OnRightButtonDown();

protected:
   vtkInteractorStyleTrackballActorConstrained();
   ~vtkInteractorStyleTrackballActorConstrained();

//BTX
   vtkstd::vector<vtkProp3D *> ActiveProps;
//ETX

private:
   vtkInteractorStyleTrackballActorConstrained(const vtkInteractorStyleTrackballActorConstrained&);  // Not implemented.
   void operator=(const vtkInteractorStyleTrackballActor&);  // Not implemented.  
};

#endif


