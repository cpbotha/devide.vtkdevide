#include "vtkInteractorStyleTrackballActorConstrained.h"
#include <vtkSystemIncludes.h> // this should give us stl thingies
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

vtkCxxRevisionMacro(vtkInteractorStyleTrackballActorConstrained, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkInteractorStyleTrackballActorConstrained);

vtkInteractorStyleTrackballActorConstrained::vtkInteractorStyleTrackballActorConstrained() : vtkInteractorStyleTrackballActor()
{
   // make sure there's nothing in our ActiveProps vector
   this->ActiveProps.clear();
}

vtkInteractorStyleTrackballActorConstrained::~vtkInteractorStyleTrackballActorConstrained()
{
   // keep it clean...
   this->ActiveProps.clear();
}

void vtkInteractorStyleTrackballActorConstrained::OnLeftButtonDown()
{
   int x = this->Interactor->GetEventPosition()[0];
   int y = this->Interactor->GetEventPosition()[1];

   this->FindPokedRenderer(x, y);
   this->FindPickedActor(x, y);
   if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
   {
      return;
   }

   if (this->Interactor->GetShiftKey())
   {
      this->StartPan();
   }
   else if (this->Interactor->GetControlKey())
   {
      this->StartSpin();
   }
   else
   {
      this->StartRotate();
   }
}

void vtkInteractorStyleTrackballActorConstrained::OnMiddleButtonDown()
{
   int x = this->Interactor->GetEventPosition()[0];
   int y = this->Interactor->GetEventPosition()[1];

   this->FindPokedRenderer(x, y);
   this->FindPickedActor(x, y);
   if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
   {
      return;
   }

   if (this->Interactor->GetControlKey())
   {
      this->StartDolly();
   }
   else
   {
      this->StartPan();
   }
}

void vtkInteractorStyleTrackballActorConstrained::OnRightButtonDown()
{
   int x = this->Interactor->GetEventPosition()[0];
   int y = this->Interactor->GetEventPosition()[1];

   this->FindPokedRenderer(x, y);
   this->FindPickedActor(x, y);
   if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
   {
      return;
   }

   this->StartUniformScale();
}

