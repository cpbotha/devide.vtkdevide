#include "vtkInteractorStyleTrackballActorConstrained.h"
#include <vtkSystemIncludes.h> // this should give us stl thingies
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkstd/algorithm>

vtkCxxRevisionMacro(vtkInteractorStyleTrackballActorConstrained, "$Revision: 1.3 $");
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

void vtkInteractorStyleTrackballActorConstrained::AddActiveProp(vtkProp3D *theProp)
{
   if (!theProp)
   {
      return;
   }

   vtkstd::vector<vtkProp3D *>::iterator ii;
   ii = vtkstd::find(this->ActiveProps.begin(), this->ActiveProps.end(), theProp);

   if (ii == this->ActiveProps.end())
   {
      // not found, so we can add it
      this->ActiveProps.push_back(theProp);
   }
}

void vtkInteractorStyleTrackballActorConstrained::RemoveActiveProp(vtkProp3D *theProp)
{
   if (!theProp)
   {
      return;
   }

   vtkstd::vector<vtkProp3D *>::iterator ii;
   ii = vtkstd::find(this->ActiveProps.begin(), this->ActiveProps.end(), theProp);

   if (ii != this->ActiveProps.end())
   {
      // found it, so we can nuke it
      this->ActiveProps.erase(ii);
   }
}

void vtkInteractorStyleTrackballActorConstrained::RemoveAllActiveProps(void)
{
   this->ActiveProps.clear();
}

void vtkInteractorStyleTrackballActorConstrained::SetActiveProp(vtkProp3D *theProp)
{
   // we have to do this either way...
   this->RemoveAllActiveProps();
   // AddActiveProp will check for NULL
   this->AddActiveProp(theProp);
}

bool vtkInteractorStyleTrackballActorConstrained::IsPropActive(vtkProp3D *theProp)
{
   vtkstd::vector<vtkProp3D *>::iterator ii;
   ii = vtkstd::find(this->ActiveProps.begin(), this->ActiveProps.end(), theProp);

   return (ii != this->ActiveProps.end());
}

void vtkInteractorStyleTrackballActorConstrained::OnLeftButtonDown()
{
   int x = this->Interactor->GetEventPosition()[0];
   int y = this->Interactor->GetEventPosition()[1];

   this->FindPokedRenderer(x, y);
   this->FindPickedActor(x, y);

   if (!this->IsPropActive(this->InteractionProp))
   {
      this->InteractionProp = NULL;
   }

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

   if (!this->IsPropActive(this->InteractionProp))
   {
      this->InteractionProp = NULL;
   }

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

   if (!this->IsPropActive(this->InteractionProp))
   {
      this->InteractionProp = NULL;
   }

   if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
   {
      return;
   }

   this->StartUniformScale();
}

