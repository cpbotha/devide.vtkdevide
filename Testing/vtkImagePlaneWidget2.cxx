/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImagePlaneWidget2.cxx,v $
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
#include "vtkImagePlaneWidget2.h"

#include "vtkActor.h"
#include "vtkAssemblyNode.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCellPicker.h"
#include "vtkDataSetMapper.h"
#include "vtkImageData.h"
#include "vtkImageMapToColors.h"
#include "vtkImageReslice.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkPlaneSource.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkTexture.h"
#include "vtkTextureMapToPlane.h"
#include "vtkTransform.h"

#include "vtkCellLocator.h"
#include "vtkGenericCell.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVoxel.h"
#include "vtkIdList.h"

vtkCxxRevisionMacro(vtkImagePlaneWidget2, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkImagePlaneWidget2);

vtkCxxSetObjectMacro(vtkImagePlaneWidget2, PlaneProperty, vtkProperty);
vtkCxxSetObjectMacro(vtkImagePlaneWidget2, SelectedPlaneProperty, vtkProperty);
vtkCxxSetObjectMacro(vtkImagePlaneWidget2, CursorProperty, vtkProperty);
vtkCxxSetObjectMacro(vtkImagePlaneWidget2, MarginProperty, vtkProperty);
vtkCxxSetObjectMacro(vtkImagePlaneWidget2, TexturePlaneProperty, vtkProperty);

#if defined(WIN32)
#define for if(0);else for
#endif

vtkImagePlaneWidget2::vtkImagePlaneWidget2() : vtkPolyDataSourceWidget()
{
  this->State = vtkImagePlaneWidget2::Start;
  this->EventCallbackCommand->SetCallback(vtkImagePlaneWidget2::ProcessEvents);

  this->Interaction              = 1;
  this->PlaneOrientation         = 0;
  this->RestrictPlaneToVolume    = 1;
  this->OriginalWindow           = 1.0;
  this->OriginalLevel            = 0.5;
  this->CurrentWindow            = 1.0;
  this->CurrentLevel             = 0.5;
  this->TextureInterpolate       = 1;
  this->ResliceInterpolate       = VTK_LINEAR_RESLICE;
  this->UserPickerEnabled        = 0;
  this->UserLookupTableEnabled   = 0;
  this->UserControlledLookupTable= 0;
  this->DisplayText              = 0;
  this->CurrentCursorPosition[0] = 0;
  this->CurrentCursorPosition[1] = 0;
  this->CurrentCursorPosition[2] = 0;
  this->CurrentImageValue        = VTK_FLOAT_MAX;
  this->MarginSelectMode         = 8;

  // attempt at cursor snapping
  this->Locator = vtkCellLocator::New();
  this->Voxel = vtkGenericCell::New();
  this->Voxel->SetCellTypeToVoxel();


  // Represent the plane's outline
  //
  this->PlaneSource = vtkPlaneSource::New();
  this->PlaneSource->SetXResolution(1);
  this->PlaneSource->SetYResolution(1);
  this->PlaneOutlinePolyData = vtkPolyData::New();
  this->PlaneOutlineMapper   = vtkPolyDataMapper::New();
  this->PlaneOutlineActor    = vtkActor::New();
  this->PlaneOutlineActor->PickableOff();
  this->PlaneOutlineActor->GetProperty()->SetRepresentationToWireframe();

  // Represent the resliced image plane
  //
  this->LookupTable        = vtkLookupTable::New();
  this->ColorMap           = vtkImageMapToColors::New();
  this->Reslice            = vtkImageReslice::New();
  this->ResliceAxes        = vtkMatrix4x4::New();
  this->Texture            = vtkTexture::New();
  this->TexturePlaneCoords = vtkTextureMapToPlane::New();
  this->TexturePlaneMapper = vtkDataSetMapper::New();
  this->TexturePlaneActor  = vtkActor::New();
  this->Transform          = vtkTransform::New();
  this->ImageData          = 0;

  // Represent the cross hair cursor
  //
  this->CursorPolyData = vtkPolyData::New();
  this->CursorMapper   = vtkPolyDataMapper::New();
  this->CursorActor    = vtkActor::New();
  this->CursorPoints   = vtkPoints::New(VTK_DOUBLE);

  // Represent the oblique positioning margins
  //
  this->MarginPolyData = vtkPolyData::New();
  this->MarginMapper   = vtkPolyDataMapper::New();
  this->MarginActor    = vtkActor::New();
  this->MarginPoints   = vtkPoints::New(VTK_DOUBLE);

  // Represent the text: annotation for cursor position and W/L
  //
  this->TextActor = vtkTextActor::New();

  this->GeneratePlaneOutline();

  // Define some default point coordinates
  //
  float bounds[6];
  bounds[0] = -0.5;
  bounds[1] =  0.5;
  bounds[2] = -0.5;
  bounds[3] =  0.5;
  bounds[4] = -0.5;
  bounds[5] =  0.5;

  // Initial creation of the widget, serves to initialize it
  //
  this->PlaceWidget(bounds);


  this->VoxelGrid = vtkUnstructuredGrid::New();
  this->VoxelMapper = vtkDataSetMapper::New();
  this->VoxelActor = vtkActor::New();


  this->GenerateTexturePlane();
  this->GenerateCursor();
  this->GenerateMargins();
  this->GenerateText();
  this->GenerateVoxel();

  // Manage the picking stuff
  //
  this->PlanePicker = vtkCellPicker::New();
  this->PlanePicker->SetTolerance(0.005); //need some fluff
  this->PlanePicker->AddPickList(this->TexturePlaneActor);
  this->PlanePicker->PickFromListOn();

  // Set up the initial properties
  //
  this->PlaneProperty         = 0;
  this->SelectedPlaneProperty = 0;
  this->TexturePlaneProperty  = 0;
  this->CursorProperty        = 0;
  this->MarginProperty        = 0;
  this->CreateDefaultProperties();
}

vtkImagePlaneWidget2::~vtkImagePlaneWidget2()
{
  this->PlaneOutlineActor->Delete();
  this->PlaneOutlineMapper->Delete();
  this->PlaneOutlinePolyData->Delete();
  this->PlaneSource->Delete();

  this->Locator->Delete();
  this->Voxel->Delete();
  this->VoxelGrid->Delete();
  this->VoxelMapper->Delete();
  this->VoxelActor->Delete();


  if ( !this->UserPickerEnabled )
    {
    this->PlanePicker->Delete();
    }
  else
    {
    this->PlanePicker = 0;
    }

  if ( this->PlaneProperty )
    {
    this->PlaneProperty->Delete();
    }

  if ( this->SelectedPlaneProperty )
    {
    this->SelectedPlaneProperty->Delete();
    }

  if ( this->CursorProperty )
    {
    this->CursorProperty->Delete();
    }

  if ( this->MarginProperty )
    {
    this->MarginProperty->Delete();
    }

  this->ResliceAxes->Delete();
  this->Transform->Delete();
  this->Reslice->Delete();

  if ( !this->UserLookupTableEnabled )
    {
    this->LookupTable->Delete();
    }
  else
    {
    this->LookupTable = 0;
    }

  this->TexturePlaneCoords->Delete();
  this->TexturePlaneMapper->Delete();
  this->TexturePlaneActor->Delete();
  this->ColorMap->Delete();
  this->Texture->Delete();

  if ( this->TexturePlaneProperty )
    {
    this->TexturePlaneProperty->Delete();
    }

  if ( this->ImageData )
    {
    this->ImageData = 0;
    }

  this->CursorActor->Delete();
  this->CursorMapper->Delete();
  this->CursorPoints->Delete();
  this->CursorPolyData->Delete();

  this->MarginActor->Delete();
  this->MarginMapper->Delete();
  this->MarginPoints->Delete();
  this->MarginPolyData->Delete();

  this->TextActor->Delete();
}

void vtkImagePlaneWidget2::SetEnabled(int enabling)
{

  if ( ! this->Interactor )
    {
    vtkErrorMacro(<<"The interactor must be set prior to enabling/disabling widget");
    return;
    }

  if ( enabling ) //----------------------------------------------------------
    {
    vtkDebugMacro(<<"Enabling plane widget");

    if ( this->Enabled ) //already enabled, just return
      {
      return;
      }
    
    this->CurrentRenderer = 
      this->Interactor->FindPokedRenderer(this->Interactor->GetLastEventPosition()[0],
                                          this->Interactor->GetLastEventPosition()[1]);
    if ( this->CurrentRenderer == 0 )
      {
      return;
      }

    this->Enabled = 1;

    // listen for the following events
    vtkRenderWindowInteractor *i = this->Interactor;
    i->AddObserver(vtkCommand::MouseMoveEvent, this->EventCallbackCommand,
                   this->Priority);
    i->AddObserver(vtkCommand::LeftButtonPressEvent,
                   this->EventCallbackCommand, this->Priority);
    i->AddObserver(vtkCommand::LeftButtonReleaseEvent,
                   this->EventCallbackCommand, this->Priority);
    i->AddObserver(vtkCommand::MiddleButtonPressEvent,
                   this->EventCallbackCommand, this->Priority);
    i->AddObserver(vtkCommand::MiddleButtonReleaseEvent,
                   this->EventCallbackCommand, this->Priority);
    i->AddObserver(vtkCommand::RightButtonPressEvent,
                   this->EventCallbackCommand, this->Priority);
    i->AddObserver(vtkCommand::RightButtonReleaseEvent,
                   this->EventCallbackCommand, this->Priority);

    // Add the plane
    this->CurrentRenderer->AddProp(this->PlaneOutlineActor);
    this->PlaneOutlineActor->SetProperty(this->PlaneProperty);

    //add the TexturePlaneActor
    this->CurrentRenderer->AddProp(this->TexturePlaneActor);
    this->TexturePlaneActor->SetProperty(this->TexturePlaneProperty);

    // Add the cross-hair cursor
    this->CurrentRenderer->AddProp(this->CursorActor);
    this->CursorActor->SetProperty(this->CursorProperty);

    // Add the margins
    this->CurrentRenderer->AddProp(this->MarginActor);
    this->MarginActor->SetProperty(this->MarginProperty);

    // Add the image data annotation
    this->CurrentRenderer->AddProp(this->TextActor);

    //add the VoxelActor
    this->CurrentRenderer->AddProp(this->VoxelActor);


    if ( this->PlanePicker )
      {
      this->TexturePlaneActor->PickableOn();
      }

    this->InvokeEvent(vtkCommand::EnableEvent,0);

    }

  else //disabling----------------------------------------------------------
    {
    vtkDebugMacro(<<"Disabling plane widget");

    if ( ! this->Enabled ) //already disabled, just return
      {
      return;
      }

    this->Enabled = 0;

    // don't listen for events any more
    this->Interactor->RemoveObserver(this->EventCallbackCommand);

    // turn off the plane
    this->CurrentRenderer->RemoveProp(this->PlaneOutlineActor);

    //turn off the texture plane
    this->CurrentRenderer->RemoveProp(this->TexturePlaneActor);

    //turn off the cursor
    this->CurrentRenderer->RemoveProp(this->CursorActor);

    //turn off the margins
    this->CurrentRenderer->RemoveProp(this->MarginActor);    

    //turn off the image data annotation
    this->CurrentRenderer->RemoveProp(this->TextActor);

    //turn off the VoxelActor
    this->CurrentRenderer->RemoveProp(this->VoxelActor);


    if ( this->PlanePicker )
      {
      this->TexturePlaneActor->PickableOff();
      }

    this->InvokeEvent(vtkCommand::DisableEvent,0);
    }

  this->Interactor->Render();
}

void vtkImagePlaneWidget2::ProcessEvents(vtkObject* vtkNotUsed(object),
                                        unsigned long event,
                                        void* clientdata,
                                        void* vtkNotUsed(calldata))
{
  vtkImagePlaneWidget2* self =
    reinterpret_cast<vtkImagePlaneWidget2 *>( clientdata );

  //okay, let's do the right thing
  switch ( event )
    {
    case vtkCommand::LeftButtonPressEvent:
      self->OnLeftButtonDown();
      break;
    case vtkCommand::LeftButtonReleaseEvent:
      self->OnLeftButtonUp();
      break;
    case vtkCommand::MiddleButtonPressEvent:
      self->OnMiddleButtonDown();
      break;
    case vtkCommand::MiddleButtonReleaseEvent:
      self->OnMiddleButtonUp();
      break;
    case vtkCommand::RightButtonPressEvent:
      self->OnRightButtonDown();
      break;
    case vtkCommand::RightButtonReleaseEvent:
      self->OnRightButtonUp();
      break;
    case vtkCommand::MouseMoveEvent:
      self->OnMouseMove();
      break;
    }
}

void vtkImagePlaneWidget2::SetInteraction(int interact)
{
  if (this->Interactor && this->Enabled)
    {
    if (this->Interaction == interact)
      {
      return;
      }
    if (interact == 0)
      {
      this->Interactor->RemoveObserver(this->EventCallbackCommand);
      }
    else
      {
      vtkRenderWindowInteractor *i = this->Interactor;
      i->AddObserver(vtkCommand::MouseMoveEvent, this->EventCallbackCommand,
                     this->Priority);
      i->AddObserver(vtkCommand::LeftButtonPressEvent,
                     this->EventCallbackCommand, this->Priority);
      i->AddObserver(vtkCommand::LeftButtonReleaseEvent,
                     this->EventCallbackCommand, this->Priority);
      i->AddObserver(vtkCommand::MiddleButtonPressEvent,
                     this->EventCallbackCommand, this->Priority);
      i->AddObserver(vtkCommand::MiddleButtonReleaseEvent,
                     this->EventCallbackCommand, this->Priority);
      i->AddObserver(vtkCommand::RightButtonPressEvent,
                     this->EventCallbackCommand, this->Priority);
      i->AddObserver(vtkCommand::RightButtonReleaseEvent,
                     this->EventCallbackCommand, this->Priority);
      }
    this->Interaction = interact;
    }
  else
    {
    vtkGenericWarningMacro(<<"set interactor and Enabled before changing interaction...");
    }
}

void vtkImagePlaneWidget2::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  if ( this->PlaneProperty )
    {
    os << indent << "Plane Property: " << this->PlaneProperty << "\n";
    }
  else
    {
    os << indent << "Plane Property: (none)\n";
    }

  if ( this->SelectedPlaneProperty )
    {
    os << indent << "Selected Plane Property: "
       << this->SelectedPlaneProperty << "\n";
    }
  else
    {
    os << indent << "Selected Plane Property: (none)\n";
    }

  if ( this->LookupTable )
    {
    os << indent << "LookupTable: "
       << this->LookupTable << "\n";
    }
  else
    {
    os << indent << "LookupTable: (none)\n";
    }

  if ( this->CursorProperty )
    {
    os << indent << "Cursor Property: "
       << this->CursorProperty << "\n";
    }
  else
    {
    os << indent << "Cursor Property: (none)\n";
    }

  if ( this->MarginProperty )
    {
    os << indent << "Margin Property: "
       << this->MarginProperty << "\n";
    }
  else
    {
    os << indent << "Margin Property: (none)\n";
    }

  if ( this->TexturePlaneProperty )
    {
    os << indent << "TexturePlane Property: "
       << this->TexturePlaneProperty << "\n";
    }
  else
    {
    os << indent << "TexturePlane Property: (none)\n";
    }

  float *o = this->PlaneSource->GetOrigin();
  float *pt1 = this->PlaneSource->GetPoint1();
  float *pt2 = this->PlaneSource->GetPoint2();

  os << indent << "Origin: (" << o[0] << ", "
     << o[1] << ", "
     << o[2] << ")\n";
  os << indent << "Point 1: (" << pt1[0] << ", "
     << pt1[1] << ", "
     << pt1[2] << ")\n";
  os << indent << "Point 2: (" << pt2[0] << ", "
     << pt2[1] << ", "
     << pt2[2] << ")\n";

  os << indent << "Plane Orientation: " << this->PlaneOrientation << "\n";
  os << indent << "Reslice Interpolate: " << this->ResliceInterpolate << "\n";
  os << indent << "Texture Interpolate: " 
     << (this->TextureInterpolate ? "On\n" : "Off\n") ;
  os << indent << "Restrict Plane To Volume: " 
     << (this->RestrictPlaneToVolume ? "On\n" : "Off\n") ;
  os << indent << "Display Text: "
     << (this->DisplayText ? "On\n" : "Off\n") ;
  os << indent << "Interaction: "
     << (this->Interaction ? "On\n" : "Off\n") ;
  os << indent << "User Controlled Lookup Table: " 
     << (this->UserControlledLookupTable ? "On\n" : "Off\n") ;
}

void vtkImagePlaneWidget2::BuildRepresentation()
{
  float *o = this->PlaneSource->GetOrigin();
  float *pt1 = this->PlaneSource->GetPoint1();
  float *pt2 = this->PlaneSource->GetPoint2();

  float x[3];
  x[0] = o[0] + (pt1[0]-o[0]) + (pt2[0]-o[0]);
  x[1] = o[1] + (pt1[1]-o[1]) + (pt2[1]-o[1]);
  x[2] = o[2] + (pt1[2]-o[2]) + (pt2[2]-o[2]);

  vtkPoints* points = this->PlaneOutlinePolyData->GetPoints();
  points->SetPoint(0,o);
  points->SetPoint(1,pt1);
  points->SetPoint(2,x);
  points->SetPoint(3,pt2);
  this->PlaneOutlineMapper->Modified();

  this->PlaneSource->GetNormal(this->Normal);
  vtkMath::Normalize(this->Normal);
}

void vtkImagePlaneWidget2::HighlightPlane(int highlight)
{
  if ( highlight )
    {
    this->PlaneOutlineActor->SetProperty(this->SelectedPlaneProperty);
    this->PlanePicker->GetPickPosition(this->LastPickPosition);
    }
  else
    {
    this->PlaneOutlineActor->SetProperty(this->PlaneProperty);
    }
}

void vtkImagePlaneWidget2::OnLeftButtonDown()
{
  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  vtkRenderer *ren = this->Interactor->FindPokedRenderer(X,Y);
  if ( ren != this->CurrentRenderer )
    {
    this->State = vtkImagePlaneWidget2::Outside;
    return;
    }

  // Okay, we can process this. If anything is picked, then we
  // can start pushing the plane.
  vtkAssemblyPath *path;
  this->PlanePicker->Pick(X,Y,0.0,this->CurrentRenderer);
  path = this->PlanePicker->GetPath();

  int found = 0;
  int i;
  if ( path != 0 )
    {
// Deal with the possibility that we may be using a shared picker
    path->InitTraversal();
    vtkAssemblyNode *node;
    for ( i = 0; i < path->GetNumberOfItems() && !found ; i++ )
      {
      node = path->GetNextNode();
      if ( node->GetProp() == vtkProp::SafeDownCast(this->TexturePlaneActor) )
        {
        found = 1;
        }
      }
    }

  if( ! found || path == 0 )
    {
    this->State = vtkImagePlaneWidget2::Outside;
    this->HighlightPlane(0);
    this->ActivateCursor(0);
    this->ActivateText(0);
    return;
    }
  else
    {
    this->State = vtkImagePlaneWidget2::Cursoring;
    this->HighlightPlane(1);
    this->ActivateCursor(1);
    this->ActivateText(1);
    this->UpdateCursor(X,Y);
    this->ManageTextDisplay();
    }

  this->EventCallbackCommand->SetAbortFlag(1);
  this->StartInteraction();
  this->InvokeEvent(vtkCommand::StartInteractionEvent,0);
  this->Interactor->Render();
}

void vtkImagePlaneWidget2::OnLeftButtonUp()
{
  if ( this->State == vtkImagePlaneWidget2::Outside ||
       this->State == vtkImagePlaneWidget2::Start )
    {
    return;
    }

  this->State = vtkImagePlaneWidget2::Start;
  this->HighlightPlane(0);
  this->ActivateCursor(0);
  this->ActivateText(0);

  this->EventCallbackCommand->SetAbortFlag(1);
  this->EndInteraction();
  this->InvokeEvent(vtkCommand::EndInteractionEvent,0);
  this->Interactor->Render();
}

void vtkImagePlaneWidget2::OnMiddleButtonDown()
{
  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  vtkRenderer *ren = this->Interactor->FindPokedRenderer(X,Y);
  if ( ren != this->CurrentRenderer )
    {
    this->State = vtkImagePlaneWidget2::Outside;
    return;
    }

  // Okay, we can process this. If anything is picked, then we
  // can start pushing or check for adjusted states.
  vtkAssemblyPath *path;
  this->PlanePicker->Pick(X,Y,0.0,this->CurrentRenderer);
  path = this->PlanePicker->GetPath();

  int found = 0;
  int i;
  if ( path != 0 )
    {
// Deal with the possibility that we may be using a shared picker
    path->InitTraversal();
    vtkAssemblyNode *node;
    for(i = 0; i< path->GetNumberOfItems() && !found ;i++)
      {
      node = path->GetNextNode();
      if(node->GetProp() == vtkProp::SafeDownCast(this->TexturePlaneActor) )
        {
        found = 1;
        }
      }
    }

  if ( !found || path == 0 )
    {
    this->State = vtkImagePlaneWidget2::Outside;
    this->HighlightPlane(0);
    this->ActivateMargins(0);
    return;
    }
  else
    {
    this->State = vtkImagePlaneWidget2::Pushing;
    this->HighlightPlane(1);
    this->ActivateMargins(1);
    this->AdjustState();
    this->UpdateMargins();
    }

  this->EventCallbackCommand->SetAbortFlag(1);
  this->StartInteraction();
  this->InvokeEvent(vtkCommand::StartInteractionEvent,0);
  this->Interactor->Render();
}

void vtkImagePlaneWidget2::OnMiddleButtonUp()
{
  if ( this->State == vtkImagePlaneWidget2::Outside ||
       this->State == vtkImagePlaneWidget2::Start )
    {
    return;
    }

  this->State = vtkImagePlaneWidget2::Start;
  this->HighlightPlane(0);
  this->ActivateMargins(0);

  this->EventCallbackCommand->SetAbortFlag(1);
  this->EndInteraction();
  this->InvokeEvent(vtkCommand::EndInteractionEvent,0);
  this->Interactor->Render();
}

void vtkImagePlaneWidget2::OnRightButtonDown()
{
  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  vtkRenderer *ren = this->Interactor->FindPokedRenderer(X,Y);
  if ( ren != this->CurrentRenderer )
    {
    this->State = vtkImagePlaneWidget2::Outside;
    return;
    }

  // Okay, we can process this. If anything is picked, then we
  // can start window-levelling.
  vtkAssemblyPath *path;
  this->PlanePicker->Pick(X,Y,0.0,this->CurrentRenderer);
  path = this->PlanePicker->GetPath();

  int found = 0;
  int i;
  if ( path != 0 )
    {
// Deal with the possibility that we may be using a shared picker
    path->InitTraversal();
    vtkAssemblyNode *node;
    for ( i = 0; i < path->GetNumberOfItems() && !found ; i++ )
      {
      node = path->GetNextNode();
      if ( node->GetProp() == vtkProp::SafeDownCast(this->TexturePlaneActor) )
        {
        found = 1;
        }
      }
    }

  if( ! found || path == 0 )
    {
    this->State = vtkImagePlaneWidget2::Outside;
    this->HighlightPlane(0);
    this->ActivateText(0);
    return;
    }
  else
    {
    this->State = vtkImagePlaneWidget2::WindowLevelling;
    this->HighlightPlane(1);
    this->ActivateText(1);
    this->WindowLevel(X,Y);
    this->ManageTextDisplay();
    }

  this->EventCallbackCommand->SetAbortFlag(1);
  this->StartInteraction();
  this->InvokeEvent(vtkCommand::StartInteractionEvent,0);
  this->Interactor->Render();
}

void vtkImagePlaneWidget2::OnRightButtonUp()
{
  if ( this->State == vtkImagePlaneWidget2::Outside ||
       this->State == vtkImagePlaneWidget2::Start )
    {
    return;
    }

  this->State = vtkImagePlaneWidget2::Start;
  this->HighlightPlane(0);
  this->ActivateText(0);

  this->EventCallbackCommand->SetAbortFlag(1);
  this->EndInteraction();
  this->InvokeEvent(vtkCommand::EndInteractionEvent,0);
  this->Interactor->Render();
}

void vtkImagePlaneWidget2::OnMouseMove()
{
  // See whether we're active
  //
  if ( this->State == vtkImagePlaneWidget2::Outside ||
       this->State == vtkImagePlaneWidget2::Start )
    {
    return;
    }

  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  // Do different things depending on state
  // Calculations everybody does
  //
  double focalPoint[4], pickPoint[4], prevPickPoint[4];
  double z, vpn[3];

  vtkRenderer *renderer = this->Interactor->FindPokedRenderer(X,Y);
  vtkCamera *camera = renderer->GetActiveCamera();
  if ( ! camera )
    {
    return;
    }

  // Compute the two points defining the motion vector
  //
  this->ComputeWorldToDisplay(this->LastPickPosition[0], this->LastPickPosition[1],
                              this->LastPickPosition[2], focalPoint);
  z = focalPoint[2];

  this->ComputeDisplayToWorld(double(this->Interactor->GetLastEventPosition()[0]),
                              double(this->Interactor->GetLastEventPosition()[1]),
                              z, prevPickPoint);

  this->ComputeDisplayToWorld(double(X), double(Y), z, pickPoint);

  if ( this->State == vtkImagePlaneWidget2::WindowLevelling )
    {
    this->WindowLevel(X,Y);
    this->ManageTextDisplay();
    }
  else if ( this->State == vtkImagePlaneWidget2::Pushing )
    {
    this->Push(prevPickPoint, pickPoint);
    this->UpdateNormal();
    this->UpdateOrigin();
    this->UpdateMargins();
    }
  else if ( this->State == vtkImagePlaneWidget2::Spinning )
    {
    this->Spin(prevPickPoint, pickPoint);
    this->UpdateNormal();
    this->UpdateOrigin();
    this->UpdateMargins();
    }
  else if ( this->State == vtkImagePlaneWidget2::Rotating )
    {
    camera->GetViewPlaneNormal(vpn);
    this->Rotate(prevPickPoint, pickPoint, vpn);
    this->UpdateNormal();
    this->UpdateOrigin();
    this->UpdateMargins();
    }
  else if ( this->State == vtkImagePlaneWidget2::Scaling )
    {
    this->Scale(prevPickPoint, pickPoint, X, Y);
    this->UpdateNormal();
    this->UpdateOrigin();
    this->UpdateMargins();
    }
  else if ( this->State == vtkImagePlaneWidget2::Moving )
    {
    this->Translate(prevPickPoint, pickPoint);
    this->UpdateNormal();
    this->UpdateOrigin();
    this->UpdateMargins();
    }
  else if ( this->State == vtkImagePlaneWidget2::Cursoring )
    {
    this->UpdateCursor(X,Y);
    this->ManageTextDisplay();
    }

  // Interact, if desired
  //
  this->EventCallbackCommand->SetAbortFlag(1);
  this->InvokeEvent(vtkCommand::InteractionEvent,0);

  this->Interactor->Render();
}

void vtkImagePlaneWidget2::WindowLevel(int X, int Y)
{
  if ( ! this->LookupTable )
    {
    return;
    }

  float range[2];
  this->LookupTable->GetTableRange(range);
  float window = range[1] - range[0];
  float level = 0.5*(range[0] + range[1]);

  float owin = this->OriginalWindow;

  level = level + (X - this->Interactor->GetLastEventPosition()[0])*owin/500.0;
  window = window + (this->Interactor->GetLastEventPosition()[1] - Y)*owin/250.0;

  if ( window == 0.0 )
    {
    window = 0.001;
    }

  float rmin = level - window*0.5;
  float rmax = level + window*0.5;

  if( rmin < rmax )
    {
    this->CurrentWindow = window;
    this->CurrentLevel = level;
    if (!this->UserControlledLookupTable)
      {
      this->LookupTable->SetTableRange(rmin,rmax);
      }
    }
}

void vtkImagePlaneWidget2::GetWindowLevel(float wl[2])
{
  float range[2];
  this->LookupTable->GetTableRange(range);
  wl[0] = range[1] - range[0];
  wl[1] = 0.5*(range[0]+range[1]);
}

int vtkImagePlaneWidget2::GetCursorData(float xyzv[4])
{
  if ( this->State != vtkImagePlaneWidget2::Cursoring  || \
    this->CurrentImageValue == VTK_FLOAT_MAX )
    {
    return 0;
    }

  xyzv[0] = this->CurrentCursorPosition[0];
  xyzv[1] = this->CurrentCursorPosition[1];
  xyzv[2] = this->CurrentCursorPosition[2];
  xyzv[3] = this->ImageData->GetScalarComponentAsFloat( \
                   this->CurrentCursorPosition[0],
                   this->CurrentCursorPosition[1],
                   this->CurrentCursorPosition[2],0);
  return 1;
}

void vtkImagePlaneWidget2::ManageTextDisplay()
{
  if ( !this->DisplayText )
    {
    return;
    }

  if ( this->State == vtkImagePlaneWidget2::WindowLevelling )
    {
    sprintf(this->TextBuff,"Window, Level: ( %g, %g )",
            this->CurrentWindow, this->CurrentLevel );
    }
  else if ( this->State == vtkImagePlaneWidget2::Cursoring )
    {
    if( this->CurrentImageValue == VTK_FLOAT_MAX )
      {
      sprintf(this->TextBuff,"Off Image");
      }
    else
      {
      float val = this->ImageData->GetScalarComponentAsFloat( \
                   this->CurrentCursorPosition[0],
                   this->CurrentCursorPosition[1],
                   this->CurrentCursorPosition[2],0);
      sprintf(this->TextBuff,"( %3d, %3d, %3d ): %g",
                   this->CurrentCursorPosition[0],
                   this->CurrentCursorPosition[1],
                   this->CurrentCursorPosition[2],val);
      }
    }

  this->TextActor->SetInput(this->TextBuff);
  this->TextActor->Modified();
}

void vtkImagePlaneWidget2::Push(double *p1, double *p2)
{
  // Get the motion vector
  //
  float v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];
  
  this->PlaneSource->Push( vtkMath::Dot(v,this->Normal) );
  this->PlaneSource->Update();
  this->BuildRepresentation();
}

void vtkImagePlaneWidget2::CreateDefaultProperties()
{
  if ( ! this->PlaneProperty )
    {
    this->PlaneProperty = vtkProperty::New();
    this->PlaneProperty->SetAmbient(1);
    this->PlaneProperty->SetColor(1,1,1);
    this->PlaneProperty->SetRepresentationToWireframe();
    this->PlaneProperty->SetInterpolationToFlat();
    }

  if ( ! this->SelectedPlaneProperty )
    {
    this->SelectedPlaneProperty = vtkProperty::New();
    this->SelectedPlaneProperty->SetAmbient(1);
    this->SelectedPlaneProperty->SetColor(0,1,0);
    this->SelectedPlaneProperty->SetRepresentationToWireframe();
    this->SelectedPlaneProperty->SetInterpolationToFlat();
    }

  if ( ! this->CursorProperty )
    {
    this->CursorProperty = vtkProperty::New();
    this->CursorProperty->SetAmbient(1);
    this->CursorProperty->SetColor(1,0,0);
    this->CursorProperty->SetRepresentationToWireframe();
    this->CursorProperty->SetInterpolationToFlat();
    }

  if ( ! this->MarginProperty )
    {
    this->MarginProperty = vtkProperty::New();
    this->MarginProperty->SetAmbient(1);
    this->MarginProperty->SetColor(0,0,1);
    this->MarginProperty->SetRepresentationToWireframe();
    this->MarginProperty->SetInterpolationToFlat();
    }

  if ( ! this->TexturePlaneProperty )
    {
    this->TexturePlaneProperty = vtkProperty::New();
    this->TexturePlaneProperty->SetAmbient(1);
    this->TexturePlaneProperty->SetInterpolationToFlat();
    }
}

void vtkImagePlaneWidget2::PlaceWidget(float bds[6])
{
  float bounds[6], center[3];

  this->AdjustBounds(bds, bounds, center);

  if ( this->PlaneOrientation == 1 )
    {
    this->PlaneSource->SetOrigin(bounds[0],center[1],bounds[4]);
    this->PlaneSource->SetPoint1(bounds[1],center[1],bounds[4]);
    this->PlaneSource->SetPoint2(bounds[0],center[1],bounds[5]);
    }
  else if ( this->PlaneOrientation == 2 )
    {
    this->PlaneSource->SetOrigin(bounds[0],bounds[2],center[2]);
    this->PlaneSource->SetPoint1(bounds[1],bounds[2],center[2]);
    this->PlaneSource->SetPoint2(bounds[0],bounds[3],center[2]);
    }
  else //default or x-normal
    {
    this->PlaneSource->SetOrigin(center[0],bounds[2],bounds[4]);
    this->PlaneSource->SetPoint1(center[0],bounds[3],bounds[4]);
    this->PlaneSource->SetPoint2(center[0],bounds[2],bounds[5]);
    }
  this->PlaneSource->Update();
  this->BuildRepresentation();
}

void vtkImagePlaneWidget2::SetPlaneOrientation(int i)
{
  // Generate a XY plane if i = 2, z-normal
  // or a YZ plane if i = 0, x-normal
  // or a ZX plane if i = 1, y-normal
  //
  this->PlaneOrientation = i;
  this->Modified();

  // This method must be called _after_ SetInput
  //
  this->ImageData = this->Reslice->GetInput();
  if ( !this->ImageData )
    {
    vtkErrorMacro(<<"SetInput() before setting plane orientation.");
    return;
    }
  this->ImageData->UpdateInformation();
  int extent[6];
  this->ImageData->GetWholeExtent(extent);
  float origin[3];
  this->ImageData->GetOrigin(origin);
  float spacing[3];
  this->ImageData->GetSpacing(spacing);

  // Prevent obscuring voxels by offsetting the plane geometry
  //
  float xbounds[] = {origin[0] + spacing[0] * (extent[0] - 0.5),
                     origin[0] + spacing[0] * (extent[1] + 0.5)};
  float ybounds[] = {origin[1] + spacing[1] * (extent[2] - 0.5),
                     origin[1] + spacing[1] * (extent[3] + 0.5)};
  float zbounds[] = {origin[2] + spacing[2] * (extent[4] - 0.5),
                     origin[2] + spacing[2] * (extent[5] + 0.5)};

  if ( spacing[0] < 0.0f )
    {
    float t = xbounds[0];
    xbounds[0] = xbounds[1];
    xbounds[1] = t;
    }
  if ( spacing[1] < 0.0f )
    {
    float t = ybounds[0];
    ybounds[0] = ybounds[1];
    ybounds[1] = t;
    }
  if ( spacing[2] < 0.0f )
    {
    float t = zbounds[0];
    zbounds[0] = zbounds[1];
    zbounds[1] = t;
    }

  if ( i == 2 ) //XY, z-normal
    {
    this->PlaneSource->SetOrigin(xbounds[0],ybounds[0],zbounds[0]);
    this->PlaneSource->SetPoint1(xbounds[1],ybounds[0],zbounds[0]);
    this->PlaneSource->SetPoint2(xbounds[0],ybounds[1],zbounds[0]);
    }
  else if ( i == 0 ) //YZ, x-normal
    {
    this->PlaneSource->SetOrigin(xbounds[0],ybounds[0],zbounds[0]);
    this->PlaneSource->SetPoint1(xbounds[0],ybounds[1],zbounds[0]);
    this->PlaneSource->SetPoint2(xbounds[0],ybounds[0],zbounds[1]);
    }
  else  //ZX, y-normal
    {
    this->PlaneSource->SetOrigin(xbounds[0],ybounds[0],zbounds[0]);
    this->PlaneSource->SetPoint1(xbounds[0],ybounds[0],zbounds[1]);
    this->PlaneSource->SetPoint2(xbounds[1],ybounds[0],zbounds[0]);
    }

  this->PlaneSource->Update();
  this->BuildRepresentation();
  this->UpdateNormal();
  this->UpdateOrigin();
}

void vtkImagePlaneWidget2::SetInput(vtkDataSet* input)
{
  this->Superclass::SetInput(input);
  this->ImageData = vtkImageData::SafeDownCast(this->GetInput());

  if( !this->ImageData )
    {
  // If NULL is passed, remove any reference that Reslice had
  // on the old ImageData
  //
    this->Reslice->SetInput(NULL);
    return;
    }

  float range[2];
  this->ImageData->GetScalarRange(range);

  if ( !this->UserControlledLookupTable )
    {
    this->LookupTable->SetTableRange(range[0],range[1]);
    this->LookupTable->Build();
    }

  this->OriginalWindow = range[1] - range[0];
  this->OriginalLevel = 0.5*(range[0] + range[1]);

  this->Reslice->SetInput(this->ImageData);
  this->SetResliceInterpolate(this->ResliceInterpolate);

  this->ColorMap->SetInput(this->Reslice->GetOutput());

  this->Texture->SetInput(this->ColorMap->GetOutput());
  this->Texture->SetInterpolate(this->TextureInterpolate);


  this->SetPlaneOrientation(this->PlaneOrientation);

  this->Locator->SetDataSet(this->ImageData);
  this->Locator->AutomaticOn();
  this->Locator->CacheCellBoundsOn();
  this->Locator->BuildLocator();
}

void vtkImagePlaneWidget2::UpdateOrigin()
{
  int i;

  if ( this->RestrictPlaneToVolume )
    {
    if ( !this->Reslice )
      {
      return;
      }
    this->ImageData = this->Reslice->GetInput();
    if ( !this->ImageData )
      {
      return;
      }
    this->ImageData->UpdateInformation();
    float origin[3];
    this->ImageData->GetOrigin(origin);
    float spacing[3];
    this->ImageData->GetSpacing(spacing);
    int extent[6];
    this->ImageData->GetWholeExtent(extent);
    float bounds[] = {origin[0] + spacing[0]*extent[0],
                      origin[0] + spacing[0]*extent[1],
                      origin[1] + spacing[1]*extent[2],
                      origin[1] + spacing[1]*extent[3],
                      origin[2] + spacing[2]*extent[4],
                      origin[2] + spacing[2]*extent[5]};

    for ( i = 0; i <= 4; i += 2 ) // reverse bounds if necessary
      {
      if ( bounds[i] > bounds[i+1] )
        {
        float t = bounds[i+1];
        bounds[i+1] = bounds[i];
        bounds[i] = t;
        }
      }

    float abs_normal[3];
    this->PlaneSource->GetNormal(abs_normal);
    float planeCenter[3];
    this->PlaneSource->GetCenter(planeCenter);
    float nmax = 0.0f;
    int k = 0;
    for ( i = 0; i < 3; i++ )
      {
      abs_normal[i] = fabs(abs_normal[i]);
      if ( abs_normal[i]>nmax )
        {
        nmax = abs_normal[i];
        k = i;
        }
      }
  // Force the plane to lie within the true image bounds along its normal
  //
    if ( planeCenter[k] > bounds[2*k+1] )
      {
      planeCenter[k] = bounds[2*k+1];
      this->PlaneSource->SetCenter(planeCenter);
      this->PlaneSource->Update();
      this->BuildRepresentation();
      }
    else if ( planeCenter[k] < bounds[2*k] )
      {
      planeCenter[k] = bounds[2*k];
      this->PlaneSource->SetCenter(planeCenter);
      this->PlaneSource->Update();
      this->BuildRepresentation();
      }
    }

  this->ResliceAxes->DeepCopy(this->Reslice->GetResliceAxes());
  this->ResliceAxes->SetElement(0,3,0);
  this->ResliceAxes->SetElement(1,3,0);
  this->ResliceAxes->SetElement(2,3,0);

  // Transpose is an exact way to invert a pure rotation matrix
  //
  this->ResliceAxes->Transpose();

  float planeOrigin[4];
  this->PlaneSource->GetOrigin(planeOrigin);
  planeOrigin[3] = 1.0;
  float originXYZW[4];
  this->ResliceAxes->MultiplyPoint(planeOrigin,originXYZW);

  this->ResliceAxes->Transpose();
  float neworiginXYZW[4];
  float point[] =  {0.0,0.0,originXYZW[2],1.0};
  this->ResliceAxes->MultiplyPoint(point,neworiginXYZW);

  this->ResliceAxes->SetElement(0,3,neworiginXYZW[0]);
  this->ResliceAxes->SetElement(1,3,neworiginXYZW[1]);
  this->ResliceAxes->SetElement(2,3,neworiginXYZW[2]);

  this->Reslice->SetResliceAxes(this->ResliceAxes);

  float spacingXYZ[3];
  this->Reslice->GetOutputSpacing(spacingXYZ);
  this->Reslice->SetOutputOrigin(0.5*spacingXYZ[0] + originXYZW[0],
                                 0.5*spacingXYZ[1] + originXYZW[1],
                                 0.0);
}

void vtkImagePlaneWidget2::UpdateNormal()
{
  float planeAxis1[3];
  float planeAxis2[3];

  this->GetVector1(planeAxis1);
  this->GetVector2(planeAxis2);

  // The x,y dimensions of the plane
  //
  float planeSizeX = vtkMath::Normalize(planeAxis1);
  float planeSizeY = vtkMath::Normalize(planeAxis2);

  this->PlaneSource->GetNormal(this->Normal);

  // Generate the slicing matrix
  //
  int i;
  this->ResliceAxes->Identity();
  for ( i = 0; i < 3; i++ )
     {
     this->ResliceAxes->SetElement(i,0,planeAxis1[i]);
     this->ResliceAxes->SetElement(i,1,planeAxis2[i]);
     this->ResliceAxes->SetElement(i,2,this->Normal[i]);
     }

  // Transpose is an exact way to invert a pure rotation matrix
  //
  this->ResliceAxes->Transpose();

  float planeOrigin[4];
  this->PlaneSource->GetOrigin(planeOrigin);
  planeOrigin[3] = 1.0;
  float originXYZW[4];
  this->ResliceAxes->MultiplyPoint(planeOrigin,originXYZW);

  this->ResliceAxes->Transpose();
  float neworiginXYZW[4];
  float point[] =  {0.0,0.0,originXYZW[2],1.0};
  this->ResliceAxes->MultiplyPoint(point,neworiginXYZW);

  this->ResliceAxes->SetElement(0,3,neworiginXYZW[0]);
  this->ResliceAxes->SetElement(1,3,neworiginXYZW[1]);
  this->ResliceAxes->SetElement(2,3,neworiginXYZW[2]);

  this->Reslice->SetResliceAxes(this->ResliceAxes);

  this->ImageData = this->Reslice->GetInput();

  // Calculate appropriate pixel spacing for the reslicing
  //
  this->ImageData->UpdateInformation();
  float spacing[3];
  this->ImageData->GetSpacing(spacing);

  float spacingX = fabs(planeAxis1[0]*spacing[0])+\
                   fabs(planeAxis1[1]*spacing[1])+\
                   fabs(planeAxis1[2]*spacing[2]);

  float spacingY = fabs(planeAxis2[0]*spacing[0])+\
                   fabs(planeAxis2[1]*spacing[1])+\
                   fabs(planeAxis2[2]*spacing[2]);


  // Pad extent up to a power of two for efficient texture mapping
  //
  int extentX = 1;
  while (extentX < planeSizeX/spacingX)
    {
    extentX = extentX << 1;
    }

  int extentY = 1;
  while (extentY < planeSizeY/spacingY)
    {
    extentY = extentY << 1;
    }

  this->Reslice->SetOutputSpacing(spacingX,spacingY,1);
  this->Reslice->SetOutputOrigin(0.5*spacingX + originXYZW[0],
                                 0.5*spacingY + originXYZW[1],
                                 0.0);

  this->Reslice->SetOutputExtent(0,extentX-1,0,extentY-1,0,0);

  // Find expansion factor to account for increasing the extent
  // to a power of two
  //
  float expand1 = extentX*spacingX;
  float expand2 = extentY*spacingY;

  // Set the texture coordinates to map the image to the plane
  //
  this->TexturePlaneCoords->SetOrigin(planeOrigin[0],
                                      planeOrigin[1],planeOrigin[2]);
  this->TexturePlaneCoords->SetPoint1(planeOrigin[0] + planeAxis1[0]*expand1,
                                      planeOrigin[1] + planeAxis1[1]*expand1,
                                      planeOrigin[2] + planeAxis1[2]*expand1);
  this->TexturePlaneCoords->SetPoint2(planeOrigin[0] + planeAxis2[0]*expand2,
                                      planeOrigin[1] + planeAxis2[1]*expand2,
                                      planeOrigin[2] + planeAxis2[2]*expand2);
}

vtkImageData* vtkImagePlaneWidget2::GetResliceOutput()
{
  if ( ! this->Reslice )
    {
    return 0;
    }
  return this->Reslice->GetOutput();
}

void vtkImagePlaneWidget2::SetResliceInterpolate(int i)
{
  if ( this->ResliceInterpolate == i )
    {
    return;
    }
  this->ResliceInterpolate = i;
  this->Modified();

  if ( !this->Reslice )
    {
    return;
    }
  
  if ( i == VTK_NEAREST_RESLICE )
    {
    this->Reslice->SetInterpolationModeToNearestNeighbor();
    }
  else if ( i == VTK_LINEAR_RESLICE)
    {
    this->Reslice->SetInterpolationModeToLinear();
    }
  else
    {
    this->Reslice->SetInterpolationModeToCubic();
    }
  this->Texture->SetInterpolate(this->TextureInterpolate);
}

void vtkImagePlaneWidget2::SetPicker(vtkCellPicker* picker)
{
  if ( this->UserPickerEnabled )
    {
    this->PlanePicker = picker;
    if (picker == 0 ) //reset and allocate an internal picker
      {
      this->PlanePicker = vtkCellPicker::New();
      this->UserPickerEnabled = 0;
      }
    }
  else
    {
    if (picker != 0 )
      {
      this->PlanePicker->Delete();
      this->PlanePicker = picker;
      this->UserPickerEnabled = 1;
      }
    else
      {
      return;
      }
    }

  this->PlanePicker->SetTolerance(0.005); //need some fluff
  this->PlanePicker->AddPickList(this->TexturePlaneActor);
  this->PlanePicker->PickFromListOn();
}

void vtkImagePlaneWidget2::SetLookupTable(vtkLookupTable* table)
{
  if ( this->UserLookupTableEnabled )
    {
    this->LookupTable = table;
    if ( table == 0 ) //reset and allocate an internal lut
      {
      this->LookupTable = vtkLookupTable::New();
      this->UserLookupTableEnabled = 0;
      }
    }
  else
    {
    if ( table != 0 )
      {
      this->LookupTable->Delete();
      this->LookupTable = table;
      this->UserLookupTableEnabled = 1;
      }
    else
      {
      return;
      }
    }

  if (!this->UserControlledLookupTable)
    {
    this->LookupTable->SetNumberOfColors( 256);
    this->LookupTable->SetHueRange( 0, 0);
    this->LookupTable->SetSaturationRange( 0, 0);
    this->LookupTable->SetValueRange( 0 ,1);
    this->LookupTable->SetAlphaRange( 1, 1);
    this->LookupTable->Build();
    }
  this->ColorMap->SetLookupTable(this->LookupTable);
  this->Texture->SetLookupTable(this->LookupTable);

  if( !this->ImageData )
    {
    return;
    }

  if (!this->UserControlledLookupTable)
    {
    float range[2];
    this->ImageData->GetScalarRange(range);

    this->LookupTable->SetTableRange(range[0],range[1]);
    this->LookupTable->Build();

    this->OriginalWindow = range[1] - range[0];
    this->OriginalLevel = 0.5*(range[0] + range[1]);
    }
}

void vtkImagePlaneWidget2::SetSlicePosition(float position)
{
  float amount = 0.0f;
  float planeOrigin[3];
  this->PlaneSource->GetOrigin(planeOrigin);

  if ( this->PlaneOrientation == 2 ) // z axis
    {
    amount = position - planeOrigin[2];
    }
  else if ( this->PlaneOrientation == 0 ) // x axis
    {
    amount = position - planeOrigin[0];
    }
  else if ( this->PlaneOrientation == 1 )  //y axis
    {
    amount = position - planeOrigin[1];
    }
  else
    {
    vtkGenericWarningMacro("only works for ortho planes: set plane orientation first");
    return;
    } 

  this->PlaneSource->Push(amount);
  this->PlaneSource->Update();
  this->BuildRepresentation();
  this->UpdateOrigin();
}

float vtkImagePlaneWidget2::GetSlicePosition()
{
  float planeOrigin[3];
  this->PlaneSource->GetOrigin(planeOrigin);

  if ( this->PlaneOrientation == 2 )
    {
    return planeOrigin[2];
    }
  else if ( this->PlaneOrientation == 1 )
    { 
    return planeOrigin[1];
    }
  else if ( this->PlaneOrientation == 0 )
    {  
    return planeOrigin[0];
    } 
  else
    {
    vtkGenericWarningMacro("only works for ortho planes: set plane orientation first");
    }

  return 0.0f;
}

void vtkImagePlaneWidget2::SetSliceIndex(int index)
{
  if ( !this->Reslice )
    {
      return;
    }
  this->ImageData = this->Reslice->GetInput();
  if ( !this->ImageData )
    {
    return;
    } 
  this->ImageData->UpdateInformation();
  float origin[3];
  this->ImageData->GetOrigin(origin);
  float spacing[3];
  this->ImageData->GetSpacing(spacing);
  float planeOrigin[3];
  this->PlaneSource->GetOrigin(planeOrigin);
  float pt1[3];
  this->PlaneSource->GetPoint1(pt1);
  float pt2[3];
  this->PlaneSource->GetPoint2(pt2);

  if ( this->PlaneOrientation == 2 )
    {
    planeOrigin[2] = origin[2] + index*spacing[2];
    pt1[2] = planeOrigin[2];
    pt2[2] = planeOrigin[2];
    }
  else if ( this->PlaneOrientation == 1 )
    {
    planeOrigin[1] = origin[1] + index*spacing[1]; 
    pt1[1] = planeOrigin[1];
    pt2[1] = planeOrigin[1];
    }
  else if ( this->PlaneOrientation == 0 )
    {
    planeOrigin[0] = origin[0] + index*spacing[0]; 
    pt1[0] = planeOrigin[0];
    pt2[0] = planeOrigin[0];
    }
  else
    {
    vtkGenericWarningMacro("only works for ortho planes: set plane orientation first");
    return; 
    }

  this->PlaneSource->SetOrigin(planeOrigin);
  this->PlaneSource->SetPoint1(pt1);
  this->PlaneSource->SetPoint2(pt2);
  this->PlaneSource->Update();
  this->BuildRepresentation();
  this->UpdateOrigin();
}

int vtkImagePlaneWidget2::GetSliceIndex()
{
  if ( ! this->Reslice )
    {
    return 0;
    }
  this->ImageData = this->Reslice->GetInput();
  if ( ! this->ImageData )
    {
    return 0;
    } 
  this->ImageData->UpdateInformation();
  float origin[3];
  this->ImageData->GetOrigin(origin);
  float spacing[3];
  this->ImageData->GetSpacing(spacing);
  float planeOrigin[3];
  this->PlaneSource->GetOrigin(planeOrigin);

  if ( this->PlaneOrientation == 2 )
    {
    return vtkMath::Round((planeOrigin[2]-origin[2])/spacing[2]);
    }
  else if ( this->PlaneOrientation == 1 )
    {
    return vtkMath::Round((planeOrigin[1]-origin[1])/spacing[1]);
    }
  else if ( this->PlaneOrientation == 0 )
    {
    return vtkMath::Round((planeOrigin[0]-origin[0])/spacing[0]);
    }
  else
    {
    vtkGenericWarningMacro("only works for ortho planes: set plane orientation first");
    }

  return 0;
}

void vtkImagePlaneWidget2::ActivateCursor(int i)
{

  if( !this->CurrentRenderer )
    {
    return;
    }

  if( i == 0 )
    {
    this->CursorActor->VisibilityOff();
    this->VoxelActor->VisibilityOff();
    }
  else
    {
    this->CursorActor->VisibilityOn();
    this->VoxelActor->VisibilityOn();
    }
}

void vtkImagePlaneWidget2::ActivateMargins(int i)
{

  if( !this->CurrentRenderer )
    {
    return;
    }

  if( i == 0 )
    {
    this->MarginActor->VisibilityOff();
    }
  else
    {
    this->MarginActor->VisibilityOn();
    }
}

void vtkImagePlaneWidget2::ActivateText(int i)
{
  if( !this->CurrentRenderer || !DisplayText)
    {
    return;
    }

  if( i == 0 )
    {
    this->TextActor->VisibilityOff();
    }
  else
    {
    this->TextActor->VisibilityOn();
    }
}

void vtkImagePlaneWidget2::UpdateCursor(int X, int Y )
{
  vtkAssemblyPath *path;
  this->PlanePicker->Pick(X,Y,0.0,this->CurrentRenderer);
  path = this->PlanePicker->GetPath();
  this->CurrentImageValue = VTK_FLOAT_MAX;

  int found = 0;
  int i;
  if ( path != 0 )
    {
  // Deal with the possibility that we may be using a shared picker
  //
    path->InitTraversal();
    vtkAssemblyNode *node;
    for ( i = 0; i< path->GetNumberOfItems() && !found ; i++ )
      {
      node = path->GetNextNode();
      if ( node->GetProp() == vtkProp::SafeDownCast(this->TexturePlaneActor) )
        {
        found = 1;
        }
      }
    }

  if( !found || path == 0 )
    {
    this->CursorActor->VisibilityOff();
    this->VoxelActor->VisibilityOff();
    return;
    }
  else
    {
    this->CursorActor->VisibilityOn();
    this->VoxelActor->VisibilityOn();
    }

  float q[3];
  this->PlanePicker->GetPickPosition(q);

  float closestPt[3];
  vtkIdType cellId;
  int subId;
  float dist2;
  this->Locator->FindClosestPoint(q,closestPt,this->Voxel,cellId,subId,dist2);

  if ( cellId == -1 )
    {
    this->CursorActor->VisibilityOff();
    this->VoxelActor->VisibilityOff();
    return;
    }

  this->Voxel->DeepCopy(this->ImageData->GetCell(cellId));
  this->VoxelGrid->SetPoints(this->Voxel->GetPoints());
  this->VoxelMapper->Update();

  float o[3];
  this->PlaneSource->GetOrigin(o);

  float origin[3];
  this->ImageData->GetOrigin(origin);
  float spacing[3];
  this->ImageData->GetSpacing(spacing);
  int extent[6];
  this->ImageData->GetExtent(extent);

  // Now query the original unsliced data
  //
  float qi[3];
  // Compute world to image coords
  for (int i = 0; i < 3; i++)
    {
    qi[i] = (closestPt[i]-origin[i])/spacing[i];
    }

  int iq[3];
  iq[0] = vtkMath::Round(qi[0]);
  iq[1] = vtkMath::Round(qi[1]);
  iq[2] = vtkMath::Round(qi[2]);

  if( iq[0] < extent[0] || iq[1] < extent[2] || iq[2] < extent[4] || \
      iq[0] > extent[1] || iq[1] > extent[3] || iq[2] > extent[5])
    {
    this->CursorActor->VisibilityOff();
    return;
    }
  else
    {
    memcpy(this->CurrentCursorPosition,iq,3*sizeof(int));
    this->CurrentImageValue = 0.0;
    }

  // Compute image to world coords
  for (i = 0; i < 3; i++)
    {
    q[i] = iq[i]*spacing[i] + origin[i];
    }

  // q relative to the plane origin
  //
  float qro[3];
  qro[0]= q[0] - o[0];
  qro[1]= q[1] - o[1];
  qro[2]= q[2] - o[2];

  float p1o[3];
  float p2o[3];

  this->GetVector1(p1o);
  this->GetVector2(p2o);

  float Lp1 = vtkMath::Dot(qro,p1o)/vtkMath::Dot(p1o,p1o);
  float Lp2 = vtkMath::Dot(qro,p2o)/vtkMath::Dot(p2o,p2o);

  float p1[3];
  this->PlaneSource->GetPoint1(p1);
  float p2[3];
  this->PlaneSource->GetPoint2(p2);

  float a[3];
  float b[3];
  float c[3];
  float d[3];

  for ( i = 0; i < 3; i++ )
    {
    a[i] = o[i]  + Lp2*p2o[i];   // left
    b[i] = p1[i] + Lp2*p2o[i];   // right
    c[i] = o[i]  + Lp1*p1o[i];   // bottom
    d[i] = p2[i] + Lp1*p1o[i];   // top
    }

  this->CursorPoints->SetPoint(0,a);
  this->CursorPoints->SetPoint(1,b);
  this->CursorPoints->SetPoint(2,c);
  this->CursorPoints->SetPoint(3,d);

  this->CursorMapper->Modified();
}

void vtkImagePlaneWidget2::ComputeWorldToImageCoords(float* in, float* out)
{
  this->ImageData = this->Reslice->GetInput();
  if( !this->ImageData )
    {
    return;
    }

  float origin[3];
  this->ImageData->GetOrigin(origin);
  float spacing[3];
  this->ImageData->GetSpacing(spacing);
  int extent[6];
  this->ImageData->GetExtent(extent);

  for (int i = 0; i < 3; i++)
    {
    out[i] = (in[i]-origin[i])/spacing[i];
    }
}

void vtkImagePlaneWidget2::ComputeImageToWorldCoords(float* in, float* out)
{
  this->ImageData = this->Reslice->GetInput();
  if( !this->ImageData )
    {
    return;
    }

  float origin[3];
  this->ImageData->GetOrigin(origin);
  float spacing[3];
  this->ImageData->GetSpacing(spacing);
  int extent[6];
  this->ImageData->GetExtent(extent);

  out[0] = in[0]*spacing[0] + origin[0];
  out[1] = in[1]*spacing[1] + origin[1];
  out[2] = in[2]*spacing[2] + origin[2];
}

void vtkImagePlaneWidget2::SetOrigin(float x, float y, float z)
{
  this->PlaneSource->SetOrigin(x,y,z);
}

void vtkImagePlaneWidget2::SetOrigin(float x[3])
{
  this->PlaneSource->SetOrigin(x);
}

float* vtkImagePlaneWidget2::GetOrigin()
{
  return this->PlaneSource->GetOrigin();
}

void vtkImagePlaneWidget2::GetOrigin(float xyz[3])
{
  this->PlaneSource->GetOrigin(xyz);
}

void vtkImagePlaneWidget2::SetPoint1(float x, float y, float z)
{
  this->PlaneSource->SetPoint1(x,y,z);
}

void vtkImagePlaneWidget2::SetPoint1(float x[3])
{
  this->PlaneSource->SetPoint1(x);
}

float* vtkImagePlaneWidget2::GetPoint1()
{
  return this->PlaneSource->GetPoint1();
}
void vtkImagePlaneWidget2::GetPoint1(float xyz[3])
{
  this->PlaneSource->GetPoint1(xyz);
}

void vtkImagePlaneWidget2::SetPoint2(float x, float y, float z)
{
  this->PlaneSource->SetPoint2(x,y,z);
}

void vtkImagePlaneWidget2::SetPoint2(float x[3])
{
  this->PlaneSource->SetPoint2(x);
}

float* vtkImagePlaneWidget2::GetPoint2()
{
  return this->PlaneSource->GetPoint2();
}

void vtkImagePlaneWidget2::GetPoint2(float xyz[3])
{
  this->PlaneSource->GetPoint2(xyz);
}

float* vtkImagePlaneWidget2::GetCenter() 
{
  return this->PlaneSource->GetCenter();
}

void vtkImagePlaneWidget2::GetCenter(float xyz[3]) 
{
  this->PlaneSource->GetCenter(xyz);
}

float* vtkImagePlaneWidget2::GetNormal() 
{
  return this->PlaneSource->GetNormal();
}

void vtkImagePlaneWidget2::GetNormal(float xyz[3])
{
  this->PlaneSource->GetNormal(xyz);
}

void vtkImagePlaneWidget2::GetPolyData(vtkPolyData *pd)
{
  pd->ShallowCopy(this->PlaneSource->GetOutput());
}

vtkPolyDataSource *vtkImagePlaneWidget2::GetPolyDataSource()
{
  return this->PlaneSource;
}

void vtkImagePlaneWidget2::UpdatePlacement(void)
{
  this->PlaneSource->Update();
  this->BuildRepresentation();
  this->UpdateNormal();
  this->UpdateOrigin();
  this->UpdateMargins();
}

void vtkImagePlaneWidget2::SetTextProperty(vtkTextProperty* tprop)
{
  this->TextActor->SetTextProperty(tprop);
}

vtkTextProperty* vtkImagePlaneWidget2::GetTextProperty()
{
  return this->TextActor->GetTextProperty();
}

vtkTexture *vtkImagePlaneWidget2::GetTexture()
{
  return this->Texture;
}

void vtkImagePlaneWidget2::GetVector1(float v1[3])
{
  float* p1 = this->PlaneSource->GetPoint1();
  float* o =  this->PlaneSource->GetOrigin();
  v1[0] = p1[0] - o[0];
  v1[1] = p1[1] - o[1];
  v1[2] = p1[2] - o[2];
}

void vtkImagePlaneWidget2::GetVector2(float v2[3])
{
  float* p2 = this->PlaneSource->GetPoint2();
  float* o =  this->PlaneSource->GetOrigin();
  v2[0] = p2[0] - o[0];
  v2[1] = p2[1] - o[1];
  v2[2] = p2[2] - o[2];
}

void vtkImagePlaneWidget2::AdjustState()
{
  if ( this->Interactor->GetShiftKey() )
    {
    this->State = vtkImagePlaneWidget2::Scaling;
    return;
    }

  float v1[3];
  this->GetVector1(v1);
  float v2[3];
  this->GetVector2(v2);
  float planeSize1 = vtkMath::Normalize(v1);
  float planeSize2 = vtkMath::Normalize(v2);
  float* planeOrigin = this->PlaneSource->GetOrigin();

  float ppo[3] = {this->LastPickPosition[0] - planeOrigin[0],
                  this->LastPickPosition[1] - planeOrigin[1],
                  this->LastPickPosition[2] - planeOrigin[2] };

  float x2D = vtkMath::Dot(ppo,v1);
  float y2D = vtkMath::Dot(ppo,v2);

  // Divide plane into three zones for different user interactions:
  // four corners -- spin around the plane's normal at its center
  // four edges   -- rotate around one of the plane's axes at its center
  // center area  -- push
  //
  float marginX = planeSize1 * 0.05;
  float marginY = planeSize2 * 0.05;

  float x0 = marginX;
  float y0 = marginY;
  float x1 = planeSize1 - marginX;
  float y1 = planeSize2 - marginY;

  if ( x2D < x0  )       // left margin
    {
    if (y2D < y0)        // bottom left corner
      {
      this->MarginSelectMode =  0;
      }
    else if (y2D > y1)   // top left corner
      {
      this->MarginSelectMode =  3;
      }
    else                 // left edge
      {
      this->MarginSelectMode =  4;
      }
    }
  else if ( x2D > x1 )   // right margin
    {
    if (y2D < y0)        // bottom right corner
      {
      this->MarginSelectMode =  1;
      }
    else if (y2D > y1)   // top right corner
      {
      this->MarginSelectMode =  2;
      }
    else                 // right edge
      {
      this->MarginSelectMode =  5;
      }
    }
  else                   // middle
    {
    if (y2D < y0)        // bottom edge
      {
      this->MarginSelectMode =  6;
      }
    else if (y2D > y1)   // top edge
      {
      this->MarginSelectMode =  7;
      }
    else                 // central area
      {
      this->MarginSelectMode =  8;
      }
    }

  if ( this->Interactor->GetControlKey() )
    {
    this->State = vtkImagePlaneWidget2::Moving;
    }
  else
    {
    if (this->MarginSelectMode >= 0 && this->MarginSelectMode < 4)
      {
      this->State = vtkImagePlaneWidget2::Spinning;
      return;
      }
    else if (this->MarginSelectMode == 8)
      {
      this->State = vtkImagePlaneWidget2::Pushing;
      return;
      }
    else
      {
      this->State = vtkImagePlaneWidget2::Rotating;
      }
    }

  float *raPtr = 0;
  float *rvPtr = 0;
  float rvfac = 1.0;
  float rafac = 1.0;

  switch ( this->MarginSelectMode )
    {
     // left bottom corner
    case 0: raPtr = v2; rvPtr = v1; rvfac = -1.0; rafac = -1.0; break;
     // right bottom corner
    case 1: raPtr = v2; rvPtr = v1;               rafac = -1.0; break;
     // right top corner
    case 2: raPtr = v2; rvPtr = v1;               break;
     // left top corner
    case 3: raPtr = v2; rvPtr = v1; rvfac = -1.0; break;
    case 4: raPtr = v2; rvPtr = v1; rvfac = -1.0; break; // left
    case 5: raPtr = v2; rvPtr = v1;               break; // right
    case 6: raPtr = v1; rvPtr = v2; rvfac = -1.0; break; // bottom
    case 7: raPtr = v1; rvPtr = v2;               break; // top
    default: raPtr = v1; rvPtr = v2; break;
    }

  for (int i = 0; i < 3; i++)
    {
    this->RotateAxis[i] = *raPtr++ * rafac;
    this->RadiusVector[i] = *rvPtr++ * rvfac;
    }
}

void vtkImagePlaneWidget2::Spin(double *p1, double *p2)
{
  // Disable cursor snap
  //
  this->PlaneOrientation = 3;
  
  // Get the motion vector, in world coords
  //
  float v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  // Plane center and normal before transform
  //
  float* wc = this->PlaneSource->GetCenter();
  float* wn = this->Normal;

  // Radius vector from center to cursor position
  //
  float rv[3] = {p2[0]-wc[0], p2[1]-wc[1], p2[2]-wc[2]};

  // Distance between center and cursor location
  //
  float rs = vtkMath::Normalize(rv);

  // Spin direction
  //
  float wn_cross_rv[3];
  vtkMath::Cross(wn,rv,wn_cross_rv);

  // Spin angle
  //
  float dw = vtkMath::RadiansToDegrees() * vtkMath::Dot(v,wn_cross_rv) / rs;

  this->Transform->Identity();
  this->Transform->Translate(wc[0],wc[1],wc[2]);
  this->Transform->RotateWXYZ(dw,wn);
  this->Transform->Translate(-wc[0],-wc[1],-wc[2]);

  float newpt[3];
  this->Transform->TransformPoint(this->PlaneSource->GetPoint1(),newpt);
  this->PlaneSource->SetPoint1(newpt);
  this->Transform->TransformPoint(this->PlaneSource->GetPoint2(),newpt);
  this->PlaneSource->SetPoint2(newpt);
  this->Transform->TransformPoint(this->PlaneSource->GetOrigin(),newpt);
  this->PlaneSource->SetOrigin(newpt);

  this->PlaneSource->Update();
  this->BuildRepresentation();
}

void vtkImagePlaneWidget2::Rotate(double *p1, double *p2, double *vpn)
{
  // Disable cursor snap
  //
  this->PlaneOrientation = 3;

  // Get the motion vector, in world coords
  //
  float v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  // Plane center and normal
  //
  float* wc = this->PlaneSource->GetCenter();

  // Radius of the rotating circle of the picked point
  //
  float radius = fabs( this->RadiusVector[0]*(p2[0]-wc[0]) +
                       this->RadiusVector[1]*(p2[1]-wc[1]) +
                       this->RadiusVector[2]*(p2[2]-wc[2]) );

  // Rotate direction ra_cross_rv
  //
  float rd[3];
  vtkMath::Cross(this->RotateAxis,this->RadiusVector,rd);

  // Direction cosin between rotating direction and view normal
  //
  float rd_dot_vpn = rd[0]*vpn[0] + rd[1]*vpn[1] + rd[2]*vpn[2];

  // 'push' plane edge when mouse moves away from plane center
  // 'pull' plane edge when mouse moves toward plane center
  //
  float dw = vtkMath::RadiansToDegrees() * (vtkMath::Dot(this->RadiusVector,v))/radius * (-rd_dot_vpn);

  this->Transform->Identity();
  this->Transform->Translate(wc[0],wc[1],wc[2]);
  this->Transform->RotateWXYZ(dw,this->RotateAxis);
  this->Transform->Translate(-wc[0],-wc[1],-wc[2]);

  float newpt[3];
  this->Transform->TransformPoint(this->PlaneSource->GetPoint1(),newpt);
  this->PlaneSource->SetPoint1(newpt);
  this->Transform->TransformPoint(this->PlaneSource->GetPoint2(),newpt);
  this->PlaneSource->SetPoint2(newpt);
  this->Transform->TransformPoint(this->PlaneSource->GetOrigin(),newpt);
  this->PlaneSource->SetOrigin(newpt);

  this->PlaneSource->Update();
  this->BuildRepresentation();
}

void vtkImagePlaneWidget2::GeneratePlaneOutline()
{
  vtkPoints* points   = vtkPoints::New(VTK_DOUBLE);
  points->SetNumberOfPoints(4);
  int i;
  for (i = 0; i < 4; i++)
    {
    points->SetPoint(i,0.0,0.0,0.0);
    }

  vtkCellArray *cells = vtkCellArray::New();
  cells->Allocate(cells->EstimateSize(4,2));
  vtkIdType pts[2];
  pts[0] = 3; pts[1] = 2;       // top edge
  cells->InsertNextCell(2,pts);
  pts[0] = 0; pts[1] = 1;       // bottom edge
  cells->InsertNextCell(2,pts);
  pts[0] = 0; pts[1] = 3;       // left edge
  cells->InsertNextCell(2,pts);
  pts[0] = 1; pts[1] = 2;       // right edge
  cells->InsertNextCell(2,pts);

  this->PlaneOutlinePolyData->SetPoints(points);
  points->Delete();
  this->PlaneOutlinePolyData->SetLines(cells);
  cells->Delete();

  this->PlaneOutlineMapper->SetInput( this->PlaneOutlinePolyData );
  this->PlaneOutlineMapper->SetResolveCoincidentTopologyToPolygonOffset();
  this->PlaneOutlineActor->SetMapper(this->PlaneOutlineMapper);
  this->PlaneOutlineActor->PickableOff();
}

void vtkImagePlaneWidget2::GenerateTexturePlane()
{
  if (!this->UserControlledLookupTable)
    {
    this->LookupTable->SetNumberOfColors( 256);
    this->LookupTable->SetHueRange( 0, 0);
    this->LookupTable->SetSaturationRange( 0, 0);
    this->LookupTable->SetValueRange( 0 ,1);
    this->LookupTable->SetAlphaRange( 1, 1);
    this->LookupTable->Build();
    }
  this->SetResliceInterpolate(this->ResliceInterpolate);

  this->ColorMap->SetLookupTable(this->LookupTable);
  this->ColorMap->SetOutputFormatToRGBA();
  this->ColorMap->PassAlphaToOutputOn();

  this->TexturePlaneCoords->SetInput(this->PlaneSource->GetOutput());
  this->TexturePlaneCoords->AutomaticPlaneGenerationOff();

  this->TexturePlaneMapper->SetInput(this->TexturePlaneCoords->GetOutput());

  this->Texture->SetQualityTo32Bit();
  this->Texture->MapColorScalarsThroughLookupTableOff();
  this->Texture->SetInterpolate(this->TextureInterpolate);
  this->Texture->RepeatOff();
  this->Texture->SetLookupTable(this->LookupTable);

  this->TexturePlaneActor->SetMapper(this->TexturePlaneMapper);
  this->TexturePlaneActor->SetTexture(this->Texture);
  this->TexturePlaneActor->PickableOn();
}

void vtkImagePlaneWidget2::GenerateMargins()
{
  // Construct initial points
  this->MarginPoints->SetNumberOfPoints(8);
  int i;
  for (i = 0; i < 8; i++)
    {
    this->MarginPoints->SetPoint(i,0.0,0.0,0.0);
    }

  vtkCellArray *cells = vtkCellArray::New();
  cells->Allocate(cells->EstimateSize(4,2));
  vtkIdType pts[2];
  pts[0] = 0; pts[1] = 1;       // top margin
  cells->InsertNextCell(2,pts);
  pts[0] = 2; pts[1] = 3;       // bottom margin
  cells->InsertNextCell(2,pts);
  pts[0] = 4; pts[1] = 5;       // left margin
  cells->InsertNextCell(2,pts);
  pts[0] = 6; pts[1] = 7;       // right margin
  cells->InsertNextCell(2,pts);

  this->MarginPolyData->SetPoints(this->MarginPoints);
  this->MarginPolyData->SetLines(cells);
  cells->Delete();

  this->MarginMapper->SetInput(this->MarginPolyData);
  this->MarginMapper->SetResolveCoincidentTopologyToPolygonOffset();
  this->MarginActor->SetMapper(this->MarginMapper);
  this->MarginActor->PickableOff();
  this->MarginActor->VisibilityOff();
}

void vtkImagePlaneWidget2::GenerateCursor()
{
  // Construct initial points
  //
  this->CursorPoints->SetNumberOfPoints(4);
  int i;
  for (i = 0; i < 4; i++)
    {
    this->CursorPoints->SetPoint(i,0.0,0.0,0.0);
    }

  vtkCellArray *cells = vtkCellArray::New();
  cells->Allocate(cells->EstimateSize(2,2));
  vtkIdType pts[2];
  pts[0] = 0; pts[1] = 1;       // horizontal segment
  cells->InsertNextCell(2,pts);
  pts[0] = 2; pts[1] = 3;       // vertical segment
  cells->InsertNextCell(2,pts);

  this->CursorPolyData->SetPoints(this->CursorPoints);
  this->CursorPolyData->SetLines(cells);
  cells->Delete();

  this->CursorMapper->SetInput(this->CursorPolyData);
  this->CursorMapper->SetResolveCoincidentTopologyToPolygonOffset();
  this->CursorActor->SetMapper(this->CursorMapper);
  this->CursorActor->PickableOff();
  this->CursorActor->VisibilityOff();
}

void vtkImagePlaneWidget2::GenerateText()
{
  sprintf(this->TextBuff,"NA");
  this->TextActor->SetInput(this->TextBuff);
  this->TextActor->ScaledTextOff();

  vtkTextProperty* textprop = this->TextActor->GetTextProperty();
  textprop->SetColor(1,1,1);
  textprop->SetFontFamilyToArial();
  textprop->SetFontSize(18);
  textprop->BoldOff();
  textprop->ItalicOff();
  textprop->ShadowOff();
  textprop->SetJustificationToLeft();
  textprop->SetVerticalJustificationToBottom();

  vtkCoordinate* coord = this->TextActor->GetPositionCoordinate();
  coord->SetCoordinateSystemToNormalizedDisplay();
  coord->SetValue(0.01, 0.01);

  this->TextActor->VisibilityOff();
}

void vtkImagePlaneWidget2::UpdateMargins()
{
  float v1[3];
  this->GetVector1(v1);
  float v2[3];
  this->GetVector2(v2);
  float o[3];
  this->PlaneSource->GetOrigin(o);
  float p1[3];
  this->PlaneSource->GetPoint1(p1);
  float p2[3];
  this->PlaneSource->GetPoint2(p2);

  float a[3];
  float b[3];
  float c[3];
  float d[3];

  float s = 0.05;
  float t = 0.05;

  int i;
  for ( i = 0; i < 3; i++)
    {
    a[i] = o[i] + v2[i]*(1-t);
    b[i] = p1[i] + v2[i]*(1-t);
    c[i] = o[i] + v2[i]*t;
    d[i] = p1[i] + v2[i]*t;
    }

  this->MarginPoints->SetPoint(0,a);
  this->MarginPoints->SetPoint(1,b);
  this->MarginPoints->SetPoint(2,c);
  this->MarginPoints->SetPoint(3,d);

  for ( i = 0; i < 3; i++)
    {
    a[i] = o[i] + v1[i]*s;
    b[i] = p2[i] + v1[i]*s;
    c[i] = o[i] + v1[i]*(1-s);
    d[i] = p2[i] + v1[i]*(1-s);
    }

  this->MarginPoints->SetPoint(4,a);
  this->MarginPoints->SetPoint(5,b);
  this->MarginPoints->SetPoint(6,c);
  this->MarginPoints->SetPoint(7,d);

  this->MarginMapper->Modified();
}

void vtkImagePlaneWidget2::Translate(double *p1, double *p2)
{
  // Get the motion vector
  //
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  float *o = this->PlaneSource->GetOrigin();
  float *pt1 = this->PlaneSource->GetPoint1();
  float *pt2 = this->PlaneSource->GetPoint2();
  float origin[3], point1[3], point2[3];

  float vdrv = this->RadiusVector[0]*v[0] + \
               this->RadiusVector[1]*v[1] + \
               this->RadiusVector[2]*v[2];
  float vdra = this->RotateAxis[0]*v[0] + \
               this->RotateAxis[1]*v[1] + \
               this->RotateAxis[2]*v[2];

  int i;
  if ( this->MarginSelectMode == 8 )       // everybody comes along
    {
    for (i=0; i<3; i++)
      {
      origin[i] = o[i] + v[i];
      point1[i] = pt1[i] + v[i];
      point2[i] = pt2[i] + v[i];
      }
    this->PlaneSource->SetOrigin(origin);
    this->PlaneSource->SetPoint1(point1);
    this->PlaneSource->SetPoint2(point2);
    }
  else if ( this->MarginSelectMode == 4 ) // left edge
    {
    for (i=0; i<3; i++)
      {
      origin[i] = o[i]   + vdrv*this->RadiusVector[i];
      point2[i] = pt2[i] + vdrv*this->RadiusVector[i];
      }
    this->PlaneSource->SetOrigin(origin);
    this->PlaneSource->SetPoint2(point2);
    }
  else if ( this->MarginSelectMode == 5 ) // right edge
    {
    for (i=0; i<3; i++)
      {
      point1[i] = pt1[i] + vdrv*this->RadiusVector[i];
      }
    this->PlaneSource->SetPoint1(point1);
    }
  else if ( this->MarginSelectMode == 6 ) // bottom edge
    {
    for (i=0; i<3; i++)
      {
      origin[i] = o[i]   + vdrv*this->RadiusVector[i];
      point1[i] = pt1[i] + vdrv*this->RadiusVector[i];
      }
    this->PlaneSource->SetOrigin(origin);
    this->PlaneSource->SetPoint1(point1);
    }
  else if ( this->MarginSelectMode == 7 ) // top edge
    {
    for (i=0; i<3; i++)
      {
      point2[i] = pt2[i] + vdrv*this->RadiusVector[i];
      }
    this->PlaneSource->SetPoint2(point2);
    }
  else if ( this->MarginSelectMode == 3 ) // top left corner
    {
    for (i=0; i<3; i++)
      {
      origin[i] = o[i]   + vdrv*this->RadiusVector[i];
      point2[i] = pt2[i] + vdrv*this->RadiusVector[i] + vdra*this->RotateAxis[i];
      }
    this->PlaneSource->SetOrigin(origin);
    this->PlaneSource->SetPoint2(point2);
    }
  else if ( this->MarginSelectMode == 0 ) // bottom left corner
    {
    for (int i=0; i<3; i++)
      {
      origin[i] = o[i]   + vdrv*this->RadiusVector[i] + vdra*this->RotateAxis[i];
      point1[i] = pt1[i] + vdra*this->RotateAxis[i];
      point2[i] = pt2[i] + vdrv*this->RadiusVector[i];
      }
    this->PlaneSource->SetOrigin(origin);
    this->PlaneSource->SetPoint1(point1);
    this->PlaneSource->SetPoint2(point2);
    }
  else if ( this->MarginSelectMode == 2 ) // top right corner
    {
    for (i=0; i<3; i++)
      {
      point1[i] = pt1[i] + vdrv*this->RadiusVector[i];
      point2[i] = pt2[i] + vdra*this->RotateAxis[i];
      }
    this->PlaneSource->SetPoint1(point1);
    this->PlaneSource->SetPoint2(point2);
    }
  else                                   // bottom right corner
    {
    for (i=0; i<3; i++)
      {
      origin[i] = o[i]   + vdra*this->RotateAxis[i];
      point1[i] = pt1[i] + vdrv*this->RadiusVector[i] + vdra*this->RotateAxis[i];
      }
    this->PlaneSource->SetPoint1(point1);
    this->PlaneSource->SetOrigin(origin);
    }

  this->PlaneSource->Update(); 
  this->BuildRepresentation();
}

void vtkImagePlaneWidget2::Scale(double *p1, double *p2, int vtkNotUsed(X), int Y)
{
  // Get the motion vector
  //
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  float *o = this->PlaneSource->GetOrigin();
  float *pt1 = this->PlaneSource->GetPoint1();
  float *pt2 = this->PlaneSource->GetPoint2();

  float center[3];
  center[0] = o[0] + (pt1[0]-o[0])/2.0 + (pt2[0]-o[0])/2.0;
  center[1] = o[1] + (pt1[1]-o[1])/2.0 + (pt2[1]-o[1])/2.0;
  center[2] = o[2] + (pt1[2]-o[2])/2.0 + (pt2[2]-o[2])/2.0;

  // Compute the scale factor
  //
  float sf = vtkMath::Norm(v) / sqrt(vtkMath::Distance2BetweenPoints(pt1,pt2));
  if ( Y > this->Interactor->GetLastEventPosition()[1] )
    {
    sf = 1.0 + sf;
    }
  else
    {
    sf = 1.0 - sf;
    }

  // Move the corner points
  //
  float origin[3], point1[3], point2[3];

  for (int i=0; i<3; i++)
    {
    origin[i] = sf * (o[i] - center[i]) + center[i];
    point1[i] = sf * (pt1[i] - center[i]) + center[i];
    point2[i] = sf * (pt2[i] - center[i]) + center[i];
    }

  this->PlaneSource->SetOrigin(origin);
  this->PlaneSource->SetPoint1(point1);
  this->PlaneSource->SetPoint2(point2);
  this->PlaneSource->Update();
  this->BuildRepresentation();
}


void vtkImagePlaneWidget2::GenerateVoxel()
{
vtkPoints* voxelPoints = vtkPoints::New();

  voxelPoints->SetNumberOfPoints(8);
  voxelPoints->InsertPoint( 0, 0, 0, 0);
  voxelPoints->InsertPoint( 1, 1, 0, 0);
  voxelPoints->InsertPoint( 2, 0, 1, 0);
  voxelPoints->InsertPoint( 3, 1, 1, 0);
  voxelPoints->InsertPoint( 4, 0, 0, 1);
  voxelPoints->InsertPoint( 5, 1, 0, 1);
  voxelPoints->InsertPoint( 6, 0, 1, 1);
  voxelPoints->InsertPoint( 7, 1, 1, 1);

  vtkVoxel* aVoxel = vtkVoxel::New();
  vtkIdList* Ids = aVoxel->GetPointIds();
  for(int i=0;i<8;i++)
    {
    Ids->SetId( i,i);
    }

  this->Voxel->DeepCopy(aVoxel);

  this->VoxelGrid->Allocate( 1, 1);
  this->VoxelGrid->InsertNextCell(this->Voxel->GetCellType(),this->Voxel->GetPointIds());
  this->VoxelGrid->SetPoints(voxelPoints);

  this->VoxelMapper->SetInput(VoxelGrid);

  this->VoxelActor->SetMapper(VoxelMapper);

  this->VoxelActor->PickableOff();
  this->VoxelActor->VisibilityOff();

  vtkProperty* prop = this->VoxelActor->GetProperty();
  prop->SetAmbient(1);
  prop->SetColor(0,0,1);
  prop->SetRepresentationToWireframe();
  prop->SetInterpolationToFlat();

  voxelPoints->Delete();
  aVoxel->Delete();
}




