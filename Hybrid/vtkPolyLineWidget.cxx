#include "vtkPolyLineWidget.h"

#include "vtkActor.h"
#include "vtkAssemblyNode.h"
#include "vtkAssemblyPath.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellPicker.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPlaneSource.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"

vtkCxxRevisionMacro(vtkPolyLineWidget, "$Revision: 1.3 $");
vtkStandardNewMacro(vtkPolyLineWidget);

vtkCxxSetObjectMacro(vtkPolyLineWidget, HandleProperty, vtkProperty);
vtkCxxSetObjectMacro(vtkPolyLineWidget, SelectedHandleProperty, vtkProperty);
vtkCxxSetObjectMacro(vtkPolyLineWidget, LineProperty, vtkProperty);
vtkCxxSetObjectMacro(vtkPolyLineWidget, SelectedLineProperty, vtkProperty);

vtkPolyLineWidget::vtkPolyLineWidget()
{
  this->State = vtkPolyLineWidget::Start;
  this->EventCallbackCommand->SetCallback(vtkPolyLineWidget::ProcessEvents);
  this->ProjectToPlane = 0;  //default off
  this->ProjectionNormal = 0;  //default YZ not used
  this->ProjectionPosition = 0.0;
  this->PlaneSource = NULL;
  this->Closed = 0;
  this->Offset = 0.0;

  // Default bounds to get started
  double bounds[6];
  bounds[0] = -0.5;
  bounds[1] = 0.5;
  bounds[2] = -0.5;
  bounds[3] = 0.5;
  bounds[4] = -0.5;
  bounds[5] = 0.5;

  // Create the handles along a straight line within the data bounds
  this->NumberOfHandles = 5;
  this->Handle         = new vtkActor* [this->NumberOfHandles];
  this->HandleMapper   = new vtkPolyDataMapper* [this->NumberOfHandles];
  this->HandleGeometry = new vtkSphereSource* [this->NumberOfHandles];
  int i;
  double position;
  double x0 = bounds[0];
  double x1 = bounds[1];
  double y0 = bounds[2];
  double y1 = bounds[3];
  double z0 = bounds[4];
  double z1 = bounds[5];
  double x;
  double y;
  double z;

  vtkPoints* points = vtkPoints::New();
  points->Allocate(this->NumberOfHandles);
  
  for (i=0; i<this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i] = vtkSphereSource::New();
    this->HandleGeometry[i]->SetThetaResolution(16);
    this->HandleGeometry[i]->SetPhiResolution(8);
    this->HandleMapper[i] = vtkPolyDataMapper::New();
    this->HandleMapper[i]->SetInput(this->HandleGeometry[i]->GetOutput());
    this->Handle[i] = vtkActor::New();
    this->Handle[i]->SetMapper(this->HandleMapper[i]);
    position = i / (this->NumberOfHandles - 1.0);
    x = (1.0-position)*x0 + position*x1;
    y = (1.0-position)*y0 + position*y1;
    z = (1.0-position)*z0 + position*z1;
    this->HandleGeometry[i]->SetCenter(x,y,z);

    points->InsertPoint(i, x, y, z);
    }

  vtkCellArray* lines = vtkCellArray::New();
  lines->Allocate(lines->EstimateSize(this->NumberOfHandles,2));
  lines->InsertNextCell(this->NumberOfHandles);

  for (i=0 ; i<this->NumberOfHandles; i++)
    {
    lines->InsertCellPoint(i);
    }

  this->LineData = vtkPolyData::New();
  this->LineData->SetPoints(points);
  points->Delete();
  this->LineData->SetLines(lines);
  lines->Delete();

  this->LineMapper = vtkPolyDataMapper::New();
  this->LineMapper->SetInput( this->LineData ) ;
  this->LineMapper->ImmediateModeRenderingOn();
  this->LineMapper->SetResolveCoincidentTopologyToPolygonOffset();

  this->LineActor = vtkActor::New();
  this->LineActor->SetMapper( this->LineMapper);

  // Initial creation of the widget, serves to initialize it
  this->PlaceFactor = 1.0;
  this->PlaceWidget(bounds);

  // Manage the picking stuff
  this->HandlePicker = vtkCellPicker::New();
  this->HandlePicker->SetTolerance(0.005);
  for (i=0; i<this->NumberOfHandles; i++)
    {
    this->HandlePicker->AddPickList(this->Handle[i]);
    }
  this->HandlePicker->PickFromListOn();

  this->LinePicker = vtkCellPicker::New();
  this->LinePicker->SetTolerance(0.01);
  this->LinePicker->AddPickList(this->LineActor);
  this->LinePicker->PickFromListOn();

  this->CurrentHandle = NULL;
  this->CurrentHandleIndex = -1;

  this->Transform = vtkTransform::New();

  // Set up the initial properties
  this->HandleProperty = NULL;
  this->SelectedHandleProperty = NULL;
  this->LineProperty = NULL;
  this->SelectedLineProperty = NULL;
  this->CreateDefaultProperties();
}

vtkPolyLineWidget::~vtkPolyLineWidget()
{
  this->LineActor->Delete();
  this->LineMapper->Delete();
  this->LineData->Delete();

  for (int i=0; i<this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i]->Delete();
    this->HandleMapper[i]->Delete();
    this->Handle[i]->Delete();
    }
  delete [] this->Handle;
  delete [] this->HandleMapper;
  delete [] this->HandleGeometry;

  this->HandlePicker->Delete();
  this->LinePicker->Delete();

  if ( this->HandleProperty )
    {
    this->HandleProperty->Delete();
    }
  if ( this->SelectedHandleProperty )
    {
    this->SelectedHandleProperty->Delete();
    }
  if ( this->LineProperty )
    {
    this->LineProperty->Delete();
    }
  if ( this->SelectedLineProperty )
    {
    this->SelectedLineProperty->Delete();
    }

  this->Transform->Delete();
}

void vtkPolyLineWidget::SetClosed(int closed)
{
  if (this->Closed == closed)
    {
    return;
    }
  this->Closed = closed;

  this->BuildRepresentation();
}


void vtkPolyLineWidget::SetHandlePosition(int handle, double x, 
                                        double y, double z)
{
  if(handle < 0 || handle >= this->NumberOfHandles)
    {
    vtkErrorMacro(<<"vtkPolyLineWidget: handle index out of range.");
    return;
    }
  this->HandleGeometry[handle]->SetCenter(x,y,z);
  this->HandleGeometry[handle]->Update();
  if ( this->ProjectToPlane )
    {
    this->ProjectPointsToPlane();
    }
  this->BuildRepresentation();
}

void vtkPolyLineWidget::SetHandlePosition(int handle, double xyz[3])
{
  this->SetHandlePosition(handle,xyz[0],xyz[1],xyz[2]);
}

void vtkPolyLineWidget::GetHandlePosition(int handle, double xyz[3])
{
  if(handle < 0 || handle >= this->NumberOfHandles)
    {
    vtkErrorMacro(<<"vtkPolyLineWidget: handle index out of range.");
    return;
    }

  this->HandleGeometry[handle]->GetCenter(xyz);
}

double* vtkPolyLineWidget::GetHandlePosition(int handle)
{
  if(handle < 0 || handle >= this->NumberOfHandles)
    {
    vtkErrorMacro(<<"vtkPolyLineWidget: handle index out of range.");
    return NULL;
    }

  return this->HandleGeometry[handle]->GetCenter();
}

void vtkPolyLineWidget::SetEnabled(int enabling)
{
  if ( ! this->Interactor )
    {
    vtkErrorMacro(<<"The interactor must be set prior to enabling/disabling widget");
    return;
    }

  if ( enabling ) //------------------------------------------------------------
    {
    vtkDebugMacro(<<"Enabling line widget");

    if ( this->Enabled ) //already enabled, just return
      {
      return;
      }

    if ( ! this->CurrentRenderer )
      {
      this->SetCurrentRenderer(this->Interactor->FindPokedRenderer(
        this->Interactor->GetLastEventPosition()[0],
        this->Interactor->GetLastEventPosition()[1]));
      if (this->CurrentRenderer == NULL)
        {
        return;
        }
      }

    this->Enabled = 1;

    // Listen for the following events
    vtkRenderWindowInteractor *i = this->Interactor;
    i->AddObserver(vtkCommand::MouseMoveEvent, this->EventCallbackCommand,
                   this->Priority);
    i->AddObserver(vtkCommand::LeftButtonPressEvent, this->EventCallbackCommand,
                   this->Priority);
    i->AddObserver(vtkCommand::LeftButtonReleaseEvent, this->EventCallbackCommand,
                   this->Priority);
    i->AddObserver(vtkCommand::MiddleButtonPressEvent, this->EventCallbackCommand,
                   this->Priority);
    i->AddObserver(vtkCommand::MiddleButtonReleaseEvent, this->EventCallbackCommand,
                   this->Priority);
    i->AddObserver(vtkCommand::RightButtonPressEvent, this->EventCallbackCommand,
                   this->Priority);
    i->AddObserver(vtkCommand::RightButtonReleaseEvent, this->EventCallbackCommand,
                   this->Priority);

    // Add the line
    this->CurrentRenderer->AddActor(this->LineActor);
    this->LineActor->SetProperty(this->LineProperty);

    // Turn on the handles
    for (int j=0; j<this->NumberOfHandles; j++)
      {
      this->CurrentRenderer->AddActor(this->Handle[j]);
      this->Handle[j]->SetProperty(this->HandleProperty);
      }
    this->BuildRepresentation();
    this->SizeHandles();

    this->InvokeEvent(vtkCommand::EnableEvent,NULL);
    }

  else //disabling----------------------------------------------------------
    {
    vtkDebugMacro(<<"Disabling line widget");

    if ( ! this->Enabled ) //already disabled, just return
      {
      return;
      }

    this->Enabled = 0;

    // Don't listen for events any more
    this->Interactor->RemoveObserver(this->EventCallbackCommand);

    // Turn off the line
    this->CurrentRenderer->RemoveActor(this->LineActor);

    // Turn off the handles
    for (int i=0; i<this->NumberOfHandles; i++)
      {
      this->CurrentRenderer->RemoveActor(this->Handle[i]);
      }

    this->CurrentHandle = NULL;
    this->InvokeEvent(vtkCommand::DisableEvent,NULL);
    this->SetCurrentRenderer(NULL);
    }

  this->Interactor->Render();
}

void vtkPolyLineWidget::ProcessEvents(vtkObject* vtkNotUsed(object),
                                  unsigned long event,
                                  void* clientdata,
                                  void* vtkNotUsed(calldata))
{
  vtkPolyLineWidget* self = reinterpret_cast<vtkPolyLineWidget *>( clientdata );

  // Okay, let's do the right thing
  switch(event)
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

void vtkPolyLineWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  if ( this->HandleProperty )
    {
    os << indent << "Handle Property: " << this->HandleProperty << "\n";
    }
  else
    {
    os << indent << "Handle Property: (none)\n";
    }
  if ( this->SelectedHandleProperty )
    {
    os << indent << "Selected Handle Property: "
       << this->SelectedHandleProperty << "\n";
    }
  else
    {
    os << indent << "Selected Handle Property: (none)\n";
    }
  if ( this->LineProperty )
    {
    os << indent << "Line Property: " << this->LineProperty << "\n";
    }
  else
    {
    os << indent << "Line Property: (none)\n";
    }
  if ( this->SelectedLineProperty )
    {
    os << indent << "Selected Line Property: "
       << this->SelectedLineProperty << "\n";
    }
  else
    {
    os << indent << "Selected Line Property: (none)\n";
    }

  os << indent << "Project To Plane: "
     << (this->ProjectToPlane ? "On" : "Off") << "\n";
  os << indent << "Projection Normal: " << this->ProjectionNormal << "\n";
  os << indent << "Projection Position: " << this->ProjectionPosition << "\n";
  os << indent << "Number Of Handles: " << this->NumberOfHandles << "\n";
  os << indent << "Closed: "
     << (this->Closed ? "On" : "Off") << "\n";
}

void vtkPolyLineWidget::ProjectPointsToPlane()
{
  if (this->ProjectionNormal == VTK_PROJECTION_OBLIQUE)
    {
    if(this->PlaneSource != NULL)
      {
      this->ProjectPointsToObliquePlane();
      }
    else
      {
      vtkGenericWarningMacro(<<"Set the plane source for oblique projections...");
      }
    }
  else
    {
    this->ProjectPointsToOrthoPlane();
    }
}

void vtkPolyLineWidget::ProjectPointsToObliquePlane()
{
  double o[3];
  double u[3];
  double v[3];

  this->PlaneSource->GetPoint1(u);
  this->PlaneSource->GetPoint2(v);
  this->PlaneSource->GetOrigin(o);

  int i;
  for(i=0;i<3;i++)
    {
    u[i]=u[i]-o[i];
    v[i]=v[i]-o[i];
    }
  vtkMath::Normalize(u);
  vtkMath::Normalize(v);

  double o_dot_u = vtkMath::Dot(o,u);
  double o_dot_v = vtkMath::Dot(o,v);
  double fac1;
  double fac2;
  double ctr[3];
  for (i=0; i<this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i]->GetCenter(ctr);
    fac1 = vtkMath::Dot(ctr,u) - o_dot_u;
    fac2 = vtkMath::Dot(ctr,v) - o_dot_v;
    ctr[0] = o[0] + fac1*u[0] + fac2*v[0];
    ctr[1] = o[1] + fac1*u[1] + fac2*v[1];
    ctr[2] = o[2] + fac1*u[2] + fac2*v[2];
    this->HandleGeometry[i]->SetCenter(ctr);
    this->HandleGeometry[i]->Update();
    }
}

void vtkPolyLineWidget::ProjectPointsToOrthoPlane()
{
  double ctr[3];
  for (int i=0; i<this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i]->GetCenter(ctr);
    ctr[  this->ProjectionNormal ] = this->ProjectionPosition;
    this->HandleGeometry[i]->SetCenter(ctr);
    this->HandleGeometry[i]->Update();
    }
}

void vtkPolyLineWidget::BuildRepresentation()
{
  double ctr[3];
  int i;
  vtkPoints* points = this->LineData->GetPoints();
  points->SetNumberOfPoints(this->NumberOfHandles);
  for (i=0; i<this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i]->GetCenter(ctr);
    points->SetPoint(i, ctr[0], ctr[1], ctr[2]);
    }

  // let's do the lines now
  vtkCellArray *lines = this->LineData->GetLines();
  // deallocate and reset
  lines->Reset();

  int npts = this->NumberOfHandles;
  if (this->Closed)
    {
    npts += 1;
    }

  lines->Allocate(lines->EstimateSize(npts,2));
  lines->InsertNextCell(npts);

  for (i = 0; i < this->NumberOfHandles; i++)
    {
    lines->InsertCellPoint(i);
    }
  if (this->Closed)
    {
    // has to point back to the beginning if closed
    lines->InsertCellPoint(0);
    }
    
}

int vtkPolyLineWidget::HighlightHandle(vtkProp *prop)
{
  // First unhighlight anything picked
  if ( this->CurrentHandle )
    {
    this->CurrentHandle->SetProperty(this->HandleProperty);
    }

  this->CurrentHandle = (vtkActor *)prop;

  if ( this->CurrentHandle )
    {
    for (int i=0; i<this->NumberOfHandles; i++) // find handle
      {
      if ( this->CurrentHandle == this->Handle[i] )
        {
        this->ValidPick = 1;
        this->HandlePicker->GetPickPosition(this->LastPickPosition);
        this->CurrentHandle->SetProperty(this->SelectedHandleProperty);
        return i;
        }
      }
    }
  return -1;
}

void vtkPolyLineWidget::HighlightLine(int highlight)
{
  if ( highlight )
    {
    this->ValidPick = 1;
    this->LinePicker->GetPickPosition(this->LastPickPosition);
    this->LineActor->SetProperty(this->SelectedLineProperty);
    }
  else
    {
    this->LineActor->SetProperty(this->LineProperty);
    }
}

void vtkPolyLineWidget::OnLeftButtonDown()
{
  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  if (!this->CurrentRenderer || !this->CurrentRenderer->IsInViewport(X, Y))
    {
    this->State = vtkPolyLineWidget::Outside;
    return;
    }

  this->State = vtkPolyLineWidget::Moving;

  // Okay, we can process this. Try to pick handles first;
  // if no handles picked, then try to pick the line.
  vtkAssemblyPath *path;
  this->HandlePicker->Pick(X,Y,0.0,this->CurrentRenderer);
  path = this->HandlePicker->GetPath();
  if ( path != NULL )
    {
    this->CurrentHandleIndex = this->HighlightHandle(path->GetFirstNode()->GetViewProp());
    }
  else
    {
    this->LinePicker->Pick(X,Y,0.0,this->CurrentRenderer);
    path = this->LinePicker->GetPath();
    if ( path != NULL )
      {
      this->HighlightLine(1);
      }
    else
      {
      this->CurrentHandleIndex = this->HighlightHandle(NULL);
      this->State = vtkPolyLineWidget::Outside;
      return;
      }
    }

  this->EventCallbackCommand->SetAbortFlag(1);
  this->StartInteraction();
  this->InvokeEvent(vtkCommand::StartInteractionEvent,NULL);
  this->Interactor->Render();
}

void vtkPolyLineWidget::OnLeftButtonUp()
{
  if ( this->State == vtkPolyLineWidget::Outside ||
       this->State == vtkPolyLineWidget::Start )
    {
    return;
    }

  this->State = vtkPolyLineWidget::Start;
  this->HighlightHandle(NULL);
  this->HighlightLine(0);

  this->SizeHandles();

  this->EventCallbackCommand->SetAbortFlag(1);
  this->EndInteraction();
  this->InvokeEvent(vtkCommand::EndInteractionEvent,NULL);
  this->Interactor->Render();
}

void vtkPolyLineWidget::OnMiddleButtonDown()
{
  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  if (!this->CurrentRenderer || !this->CurrentRenderer->IsInViewport(X, Y))
    {
    this->State = vtkPolyLineWidget::Outside;
    return;
    }

  if ( this->Interactor->GetControlKey() )
    {
    this->State = vtkPolyLineWidget::Spinning;
    this->CalculateCentroid();
    }
  else
    {
    this->State = vtkPolyLineWidget::Moving;
    }

  // Okay, we can process this. Try to pick handles first;
  // if no handles picked, then try to pick the line.
  vtkAssemblyPath *path;
  this->HandlePicker->Pick(X,Y,0.0,this->CurrentRenderer);
  path = this->HandlePicker->GetPath();
  if ( path == NULL )
    {
    this->LinePicker->Pick(X,Y,0.0,this->CurrentRenderer);
    path = this->LinePicker->GetPath();
    if ( path == NULL )
      {
      this->State = vtkPolyLineWidget::Outside;
      this->HighlightLine(0);
      return;
      }
    else
      {
      this->HighlightLine(1);
      }
    }
  else  //we picked a handle but lets make it look like the line is picked
    {
    this->HighlightLine(1);
    }

  this->EventCallbackCommand->SetAbortFlag(1);
  this->StartInteraction();
  this->InvokeEvent(vtkCommand::StartInteractionEvent,NULL);
  this->Interactor->Render();
}

void vtkPolyLineWidget::OnMiddleButtonUp()
{
  if ( this->State == vtkPolyLineWidget::Outside ||
       this->State == vtkPolyLineWidget::Start )
    {
    return;
    }

  this->State = vtkPolyLineWidget::Start;
  this->HighlightLine(0);

  this->SizeHandles();

  this->EventCallbackCommand->SetAbortFlag(1);
  this->EndInteraction();
  this->InvokeEvent(vtkCommand::EndInteractionEvent,NULL);
  this->Interactor->Render();
}

void vtkPolyLineWidget::OnRightButtonDown()
{
  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  if (!this->CurrentRenderer || !this->CurrentRenderer->IsInViewport(X, Y))
    {
    this->State = vtkPolyLineWidget::Outside;
    return;
    }

  this->State = vtkPolyLineWidget::Scaling;

  // Okay, we can process this. Try to pick handles first;
  // if no handles picked, then pick the bounding box.
  vtkAssemblyPath *path;
  this->HandlePicker->Pick(X,Y,0.0,this->CurrentRenderer);
  path = this->HandlePicker->GetPath();
  if ( path == NULL )
    {
    this->LinePicker->Pick(X,Y,0.0,this->CurrentRenderer);
    path = this->LinePicker->GetPath();
    if ( path == NULL )
      {
      this->State = vtkPolyLineWidget::Outside;
      this->HighlightLine(0);
      return;
      }
    else
      {
      this->HighlightLine(1);
      }
    }
  else  //we picked a handle but lets make it look like the line is picked
    {
    this->HighlightLine(1);
    }

  this->EventCallbackCommand->SetAbortFlag(1);
  this->StartInteraction();
  this->InvokeEvent(vtkCommand::StartInteractionEvent,NULL);
  this->Interactor->Render();
}

void vtkPolyLineWidget::OnRightButtonUp()
{
  if ( this->State == vtkPolyLineWidget::Outside ||
       this->State == vtkPolyLineWidget::Start )
    {
    return;
    }

  this->State = vtkPolyLineWidget::Start;
  this->HighlightLine(0);

  this->SizeHandles();

  this->EventCallbackCommand->SetAbortFlag(1);
  this->EndInteraction();
  this->InvokeEvent(vtkCommand::EndInteractionEvent,NULL);
  this->Interactor->Render();
}

void vtkPolyLineWidget::OnMouseMove()
{
  // See whether we're active
  if ( this->State == vtkPolyLineWidget::Outside ||
       this->State == vtkPolyLineWidget::Start )
    {
    return;
    }

  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  // Do different things depending on state
  // Calculations everybody does
  double focalPoint[4], pickPoint[4], prevPickPoint[4];
  double z, vpn[3];

  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  if ( !camera )
    {
    return;
    }

  // Compute the two points defining the motion vector
  this->ComputeWorldToDisplay(this->LastPickPosition[0], this->LastPickPosition[1],
                              this->LastPickPosition[2], focalPoint);
  z = focalPoint[2];
  this->ComputeDisplayToWorld(double(this->Interactor->GetLastEventPosition()[0]),
                              double(this->Interactor->GetLastEventPosition()[1]),
                              z, prevPickPoint);
  this->ComputeDisplayToWorld(double(X), double(Y), z, pickPoint);

  // Process the motion
  if ( this->State == vtkPolyLineWidget::Moving )
    {
    // Okay to process
    if ( this->CurrentHandle )
      {
      this->MovePoint(prevPickPoint, pickPoint);
      }
    else // Must be moving the spline
      {
      this->Translate(prevPickPoint, pickPoint);
      }
    }
  else if ( this->State == vtkPolyLineWidget::Scaling )
    {
    this->Scale(prevPickPoint, pickPoint, X, Y);
    }
  else if ( this->State == vtkPolyLineWidget::Spinning )
    {
    camera->GetViewPlaneNormal(vpn);
    this->Spin(prevPickPoint, pickPoint, vpn);
    }

  if ( this->ProjectToPlane )
    {
    this->ProjectPointsToPlane();
    }
  this->BuildRepresentation();

  // Interact, if desired
  this->EventCallbackCommand->SetAbortFlag(1);
  this->InvokeEvent(vtkCommand::InteractionEvent,NULL);
  this->Interactor->Render();
}

void vtkPolyLineWidget::MovePoint(double *p1, double *p2)
{
  if ( this->CurrentHandleIndex < 0 || this->CurrentHandleIndex >= this->NumberOfHandles )
    {
    vtkGenericWarningMacro(<<"Spline handle index out of range.");
    return;
    }
  // Get the motion vector
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  double *ctr = this->HandleGeometry[ this->CurrentHandleIndex ]->GetCenter();

  double newCtr[3];
  newCtr[0] = ctr[0] + v[0];
  newCtr[1] = ctr[1] + v[1];
  newCtr[2] = ctr[2] + v[2];

  this->HandleGeometry[this->CurrentHandleIndex]->SetCenter(newCtr);
  this->HandleGeometry[this->CurrentHandleIndex]->Update();
}

void vtkPolyLineWidget::Translate(double *p1, double *p2)
{
  // Get the motion vector
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  double newCtr[3];
  for (int i = 0; i< this->NumberOfHandles; i++)
    {
    double* ctr =  this->HandleGeometry[i]->GetCenter();
    for (int j=0; j<3; j++)
      {
      newCtr[j] = ctr[j] + v[j];
      }
     this->HandleGeometry[i]->SetCenter(newCtr);
     this->HandleGeometry[i]->Update();
    }
}

void vtkPolyLineWidget::Scale(double *p1, double *p2, int vtkNotUsed(X), int Y)
{
  // Get the motion vector
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  double center[3] = {0.0,0.0,0.0};
  double avgdist = 0.0;
  double *prevctr = this->HandleGeometry[0]->GetCenter();
  double *ctr;

  center[0] += prevctr[0];
  center[1] += prevctr[1];
  center[2] += prevctr[2];

  int i;
  for (i = 1; i<this->NumberOfHandles; i++)
    {
    ctr = this->HandleGeometry[i]->GetCenter();
    center[0] += ctr[0];
    center[1] += ctr[1];
    center[2] += ctr[2];
    avgdist += sqrt(vtkMath::Distance2BetweenPoints(ctr,prevctr));
    prevctr = ctr;
    }

  avgdist /= this->NumberOfHandles;

  center[0] /= this->NumberOfHandles;
  center[1] /= this->NumberOfHandles;
  center[2] /= this->NumberOfHandles;

  // Compute the scale factor
  double sf = vtkMath::Norm(v) / avgdist;
  if ( Y > this->Interactor->GetLastEventPosition()[1] )
    {
    sf = 1.0 + sf;
    }
  else
    {
    sf = 1.0 - sf;
    }

  // Move the handle points
  double newCtr[3];
  for (i = 0; i< this->NumberOfHandles; i++)
    {
    ctr = this->HandleGeometry[i]->GetCenter();
    for (int j=0; j<3; j++)
      {
      newCtr[j] = sf * (ctr[j] - center[j]) + center[j];
      }
    this->HandleGeometry[i]->SetCenter(newCtr);
    this->HandleGeometry[i]->Update();
    }
}

void vtkPolyLineWidget::Spin(double *p1, double *p2, double *vpn)
{
  // Mouse motion vector in world space
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  // Axis of rotation
  double axis[3] = {0.0,0.0,0.0};

  if ( this->ProjectToPlane )
    {
    if(this->ProjectionNormal == VTK_PROJECTION_OBLIQUE && this->PlaneSource != NULL)
      {
      double* normal = this->PlaneSource->GetNormal();
      axis[0] = normal[0];
      axis[1] = normal[1];
      axis[2] = normal[2];
      vtkMath::Normalize(axis);
      }
    else
      {
      axis[ this->ProjectionNormal ] = 1.0;
      }
    }
  else
    {
  // Create axis of rotation and angle of rotation
    vtkMath::Cross(vpn,v,axis);
    if ( vtkMath::Normalize(axis) == 0.0 )
      {
      return;
      }
    }

  // Radius vector (from mean center to cursor position)
  double rv[3] = {p2[0] - this->Centroid[0],
                  p2[1] - this->Centroid[1],
                  p2[2] - this->Centroid[2]};

  // Distance between center and cursor location
  double rs = vtkMath::Normalize(rv);

  // Spin direction
  double ax_cross_rv[3];
  vtkMath::Cross(axis,rv,ax_cross_rv);

  // Spin angle
  double theta = 360.0 * vtkMath::Dot(v,ax_cross_rv) / rs;

  // Manipulate the transform to reflect the rotation
  this->Transform->Identity();
  this->Transform->Translate(this->Centroid[0],this->Centroid[1],this->Centroid[2]);
  this->Transform->RotateWXYZ(theta,axis);
  this->Transform->Translate(-this->Centroid[0],-this->Centroid[1],-this->Centroid[2]);

  // Set the handle points
  double newCtr[3];
  double ctr[3];
  for (int i=0; i<this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i]->GetCenter(ctr);
    this->Transform->TransformPoint(ctr,newCtr);
    this->HandleGeometry[i]->SetCenter(newCtr);
    this->HandleGeometry[i]->Update();
    }
}

void vtkPolyLineWidget::CreateDefaultProperties()
{
  if ( ! this->HandleProperty )
    {
    this->HandleProperty = vtkProperty::New();
    this->HandleProperty->SetColor(1,1,1);
    }
  if ( ! this->SelectedHandleProperty )
    {
    this->SelectedHandleProperty = vtkProperty::New();
    this->SelectedHandleProperty->SetColor(1,0,0);
    }

  if ( ! this->LineProperty )
    {
    this->LineProperty = vtkProperty::New();
    this->LineProperty->SetRepresentationToWireframe();
    this->LineProperty->SetAmbient(1.0);
    this->LineProperty->SetColor(1.0,1.0,0.0);
    this->LineProperty->SetLineWidth(2.0);
    }
  if ( ! this->SelectedLineProperty )
    {
    this->SelectedLineProperty = vtkProperty::New();
    this->SelectedLineProperty->SetRepresentationToWireframe();
    this->SelectedLineProperty->SetAmbient(1.0);
    this->SelectedLineProperty->SetAmbientColor(0.0,1.0,0.0);
    this->SelectedLineProperty->SetLineWidth(2.0);
    }
}

void vtkPolyLineWidget::PlaceWidget(double bds[6])
{
  int i;
  double bounds[6], center[3];
  this->AdjustBounds(bds, bounds, center);

  if ( this->ProjectToPlane )
    {
    this->ProjectPointsToPlane();
    }
  else  //place the center
    {
    // Create a default straight line within the data bounds
    double x0 = bounds[0];
    double x1 = bounds[1];
    double y0 = bounds[2];
    double y1 = bounds[3];
    double z0 = bounds[4];
    double z1 = bounds[5];
    double x;
    double y;
    double z;
    double position;
    for (i=0; i<this->NumberOfHandles; i++)
      {
      position = i / (this->NumberOfHandles - 1.0);
      x = (1.0-position)*x0 + position*x1;
      y = (1.0-position)*y0 + position*y1;
      z = (1.0-position)*z0 + position*z1;
      this->HandleGeometry[i]->SetCenter(x,y,z);
      }
    }

  for (i=0; i<6; i++)
    {
    this->InitialBounds[i] = bounds[i];
    }
  this->InitialLength = sqrt((bounds[1]-bounds[0])*(bounds[1]-bounds[0]) +
                             (bounds[3]-bounds[2])*(bounds[3]-bounds[2]) +
                             (bounds[5]-bounds[4])*(bounds[5]-bounds[4]));

  // Re-compute the spline coeffs
  this->BuildRepresentation();
  this->SizeHandles();
}

void vtkPolyLineWidget::SetProjectionPosition(double position)
{
  this->ProjectionPosition = position; 
  if ( this->ProjectToPlane )
    {
    this->ProjectPointsToPlane();
    }
  this->BuildRepresentation();
}

void vtkPolyLineWidget::SetPlaneSource(vtkPlaneSource* plane)
{
  if (this->PlaneSource == plane)
    {
    return;
    }
  this->PlaneSource = plane;
}

void vtkPolyLineWidget::SetNumberOfHandles(int npts)
{

  if (this->NumberOfHandles == npts)
    {
    return;
    }
  if (npts < 2)
    {
    vtkGenericWarningMacro(
      <<"vtkPolyLineWidget: minimum of 2 points required.");
    return;
    }

  // before we NUKE the whole thing, let's store some useful information
  double radius = this->HandleGeometry[0]->GetRadius();
  double factor = (this->NumberOfHandles - 1.0)/(npts - 1.0);
  // we're going to store x, y and z
  int NumberOfStoredHandles = this->NumberOfHandles;
  double *StoredHandles = new double [this->NumberOfHandles * 3];
  int i;
  for (i = 0; i < this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i]->GetCenter(StoredHandles + i * 3);
    }
  
  
  this->Initialize();

  this->NumberOfHandles = npts;

  // Create the handles
  this->Handle         = new vtkActor* [this->NumberOfHandles];
  this->HandleMapper   = new vtkPolyDataMapper* [this->NumberOfHandles];
  this->HandleGeometry = new vtkSphereSource* [this->NumberOfHandles];

  double x,y,z;
  double px,py,pz;
  double nx,ny,nz;

  int numIntervals = NumberOfStoredHandles;
//  if (!this->Closed)
//    {
//    numIntervals -= 1;
//    }
  int numHandlesPerInterval =
    (this->NumberOfHandles) / numIntervals;

//  int numNewHandlesPerInterval =
//    (this->NumberOfHandles - NumberOfStoredHandles) / numIntervals;
  
  // insert and interpolate new points
  for (i=0; i<this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i] = vtkSphereSource::New();
    this->HandleGeometry[i]->SetThetaResolution(16);
    this->HandleGeometry[i]->SetPhiResolution(8);
    this->HandleMapper[i] = vtkPolyDataMapper::New();
    this->HandleMapper[i]->SetInput(this->HandleGeometry[i]->GetOutput());
    this->Handle[i] = vtkActor::New();
    this->Handle[i]->SetMapper(this->HandleMapper[i]);
    this->Handle[i]->SetProperty(this->HandleProperty);

    // let's find out what the coord should be
    if (this->NumberOfHandles <= NumberOfStoredHandles)
      {
      // numIntervals is negative, that means there are now less
      // points than their are stored, so we can just use the old ones
      // up to the number of new ones required
      x = StoredHandles[i*3];
      y = StoredHandles[i*3 + 1];
      z = StoredHandles[i*3 + 2];
      }
    else
      {
      int oldIndex = i / numHandlesPerInterval;      
      if (i % numHandlesPerInterval == 0 && oldIndex < NumberOfStoredHandles)
        {
        // we are on one of the old points (i.e. don't interpolate)

        x = StoredHandles[oldIndex*3];
        y = StoredHandles[oldIndex*3 + 1];
        z = StoredHandles[oldIndex*3 + 2];
        }
      else
        {
        // we are between old nodes, so we have to interpolate
        // find indices in stored nodes
        int pi = i / numHandlesPerInterval;
        int ni = pi + 1;

        if (ni >= NumberOfStoredHandles - 1)
          {
          // this means we've gone past the end, so:
          if (this->Closed)
            {
            // we're closed, so the next point is the beginning
            pi = NumberOfStoredHandles - 1;
            ni = 0;
            }
          else
            {
            // we're open, so we move one interval back
            
            ni = NumberOfStoredHandles -1;
            pi = ni -1;
            }
          }

        // distance along
        double d = i / (double)numHandlesPerInterval - (double)pi;
        //
        px = StoredHandles[pi*3];
        py = StoredHandles[pi*3 + 1];
        pz = StoredHandles[pi*3 + 2];
        nx = StoredHandles[ni*3];
        ny = StoredHandles[ni*3 + 1];
        nz = StoredHandles[ni*3 + 2];
        x = px + d * (nx - px);
        y = py + d * (ny - py);
        z = pz + d * (nz - pz);        
        }
      } // had to add more points
    
    //cout << "xyz: " << x << "," << y << "," << z << endl;

    this->HandleGeometry[i]->SetCenter(x,y,z);
    this->HandleGeometry[i]->SetRadius(radius);
    this->HandlePicker->AddPickList(this->Handle[i]);
    }
  
  //cout << this->NumberOfHandles << endl;

  // take care of this mem on the heap
  delete[] StoredHandles;
      
  this->BuildRepresentation();

  if ( this->Interactor )
    {
    if (!this->CurrentRenderer)
      {
      this->SetCurrentRenderer(this->Interactor->FindPokedRenderer(
        this->Interactor->GetLastEventPosition()[0],
        this->Interactor->GetLastEventPosition()[1]));
      }
    if (this->CurrentRenderer != NULL)
      {
      for (i=0; i<this->NumberOfHandles; i++)
        {
        this->CurrentRenderer->AddViewProp(this->Handle[i]);
        }
      }
      this->Interactor->Render();
    }
}

void vtkPolyLineWidget::Initialize(void)
{
  int i;
  if ( this->Interactor )
    {
    if (!this->CurrentRenderer)
      {
      this->SetCurrentRenderer(this->Interactor->FindPokedRenderer(
        this->Interactor->GetLastEventPosition()[0],
        this->Interactor->GetLastEventPosition()[1]));
      }
    if ( this->CurrentRenderer != NULL)
      {
      for (i=0; i<this->NumberOfHandles; i++)
        {
        this->CurrentRenderer->RemoveViewProp(this->Handle[i]);
        }
      }
    }

  for (i=0; i<this->NumberOfHandles; i++)
    {
    this->HandlePicker->DeletePickList(this->Handle[i]);
    this->HandleGeometry[i]->Delete();
    this->HandleMapper[i]->Delete();
    this->Handle[i]->Delete();
    }

  this->NumberOfHandles = 0;

  delete [] this->Handle;
  delete [] this->HandleMapper;
  delete [] this->HandleGeometry;
}

void vtkPolyLineWidget::GetPolyData(vtkPolyData *pd)
{
  pd->ShallowCopy(this->LineData);
}

void vtkPolyLineWidget::SizeHandles()
{
  double radius = this->vtk3DWidget::SizeHandles(1.0);
  for (int i=0; i<this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i]->SetRadius(radius);
    }
}

double vtkPolyLineWidget::GetSummedLength()
{
  vtkPoints* points = this->LineData->GetPoints();
  int npts = points->GetNumberOfPoints();

  if (npts < 2) { return 0.0; }

  double a[3];
  double b[3];
  double sum = 0.0;
  int i = 0;
  points->GetPoint(i,a);
  int imax = (npts%2 == 0)?npts-2:npts-1;

  while(i<imax)
    {
    points->GetPoint(i+1,b);
    sum = sum + sqrt(vtkMath::Distance2BetweenPoints(a,b));
    i = i + 2;
    points->GetPoint(i,a);
    sum = sum + sqrt(vtkMath::Distance2BetweenPoints(a,b));
    }

  if(npts%2 == 0)
    {
    points->GetPoint(i+1,b);
    sum = sum + sqrt(vtkMath::Distance2BetweenPoints(a,b));
    }

  return sum;
}

void vtkPolyLineWidget::CalculateCentroid()
{
  this->Centroid[0] = 0.0;
  this->Centroid[1] = 0.0;
  this->Centroid[2] = 0.0;

  double ctr[3];
  for (int i = 0; i<this->NumberOfHandles; i++)
    {
    this->HandleGeometry[i]->GetCenter(ctr);
    this->Centroid[0] += ctr[0];
    this->Centroid[1] += ctr[1];
    this->Centroid[2] += ctr[2];
    }

  this->Centroid[0] /= this->NumberOfHandles;
  this->Centroid[1] /= this->NumberOfHandles;
  this->Centroid[2] /= this->NumberOfHandles;
}

