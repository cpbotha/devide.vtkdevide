#include "vtkBoxWidgetConstrained.h"

#include "vtkActor.h"
#include "vtkAssemblyNode.h"
#include "vtkAssemblyPath.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellPicker.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPlanes.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"

vtkCxxRevisionMacro(vtkBoxWidgetConstrained, "$Revision: 1.7 $");
vtkStandardNewMacro(vtkBoxWidgetConstrained);

vtkBoxWidgetConstrained::vtkBoxWidgetConstrained() : vtkBoxWidget()
{
  // default: don't constrain
  this->SetConstraintType(0);
  // vector is set to some default value
  this->SetConstraintVector(0.0, 0.0, 0.0);
}

// Loop through all points and translate them
void vtkBoxWidgetConstrained::Translate(double *p1, double *p2)
{
  double *pts = ((vtkDoubleArray *)this->Points->GetData())->GetPointer(0);
  double v[3];

  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  // constrained to line  
  if (this->GetConstraintType() == 1)
    {
    // calculate the component of the motion vector that is collinear
    // with the constraintvector
    double goodMagnitude = vtkMath::Dot(v, this->ConstraintVector);
    for (int i = 0; i < 3; i++)
      {
      // that will be the new motion vector!
      v[i] = goodMagnitude * this->ConstraintVector[i];
      }

    }
  if (this->GetConstraintType() == 2)
    {
    // calculate the component of the motion vector that is collinear
    // with the constraintvector
    double bdp = vtkMath::Dot(v, this->ConstraintVector);
    double badvector[3];
    for (int i=0; i<3; i++)
      {
      badvector[i] = bdp * this->ConstraintVector[i];
      // this has to be substracted from the motion vector, as we
      // only tolerate motion parallel to the constraintplane
      v[i] -= badvector[i];
      }
    }

  // Move the corners
  for (int i=0; i<8; i++)
    {
    *pts++ += v[0];
    *pts++ += v[1];
    *pts++ += v[2];
    }
  this->PositionHandles();
}

void vtkBoxWidgetConstrained::Rotate(int X, int Y, double *p1, double *p2, double *vpn)
{
  double *pts = ((vtkDoubleArray *)this->Points->GetData())->GetPointer(0);
  double *center = ((vtkDoubleArray *)this->Points->GetData())->GetPointer(3*14);
  double v[3]; //vector of motion
  double axis[3]; //axis of rotation
  double theta; //rotation angle
  int i;

  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  // if we're constrained, constrain man!
  if (this->ConstraintType)
    {
    double blaat[3];
    vtkMath::Cross(vpn, v, blaat);

    double bdp = vtkMath::Dot(v, this->ConstraintVector);
    double badvector[3];
    for (int i=0; i<3; i++)
      {
      badvector[i] = bdp * this->ConstraintVector[i];
      v[i] -= badvector[i];
      axis[i] = this->ConstraintVector[i];
      }

	
    if (vtkMath::Dot(axis, blaat) < 0)
      {
      for (int i=0; i<3; i++)
        {
        axis[i] *= -1;
        }
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

  int *size = this->CurrentRenderer->GetSize();
  double l2 = (X-this->Interactor->GetLastEventPosition()[0])*(X-this->Interactor->GetLastEventPosition()[0]) + (Y-this->Interactor->GetLastEventPosition()[1])*(Y-this->Interactor->GetLastEventPosition()[1]);
  theta = 360.0 * sqrt(l2/((double)size[0]*size[0]+size[1]*size[1]));

  //Manipulate the transform to reflect the rotation
  this->Transform->Identity();
  this->Transform->Translate(center[0],center[1],center[2]);
  this->Transform->RotateWXYZ(theta,axis);
  this->Transform->Translate(-center[0],-center[1],-center[2]);

  //Set the corners
  vtkPoints *newPts = vtkPoints::New(VTK_DOUBLE);
  this->Transform->TransformPoints(this->Points,newPts);

  for (i=0; i<8; i++, pts+=3)
    {
    this->Points->SetPoint(i, newPts->GetPoint(i));
    }

  newPts->Delete();
  this->PositionHandles();
}
  
