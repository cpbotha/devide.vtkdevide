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

vtkCxxRevisionMacro(vtkBoxWidgetConstrained, "$Revision: 1.11 $");
vtkStandardNewMacro(vtkBoxWidgetConstrained);

vtkBoxWidgetConstrained::vtkBoxWidgetConstrained() : vtkBoxWidget()
{
  // default: don't constrain
  this->SetConstraintType(0);
  // vector is set to some default value
  this->SetConstraintVector(0.0, 0.0, 0.0);
}

void vtkBoxWidgetConstrained::SetConstraintVector(double v0, double v1, double v2)
{
  this->ConstraintVector[0] = v0;
  this->ConstraintVector[1] = v1;
  this->ConstraintVector[2] = v2;

  // make sure it's normalised
  vtkMath::Normalize(this->ConstraintVector);
  this->Modified();
}

void vtkBoxWidgetConstrained::SetConstraintVector(double data[])
{
  this->SetConstraintVector(data[0], data[1], data[2]);
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
  else if (this->GetConstraintType() == 2)
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

  // constrained to plane OR line
  if (this->ConstraintType)
    {
    // by definition, the axis is the ConstraintVector
    // we also flatten v to be in the plane
    // in the case of the line constraint, the line is the axis of rotation
    // in the case of the plane constraint, the plane normal is rot. axis
    // i.e. in both cases the ConstraintVector
    double bdp = vtkMath::Dot(v, this->ConstraintVector);
    for (i=0; i<3; i++)
      {
      v[i] -= this->ConstraintVector[i] * bdp;
      axis[i] = this->ConstraintVector[i];
      }

    // now calculate cross-product of our axis and v
    double avc[3];
    vtkMath::Cross(axis, v, avc);
    vtkMath::Normalize(avc);
    // this cross-product (along with center) defines two hemispheres
    // the hemisphere determines the correct rotation direction
    // so let's determine the hemisphere of p2
    double p2v[3];
    for (i = 0; i < 3; i++)
      {
      p2v[i] = p2[i] - center[i];
      }
    // if it's incorrect, flip the rotation axis and all is well!
    // this is simple code, but it did take some thinking :)
    double hdp = vtkMath::Dot(p2v, avc);
    if (hdp >= 0)
      {
      for (i = 0; i < 3; i++)
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
  double l2 = (X-this->Interactor->GetLastEventPosition()[0])*(X-this->Interactor->GetLastEventPosition()[0]) + 
    (Y-this->Interactor->GetLastEventPosition()[1])*(Y-this->Interactor->GetLastEventPosition()[1]);
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
  
