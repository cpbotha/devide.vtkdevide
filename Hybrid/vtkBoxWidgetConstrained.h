#ifndef __vtkBoxWidgetConstrained_h
#define __vtkBoxWidgetConstrained_h

#include "vtkBoxWidget.h"
#include "vtkdscasHybridWin32Header.h"

class vtkActor;
class vtkCellPicker;
class vtkPlanes;
class vtkPoints;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkProp;
class vtkProperty;
class vtkSphereSource;
class vtkTransform;

class VTK_DSCAS_HYBRID_EXPORT vtkBoxWidgetConstrained : public vtkBoxWidget
{
public:
  // Description:
  // Instantiate the object.
  static vtkBoxWidgetConstrained *New();
  vtkTypeRevisionMacro(vtkBoxWidgetConstrained,vtkBoxWidget);

  // BEGIN cpbotha
  // are we constrained to a plane?
  vtkSetMacro(ConstrainToPlane, int);
  vtkGetMacro(ConstrainToPlane, int);
  vtkBooleanMacro(ConstrainToPlane, int);

  vtkSetVector3Macro(ConstrainPlaneNormal, double);
  vtkGetVector3Macro(ConstrainPlaneNormal, double);
  // END cpbotha

protected:
  vtkBoxWidgetConstrained();

  // Methods to manipulate the hexahedron.
  virtual void Translate(double *p1, double *p2);
  virtual void Rotate(int X, int Y, double *p1, double *p2, double *vpn);
  // BEGIN cpbotha
  // Set this if you want motion to be constrained to the given plane
  int ConstrainToPlane;
  double ConstrainPlaneNormal[3];
  // END cpbotha

private:
  vtkBoxWidgetConstrained(const vtkBoxWidgetConstrained&);  //Not implemented
  void operator=(const vtkBoxWidgetConstrained&);  //Not implemented
};

#endif
