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

  vtkSetMacro(ConstraintType, int);
  vtkGetMacro(ConstraintType, int);

  vtkSetVector3Macro(ConstraintVector, double);
  vtkGetVector3Macro(ConstraintVector, double);

protected:
  vtkBoxWidgetConstrained();

  // Methods to manipulate the hexahedron.
  virtual void Translate(double *p1, double *p2);
  virtual void Rotate(int X, int Y, double *p1, double *p2, double *vpn);

  /// ConstraintType determines how exactly the box motion will be constrained.
  /// 0 - no constraint
  /// 1 - constrain to line
  /// 2 - constrain to plane
  int ConstraintType;
  /// ConstraintVector defines the constraint geometry
  /// ConstraintType == 1: ConstraintVector is a line
  /// ConstraintType == 2: ConstraintVector is a plane normal
  double ConstraintVector[3];

private:
  vtkBoxWidgetConstrained(const vtkBoxWidgetConstrained&);  //Not implemented
  void operator=(const vtkBoxWidgetConstrained&);  //Not implemented
};

#endif
