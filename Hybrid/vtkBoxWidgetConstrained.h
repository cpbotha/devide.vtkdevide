#ifndef __vtkBoxWidgetConstrained_h
#define __vtkBoxWidgetConstrained_h

#include "vtk3DWidget.h"

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

class VTK_HYBRID_EXPORT vtkBoxWidgetConstrained : public vtk3DWidget
{
public:
  // Description:
  // Instantiate the object.
  static vtkBoxWidgetConstrained *New();

  vtkTypeRevisionMacro(vtkBoxWidgetConstrained,vtk3DWidget);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Methods that satisfy the superclass' API.
  virtual void SetEnabled(int);
  virtual void PlaceWidget(float bounds[6]);
  void PlaceWidget()
    {this->Superclass::PlaceWidget();}
  void PlaceWidget(float xmin, float xmax, float ymin, float ymax, 
                   float zmin, float zmax)
    {this->Superclass::PlaceWidget(xmin,xmax,ymin,ymax,zmin,zmax);}


  // Description:
  // Get the planes describing the implicit function defined by the box
  // widget. The user must provide the instance of the class vtkPlanes. Note
  // that vtkPlanes is a subclass of vtkImplicitFunction, meaning that it can
  // be used by a variety of filters to perform clipping, cutting, and
  // selection of data.  (The direction of the normals of the planes can be
  // reversed enabling the InsideOut flag.)
  void GetPlanes(vtkPlanes *planes);

  // Description:
  // Set/Get the InsideOut flag. When off, the normals point out of the
  // box. When on, the normals point into the hexahedron.  InsideOut
  // is off by default.
  vtkSetMacro(InsideOut,int);
  vtkGetMacro(InsideOut,int);
  vtkBooleanMacro(InsideOut,int);

  // Description:
  // Retrieve a linear transform characterizing the transformation of the
  // box. Note that the transformation is relative to where PlaceWidget
  // was initially called. This method modifies the transform provided. The
  // transform can be used to control the position of vtkProp3D's, as well as
  // other transformation operations (e.g., vtkTranformPolyData).
  void GetTransform(vtkTransform *t);

  // Description:
  // Grab the polydata (including points) that define the box widget. The
  // polydata consists of 6 quadrilateral faces and 15 points. The first
  // eight points define the eight corner vertices; the next six define the
  // -x,+x, -y,+y, -z,+z face points; and the final point (the 15th out of 15
  // points) defines the center of the hexahedron. These point values are
  // guaranteed to be up-to-date when either the InteractionEvent or
  // EndInteractionEvent events are invoked. The user provides the
  // vtkPolyData and the points and cells are added to it.
  void GetPolyData(vtkPolyData *pd);

  // Description:
  // Get the handle properties (the little balls are the handles). The 
  // properties of the handles when selected and normal can be 
  // set.
  vtkGetObjectMacro(HandleProperty,vtkProperty);
  vtkGetObjectMacro(SelectedHandleProperty,vtkProperty);
  
  // Description:
  // Get the face properties (the faces of the box). The 
  // properties of the face when selected and normal can be 
  // set.
  vtkGetObjectMacro(FaceProperty,vtkProperty);
  vtkGetObjectMacro(SelectedFaceProperty,vtkProperty);
  
  // Description:
  // Get the outline properties (the outline of the box). The 
  // properties of the outline when selected and normal can be 
  // set.
  vtkGetObjectMacro(OutlineProperty,vtkProperty);
  vtkGetObjectMacro(SelectedOutlineProperty,vtkProperty);
  
  // Description:
  // Control the representation of the outline. This flag enables
  // face wires. By default face wires are off.
  vtkSetMacro(OutlineFaceWires,int);
  vtkGetMacro(OutlineFaceWires,int);
  vtkBooleanMacro(OutlineFaceWires,int);

  // Description:
  // Control the representation of the outline. This flag enables
  // the cursor lines running between the handles. By default cursor
  // wires are on.
  vtkSetMacro(OutlineCursorWires,int);
  vtkGetMacro(OutlineCursorWires,int);
  vtkBooleanMacro(OutlineCursorWires,int);

  // Description:
  // Control the behavior of the widget. Translation, rotation, and
  // scaling can all be enabled and disabled.
  vtkSetMacro(TranslationEnabled,int);
  vtkGetMacro(TranslationEnabled,int);
  vtkBooleanMacro(TranslationEnabled,int);
  vtkSetMacro(ScalingEnabled,int);
  vtkGetMacro(ScalingEnabled,int);
  vtkBooleanMacro(ScalingEnabled,int);
  vtkSetMacro(RotationEnabled,int);
  vtkGetMacro(RotationEnabled,int);
  vtkBooleanMacro(RotationEnabled,int);

protected:
  vtkBoxWidgetConstrained();
  ~vtkBoxWidgetConstrained();

//BTX - manage the state of the widget
  int State;
  enum WidgetState
  {
    Start=0,
    Moving,
    Scaling,
    Outside
  };
//ETX
    
  // Handles the events
  static void ProcessEvents(vtkObject* object, 
                            unsigned long event,
                            void* clientdata, 
                            void* calldata);

  // ProcessEvents() dispatches to these methods.
  void OnMouseMove();
  void OnLeftButtonDown();
  void OnLeftButtonUp();
  void OnMiddleButtonDown();
  void OnMiddleButtonUp();
  void OnRightButtonDown();
  void OnRightButtonUp();
  
  // the hexahedron (6 faces)
  vtkActor          *HexActor;
  vtkPolyDataMapper *HexMapper;
  vtkPolyData       *HexPolyData;
  vtkPoints         *Points;  //used by others as well
  double             N[6][3]; //the normals of the faces

  // A face of the hexahedron
  vtkActor          *HexFace;
  vtkPolyDataMapper *HexFaceMapper;
  vtkPolyData       *HexFacePolyData;

  // glyphs representing hot spots (e.g., handles)
  vtkActor          **Handle;
  vtkPolyDataMapper **HandleMapper;
  vtkSphereSource   **HandleGeometry;
  void PositionHandles();
  void HandlesOn(double length);
  void HandlesOff();
  int HighlightHandle(vtkProp *prop); //returns cell id
  void HighlightFace(int cellId);
  void HighlightOutline(int highlight);
  void ComputeNormals();
  virtual void SizeHandles();
  
  // wireframe outline
  vtkActor          *HexOutline;
  vtkPolyDataMapper *OutlineMapper;
  vtkPolyData       *OutlinePolyData;

  // Do the picking
  vtkCellPicker *HandlePicker;
  vtkCellPicker *HexPicker;
  vtkActor *CurrentHandle;
  int      CurrentHexFace;
  
  // Methods to manipulate the hexahedron.
  void Translate(double *p1, double *p2);
  void Scale(double *p1, double *p2, int X, int Y);
  void Rotate(int X, int Y, double *p1, double *p2, double *vpn);
  void MovePlusXFace(double *p1, double *p2);
  void MoveMinusXFace(double *p1, double *p2);
  void MovePlusYFace(double *p1, double *p2);
  void MoveMinusYFace(double *p1, double *p2);
  void MovePlusZFace(double *p1, double *p2);
  void MoveMinusZFace(double *p1, double *p2);

  void MoveFace(double *p1, double *p2, double *h1, double *h2, 
                double *x1, double *x2, double *x3, double *x4,
                double *x5);
  
  // Transform the hexahedral points (used for rotations)
  vtkTransform *Transform;
  
  // Properties used to control the appearance of selected objects and
  // the manipulator in general.
  vtkProperty *HandleProperty;
  vtkProperty *SelectedHandleProperty;
  vtkProperty *FaceProperty;
  vtkProperty *SelectedFaceProperty;
  vtkProperty *OutlineProperty;
  vtkProperty *SelectedOutlineProperty;
  void CreateDefaultProperties();
  
  // Control the orientation of the normals
  int InsideOut;
  int OutlineFaceWires;
  int OutlineCursorWires;
  void GenerateOutline();
  
  // Control whether scaling, rotation, and translation are supported
  int TranslationEnabled;
  int ScalingEnabled;
  int RotationEnabled;

private:
  vtkBoxWidgetConstrained(const vtkBoxWidgetConstrained&);  //Not implemented
  void operator=(const vtkBoxWidgetConstrained&);  //Not implemented
};

#endif
