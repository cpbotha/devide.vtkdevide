// .NAME vtkImageBorderMask - Simple example of an image-image filter.
// .SECTION Description
// This is an example of a simple image-image filter. It copies it's input
// to it's output (point by point). It shows how templates can be used
// to support various data types.
// .SECTION See also
// vtkImageBorderMask

#ifndef __vtkImageBorderMask_h
#define __vtkImageBorderMask_h

#include "vtkSimpleImageToImageFilter.h"
#include "vtkdevideImagingWin32Header.h"

class VTK_DEVIDE_IMAGING_EXPORT vtkImageBorderMask :
  public vtkSimpleImageToImageFilter
{
public:
  static vtkImageBorderMask *New();
  vtkTypeRevisionMacro(vtkImageBorderMask,vtkSimpleImageToImageFilter);

  vtkSetMacro(BorderValue, int);
  vtkGetMacro(BorderValue, int);

  vtkSetMacro(InteriorValue, int);
  vtkGetMacro(InteriorValue, int);

  // 0 - use BorderValue for Border
  // 1 - use input for Border
  // 2 - use min(input) for Border
  // 3 - use max(input) for Border
  vtkSetMacro(BorderMode, int);
  vtkGetMacro(BorderMode, int);

  // 0 - use InteriorValue for Interior
  // 1 - use input for Border
  // 2 - use min(input) for Border
  // 3 - use max(input) for Border
  vtkSetMacro(InteriorMode, int);
  vtkGetMacro(InteriorMode, int);

  vtkSetVector3Macro(Borders, int);
  vtkGetVector3Macro(Borders, int)

protected:
  vtkImageBorderMask();
  ~vtkImageBorderMask() {};

  int BorderValue;
  int BorderMode;
  int Borders[3];
  
  int InteriorValue;
  int InteriorMode;
  

  virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);
  
private:
  vtkImageBorderMask(const vtkImageBorderMask&);  // Not implemented.
  void operator=(const vtkImageBorderMask&);  // Not implemented.
};

#endif







