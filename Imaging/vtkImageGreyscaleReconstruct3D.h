
#ifndef __vtkImageGreyscaleReconstruct3D_h
#define __vtkImageGreyscaleReconstruct3D_h

#include "vtkImageTwoInputFilter.h"
#include "vtkdevideImagingWin32Header.h"

/**
 * Perform greyscale reconstruction of mask I from marker (seed) J.  I
 * is first input, J is the second.
 */

class VTK_DEVIDE_IMAGING_EXPORT vtkImageGreyscaleReconstruct3D :
public vtkImageTwoInputFilter
{
public:
  static vtkImageGreyscaleReconstruct3D *New();
  vtkTypeRevisionMacro(vtkImageGreyscaleReconstruct3D,vtkImageTwoInputFilter);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  vtkSetMacro(Dual, int);
  vtkGetMacro(Dual, int);
  
protected:
  vtkImageGreyscaleReconstruct3D();
  ~vtkImageGreyscaleReconstruct3D() {};

  int Dual;
  
  void ExecuteInformation(vtkImageData **inDatas, vtkImageData *outData);
  virtual void ComputeInputUpdateExtent(int inExt[6], int outExt[6],
                                        int whichInput);
  // this is called by the SuperClass - we just chain to the other EI
  void ExecuteInformation(){this->vtkImageTwoInputFilter::ExecuteInformation();};
  
private:
  void ExecuteData(vtkDataObject *out);
  // Not implemented.  
  vtkImageGreyscaleReconstruct3D(const vtkImageGreyscaleReconstruct3D&);  
  void operator=(const vtkImageGreyscaleReconstruct3D&);  // Not implemented.
};

#endif



