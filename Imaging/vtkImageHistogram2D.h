
#ifndef __vtkImageHistogram2D_h
#define __vtkImageHistogram2D_h

#include "vtkImageTwoInputFilter.h"
#include "vtkdevideImagingWin32Header.h"

class VTK_DEVIDE_IMAGING_EXPORT vtkImageHistogram2D :
public vtkImageTwoInputFilter
{
public:
  static vtkImageHistogram2D *New();
  vtkTypeRevisionMacro(vtkImageHistogram2D,vtkImageTwoInputFilter);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  vtkSetMacro(Input1Bins,int);
  vtkGetMacro(Input1Bins,int);

  vtkSetMacro(Input2Bins,int);
  vtkGetMacro(Input2Bins,int);
  
  
protected:
  vtkImageHistogram2D();
  ~vtkImageHistogram2D() {};

  int Input1Bins;
  int Input2Bins;
  
  void ExecuteInformation(vtkImageData **inDatas, vtkImageData *outData);
  virtual void ComputeInputUpdateExtent(int inExt[6], int outExt[6],
                                        int whichInput);
  // this is called by the SuperClass - we just chain to the other EI
  void ExecuteInformation(){this->vtkImageTwoInputFilter::ExecuteInformation();};
  
private:
  void ExecuteData(vtkDataObject *out);
  
  vtkImageHistogram2D(const vtkImageHistogram2D&);  // Not implemented.
  void operator=(const vtkImageHistogram2D&);  // Not implemented.
};

#endif



